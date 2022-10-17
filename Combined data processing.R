#author: Joris Jean van der Lugt
#date: 05-08-2022
#Combined data processing of Human Cinnamaldehyde pbk models
#Before running code first run all model associated with this data processing
#Human desolve model
#human rxode model
#Human population model

library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(shiny)
library(truncnorm)
library(reshape2)
library(plotly)
library(PKNCA)
library(MESS)

#Rxode
#cinnamaldehyde model Human
#Mass balance calculation rxode inhalation complete
mass_df <- solve.pbk_nonpop/BW * MW /1e+3
mass_df <- mass_df[,c(67:80,82,83,86:91,95:99)]
mass_at_t <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])


#cinnamaldehyde model Rat
#Mass balance calculation rxode inhalation complete
mass_df <- solve.pbk_rat/BW * MW /1e+3
mass_df <- mass_df[,c(66:84,86,88:92,94,96.97)]
mass_at_t <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])


#Rxode 
#population mass balance
mass_df <- solve.pbk_popgen/phys[1,3] * MW /1e+3
mass_df <- mass_df[1:81,c(71:85,87,88,91:98,102:106)]
mass_at_t <- data.frame(mass=as.numeric())


for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])


#Non inhalation mass balance 
mass_df <- solve.pbk/BW * MW /1e+3
mass_df <- mass_df[,c(59:68,70,71,74:81,85:89)]
mass_at_t <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])

#non inhalation Rat
mass_df <- solve.pbk_rat/BW * MW /1e+3
mass_df <- mass_df[,c(57:71,73,75,77:80,82,84)]
mass_at_t <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])



#Rxode
#AUC calculations
#PKNCA
#After running model
#Extracting organ concentration, time and sim-id from simulation results

#Lung compartment concentration
sub_set <- solve.pbk_popgen[1:162000,c(1,2,12)]

#Blood compartment concentration
sub_set <- solve.pbk_popgen[1:162000,c(1,2,4,5)]

#Combining Arterial and venous blood into one general blood compartment
sub_set[3]<- sub_set[3]+sub_set[4]

#Dropping column as data is already added to column 3
sub_set[4] <- NULL
#Renamning the columns
colnames(sub_set) <- c("id","time","C_B")

#Liver compartment concentration
sub_set <- solve.pbk_popgen[1:162000,c(1,2,50)]

#Slowly perfused compartment
sub_set <- solve.pbk_popgen[1:162000,c(1,2,30)]

#richly perfused compartment
sub_set <- solve.pbk_popgen[1:162000,c(1,2,24)]

#Fat compartment
sub_set <- solve.pbk_popgen[1:162000,c(1,2,18)]

#sI compartment
sub_set <- solve.pbk_popgen[1:162000,c(1,2,36)]




conc_C <- PKNCAconc(sub_set, C_SI~time|id)


#whole dataset
#AUC_data <-PKNCAconc(solve.pbk_popgen, C_L~time|id)

#Dosing data per subject is part of the parameter file but it is missing sim id and the time variable 
#these will be added.
dose_extraction <- as.data.frame(parameters[,64])
sim_extraction <- unique(solve.pbk_popgen[solve.pbk_popgen$time == 0,c("time", "id")])

#Combining into 1 file that can be used
d_dose <- cbind(sim_extraction,dose_extraction$`parameters[, 64]`)                       
d_dose <- set_names(d_dose, c("time","id","dose"))                        

#d_dose for 2000 results is to big for laptop so to see if it works smaller sample will be used
#d_dose <- d_dose[1:400,]
dose_obj <- PKNCAdose(d_dose, dose~time|id)

#Setting the end of the auc calculation at 8 hours
intervals_manual <- data.frame(start=0,
                               end=8,
                               cmax=TRUE,
                               tmax=TRUE,
                               aucinf.obs=FALSE,
                               auclast=TRUE)

#Use this when evaluating the whole dataset 
#data_obj_manual <- PKNCAdata(AUC_data, dose_obj,
#intervals=intervals_manual)

data_obj_manual<- PKNCAdata(conc_C, dose_obj,
                            intervals=intervals_manual)

#letting pknc chose the end time of the auc calc
#data_obj_automatic <- PKNCAdata(conc_C_Pu, dose_obj)

#Computing the data both manual and automatic
#results_obj_automatic <- pk.nca(data_obj_automatic)
results_obj_manual <- pk.nca(data_obj_manual)

#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_manual)))

summary(results_obj_manual)



#Wrting the results into a CSV file 
write.csv(results_obj_manual$result, "results_250_inhalation_C_SI.csv")


#C_Pu data extraction
#Loading data back in
results_250_oral_C_Pu <- read_csv("results_250_oral_C_Pu.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_Oral_C_Pu <- unique(results_250_oral_C_Pu[results_250_oral_C_Pu$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_Pu <- unique(results_250_oral_C_Pu[results_250_oral_C_Pu$PPTESTCD == "cmax",c("id","PPORRES")])


#Lung AUC
#Loading data back in
results_250_oral_C_Pu <- read_csv("results_250_oral_C_Pu.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_Oral__C_Pu <- unique(results_250_oral_C_Pu[results_250_oral_C_Pu$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_Pu <- unique(results_250_oral_C_Pu[results_250_oral_C_Pu$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_Oral__C_Pu$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_Oral__C_Pu$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_Oral__C_Pu$PPORRES[1001:2000])




#---------------BLOOD------------------#
#Loading data back in
results_250_oral_C_B <- read_csv("results_250_oral_C_B.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_Oral__C_B <- unique(results_250_oral_C_B[results_250_oral_C_B$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_B <- unique(results_250_oral_C_B[results_250_oral_C_B$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_Oral__C_B$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_Oral__C_B$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_Oral__C_B$PPORRES[1001:2000])




#--------------------FAT--------------------------#
#Loading data back in
results_250_oral_C_F <- read_csv("results_250_oral_C_F.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_Oral__C_F <- unique(results_250_oral_C_F[results_250_oral_C_F$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_F <- unique(results_250_oral_C_F[results_250_oral_C_F$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_Oral__C_F$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_Oral__C_F$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_Oral__C_F$PPORRES[1001:2000])




#----------------------Richly perfused--------------------------#
#Loading data back in
results_250_oral_C_RP <- read_csv("results_250_oral_C_RP.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_Oral__C_RP <- unique(results_250_oral_C_RP[results_250_oral_C_RP$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_RP <- unique(results_250_oral_C_RP[results_250_oral_C_RP$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_Oral__C_RP$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_Oral__C_RP$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_Oral__C_RP$PPORRES[1001:2000])


#-----------------------------Slowly perfused----------------------#
#Loading data back in
results_250_oral_C_SP <- read_csv("results_250_oral_C_SP.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_Oral__C_SP <- unique(results_250_oral_C_SP[results_250_oral_C_SP$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_SP <- unique(results_250_oral_C_SP[results_250_oral_C_SP$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_Oral__C_SP$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_Oral__C_SP$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_Oral__C_SP$PPORRES[1001:2000])



#-----------------------------Small intestine----------------------#
#Loading data back in
results_250_oral_C_SI <- read_csv("results_250_oral_C_SI.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_Oral__C_SI <- unique(results_250_oral_C_SI[results_250_oral_C_SI$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_SI <- unique(results_250_oral_C_SI[results_250_oral_C_SI$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_Oral__C_SI$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_Oral__C_SI$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_Oral__C_SI$PPORRES[1001:2000])



#-----------------------------Liver----------------------#
#Loading data back in
results_250_oral_C_L <- read_csv("results_250_oral_C_L.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_Oral__C_L <- unique(results_250_oral_C_L[results_250_oral_C_L$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_L <- unique(results_250_oral_C_L[results_250_oral_C_L$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_Oral__C_L$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_Oral__C_L$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_Oral__C_L$PPORRES[1001:2000])


boxplot(AUC_extraction_250_Oral__C_Pu$PPORRES[1:1000], AUC_extraction_250_Oral__C_Pu$PPORRES[1001:2000],AUC_extraction_250_Oral__C_B$PPORRES[1:1000],AUC_extraction_250_Oral__C_B$PPORRES[1001:2000],AUC_extraction_250_Oral__C_F$PPORRES[1:1000],AUC_extraction_250_Oral__C_F$PPORRES[1001:2000],
        AUC_extraction_250_Oral__C_RP$PPORRES[1:1000],AUC_extraction_250_Oral__C_RP$PPORRES[1001:2000],AUC_extraction_250_Oral__C_SP$PPORRES[1:1000],AUC_extraction_250_Oral__C_SP$PPORRES[1001:2000],AUC_extraction_250_Oral__C_SI$PPORRES[1:1000],AUC_extraction_250_Oral__C_SI$PPORRES[1001:2000],
        AUC_extraction_250_Oral__C_L$PPORRES[1:1000],AUC_extraction_250_Oral__C_L$PPORRES[1001:2000],
        Main= "Area under the curve Concentration of Cinnamaldehyde",
        ylab= "umol/l",
        names= c("Male Lung", "Female Lung", "Male Blood", "Female Blood", "male Fat","Female Fat", "Male Richly perfused", "Femal Richly perfused","Male slowly perfused", "Femal slowly perfused","Male SI","Female SI",
                 "Male Liver","Female Liver"),
        col="orange",
        las=2)

#Combinded box plot 
boxplot(AUC_extraction_250_Oral__C_Pu$PPORRES, AUC_extraction_250_Oral__C_B$PPORRES,AUC_extraction_250_Oral__C_F$PPORRES,
        AUC_extraction_250_Oral__C_RP$PPORRES,AUC_extraction_250_Oral__C_SP$PPORRES,AUC_extraction_250_Oral__C_SI$PPORRES,
        AUC_extraction_250_Oral__C_L$PPORRES,
        Main= "Area under the curve Concentration of Cinnamaldehyde",
        ylab= "umol/l",
        names= c("Lung","Blood", "Fat", "Richly perfused","slowly perfused","SI",
                 "Liver"),
        col="orange",
        las=2)

#AUC and C max box plot
boxplot(AUC_extraction_250_Oral__C_Pu$PPORRES,Cmax_extraction_C_Pu$PPORRES, AUC_extraction_250_Oral__C_B$PPORRES,Cmax_extraction_C_B$PPORRES, AUC_extraction_250_Oral__C_F$PPORRES,Cmax_extraction_C_F$PPORRES,
        AUC_extraction_250_Oral__C_RP$PPORRES,Cmax_extraction_C_RP$PPORRES, AUC_extraction_250_Oral__C_SP$PPORRES,Cmax_extraction_C_SP$PPORRES, AUC_extraction_250_Oral__C_SI$PPORRES,Cmax_extraction_C_SI$PPORRES,
        AUC_extraction_250_Oral__C_L$PPORRES,Cmax_extraction_C_L$PPORRES,
        Main= "Area under the curve Concentration of Cinnamaldehyde",
        ylab= "umol/l",
        names= c("Lung AUC","Lung Cmax", "Blood AUC","Blood Cmax", "Fat AUC", "Fat Cmax", "Richly perfused AUC","RP Cmax", "slowly perfused AUC","SP cmax", "SI AUC", "SI Cmax",
                 "Liver AUC","liver Cmax"),
        col="orange",
        las=2)



#Lung AUC
#Loading data back in
results_250_inhalation_C_Pu <- read_csv("results_250_inhalation_C_Pu.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_inhalation__C_Pu <- unique(results_250_inhalation_C_Pu[results_250_inhalation_C_Pu$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_Pu <- unique(results_250_inhalation_C_Pu[results_250_inhalation_C_Pu$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_inhalation__C_Pu$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_inhalation__C_Pu$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_inhalation__C_Pu$PPORRES[1001:2000])



#---------------BLOOD------------------#
#Loading data back in
results_250_inhalation_C_B <- read_csv("results_250_inhalation_C_B.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_inhalation__C_B <- unique(results_250_inhalation_C_B[results_250_inhalation_C_B$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_B <- unique(results_250_inhalation_C_B[results_250_inhalation_C_B$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_inhalation__C_B$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_inhalation__C_B$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_inhalation__C_B$PPORRES[1001:2000])



#-----------------------------Liver----------------------#
#Loading data back in
results_250_inhalation_C_L <- read_csv("results_250_inhalation_C_L.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_inhalation__C_L <- unique(results_250_inhalation_C_L[results_250_inhalation_C_L$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_L <- unique(results_250_inhalation_C_L[results_250_inhalation_C_L$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_inhalation__C_L$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_inhalation__C_L$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_inhalation__C_L$PPORRES[1001:2000])


#-----------------------------Slowly perfused----------------------#
#Loading data back in
results_250_inhalation_C_SP <- read_csv("results_250_inhalation_C_SP.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_inhalation__C_SP <- unique(results_250_inhalation_C_SP[results_250_inhalation_C_SP$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_SP <- unique(results_250_inhalation_C_SP[results_250_inhalation_C_SP$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_inhalation__C_SP$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_inhalation__C_SP$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_inhalation__C_SP$PPORRES[1001:2000])




#----------------------Richly perfused--------------------------#
#Loading data back in
results_250_inhalation_C_RP <- read_csv("results_250_inhalation_C_RP.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_inhalation__C_RP <- unique(results_250_inhalation_C_RP[results_250_inhalation_C_RP$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_RP <- unique(results_250_inhalation_C_RP[results_250_inhalation_C_RP$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_inhalation__C_RP$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_inhalation__C_RP$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_inhalation__C_RP$PPORRES[1001:2000])


#--------------------FAT--------------------------#
#Loading data back in
results_250_inhalation_C_F <- read_csv("results_250_inhalation_C_F.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_inhalation__C_F <- unique(results_250_inhalation_C_F[results_250_inhalation_C_F$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_F <- unique(results_250_inhalation_C_F[results_250_inhalation_C_F$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_inhalation__C_F$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_inhalation__C_F$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_inhalation__C_F$PPORRES[1001:2000])


#-----------------------------Small intestine----------------------#
#Loading data back in
results_250_inhalation_C_SI <- read_csv("results_250_inhalation_C_SI.csv")

#Extracting AUC data for data visualization
AUC_extraction_250_inhalation__C_SI <- unique(results_250_inhalation_C_SI[results_250_inhalation_C_SI$PPTESTCD == "auclast",c("id","PPORRES")])

#Extracting Cmax data for data visualization
Cmax_extraction_C_SI <- unique(results_250_inhalation_C_SI[results_250_inhalation_C_SI$PPTESTCD == "cmax",c("id","PPORRES")])

#Histogram of AUC values
hist(AUC_extraction_250_inhalation__C_SI$PPORRES)

#Histogram of male AUC values
hist(AUC_extraction_250_inhalation__C_SI$PPORRES[1:1000])

#Histogram of female AUC values
hist(AUC_extraction_250_inhalation__C_SI$PPORRES[1001:2000])


#All compartments plot
boxplot(AUC_extraction_250_inhalation__C_Pu$PPORRES[1:1000], AUC_extraction_250_inhalation__C_Pu$PPORRES[1001:2000],AUC_extraction_250_inhalation__C_B$PPORRES[1:1000],AUC_extraction_250_inhalation__C_B$PPORRES[1001:2000],AUC_extraction_250_inhalation__C_F$PPORRES[1:1000],AUC_extraction_250_inhalation__C_F$PPORRES[1001:2000],
        AUC_extraction_250_inhalation__C_RP$PPORRES[1:1000],AUC_extraction_250_inhalation__C_RP$PPORRES[1001:2000],AUC_extraction_250_inhalation__C_SP$PPORRES[1:1000],AUC_extraction_250_inhalation__C_SP$PPORRES[1001:2000],AUC_extraction_250_inhalation__C_SI$PPORRES[1:1000],AUC_extraction_250_inhalation__C_SI$PPORRES[1001:2000],
        AUC_extraction_250_inhalation__C_L$PPORRES[1:1000],AUC_extraction_250_inhalation__C_L$PPORRES[1001:2000],
        Main= "Area under the curve Concentration of Cinnamaldehyde",
        ylab= "umol/l",
        names= c("Male Lung", "Female Lung", "Male Blood", "Female Blood", "male Fat","Female Fat", "Male Richly perfused", "Femal Richly perfused","Male slowly perfused", "Femal slowly perfused","Male SI","Female SI",
                 "Male Liver","Female Liver"),
        col="orange",
        las=2)



#Combinded box plot 
boxplot(AUC_extraction_250_inhalation__C_Pu$PPORRES, AUC_extraction_250_inhalation__C_B$PPORRES,AUC_extraction_250_inhalation__C_F$PPORRES,
        AUC_extraction_250_inhalation__C_RP$PPORRES,AUC_extraction_250_inhalation__C_SP$PPORRES,AUC_extraction_250_inhalation__C_SI$PPORRES,
        AUC_extraction_250_inhalation__C_L$PPORRES,
        Main= "Area under the curve Concentration of Cinnamaldehyde",
        ylab= "umol/l",
        names= c("Lung","Blood", "Fat", "Richly perfused","slowly perfused","SI",
                 "Liver"),
        col="orange",
        las=2)

#AUC and C max box plot
boxplot(AUC_extraction_250_inhalation__C_Pu$PPORRES,Cmax_extraction_C_Pu$PPORRES, AUC_extraction_250_inhalation__C_B$PPORRES,Cmax_extraction_C_B$PPORRES, AUC_extraction_250_inhalation__C_F$PPORRES,Cmax_extraction_C_F$PPORRES,
        AUC_extraction_250_inhalation__C_RP$PPORRES,Cmax_extraction_C_RP$PPORRES, AUC_extraction_250_inhalation__C_SP$PPORRES,Cmax_extraction_C_SP$PPORRES, AUC_extraction_250_inhalation__C_SI$PPORRES,Cmax_extraction_C_SI$PPORRES,
        AUC_extraction_250_inhalation__C_L$PPORRES,Cmax_extraction_C_L$PPORRES,
        Main= "Area under the curve Concentration of Cinnamaldehyde",
        ylab= "umol/l",
        names= c("Lung AUC","Lung Cmax", "Blood AUC","Blood Cmax", "Fat AUC", "Fat Cmax", "Richly perfused AUC","RP Cmax", "slowly perfused AUC","SP cmax", "SI AUC", "SI Cmax",
                 "Liver AUC","liver Cmax"),
        col="orange",
        las=2)




#Combined oral and inhalation box plot 
boxplot(AUC_extraction_250_Oral__C_Pu$PPORRES,AUC_extraction_250_inhalation__C_Pu$PPORRES, AUC_extraction_250_Oral__C_B$PPORRES, AUC_extraction_250_inhalation__C_B$PPORRES,
        AUC_extraction_250_Oral__C_F$PPORRES, AUC_extraction_250_inhalation__C_F$PPORRES, AUC_extraction_250_Oral__C_RP$PPORRES, AUC_extraction_250_inhalation__C_RP$PPORRES,
        AUC_extraction_250_Oral__C_SP$PPORRES, AUC_extraction_250_inhalation__C_SP$PPORRES, AUC_extraction_250_Oral__C_SI$PPORRES, AUC_extraction_250_inhalation__C_SI$PPORRES,
        AUC_extraction_250_Oral__C_L$PPORRES, AUC_extraction_250_inhalation__C_L$PPORRES,
        Main= "Area under the curve Concentration of Cinnamaldehyde",
        ylab= "umol/l",
        names= c("Lung oral", "Lung inhalation","Blood oral","Blood inhalation", "Fat oral","Fat inhalation", "RP oral","RP inhalation", "SP oral", "SP inhalation", "SI oral", "SI inhalation",
                 "Liver oral", "SI inhalation"),
        col="orange",
        las=2)




#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]


#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)



write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")




Combined_data_file_for_graph_500mg <- read_delim("Combined data file for graph 500mg.csv", 
                                                 delim = ";", escape_double = FALSE, trim_ws = TRUE)
colnames(Combined_data_file_for_graph_500mg) <- c("time","umol.l","ID")

p1 <- plot_ly(Combined_data_file_for_graph_500mg, x=~time, y=Combined_data_file_for_graph_500mg$'umol.l', 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 500mg/kg-bw',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in μmol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1

g <- ggplot(Combined_data_file_for_graph_500mg,aes(time,Combined_data_file_for_graph_500mg$'umol.l',color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+
  labs(subtitle="Oral dose 500mg/kg/bw", 
       y="Cinnamaldehyde concentration in μmol/l", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")


#Calculating AUC values vor the comparison graphs made above  

AUC_data <-PKNCAconc(Combined_data_file_for_graph_500mg, umol.l~time|ID)

#Dosing data per Paper/simulation run 
#Loading dosing file
d_dose <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/500m_comparison_doses.csv", sep=";")                       

dose_obj <- PKNCAdose(d_dose, dose~time|ID)




#letting pknc chose the end time of the auc calc
data_obj_automatic <- PKNCAdata(AUC_data, dose_obj)



#Computing the data both manual and automatic
results_obj_automatic <- pk.nca(data_obj_automatic)



#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_automatic)))
summary(results_obj_automatic)

#saving data 
write.csv(results_obj_automatic$result,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//500mg_oral_results.csv")


#RAT_1 calculation
x <- c(Combined_data_file_for_graph_500mg[1:15,1])
y <- c(Combined_data_file_for_graph_500mg[1:15,2])
RAT_1_AUC <- auc(x[["time"]],y[["umol.l"]], type=c("spline"))


#RAT_2
x <- c(Combined_data_file_for_graph_500mg[16:30,1])
y <- c(Combined_data_file_for_graph_500mg[16:30,2])
RAT_2_AUC <- auc(x[["time"]],y[["umol.l"]], type=c("spline"))

#RAT_3
x <- c(Combined_data_file_for_graph_500mg[31:46,1])
y <- c(Combined_data_file_for_graph_500mg[31:46,2])
RAT_3_AUC <- auc(x[["time"]],y[["umol.l"]], type=c("spline"))

#Zhao
x <- c(Combined_data_file_for_graph_500mg[47:57,1])
y <- c(Combined_data_file_for_graph_500mg[47:57,2])
ZHAO_AUC <- auc(x[["time"]],y[["umol.l"]], type=c("spline"))

#Generating a file for the comparison with in vivo data for 250 mg exposure  
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]

#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


Combined_data_file_for_graph_250mg <- read_delim("Combined data file for graph 250mg.csv", 
                                                 delim = ";", escape_double = FALSE, trim_ws = TRUE)

p1 <- plot_ly(Combined_data_file_for_graph_250mg, x=~time, y=Combined_data_file_for_graph_250mg$`ug/ml`, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 250mg/kg-bw',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in umol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1

g <- ggplot(Combined_data_file_for_graph_250mg,aes(time,Combined_data_file_for_graph_250mg$`ug/ml`,color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+
  labs(subtitle="Oral dose 250mg/kg/bw", 
       y="Cinnamaldehyde concentration in umol/l", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="Cinnamaldehyde PBK model")
#Calculating AUC values for the comparison graphs made above  

AUC_data <-PKNCAconc(Combined_data_file_for_graph_250mg, ug.ml~time|ID)

#Dosing data per Paper/simulation run 
#Loading dosing file
d_dose <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/250mg_comparison_doses.csv", sep=";")                       

dose_obj <- PKNCAdose(d_dose, dose~time|ID)




#letting pknc chose the end time of the auc calc
data_obj_automatic <- PKNCAdata(AUC_data, dose_obj)



#Computing the data both manual and automatic
results_obj_automatic <- pk.nca(data_obj_automatic)



#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_automatic)))
summary(results_obj_automatic)

#saving data 
write.csv(results_obj_automatic$result,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//250mg_oral_results.csv")


#RAT_1 calculation
x <- c(Combined_data_file_for_graph_250mg[21:36,1])
y <- c(Combined_data_file_for_graph_250mg[21:36,2])
RAT_1_AUC <- auc(x[["time"]],y[["ug.ml"]], type=c("spline"))


#RAT_2
x <- c(Combined_data_file_for_graph_250mg[37:52,1])
y <- c(Combined_data_file_for_graph_250mg[37:52,2])
RAT_2_AUC <- auc(x[["time"]],y[["ug.ml"]], type=c("spline"))

#RAT_3
x <- c(Combined_data_file_for_graph_250mg[53:68,1])
y <- c(Combined_data_file_for_graph_250mg[53:68,2])
RAT_3_AUC <- auc(x[["time"]],y[["ug.ml"]], type=c("spline"))

#Zhao
x <- c(Combined_data_file_for_graph_250mg[1:8,1])
y <- c(Combined_data_file_for_graph_250mg[1:8,2])
ZHAO_AUC <- auc(x[["time"]],y[["ug.ml"]], type=c("spline"))

#Kiwa
x <- c(Combined_data_file_for_graph_250mg[9:20,1])
y <- c(Combined_data_file_for_graph_250mg[9:20,2])
Kiwa_AUC <- auc(x[["time"]],y[["ug.ml"]], type=c("spline"))


#iv concentration
#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]

#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


Combined_data_file_for_graph <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/IV_blood_Concentration.csv", sep=";")

p1 <- plot_ly(Combined_data_file_for_graph, x=~Time, y=Combined_data_file_for_graph$umol.L, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 20mg/kg-bw IV',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in umol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1


g <- ggplot(Combined_data_file_for_graph,aes(Time,umol.L,color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+
  labs(subtitle="IV dose 20mg/kg/bw", 
       y="Cinnamaldehyde concentration in umol/l", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")



#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]


#going from umol to ug
blood_data[2]<- blood_data[2]*136.12


#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)

#going back an amount to an amount per ml
blood_data[2]<- blood_data[2]/ 1000


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")

Combined_data_file_for_graph <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/IV_10mg_comparison.csv", sep=";")

p1 <- plot_ly(Combined_data_file_for_graph, x=~Time, y=Combined_data_file_for_graph$ug.ml, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 10mg/kg-bw IV',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in ug/ml', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1


g <- ggplot(Combined_data_file_for_graph,aes(Time,ug.ml,color=ID))

g + geom_point()+ 
  labs(subtitle="IV dose 10mg/kg/bw", 
       y="Cinnamaldehyde concentration in ug/ml", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")



####-------------DONG et al 2022 data comparison--------------####
#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]

#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")



Combined_data_file_for_graph <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/375 oral dose comparison DONG et al 2022.csv", sep=";")

p1 <- plot_ly(Combined_data_file_for_graph, x=~Time, y=Combined_data_file_for_graph$umol.l, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 375mg/kg/bw',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in umol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1


g <- ggplot(Combined_data_file_for_graph,aes(Time,umol.l,color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+ 
  labs(subtitle="oral375mg/kg/bw dose", 
       y="Cinnamaldehyde concentration in umol/l", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")

####-------------JI et al 2013 data comparison--------------####
#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]

#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")



Combined_data_file_for_graph <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/15 mg oral dose .csv", sep=";")

p1 <- plot_ly(Combined_data_file_for_graph, x=~Time, y=Combined_data_file_for_graph$umol.l, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 12.5mg/kg-bw IV',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in umol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1


g <- ggplot(Combined_data_file_for_graph,aes(Time,umol.l,color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+ 
  labs(subtitle="IV dose 12.5mg/kg/bw", 
       y="Cinnamaldehyde concentration in umol/l", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")


####-------------YOUNG et al 2020 data comparison--------------####
#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]

#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")



Combined_data_file_for_graph <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/Yong et al 2020 50mg-kg-bw.csv", sep=";")

p1 <- plot_ly(Combined_data_file_for_graph, x=~Time, y=Combined_data_file_for_graph$umol.l, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 50mg/kg/bw',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in umol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1


g <- ggplot(Combined_data_file_for_graph,aes(Time,umol.l,color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+ 
  labs(subtitle="oral 50 mg/kg/bw dose", 
       y="Cinnamaldehyde concentration in umol/l", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")







#Metabolic output checking
#human model
#Checking proportions of urinary metabolites
#Urinary metabolites are all ultimately derived from the oxidation of Cinnamaldehyde.

data_urinary_metabolites <- solve.pbk_popgen[,c(1,79,88)]

total_dose <- Oral_Dose

total_urinary_metabolites_24h <- data_urinary_metabolites[241,2] + data_urinary_metabolites[241,3]

percentage_metabolised <- total_urinary_metabolites_24h/total_dose *100


#Similar data extraction but now for the population model

meta_extraction <- unique(solve.pbk_popgen[solve.pbk_popgen$time == 24,c("time", "id","AM_L_CA","AM_SI_CA")])

dose_extraction <- unique(ex[ex$cmt == "A_GI",c("id","amt")])

for (id in meta_extraction){
  total <- meta_extraction[id,2] + meta_extraction[id,3]
  dose  <- dose_extraction[id,2]
  metabolised <- total/dose *100
  percentage_metabolised_popgen <- metabolised
}

metabolised_data_frame <- as.data.frame(percentage_metabolised_popgen, row.names="percentage")

male_metabolised_data<- metabolised_data_frame[1:1000,]
female_metabolised_data<- metabolised_data_frame[1:2000,]

boxplot(male_metabolised_data,female_metabolised_data,
        Main= "percentage metabolised into urinary metabolites",
        ylab= "%",
        names= c("Male", "Female"),
        col="orange",
        las=2)
hist(male_metabolised_data)
hist(female_metabolised_data)

set.seed(101)
dd <- data.frame(x=rnorm(100),y=rnorm(100),
                 z=rnorm(100))
dd$w <- with(dd,
             rnorm(100,mean=x+2*y+z,sd=0.5))

m <- lm(w~x+y+z,dd)
plot(predict(m),dd$w,
     xlab="predicted",ylab="actual")
abline(a=0,b=1)