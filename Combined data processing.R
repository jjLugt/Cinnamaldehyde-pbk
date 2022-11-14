#author: Joris Jean van der Lugt
#date: 05-08-2022
#Combined data processing of Human Cinnamaldehyde pbk models
#Before running code first run all model associated with this data processing
#human rxode model


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
mass_df <- solve.pbk/BW * MW /1e+3
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
mass_df <- solve.pbk/phys[1,3] * MW /1e+3
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



#AUC calculations
#single Human model
#Making a dataframe for the calculations
AUC_single_human<- solve.pbk[,c(1,11,17,23,29,35,49)]

#Combining Arterial and venous blood into one general blood compartment
AUC_single_human[,8]<- solve.pbk[,3]+ solve.pbk[,4]

colnames(AUC_single_human)<-c("time","C_Pu","C_F","C_RP","C_SP","C_SI","C_L","C_B")

#passing the concentration data frame to PKNCA
single_human_conc <- PKNCAconc(AUC_single_human, C_B~time)

#Generating a Dose data frame which includes time and dose
single_human_dose <-as.data.frame(solve.pbk[,1])
#Grabbing an easy empty colum of the correct length
single_human_dose[2] <-solve.pbk[,7]
#Adding the dose to the first time point in the dose data frame 
single_human_dose[1,2] <-ex[2,6]
colnames(single_human_dose)<-c("time","dose")


#passing the dose data frame to PKNCA
dose_obj <- PKNCAdose(single_human_dose, dose~time)

#Running the calculations 
#Setting the end of the auc calculation at 24 hours
intervals_manual <- data.frame(start=0,
                               end=24,
                               cmax=TRUE,
                               tmax=TRUE,
                               aucinf.obs=TRUE,
                               auclast=TRUE)

#Creating the combined data obj with bothe dose and concentration data 
data_obj_manual<- PKNCAdata(single_human_conc, dose_obj,
                            intervals=intervals_manual)
#Running the PKNCA analysis 
results_obj_manual <- pk.nca(data_obj_manual)

#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_manual)))

summary(results_obj_manual)
results_obj_manual$result

#Writing the results into a CSV file 
write.csv(results_obj_manual$result, "results_250mg_oral_single_human_C_B_AUC")


#250mg dose Rat AUC values 
#AUC calculations
#single Human model
#Making a dataframe for the calculations
AUC_single_rat<- solve.pbk_rat[,c(1,11,17,23,29,35,48)]

#Combining Arterial and venous blood into one general blood compartment
AUC_single_rat[,8]<- solve.pbk_rat[,3]+ solve.pbk_rat[,4]

colnames(AUC_single_rat)<-c("time","C_Pu","C_F","C_RP","C_SP","C_SI","C_L","C_B")

#Generating a Dose data frame which includes time and dose
single_rat_dose <-as.data.frame(solve.pbk_rat[,1])
#Grabbing an easy empty colum of the correct length
single_rat_dose[2] <-solve.pbk_rat[,7]
#Adding the dose to the first time point in the dose data frame 
single_rat_dose[1,2] <-ex[2,6]
colnames(single_rat_dose)<-c("time","dose")


#---------Calculating Blood auc values--------------#
#passing the concentration data frame to PKNCA
single_rat_conc <- PKNCAconc(AUC_single_rat, C_B~time)

#passing the dose data frame to PKNCA
dose_obj <- PKNCAdose(single_rat_dose, dose~time)

#Running the calculations 
#Setting the end of the auc calculation at 24 hours
intervals_manual <- data.frame(start=0,
                               end=24,
                               cmax=TRUE,
                               tmax=TRUE,
                               aucinf.obs=TRUE,
                               auclast=TRUE)

#Creating the combined data obj with bothe dose and concentration data 
data_obj_manual<- PKNCAdata(single_rat_conc, dose_obj,
                            intervals=intervals_manual)
#Running the PKNCA analysis 
results_obj_manual <- pk.nca(data_obj_manual)

#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_manual)))

summary(results_obj_manual)
results_obj_manual$result

#Writing the results into a CSV file 
write.csv(results_obj_manual$result, "results_250mg_oral_single_rat_C_B_AUC")

#------------------calculating Lung auc values----------#
#passing the concentration data frame to PKNCA
single_rat_conc <- PKNCAconc(AUC_single_rat, C_Pu~time)

#passing the dose data frame to PKNCA
dose_obj <- PKNCAdose(single_rat_dose, dose~time)

#Running the calculations 
#Setting the end of the auc calculation at 24 hours
intervals_manual <- data.frame(start=0,
                               end=24,
                               cmax=TRUE,
                               tmax=TRUE,
                               aucinf.obs=TRUE,
                               auclast=TRUE)

#Creating the combined data obj with bothe dose and concentration data 
data_obj_manual<- PKNCAdata(single_rat_conc, dose_obj,
                            intervals=intervals_manual)
#Running the PKNCA analysis 
results_obj_manual <- pk.nca(data_obj_manual)

#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_manual)))

summary(results_obj_manual)
results_obj_manual$result

#Writing the results into a CSV file 
write.csv(results_obj_manual$result, "results_250mg_oral_single_rat_C_Pu_AUC")

#------------------calculating SP tissue auc values----------#
#passing the concentration data frame to PKNCA
single_rat_conc <- PKNCAconc(AUC_single_rat, C_SP~time)

#passing the dose data frame to PKNCA
dose_obj <- PKNCAdose(single_rat_dose, dose~time)

#Running the calculations 
#Setting the end of the auc calculation at 24 hours
intervals_manual <- data.frame(start=0,
                               end=24,
                               cmax=TRUE,
                               tmax=TRUE,
                               aucinf.obs=TRUE,
                               auclast=TRUE)

#Creating the combined data obj with bothe dose and concentration data 
data_obj_manual<- PKNCAdata(single_rat_conc, dose_obj,
                            intervals=intervals_manual)
#Running the PKNCA analysis 
results_obj_manual <- pk.nca(data_obj_manual)

#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_manual)))

summary(results_obj_manual)
results_obj_manual$result

#Writing the results into a CSV file 
write.csv(results_obj_manual$result, "results_250mg_oral_single_rat_C_SP_AUC")

#------------------calculating RP auc values----------#
#passing the concentration data frame to PKNCA
single_rat_conc <- PKNCAconc(AUC_single_rat, C_RP~time)

#passing the dose data frame to PKNCA
dose_obj <- PKNCAdose(single_rat_dose, dose~time)

#Running the calculations 
#Setting the end of the auc calculation at 24 hours
intervals_manual <- data.frame(start=0,
                               end=24,
                               cmax=TRUE,
                               tmax=TRUE,
                               aucinf.obs=TRUE,
                               auclast=TRUE)

#Creating the combined data obj with bothe dose and concentration data 
data_obj_manual<- PKNCAdata(single_rat_conc, dose_obj,
                            intervals=intervals_manual)
#Running the PKNCA analysis 
results_obj_manual <- pk.nca(data_obj_manual)

#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_manual)))

summary(results_obj_manual)
results_obj_manual$result

#Writing the results into a CSV file 
write.csv(results_obj_manual$result, "results_250mg_oral_single_rat_C_RP_AUC")


#------------------calculating si auc values----------#
#passing the concentration data frame to PKNCA
single_rat_conc <- PKNCAconc(AUC_single_rat, C_SI~time)

#passing the dose data frame to PKNCA
dose_obj <- PKNCAdose(single_rat_dose, dose~time)

#Running the calculations 
#Setting the end of the auc calculation at 24 hours
intervals_manual <- data.frame(start=0,
                               end=24,
                               cmax=TRUE,
                               tmax=TRUE,
                               aucinf.obs=TRUE,
                               auclast=TRUE)

#Creating the combined data obj with bothe dose and concentration data 
data_obj_manual<- PKNCAdata(single_rat_conc, dose_obj,
                            intervals=intervals_manual)
#Running the PKNCA analysis 
results_obj_manual <- pk.nca(data_obj_manual)

#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_manual)))

summary(results_obj_manual)
results_obj_manual$result

#Writing the results into a CSV file 
write.csv(results_obj_manual$result, "results_250mg_oral_single_rat_C_SI_AUC")

#------------------calculating Liver auc values----------#
#passing the concentration data frame to PKNCA
single_rat_conc <- PKNCAconc(AUC_single_rat, C_L~time)

#passing the dose data frame to PKNCA
dose_obj <- PKNCAdose(single_rat_dose, dose~time)

#Running the calculations 
#Setting the end of the auc calculation at 24 hours
intervals_manual <- data.frame(start=0,
                               end=24,
                               cmax=TRUE,
                               tmax=TRUE,
                               aucinf.obs=TRUE,
                               auclast=TRUE)

#Creating the combined data obj with bothe dose and concentration data 
data_obj_manual<- PKNCAdata(single_rat_conc, dose_obj,
                            intervals=intervals_manual)
#Running the PKNCA analysis 
results_obj_manual <- pk.nca(data_obj_manual)

#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_manual)))

summary(results_obj_manual)
results_obj_manual$result

#Writing the results into a CSV file 
write.csv(results_obj_manual$result, "results_250mg_oral_single_rat_C_L_AUC")

#------------------calculating Lung FAT values----------#
#passing the concentration data frame to PKNCA
single_rat_conc <- PKNCAconc(AUC_single_rat, C_Pu~time)

#passing the dose data frame to PKNCA
dose_obj <- PKNCAdose(single_rat_dose, dose~time)

#Running the calculations 
#Setting the end of the auc calculation at 24 hours
intervals_manual <- data.frame(start=0,
                               end=24,
                               cmax=TRUE,
                               tmax=TRUE,
                               aucinf.obs=TRUE,
                               auclast=TRUE)

#Creating the combined data obj with bothe dose and concentration data 
data_obj_manual<- PKNCAdata(single_rat_conc, dose_obj,
                            intervals=intervals_manual)
#Running the PKNCA analysis 
results_obj_manual <- pk.nca(data_obj_manual)

#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_manual)))

summary(results_obj_manual)
results_obj_manual$result

#Writing the results into a CSV file 
write.csv(results_obj_manual$result, "results_250mg_oral_single_rat_C_F_AUC")


#Plotting rat and human results of the auc calculations
results_250mg_oral_single_human_C_B_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_human_C_B_AUC")
results_250mg_oral_single_human_C_L_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_human_C_L_AUC")
results_250mg_oral_single_human_C_SI_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_human_C_SI_AUC")
results_250mg_oral_single_human_C_SP_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_human_C_SP_AUC")
results_250mg_oral_single_human_C_RP_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_human_C_RP_AUC")
results_250mg_oral_single_human_C_F_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_human_C_F_AUC")
results_250mg_oral_single_human_C_Pu_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_human_C_Pu_AUC")

#Plotting rat and human results of the auc calculations
results_250mg_oral_single_rat_C_B_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_B_AUC")
results_250mg_oral_single_rat_C_L_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_L_AUC")
results_250mg_oral_single_rat_C_SI_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_SI_AUC")
results_250mg_oral_single_rat_C_SP_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_SP_AUC")
results_250mg_oral_single_rat_C_RP_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_RP_AUC")
results_250mg_oral_single_rat_C_F_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_F_AUC")
results_250mg_oral_single_rat_C_Pu_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_Pu_AUC")







#Rxode
#After running population model
#Extracting organ concentration, time and sim-id from simulation results

#Lung compartment concentration
sub_set <- solve.pbk[1:482000,c(1,2,12)]

conc_C <- PKNCAconc(sub_set, C_Pu~time|id)

#Blood compartment concentration
sub_set <- solve.pbk[1:482000,c(1,2,4,5)]

#Combining Arterial and venous blood into one general blood compartment
sub_set[3]<- sub_set[3]+sub_set[4]

#Dropping column as data is already added to column 3
sub_set[4] <- NULL
#Renamning the columns
colnames(sub_set) <- c("id","time","C_B")
conc_C <- PKNCAconc(sub_set, C_B~time|id)

#Liver compartment concentration
sub_set <- solve.pbk[1:482000,c(1,2,50)]
onc_C <- PKNCAconc(sub_set, C_L~time|id)
#Slowly perfused compartment
sub_set <- solve.pbk[1:482000,c(1,2,30)]
onc_C <- PKNCAconc(sub_set, C_SP~time|id)
#richly perfused compartment
sub_set <- solve.pbk[1:482000,c(1,2,24)]
onc_C <- PKNCAconc(sub_set, C_RP~time|id)
#Fat compartment
sub_set <- solve.pbk[1:482000,c(1,2,18)]
onc_C <- PKNCAconc(sub_set, C_F~time|id)
#sI compartment
sub_set <- solve.pbk[1:482000,c(1,2,36)]
onc_C <- PKNCAconc(sub_set, C_SI~time|id)




#Dosing data per subject is part of the parameter file but it is missing sim id and the time variable 
#these will be added.


#Oral dose extraction
dose_extraction <- as.data.frame(parameters[,61])
sim_extraction <- unique(solve.pbk[solve.pbk$time == 0,c("time", "id")])

#Combining into 1 file that can be used
d_dose <- cbind(sim_extraction,dose_extraction$`parameters[, 61]`)                       
d_dose <- set_names(d_dose, c("time","id","dose"))                        

dose_obj <- PKNCAdose(d_dose, dose~time|id)

#Setting the end of the auc calculation at 8 hours
intervals_manual <- data.frame(start=0,
                               end=24,
                               cmax=TRUE,
                               tmax=TRUE,
                               aucinf.obs=TRUE,
                               auclast=TRUE)


data_obj_manual<- PKNCAdata(conc_C, dose_obj,
                            intervals=intervals_manual)

results_obj_manual <- pk.nca(data_obj_manual)

#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_manual)))

summary(results_obj_manual)



#Writing the results into a CSV file 
write.csv(results_obj_manual$result, "results_250_oral_C_B.csv")


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

#boxplot
boxplot(AUC_extraction_250_Oral__C_B$PPORRES[1:1000])



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

#combining C_V + C_A to create a combined amount of cinnamaldehyde
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


RAT_data_obs_1 <- Combined_data_file_for_graph_500mg[1:15,]
RAT_data_obs_2 <- Combined_data_file_for_graph_500mg[16:30,]
RAT_data_obs_3 <- Combined_data_file_for_graph_500mg[31:46,]
RAT_data_Zao   <- Combined_data_file_for_graph_500mg[47:57,]
SIM_data_pred    <-Combined_data_file_for_graph_500mg[c(79,89,100,114,124,132,150,169,190,210,231,252,273,302,322),]
SIM_data_pred_kiwa <- Combined_data_file_for_graph_500mg[c(380,390,401,415,425,433,451,470,491,511,532,553,574,603,623),]
SIM_data_ka    <-Combined_data_file_for_graph_500mg[c(681,691,702,716,726,734,752,771,792,812,833,854,875,904,924),]

SIM_data_pred[,4]<- as.data.frame(RAT_data_obs_1[,2])
SIM_data_pred[,5]<- as.data.frame(RAT_data_obs_2[,2])
SIM_data_pred[,6]<- as.data.frame(RAT_data_obs_3[1:15,2])
SIM_data_pred[,7]<- as.data.frame(SIM_data_pred_kiwa[,2])
SIM_data_pred[,8]<- as.data.frame(SIM_data_ka[,2])
colnames(SIM_data_pred)<- c("Time","sim","ID","rat_1","rat_2","rat_3","SIM-Kiwa","ka")

p<-ggplot(SIM_data_pred, aes(x=SIM_data_pred$sim, y=SIM_data_pred$rat_1)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Predicted Values', y='Actual Values', title='Predicted vs. Actual Values 500mg oral dose Yuan data ')+
  ylim(0,50)+
  geom_point(aes(x=SIM_data_pred$sim,y=SIM_data_pred$rat_2))+
  geom_point(aes(x=SIM_data_pred$sim,y=SIM_data_pred$rat_3))+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
ggsave(plot=p,"Pred vs actual 500mg oral Yuan data.png",
       width= 11.69, height= 8.3, dpi= 250)

#plotting new plot with adjusted parameters for ka and k_l_ca
p<-ggplot(SIM_data_pred, aes(x=SIM_data_pred$ka, y=SIM_data_pred$rat_1)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Predicted Values', y='Actual Values', title='Predicted vs. Actual Values 500mg oral dose Yuan data ka ')+
  xlim(0,20)+
  ylim(0,20)+
  geom_point(aes(x=SIM_data_pred$ka,y=SIM_data_pred$rat_2))+
  geom_point(aes(x=SIM_data_pred$ka,y=SIM_data_pred$rat_3))+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
ggsave(plot=p,"Pred vs actual 500mg oral Yuan data ka .png",
       width= 11.69, height= 8.3, dpi= 250)


#Making a pred vs pred plot

Kiwa_data      <-Combined_data_file_for_graph_500mg[c(58,60:72),]
SIM_data_pred  <-Combined_data_file_for_graph_500mg[c(374,375,376,379,384,401,427,468,505,538,574,607,640,667),]
SIM_data_pred[2,2]  <-blood_data[6,2]
SIM_data_pred[3,2]  <-blood_data[16,2]
SIM_data_pred[4,3]<- blood_data[47,2]

colnames(SIM_data_pred)<- c("Time","sim","kiwa")
p<-ggplot(SIM_data_pred, aes(x=SIM_data_pred$sim, y=SIM_data_pred$kiwa)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Predicted Values R', y='Predicted values Kiwa', title='Predicted vs. Predicted Values 500mg oral dose')+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
ggsave(plot=p,"Pred vs pred 500mg oral.png",
       width= 11.69, height= 8.3, dpi= 250)

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

p1 <- plot_ly(Combined_data_file_for_graph_250mg, x=~time, y=Combined_data_file_for_graph_250mg$`ug.ml`, 
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

RAT_data_obs_1 <- Combined_data_file_for_graph_250mg[21:36,]
RAT_data_obs_2 <- Combined_data_file_for_graph_250mg[37:52,]
RAT_data_obs_3 <- Combined_data_file_for_graph_250mg[53:68,]
RAT_data_Zao   <- Combined_data_file_for_graph_250mg[1:8,]
SIM_data_pred  <-Combined_data_file_for_graph_250mg[c(70,71,74,80,84,89,94,104,114,124,134,144,154,164,179,189),]
SIM_data_Ka    <-Combined_data_file_for_graph_250mg[c(312,313,316,322,326,331,336,346,356,366,376,387,396,411,421,431),]


SIM_data_pred[,4]<- as.data.frame(RAT_data_obs_1[,2])
SIM_data_pred[,5]<- as.data.frame(RAT_data_obs_2[,2])
SIM_data_pred[,6]<- as.data.frame(RAT_data_obs_3[,2])
SIM_data_pred[,7]<- as.data.frame(SIM_data_Ka[,2])


#Calculating residual values
SIM_data_pred[,8] <- SIM_data_pred[,2]-SIM_data_pred[,4]
SIM_data_pred[,9] <- SIM_data_pred[,4]-SIM_data_pred[,7]
SIM_data_pred[,10] <- SIM_data_pred[,5]-SIM_data_pred[,7]


colnames(SIM_data_pred)<- c("Time","sim","ID","rat_1", "rat_2","rat_3","SIM_ka","Residual_Rat_1","Residual_Rat1_ka","Residual_Rat2_ka")

predvsout_250mg <- ggplot(SIM_data_pred, aes(x=SIM_data_pred$sim, y=SIM_data_pred$rat_1)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Predicted Values', y='Actual Values', title='Predicted vs. Actual Values 250mg oral Yuan data')+
  ylim(0,20)+
  xlim(0,20)+
  geom_point(aes(x=SIM_data_pred$sim,y=SIM_data_pred$rat_2))
  geom_point(aes(x=SIM_data_pred$sim,y=SIM_data_pred$rat_3))

predvsout_250mg

predvsout_250mg_ka <- ggplot(SIM_data_pred, aes(x=SIM_data_pred$SIM_ka, y=SIM_data_pred$rat_1)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Predicted Values', y='Actual Values', title='Predicted vs. Actual Values 250mg oral Yuan data')+
  ylim(0,15)+
  xlim(0,10)+
  geom_point(aes(x=SIM_data_pred$SIM_ka,y=SIM_data_pred$rat_2))+
  geom_point(aes(x=SIM_data_pred$SIM_ka,y=SIM_data_pred$rat_3))+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
ggsave(plot=p,"Pred vs actual 250mg oral Yuan data.png",
       width= 11.69, height= 8.3, dpi= 250)

predvsout_250mg_ka
  
  
#plotting residual rat 1
ggplot(SIM_data_pred, aes(x=SIM_data_pred$Time, y=SIM_data_pred$Residual_Rat_1)) +
  geom_point() +
  geom_abline(slope = 0,intercept = 0)+ 
labs(x= "Time", y='Residual values', title='Residual Values 250mg oral Yuan data')

#plotting residual rat 1  vs ka
ggplot(SIM_data_pred, aes(x=SIM_data_pred$Time, y=SIM_data_pred$Residual_Rat1_ka)) +
  geom_point() +
  geom_abline(slope = 0,intercept = 0)+ 
  labs(x= "Time", y='Residual values', title='Residual Values 250mg oral Yuan data ka')

#plotting residual rat 2  vs ka
ggplot(SIM_data_pred, aes(x=SIM_data_pred$Time, y=SIM_data_pred$Residual_Rat2_ka)) +
  geom_point() +
  geom_abline(slope = 0,intercept = 0)+ 
  labs(x= "Time", y='Residual values', title='Residual Values 250mg oral Yuan data')


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


KIWA_data <- Combined_data_file_for_graph[9:27,]
SIM_data_pred <-Combined_data_file_for_graph[c(29,30,31,32,33,34,35,36,37,38,40,43,46,49,52,55,58,61,63),]


SIM_data_pred[,4]<- as.data.frame(KIWA_data[,2])


colnames(SIM_data_pred)<- c("Time","sim","ID","KIWA")

p <-ggplot(SIM_data_pred, aes(x=SIM_data_pred$sim, y=SIM_data_pred$KIWA)) +
  geom_point(size=3) +
  geom_abline(intercept=0, slope=1) +
  labs(x='Predicted Values R', y='Predicted values Kiwa', title='Predicted vs. predicted Values 20mg IV dose Kiwa')+
  xlim(0,45)+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
ggsave(plot=p,"Pred vs pred 20mg iv.png",
       width= 11.69, height= 8.3, dpi= 250)


ZAO_data_2014 <- Combined_data_file_for_graph[1:8,]
SIM_data_kiwa    <-Combined_data_file_for_graph[c(29,30,31,33,38,43,48,58),]

SIM_data_kiwa[,4]<- as.data.frame(ZAO_data_2014[,2])


colnames(SIM_data_kiwa)<- c("Time","sim","ID","ZAO_2014")

p <-ggplot(SIM_data_kiwa, aes(x=SIM_data_kiwa$sim, y=SIM_data_kiwa$ZAO_2014)) +
  geom_point(size=3) +
  geom_abline(intercept=0, slope=1) +
  labs(x='Predicted Values', y='Actual Values', title='Predicted vs. Actual Values 20mg IV dose ZAO data')+
  ylim(0,10)+
  xlim(x,10)+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
ggsave(plot=p,"Pred vs actual 20mg iv.png",
       width= 11.69, height= 8.3, dpi= 250)


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

data_urinary_metabolites <- solve.pbk[,c(1,79,88)]

total_dose <- Oral_Dose

total_urinary_metabolites_24h <- data_urinary_metabolites[241,2] + data_urinary_metabolites[241,3]

percentage_metabolised <- total_urinary_metabolites_24h/total_dose *100


#Similar data extraction but now for the population model

meta_extraction <- unique(solve.pbk[solve.pbk$time == 24,c("time", "id","AM_L_CA","AM_SI_CA")])

dose_extraction <- unique(ex[ex$cmt == "A_GI",c("id","amt")])

for (id in meta_extraction){
  total <- meta_extraction[id,3] + meta_extraction[id,4]
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
        las=2,
        title="Amount metabolized to carboxylic acid")
hist(male_metabolised_data)
hist(female_metabolised_data)

