#AUC calculations
#single Human model
#Making a dataframe for the calculations
AUC_single_human<- solve.pbk[,c(1,11,17,23,29,35,49)]

#Combining Arterial and venous blood into one general blood compartment
AUC_single_human[,8]<- solve.pbk[,3]+ solve.pbk[,4]

colnames(AUC_single_human)<-c("time","C_Pu","C_F","C_RP","C_SP","C_SI","C_L","C_B")



#Generating a Dose data frame which includes time and dose
single_human_dose <-as.data.frame(solve.pbk[,1])
#Grabbing an easy empty colum of the correct length
single_human_dose[2] <-solve.pbk[,7]
#Adding the dose to the first time point in the dose data frame 

#Oral Dose
single_human_dose[1,2] <-ex[2,6]
colnames(single_human_dose)<-c("time","dose")

#inhalation dose 
single_human_dose[1,2] <-ex[3,6]
colnames(single_human_dose)<-c("time","dose")

#IV dose
single_human_dose[1,2] <-ex[6,6]
colnames(single_human_dose)<-c("time","dose")


#-------------Liver-----------#
#passing the concentration data frame to PKNCA
single_human_conc <- PKNCAconc(AUC_single_human, C_L~time)

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
write.csv(results_obj_manual$result, "results_250mg_inhalation_single_human_C_L_AUC")


#-------------Lung tissue-----------#
#passing the concentration data frame to PKNCA
single_human_conc <- PKNCAconc(AUC_single_human, C_Pu~time)

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
write.csv(results_obj_manual$result, "results_250mg_inhalation_single_human_C_Pu_AUC")


#-------------Fat-----------#
#passing the concentration data frame to PKNCA
single_human_conc <- PKNCAconc(AUC_single_human, C_F~time)

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
write.csv(results_obj_manual$result, "results_250mg_inhalation_single_human_C_F_AUC")


#-------------Richly perfused-----------#
#passing the concentration data frame to PKNCA
single_human_conc <- PKNCAconc(AUC_single_human, C_RP~time)

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
write.csv(results_obj_manual$result, "results_250mg_inhalation_single_human_C_RP_AUC")


#-------------Slowly perfused-----------#
#passing the concentration data frame to PKNCA
single_human_conc <- PKNCAconc(AUC_single_human, C_SP~time)

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
write.csv(results_obj_manual$result, "results_250mg_inhalation_single_human_C_SP_AUC")


#-------------Blood-----------#
#passing the concentration data frame to PKNCA
single_human_conc <- PKNCAconc(AUC_single_human, C_B~time)

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
write.csv(results_obj_manual$result, "results_250mg_inhalation_single_human_C_B_AUC")



#-------------Small intestine-----------#
#passing the concentration data frame to PKNCA
single_human_conc <- PKNCAconc(AUC_single_human, C_SI~time)

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
write.csv(results_obj_manual$result, "results_250mg_inhalation_single_human_C_SI_AUC")



#------------------250mg dose Rat AUC values------------------- 
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


#Plotting rat and human results of the auc calculations------
results_250mg_inhalation_single_human_C_B_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_inhalation_single_human_C_B_AUC")
results_250mg_inhalation_single_human_C_L_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_inhalation_single_human_C_L_AUC")
results_250mg_inhalation_single_human_C_SI_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_inhalation_single_human_C_SI_AUC")
results_250mg_inhalation_single_human_C_SP_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_inhalation_single_human_C_SP_AUC")
results_250mg_inhalation_single_human_C_RP_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_inhalation_single_human_C_RP_AUC")
results_250mg_inhalation_single_human_C_F_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_inhalation_single_human_C_F_AUC")
results_250mg_inhalation_single_human_C_Pu_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_inhalation_single_human_C_Pu_AUC")

Plot_AUC_Human_250mg<- as.data.frame(c())
Plot_AUC_Human_250mg[1,1]<-as.data.frame(results_250mg_inhalation_single_human_C_B_AUC$PPORRES[1])
Plot_AUC_Human_250mg[1,2]<-as.data.frame(results_250mg_inhalation_single_human_C_B_AUC$PPORRES[2])

Plot_AUC_Human_250mg[2,1]<-as.data.frame(results_250mg_inhalation_single_human_C_L_AUC$PPORRES[1])
Plot_AUC_Human_250mg[2,2]<-as.data.frame(results_250mg_inhalation_single_human_C_L_AUC$PPORRES[2])

Plot_AUC_Human_250mg[3,1]<-as.data.frame(results_250mg_inhalation_single_human_C_SI_AUC$PPORRES[1])
Plot_AUC_Human_250mg[3,2]<-as.data.frame(results_250mg_inhalation_single_human_C_SI_AUC$PPORRES[2])

Plot_AUC_Human_250mg[4,1]<-as.data.frame(results_250mg_inhalation_single_human_C_SP_AUC$PPORRES[2])
Plot_AUC_Human_250mg[4,2]<-as.data.frame(results_250mg_inhalation_single_human_C_SP_AUC$PPORRES[2])

Plot_AUC_Human_250mg[5,1]<-as.data.frame(results_250mg_inhalation_single_human_C_RP_AUC$PPORRES[1])
Plot_AUC_Human_250mg[5,2]<-as.data.frame(results_250mg_inhalation_single_human_C_RP_AUC$PPORRES[2])

Plot_AUC_Human_250mg[6,1]<-as.data.frame(results_250mg_inhalation_single_human_C_F_AUC$PPORRES[1])
Plot_AUC_Human_250mg[6,2]<-as.data.frame(results_250mg_inhalation_single_human_C_F_AUC$PPORRES[2])

Plot_AUC_Human_250mg[7,1]<-as.data.frame(results_250mg_inhalation_single_human_C_Pu_AUC$PPORRES[1])
Plot_AUC_Human_250mg[7,2]<-as.data.frame(results_250mg_inhalation_single_human_C_Pu_AUC$PPORRES[2])

colnames(Plot_AUC_Human_250mg)<- c("AUC","Cmax") 
rownames(Plot_AUC_Human_250mg)<- c("C_B","C_L","C_SI","C_SP","C_RP","C_F","C_Pu")

plot(Plot_AUC_Human_250mg$AUC)
plot(Plot_AUC_Human_250mg$Cmax)

#Plotting rat and human results of the auc calculations
results_250mg_oral_single_rat_C_B_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_B_AUC")
results_250mg_oral_single_rat_C_L_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_L_AUC")
results_250mg_oral_single_rat_C_SI_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_SI_AUC")
results_250mg_oral_single_rat_C_SP_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_SP_AUC")
results_250mg_oral_single_rat_C_RP_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_RP_AUC")
results_250mg_oral_single_rat_C_F_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_F_AUC")
results_250mg_oral_single_rat_C_Pu_AUC <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/results_250mg_oral_single_rat_C_Pu_AUC")


#--------------------Population model-----------------

#---------------------Lung compartment-------------------#
sub_set <- solve.pbk[1:482000,c(1,2,12)]

conc_C <- PKNCAconc(sub_set, C_Pu~time|id)

#Dosing data per subject is part of the parameter file but it is missing sim id and the time variable 
#these will be added.
#Oral dose extraction
dose_extraction <- as.data.frame(parameters[,58])
sim_extraction <- unique(solve.pbk[solve.pbk$time == 0,c("time", "id")])

#Combining into 1 file that can be used
d_dose <- cbind(sim_extraction,dose_extraction$`parameters[, 58]`)                       
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
write.csv(results_obj_manual$result, "results_250_oral_C_Pu.csv")




#---------------------Liver compartment-------------------#
#Liver compartment concentration
sub_set <- solve.pbk[1:482000,c(1,2,49)]
conc_C <- PKNCAconc(sub_set, C_L~time|id)

#Dosing data per subject is part of the parameter file but it is missing sim id and the time variable 
#these will be added.
#Oral dose extraction
dose_extraction <- as.data.frame(parameters[,58])
sim_extraction <- unique(solve.pbk[solve.pbk$time == 0,c("time", "id")])

#Combining into 1 file that can be used
d_dose <- cbind(sim_extraction,dose_extraction$`parameters[, 58]`)                       
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
write.csv(results_obj_manual$result, "results_250_oral_C_L.csv")


#---------------------Blood compartment-------------------#
#Blood compartment concentration
sub_set <- solve.pbk[1:482000,c(1,2,4,5)]

#Combining Arterial and venous blood into one general blood compartment
sub_set[3]<- sub_set[3]+sub_set[4]

#Dropping column as data is already added to column 3
sub_set[4] <- NULL
#Renamning the columns
colnames(sub_set) <- c("id","time","C_B")
conc_C <- PKNCAconc(sub_set, C_B~time|id)

#Dosing data per subject is part of the parameter file but it is missing sim id and the time variable 
#these will be added.
#Oral dose extraction
dose_extraction <- as.data.frame(parameters[,58])
sim_extraction <- unique(solve.pbk[solve.pbk$time == 0,c("time", "id")])

#Combining into 1 file that can be used
d_dose <- cbind(sim_extraction,dose_extraction$`parameters[, 58]`)                       
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

#---------------------Slowly compartment-------------------#
#Slowly perfused compartment
sub_set <- solve.pbk[1:482000,c(1,2,30)]
conc_C <- PKNCAconc(sub_set, C_SP~time|id)

#Dosing data per subject is part of the parameter file but it is missing sim id and the time variable 
#these will be added.
#Oral dose extraction
dose_extraction <- as.data.frame(parameters[,58])
sim_extraction <- unique(solve.pbk[solve.pbk$time == 0,c("time", "id")])

#Combining into 1 file that can be used
d_dose <- cbind(sim_extraction,dose_extraction$`parameters[, 58]`)                       
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
write.csv(results_obj_manual$result, "results_250_oral_C_SP.csv")

#---------------------Richly  compartment-------------------#
#richly perfused compartment
sub_set <- solve.pbk[1:482000,c(1,2,24)]
conc_C <- PKNCAconc(sub_set, C_RP~time|id)

#Dosing data per subject is part of the parameter file but it is missing sim id and the time variable 
#these will be added.
#Oral dose extraction
dose_extraction <- as.data.frame(parameters[,58])
sim_extraction <- unique(solve.pbk[solve.pbk$time == 0,c("time", "id")])

#Combining into 1 file that can be used
d_dose <- cbind(sim_extraction,dose_extraction$`parameters[, 58]`)                       
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
write.csv(results_obj_manual$result, "results_250_oral_C_RP.csv")

#---------------------Fat compartment-------------------#
#Fat compartment
sub_set <- solve.pbk[1:482000,c(1,2,18)]
conc_C <- PKNCAconc(sub_set, C_F~time|id)

#Dosing data per subject is part of the parameter file but it is missing sim id and the time variable 
#these will be added.
#Oral dose extraction
dose_extraction <- as.data.frame(parameters[,58])
sim_extraction <- unique(solve.pbk[solve.pbk$time == 0,c("time", "id")])

#Combining into 1 file that can be used
d_dose <- cbind(sim_extraction,dose_extraction$`parameters[, 58]`)                       
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
write.csv(results_obj_manual$result, "results_250_oral_C_F.csv")

#---------------------Small Intestine compartment-------------------#

#SI compartment
sub_set <- solve.pbk[1:482000,c(1,2,36)]
conc_C <- PKNCAconc(sub_set, C_SI~time|id)

#Dosing data per subject is part of the parameter file but it is missing sim id and the time variable 
#these will be added.
#Oral dose extraction
dose_extraction <- as.data.frame(parameters[,58])
sim_extraction <- unique(solve.pbk[solve.pbk$time == 0,c("time", "id")])

#Combining into 1 file that can be used
d_dose <- cbind(sim_extraction,dose_extraction$`parameters[, 58]`)                       
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
write.csv(results_obj_manual$result, "results_250_oral_C_SI.csv")




#---------------------------Extracting the calculated PKNCA data----------
#------C_Pu data extraction---#
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
        ylab= "umol/l-hr",
        names= c("Lung oral", "Lung inhalation","Blood oral","Blood inhalation", "Fat oral","Fat inhalation", "RP oral","RP inhalation", "SP oral", "SP inhalation", "SI oral", "SI inhalation",
                 "Liver oral", "SI inhalation"),
        col="orange",
        las=2)




