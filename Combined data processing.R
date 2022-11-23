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


#Rxode
#cinnamaldehyde model Human
#Mass balance calculation rxode inhalation complete
mass_df_single <- solve.pbk/BW * MW /1e+3
mass_df_single <- mass_df_single[,c(66:79,81,82,85:89,93:97)]
mass_at_t_single <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df_single)){
  mass_at_t_single[nrow(mass_at_t_single) + 1,] <- rowSums(mass_df_single[i,])
}
plot(mass_at_t_single[,1])


#cinnamaldehyde model Rat
#Mass balance calculation rxode inhalation complete
mass_df <- solve.pbk_rat/BW * MW /1e+3
mass_df <- mass_df[,c(66:84,86,88:92,94,96,97)]
mass_at_t <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])


#Rxode 
#population mass balance
mass_df <- solve.pbk/phys[2,3] * MW /1e+3
mass_df <- mass_df[242:482,c(67:80,82,83,86:90,94:98)]
mass_at_t <- data.frame(mass=as.numeric())


for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])

results_human<-solve.pbk[c(1:241),c(1:98)]
plot(results_human$AM_L_CA)



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


#Metabolic output checking
#human model
#Checking proportions of urinary metabolites
#Urinary metabolites are all ultimately derived from the oxidation of Cinnamaldehyde.

data_urinary_metabolites <- solve.pbk[,c(1,78,87)]

total_dose <- Oral_Dose
total_dose<- Inhalation_Dose

total_urinary_metabolites_24h <- data_urinary_metabolites[241,2] + data_urinary_metabolites[241,3]

percentage_metabolised <- total_urinary_metabolites_24h/total_dose *100

#-------------population model
# data extraction  

meta_extraction <- unique(solve.pbk[solve.pbk$time == 24,c("time", "id","AM_L_CA","AM_SI_CA")])

dose_extraction <- unique(ex[ex$cmt == "A_GI",c("id","amt")])

percentage_metabolised_popgen <-c()
for (i in 1:2000){
  total <- meta_extraction[i,3] + meta_extraction[i,4]
  dose  <- dose_extraction[i,2]
  metabolised <- total/dose *100
  percentage_metabolised_popgen <-append(percentage_metabolised_popgen, metabolised)
}

metabolised_data_frame <- as.data.frame(melt(percentage_metabolised_popgen))

male_metabolised_data<- metabolised_data_frame[1:1000,]
female_metabolised_data<- metabolised_data_frame[1:2000,]

boxplot(male_metabolised_data,female_metabolised_data,
        ylab= "percentage metabolised into carboxylic acid metabolites",
        names= c("Male", "Female"),
        col="orange",
        las=2,
        title="Amount metabolized to carboxylic acid")
hist(male_metabolised_data)
hist(female_metabolised_data)




