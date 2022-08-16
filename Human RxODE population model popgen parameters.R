#author: Joris Jean van der Lugt
#date: 27-01-2021
#Human cinnamaldehyde pbk Model adapted from:  "Dose-dependent DNA adduct formation by cinnamaldehyde and other food-borne α,β-unsaturated aldehydes predicted by physiologically based in silico modelling"

library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(shiny)
library(truncnorm)
library(reshape2)
library(PKNCA)

#Simulations
set.seed(15204)         #to ensure a reproducible output
amount.units            <-"umol"
time.units              <-"h"
nbr.doses               <-1        #number of doses
time.0                  <-0        #time start dosing
time.end                <-8        #time end of simulation
time.frame              <-0.1     #time steps of simulation
N                       <-1000     #Number of males
NF                      <-1000     #Number of females
Oral_Dose_in_mg_bw      <-250      #Dose in mg/kg-bw
Inhalation_Dose_in_mg_bw<-5        #The inhaled dose in mg/kg
Volume_exposure_chamber <-10       #volume exposure chamber in L
MW                      <-132.16   #The molecular weight of Cinnamaldehyde



#Importing the popgen data set 
pg_m_par <- read_csv("popgen male parameters version 06-04-2022.csv", 
                                                      col_types = cols(Age = col_number(), 
                                                                       `Body Mass` = col_number(), Height = col_number(), 
                                                                       `Cardiac Output` = col_number(), 
                                                                       `Lung mass` = col_number(), `Lung flow` = col_number(), 
                                                                       `Liver mass` = col_number(), `Liver flow` = col_number(), 
                                                                       `Small intestine mass` = col_number(), 
                                                                       `Small intestine flow` = col_number(), 
                                                                       `Adipose mass` = col_number(), `Adipose flow` = col_number(), 
                                                                       `Liver Total flow` = col_number(), 
                                                                       `Slowly Perfused mass` = col_number(), 
                                                                       `Slowly Perfused flow` = col_number(), 
                                                                       `Richly Perfused mass` = col_number(), 
                                                                       `Richly Perfused flow` = col_number(), 
                                                                       `Lung Bronchial flow` = col_number()))

pg_f_par <-read_csv("popgen female parameters version 09-08-2022.csv", 
                            col_types = cols(Age = col_number(), 
                                     `Body Mass` = col_number(), Height = col_number(), 
                                     `Cardiac Output` = col_number(), 
                                     `Lung mass` = col_number(), `Lung flow` = col_number(), 
                                     `Liver mass` = col_number(), `Liver flow` = col_number(), 
                                     `Small intestine mass` = col_number(), 
                                     `Small intestine flow` = col_number(), 
                                     `Adipose mass` = col_number(), `Adipose flow` = col_number(), 
                                     `Liver Total flow` = col_number(), 
                                     `Slowly Perfused mass` = col_number(), 
                                     `Slowly Perfused flow` = col_number(), 
                                     `Richly Perfused mass` = col_number(), 
                                     `Richly Perfused flow` = col_number(), 
                                     `Lung Bronchial flow` = col_number()))


colnames <-c("Age","Height","BW","V_L","V_F",
             "V_B","V_A","V_V","V_SI","V_Pu","V_RP","V_SP","Q_C","Q_F","Q_L","Q_SI","Q_RP","Q_SP")

par_var_m_pop <- length(colnames)
par_var_f_pop <- length(colnames)

#create data frames for population males
var_m_pop <- matrix(NA, nrow = N, ncol = par_var_m_pop)
colnames(var_m_pop) <- colnames
var_m_pop <- as.data.frame(var_m_pop)
#create data frames for the population of females
var_f_pop <- matrix(NA, nrow = NF, ncol = par_var_f_pop)
colnames(var_f_pop) <- colnames
var_f_pop <- as.data.frame(var_f_pop)


#--Physico-chemical parameters--#
#-Cinnamaldehyde-#
var_m_pop$P_F      <-  1.62    #Fat/Blood partition coefficient
var_m_pop$P_L      <-  0.59    #Fat/Blood partition coefficient
var_m_pop$P_SI     <-  0.59    #Small intestine/Blood partition coefficients
var_m_pop$P_RP     <-  0.59    #Richly perfused tissues/Blood partition coefficients
var_m_pop$P_B      <-  1.25E5  #Blood/Air Partition Coefficient 
var_m_pop$P_SP     <-  0.78    #Slowly perfused tissues/Blood partition coefficients
var_m_pop$P_Pu     <-  0.59    #Lung/Blood partition coefficients

#-Cinnamyl Alcohol-#
var_m_pop$P_OH_F    <-  1.64   #Fat/Blood partition coefficient
var_m_pop$P_OH_L    <-  0.59   #Fat/Blood partition coefficient
var_m_pop$P_OH_SI   <-  0.59   #Small intestine/Blood partition coefficients
var_m_pop$P_OH_RP   <-  0.59   #Richly perfused tissues/Blood partition coefficients
var_m_pop$P_OH_SP   <-  0.78   #Slowly perfused tissues/Blood partition coefficients
var_m_pop$P_OH_Pu   <-  0.59   #Lung/Blood partition coefficients

#-----Cinnamaldehyde----#
var_f_pop$P_F       <-  1.62   #Fat/Blood partition coefficient
var_f_pop$P_L       <-  0.59   #Fat/Blood partition coefficient
var_f_pop$P_SI      <-  0.59   #Small intestine/Blood partition coefficients
var_f_pop$P_RP      <-  0.59   #Richly perfused tissues/Blood partition coefficients
var_f_pop$P_B       <-  1.25E5 #Blood/Air Partition Coefficient 
var_f_pop$P_SP      <-  0.78   #Slowly perfused tissues/Blood partition coefficients
var_f_pop$P_Pu      <-  0.59   #Lung/Blood partition coefficients

#-Cinnamyl Alcohol-#
var_f_pop$P_OH_F    <-  1.64   #Fat/Blood partition coefficient
var_f_pop$P_OH_L    <-  0.59   #Fat/Blood partition coefficient
var_f_pop$P_OH_SI   <-  0.59   #Small intestine/Blood partition coefficients
var_f_pop$P_OH_RP   <-  0.59   #Richly perfused tissues/Blood partition coefficients
var_f_pop$P_OH_SP   <-  0.78   #Slowly perfused tissues/Blood partition coefficients
var_f_pop$P_OH_Pu   <-  0.59   #Lung/Blood partition coefficients

#--Physiological Parameters--#

#Population specific parameters (Male)
var_m_pop$Age              <- pg_m_par$Age                            #Age (years)           
var_m_pop$Height           <- pg_m_par$Height                         #Body height (cm)                                     
var_m_pop$BW               <- pg_m_par$'Body Mass'                    #Body weight (kg)

#-Tissues volumes in % body weight-#

var_m_pop$V_L      <-pg_m_par$`Liver mass`              #Volume liver tissue (l)
var_m_pop$V_F      <-pg_m_par$`Adipose mass`            #Volume adipose tissue (L)
var_m_pop$V_B      <-(((13.1 * var_m_pop$Height + 18.05 * var_m_pop$BW - 480) / 0.5723) / 1000)#Volume blood (L)
var_m_pop$V_A      <-var_m_pop$V_B / 3                  #Volume arterial blood (L)
var_m_pop$V_V      <-var_m_pop$V_B * (2/3)              #Volume venous blood (L) 
var_m_pop$V_SI     <-pg_m_par$`Small intestine mass`    #Volume Small intestine (L)
var_m_pop$V_Pu     <-pg_m_par$`Lung mass`
var_m_pop$V_RP     <-pg_m_par$`Richly Perfused mass`    #Volume richly perfused tissue (L)
var_m_pop$V_SP     <-pg_m_par$`Slowly Perfused mass`    #Volume slowly perfused tissue (L)

#-Cardiac parameters-#
var_m_pop$Q_C           <- pg_m_par$`Cardiac Output`           #Cardiac output (L/h)
var_m_pop$Q_SI          <- pg_m_par$`Small intestine flow`                                         #Blood flow to the gut (L/h)
var_m_pop$Q_F           <- pg_m_par$`Adipose flow`                                         #Blood flow to adipose tissue (L/h)
var_m_pop$Q_L           <- pg_m_par$`Liver flow`                                        #Blood flow to liver(L/h)
var_m_pop$Q_Pu          <- pg_m_par$`Lung flow`
var_m_pop$Q_RP          <- pg_m_par$`Richly Perfused flow`             #Blood flow to richly perfused tissue (L/h)
var_m_pop$Q_SP          <- pg_m_par$`Slowly Perfused flow`                            #Blood flow to slowly perfused tissue (L/h)
var_m_pop$check         <- var_m_pop$Q_SI + var_m_pop$Q_F + var_m_pop$Q_L + var_m_pop$Q_RP + var_m_pop$Q_SP

#Pulmonary ventilation
var_m_pop$P_V           <-rnorm(N,mean=540,sd=3)

#----GSH parameters----#
#--GSH synthesis in umol/kg tissue/h--#

var_m_pop$G_SYN_L     <- 1122  #Liver 
var_m_pop$G_SYN_SI    <- 27    #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
var_m_pop$k_L_GLOS    <- 0.142 #Liver
var_m_pop$k_SI_GLOS   <- 0.044 #Small intestine

#--Initial GSH concentration--#
var_m_pop$init_GSH_L  <- 5639 * var_m_pop$V_L  #initial GSH concentration in the liver in umol
var_m_pop$init_GSH_SI <- 1250 * var_m_pop$V_SI  #initial GSH concentration in the small intestine in umol

var_m_pop$k_GSH <- 6.6 * 10^(-4) #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h
var_m_pop$k_DNA <- 1.6 * 10^(-8) #The second-order rate constant of the reaction between cinnamaldehyde and 2ʹ-dG in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
var_m_pop$C_PRO_L     <- 5319  #Liver
var_m_pop$C_PRO_SI    <- 245   #Small intestine

#----DNA parameters----#
var_m_pop$C_L_dG     <-  1.36 #Concentration of 2ʹ-dG in the liver μmol/kg liver
var_m_pop$T_0.5      <-  38.5   #Half-life of DNA adduct in the liver in hours


#--Chemical parameters--#
var_m_pop$Ka <- 5.0  #Absorption rate constant for uptake in the Small intestine in per H

#----Liver----#

#-first rate order constants-#
var_m_pop$k_L_OH  <- 4.2e-02   #Scaled first rate order constant for the enzymatic oxidation of cinnamyl alcohol in the liver in umol/h

#--Michaelis menten constants--#
var_m_pop$Km_L_CA     <-  8.5  #Km for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the liver in μM
var_m_pop$Km_L_AO     <-  330  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the liver in μM
var_m_pop$Km_L_GST    <-  100  #Km for enzymatic conjugation of cinnamaldehyde with GST in the liver in μM  
var_m_pop$Km_L_GST_G  <-  1.7*10^3  #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the small intestine μM

#--Vmax values--#
var_m_pop$Vsmax_L_CA    <-  9.7  #Scaled Vmax for enzymatic oxidation of cinnamaldehyde in the liver in μmol/h 
var_m_pop$Vsmax_L_AO    <-  73   #Scaled Vmax for enzymatic reduction of cinnamaldehyde in the liver in μmol/h
var_m_pop$Vsmax_L_GST   <-  37   #Scaled Vmax for enzymatic conjugation of cinnamaldehyde with GSH in the liver in μmol/h

#----Small intestines----#
#--Michaelis menten constants--#
var_m_pop$Km_SI_CA    <- 70  #Km for enzymatic oxidation of cinnamaldehyde into cinnamic acid in the Small Intestine in μM
var_m_pop$Km_SI_AO    <- 90  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the Small Intestine in μM
var_m_pop$Km_SI_OH    <- 290 #Km for enzymatic oxidation of cinnamly alcOHol into cinnamaldehyde in the Small Intestine in μM
var_m_pop$Km_SI_GST   <- 600 #Km for enzymatic conjugation of cinnamaldehye with GST in the Small Intestine in μM RAT value
var_m_pop$Km_SI_GST_G <- 100 #Km toward cinnamaldehyde for enzymatic conjugation of cinnamaldehyde in the small intestine μM

#-Vmax values-#
var_m_pop$Vsmax_SI_CA    <- 21 #Scaled Vmax for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the Small Intestine in μmol/h 
var_m_pop$Vsmax_SI_AO    <- 30 #Scaled Vmax for enzymatic reduction of cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
var_m_pop$Vsmax_SI_OH    <- 5.0 #Scaled Vmax for enzymatic Oxidation of cinnamyl alcohol into cinnamaldehyde in the Small Intestine in μmol/h 
var_m_pop$Vsmax_SI_GST   <- 63 #Scaled Vmax for enzymatic Conjugation of cinnamaldehyde with GSH in the in the small intestine in μmol/h RAT value

#---Dose male---#
var_m_pop$Oral_Dose       <- (Oral_Dose_in_mg_bw * var_m_pop$BW)/ MW  * 1e+3     #The administered dose in umol 
var_m_pop$Inhalation_Dose <- (Inhalation_Dose_in_mg_bw * var_m_pop$BW)/ MW  * 1e+3 #The inhaled dose in μmol

var_m_pop$Volume_exposure_chamber <-Volume_exposure_chamber
##------------Population specific parameters (female)--------------------##                                     
var_f_pop$Age              <- pg_f_par$Age                                        #Age (years)   
var_f_pop$Height           <- pg_f_par$Height                                     #Body height baseline (cm)
var_f_pop$BW               <- pg_f_par$`Body Mass`                                #Body weight (kg)

#-Tissues volumes in % body weight-#

var_f_pop$V_L       <-pg_f_par$`Liver mass`             #Volume liver tissue (l)
var_f_pop$V_F       <-pg_f_par$`Adipose mass`           #Volume adipose tissue (L)
var_f_pop$V_B       <-(((35.5 * var_f_pop$Height + 2.27 * var_f_pop$BW - 3382)/ 0.6178 )/ 1000) #Volume blood (L)
var_f_pop$V_A       <-var_f_pop$V_B / 3                 #Volume arterial blood (L)
var_f_pop$V_V       <-var_f_pop$V_B * (2/3)             #Volume venous blood (L) 
var_f_pop$V_SI      <-pg_f_par$`Small intestine mass`   #Volume small intestine (L)
var_f_pop$V_Pu      <-pg_f_par$`Lung mass`              #Volume lung (L)
var_f_pop$V_RP      <-pg_f_par$`Richly Perfused flow`   #Volume richly perfused tissue (L)
var_f_pop$V_SP      <-pg_f_par$`Slowly Perfused mass`   #Volume slowly perfused tissue (L)

#-Cardiac parameters-#

var_f_pop$Q_C           <- pg_f_par$`Cardiac Output`                     #Cardiac output (L/h)
var_f_pop$Q_SI          <- pg_f_par$`Small intestine flow`               #Blood flow to the gut (L/h)
var_f_pop$Q_F           <- pg_f_par$`Adipose flow`                       #Blood flow to adipose tissue (L/h)
var_f_pop$Q_L           <- pg_f_par$`Liver flow`                         #Blood flow to liver via hepatic artery (L/h)
var_f_pop$Q_RP          <- pg_f_par$`Richly Perfused flow`               #Blood flow to richly perfused tissue (L/h)
var_f_pop$Q_SP          <- pg_f_par$`Slowly Perfused flow`               #Blood flow to Slowly perfused tissue (L/h)
var_f_pop$Q_Pu          <- pg_f_par$`Lung flow`                          #Blood flow to Lung tissue (L/h)
var_f_pop$check         <- var_f_pop$Q_SI + var_f_pop$Q_F + var_f_pop$Q_L + var_f_pop$Q_RP + var_f_pop$Q_SP + var_f_pop$Q_Pu  #check if total adds up to Q_C

#Pulmonary ventilation
var_f_pop$P_V           <-rnorm(N,mean=390,sd=3) #Pulmonary ventilation in L/h

#----GSH parameters female----#
#--GSH synthesis in umol/kg tissue/h--#

var_f_pop$G_SYN_L     <- 1122  #Liver 
var_f_pop$G_SYN_SI    <- 27    #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
var_f_pop$k_L_GLOS    <- 0.142 #Liver
var_f_pop$k_SI_GLOS   <- 0.044 #Small intestine

#--Initial GSH concentration--#
var_f_pop$init_GSH_L  <- 5639 * var_f_pop$V_L  #initial GSH concentration in the liver in umol/kg
var_f_pop$init_GSH_SI <- 1250 * var_f_pop$V_SI  #initial GSH concentration in the small intestine in umol/kg

var_f_pop$k_GSH <- 6.6 * 10^(-4) #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h
var_f_pop$k_DNA <- 1.6 * 10^(-8) #The second-order rate constant of the reaction between cinnamaldehyde and 2ʹ-dG in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
var_f_pop$C_PRO_L     <- 5319  #Liver
var_f_pop$C_PRO_SI    <- 245   #Small intestine

#----DNA parameters----#
var_f_pop$C_L_dG     <-  1.36 #Concentration of 2ʹ-dG in the liver μmol/kg liver
var_f_pop$T_0.5      <-  38.5   #Half-life of DNA adduct in the liver in hours


#--Chemical parameters--#
var_f_pop$Ka <- 5.0  #Absorption rate constant for uptake in the Small intestine in per H

#----Liver----#

#-first rate order constants-#
var_f_pop$k_L_OH  <- 4.2e-02   #Scaled first rate order constant for the enzymatic oxidation of cinnamyl alcohol in the liver in umol/h

#--Michaelis menten constants--#
var_f_pop$Km_L_CA     <-  8.5  #Km for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the liver in μM
var_f_pop$Km_L_AO     <-  330  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the liver in μM
var_f_pop$Km_L_GST    <-  100 #Km for enzymatic conjugation of cinnamaldehyde with GST in the liver in μM  
var_f_pop$Km_L_GST_G  <-  1.7*10^3 #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the liver (μM)

#--Vmax values--#
var_f_pop$Vsmax_L_CA    <-  9.7  #Scaled Vmax for enzymatic oxidation of cinnamaldehyde in the liver in μmol/h 
var_f_pop$Vsmax_L_AO    <-  73   #Scaled Vmax for enzymatic reduction of cinnamaldehyde in the liver in μmol/h
var_f_pop$Vsmax_L_GST   <-  37   #Scaled Vmax for enzymatic conjugation of cinnamaldehyde with GSH in the liver in μmol/h

#----Small intestines----#
#--Michaelis menten constants--#
var_f_pop$Km_SI_CA    <- 70  #Km for enzymatic oxidation of cinnamaldehyde into cinnamic acid in the Small Intestine in μM
var_f_pop$Km_SI_AO    <- 90  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the Small Intestine in μM
var_f_pop$Km_SI_OH    <- 290 #Km for enzymatic oxidation of cinnamly alcOHol into cinnamaldehyde in the Small Intestine in μM
var_f_pop$Km_SI_GST   <- 600 #Km for enzymatic conjugation of cinnamaldehye with GST in the Small Intestine in μM (RAT value)
var_f_pop$Km_SI_GST_G <- 100  #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the small intestine (μM)

#-Vmax values-#
var_f_pop$Vsmax_SI_CA    <- 21 #Scaled Vmax for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the Small Intestine in μmol/h 
var_f_pop$Vsmax_SI_AO    <- 30 #Scaled Vmax for enzymatic reduction of cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
var_f_pop$Vsmax_SI_OH    <- 5.0 #Scaled Vmax for enzymatic Oxidation of cinnamyl alcohol into cinnamaldehyde in the Small Intestine in μmol/h 
var_f_pop$Vsmax_SI_GST   <- 63 #Scaled Vmax for enzymatic Conjugation of cinnamaldehyde with GSH in the in the small intestine in μmol/h (RAT value)

#---Dose female---#
var_f_pop$Oral_Dose       <-(Oral_Dose_in_mg_bw * var_f_pop$BW)/ MW  * 1e+3     #The administered dose in umol 
var_f_pop$Inhalation_Dose <- (Inhalation_Dose_in_mg_bw * var_f_pop$BW)/ MW  * 1e+3 #The inhaled dose in μmol

var_f_pop$Volume_exposure_chamber <-Volume_exposure_chamber
#Combine data sets Male and Female for PBK model
phys <- rbind(var_m_pop,var_f_pop)

#ONLY MALE
#phys <- var_m_pop

#ONLY FEMALE
#  phys <- var_f_pop

P_F<-phys$P_F
P_L<-phys$P_L
P_SI<-phys$P_SI
P_RP<-phys$P_RP
P_B<-phys$P_B
P_SP<-phys$P_SP
P_Pu<-phys$P_Pu
P_OH_F<-phys$P_OH_F
P_OH_L<-phys$P_OH_L
P_OH_SI<-phys$P_OH_SI
P_OH_RP<-phys$P_OH_RP
P_OH_SP<-phys$P_OH_SP
P_OH_Pu<-phys$P_OH_Pu
Age <- phys$Age
Height <- phys$Height
BW<-phys$BW
V_L<-phys$V_L
V_F<-phys$V_F
V_B <- phys$V_B
V_A<-phys$V_A
V_V<-phys$V_V
V_SI<-phys$V_SI
V_Pu<-phys$V_Pu
V_RP<-phys$V_RP
V_SP<-phys$V_SP
Q_C<-phys$Q_C
Q_SI<-phys$Q_SI
Q_F<-phys$Q_F
Q_L<-phys$Q_L
Q_Pu<-phys$Q_Pu
Q_RP<-phys$Q_RP
Q_SP<-phys$Q_SP
P_V<-phys$P_V
G_SYN_L<-phys$G_SYN_L
G_SYN_SI<-phys$G_SYN_SI
k_L_GLOS<-phys$k_L_GLOS
k_SI_GLOS<-phys$k_SI_GLOS
init_GSH_L<-phys$init_GSH_L
init_GSH_SI<-phys$init_GSH_SI
k_GSH<-phys$k_GSH
k_DNA<-phys$k_DNA
C_PRO_L<-phys$C_PRO_L
C_PRO_SI<-phys$C_PRO_SI
C_L_dG<-phys$C_L_dG
T_0.5<-phys$T_0.5
Ka<-phys$Ka
k_L_OH <- phys$k_L_OH
Km_L_CA<-phys$Km_L_CA
Km_L_AO<-phys$Km_L_AO
Km_L_GST<-phys$Km_L_GST
Km_L_GST_G<-phys$Km_L_GST_G
Vsmax_L_CA<-phys$Vsmax_L_CA
Vsmax_L_AO<-phys$Vsmax_L_AO
Vsmax_L_GST<-phys$Vsmax_L_GST
Km_SI_CA<-phys$Km_SI_CA
Km_SI_AO<-phys$Km_SI_AO
Km_SI_OH<-phys$Km_SI_OH
Km_SI_GST<-phys$Km_SI_GST
Km_SI_GST_G<-phys$Km_SI_GST_G
Vsmax_SI_CA<-phys$Vsmax_SI_CA
Vsmax_SI_AO<-phys$Vsmax_SI_AO
Vsmax_SI_OH<-phys$Vsmax_SI_OH
Vsmax_SI_GST<-phys$Vsmax_SI_GST
Oral_Dose<- phys$Oral_Dose 
Inhalation_Dose<-phys$Inhalation_Dose
Volume_exposure_chamber<-phys$Volume_exposure_chamber


parameters=cbind(P_F,
                 P_L,
                 P_SI,
                 P_RP,
                 P_SP,
                 P_B,
                 P_Pu,
                 P_OH_F,
                 P_OH_L,
                 P_OH_SI,
                 P_OH_RP,
                 P_OH_SP,
                 P_OH_Pu,
                 Age,
                 Height,
                 BW,
                 V_F,
                 V_L,
                 V_SI,
                 V_B,
                 V_A,
                 V_V,
                 V_RP,
                 V_SP,
                 V_Pu,
                 Q_C,
                 Q_Pu,
                 Q_F,
                 Q_L,
                 Q_SI,
                 Q_RP,
                 Q_SP,
                 P_V,
                 G_SYN_L,
                 G_SYN_SI,
                 k_L_GLOS,
                 k_SI_GLOS,
                 init_GSH_L,
                 init_GSH_SI,
                 k_GSH,
                 k_DNA,
                 C_PRO_L,
                 C_PRO_SI,
                 C_L_dG,
                 T_0.5,
                 Ka,
                 k_L_OH,
                 Km_L_CA,
                 Km_L_AO,
                 Km_L_GST,
                 Km_L_GST_G,
                 Vsmax_L_CA,
                 Vsmax_L_AO,
                 Vsmax_L_GST,
                 Km_SI_CA,
                 Km_SI_AO,
                 Km_SI_OH,
                 Km_SI_GST,
                 Km_SI_GST_G,
                 Vsmax_SI_CA,
                 Vsmax_SI_AO,
                 Vsmax_SI_OH,
                 Vsmax_SI_GST,
                 Oral_Dose,
                 Inhalation_Dose,
                 Volume_exposure_chamber)

#defining the begin situation of the model Inhalation variation 
inits <- c("A_GI"         =0,
           "A_P_Art"      =0,
           "A_Inhalation" =0,
           "A_Exhalation" =0,
           "A_Pu"         =0,
           "A_OH_Pu"      =0,
           "A_V"          =0,
           "A_OH_V"       =0,
           "A_F"          =0,
           "A_OH_F"       =0,
           "AM_L_CA"      =0,
           "AM_L_AO"      =0,
           "AM_L_AG_GST"  =0,
           "AM_L_AG_CHEM" =0,
           "AM_L_AP"      =0,
           "AM_L_DA_FORM" =0,
           "AM_L_DA"      =0,
           "A_OH_M_L_C_A" =0,
           "A_OH_L"       =0,
           "A_L"          =0,
           "AM_Lc_GSH"    =init_GSH_L, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =init_GSH_SI,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#exposure
ex <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=1:6500, amt=(Oral_Dose_in_mg_bw) * phys$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=1:6500, amt=(Inhalation_Dose_in_mg_bw) * phys$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(seq(from = time.0, to = time.end, by = time.frame)) 
