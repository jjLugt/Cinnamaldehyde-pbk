#author: Joris Jean van der Lugt
#date: 05-08-2022
#Population based parameters adapted for Kasteel et al

library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(shiny)
library(truncnorm)
library(reshape2)

#Simulations
set.seed(15204)                  #to ensure a reproducible output by setting a seed for random components of the parameters
amount.units             <-"umol"
time.units               <-"h"
nbr.doses                <-1        #number of doses
time.0                   <-0        #time start dosing
time.end                 <-8        #time end of simulation
time.frame               <-0.01     #time steps of simulation
N                        <-1000     #Number of males
NF                       <-1000     #Number of females
Oral_dose_in_mg_bw       <-250      #Dose in mg/kg-bw
Inhalation_dose_in_mg_bw <-0        #The inhaled dose in mg/kg
MW                       <-132.16   #The molecular weight of Cinnamaldehyde


colnames <-c("Age","Height_start","Height_cv","Height","BW_start","BW_cv","BW","BSA","V_L","V_F","V_F_min",
             "V_B","V_A","V_V","V_SI","V_Pu","V_RP","V_SP","Q_C","Q_Pu","Q_F","Q_L","Q_SI","Q_RP","Q_SP")

par_var_m <- length(colnames)
par_var_f <- length(colnames)

#create data frames for population males
var_m <- matrix(NA, nrow = N, ncol = par_var_m)
colnames(var_m) <- colnames
var_m <- as.data.frame(var_m)
#create data frames for the population of females
var_f <- matrix(NA, nrow = NF, ncol = par_var_f)
colnames(var_f) <- colnames
var_f <- as.data.frame(var_f)


#--Physico-chemical parameters--#
#-Cinnamaldehyde-#
var_m$P_F      <-  1.62    #Fat/Blood partition coefficient
var_m$P_L      <-  0.59    #Fat/Blood partition coefficient
var_m$P_SI     <-  0.59    #Small intestine/Blood partition coefficients
var_m$P_RP     <-  0.59    #Richly perfused tissues/Blood partition coefficients
var_m$P_B      <-  1.25E5  #Blood/Air Partition Coefficient 
var_m$P_SP     <-  0.78    #Slowly perfused tissues/Blood partition coefficients
var_m$P_Pu     <-  0.59    #Lung/Blood partition coefficients

#-Cinnamyl Alcohol-#
var_m$P_OH_F    <-  1.64   #Fat/Blood partition coefficient
var_m$P_OH_L    <-  0.59   #Fat/Blood partition coefficient
var_m$P_OH_SI   <-  0.59   #Small intestine/Blood partition coefficients
var_m$P_OH_RP   <-  0.59   #Richly perfused tissues/Blood partition coefficients
var_m$P_OH_SP   <-  0.78   #Slowly perfused tissues/Blood partition coefficients
var_m$P_OH_Pu   <-  0.59   #Lung/Blood partition coefficients
#---Cinnamaldehyde----#
var_f$P_F       <-  1.62   #Fat/Blood partition coefficient
var_f$P_L       <-  0.59   #Fat/Blood partition coefficient
var_f$P_SI      <-  0.59   #Small intestine/Blood partition coefficients
var_f$P_RP      <-  0.59   #Richly perfused tissues/Blood partition coefficients
var_f$P_B       <-  1.25E5 #Blood/Air Partition Coefficient 
var_f$P_SP      <-  0.78   #Slowly perfused tissues/Blood partition coefficients
var_f$P_Pu      <-  0.59   #Lung/Blood partition coefficients

#-Cinnamyl Alcohol-#
var_f$P_OH_F    <-  40.5   #Fat/Blood partition coefficient
var_f$P_OH_L    <-  2.09   #Fat/Blood partition coefficient
var_f$P_OH_SI   <-  2.09   #Small intestine/Blood partition coefficients
var_f$P_OH_RP   <-  2.09   #Richly perfused tissues/Blood partition coefficients
var_f$P_OH_SP   <-  1.60   #Slowly perfused tissues/Blood partition coefficients
var_f$P_OH_Pu   <-  0.59   #Lung/Blood partition coefficients

#--Pyshiological Parameters--#
#Population specific parameters (Male)
Age                    <- runif(N,18,50)                                       #Age (years)
var_m$Age              <- Age
var_m$Height_start     <- 175.32 + 0.1113 * var_m$Age-0.0025 * var_m$Age^2     #Body height baseline (cm)
var_m$Height_cv        <- rnorm(N,0,0.039)                                     #Variation in body height
var_m$Height           <- var_m$Height_start * exp(var_m$Height_cv)            #Body height (cm)
var_m$BW_start         <- exp(2.643+0.0099 * var_m$Height)                     #Body weight baseline (kg)
var_m$BW_cv            <- rnorm(N,0,0.15)                                      #Variation in body weight
var_m$BW               <- var_m$BW_start * exp(var_m$BW_cv)                    #Body weight (kg)
var_m$BSA              <- 0.007184 * var_m$Height^0.725 * var_m$BW^0.425       #Body surface area (m2)

#-Tissues volumes in % body weight-#

var_m$V_L      <-(1072.8 * (var_m$BSA)-345.7) / 1000                               #Volume liver tissue (l)
var_m$V_F      <-(1.36 * var_m$BW)/(var_m$Height/100)-42                           #Volume adipose tissue (L)
var_m$V_F_min  <-0.05 * var_m$BW                                                   #Minimum of adipose tissue should be at least 5% of body weight
var_m$V_F      <-ifelse(var_m$V_F < var_m$V_F_min, var_m$V_F_min, var_m$V_F)       #To ensure that adipose tissue is at least 5% of body weight
var_m$V_B      <-(((13.1 * var_m$Height + 18.05 * var_m$BW - 480) / 0.5723) / 1000)#Volume blood (L)
var_m$V_A      <-var_m$V_B / 3                                                     #Volume arterial blood (L)
var_m$V_V      <-var_m$V_B * (2/3)                                                 #Volume venous blood (L) 
var_m$V_SI     <-0.021 * (var_m$BW - var_m$V_F * 0.92) / 1.05                      #Volume small intestine (L)
var_m$V_Pu     <-
var_m$V_RP     <-(2.331 * 10^-3 * var_m$Age + 0.1253 * var_m$BW^0.8477 + var_m$Height^0.3821 - 4.725) - var_m$V_SI - var_m$V_L   #Volume richly perfused tissue (L)
var_m$V_SP     <-var_m$BW - var_m$V_B - var_m$V_RP -var_m$V_SI - var_m$V_L - var_m$V_F  #Volume slowly perfused tissue (L)

#-Cardiac parameters-#

var_m$Q_C      <-var_m$BSA * 60 * (3 - 0.01 * (var_m$Age - 20))           #Cardiac output (L/h)
var_m$Q_Pu     <-var_m$Q_C                                                #Blood flow to the Lungs (L/h)
var_m$Q_SI     <-var_m$Q_C * 0.15                                         #Blood flow to the small intestine (L/h)
var_m$Q_F      <-var_m$Q_C * 0.05                                         #Blood flow to adipose tissue (L/h)
var_m$Q_L      <-var_m$Q_C * 0.065                                        #Blood flow to liver (L/h)
var_m$Q_RP     <-0.626 * var_m$Q_C - var_m$Q_SI - var_m$Q_L               #Blood flow to richly perfused tissue (L/h)
var_m$Q_SP     <-0.374 * var_m$Q_C - var_m$Q_F                            #Blood flow to slowly perfused tissue (L/h)

var_m$P_V      <-rnorm(N,mean=540,sd=3)                                   #Pulmonary ventilation (L/h)
#----GSH parameters----#
#--GSH synthesis in umol/kg tissue/h--#

var_m$G_SYN_L     <- 1122  #Liver 
var_m$G_SYN_SI    <- 27    #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
var_m$k_L_GLOS    <- 0.142 #Liver
var_m$k_SI_GLOS   <- 0.044 #Small intestine

#--Initial GSH concentration--#
var_m$init_GSH_L  <- 5639 * var_m$V_L  #initial GSH concentration in the liver in umol
var_m$init_GSH_SI <- 1250 * var_m$V_SI  #initial GSH concentration in the small intestine in umol

var_m$k_GSH <- 6.6 * 10^(-4) #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h
var_m$k_DNA <- 1.6 * 10^(-8) #The second-order rate constant of the reaction between cinnamaldehyde and 2ʹ-dG in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
var_m$C_PRO_L     <- 5319  #Liver
var_m$C_PRO_SI    <- 245   #Small intestine

#----DNA parameters----#
var_m$C_L_dG     <-  1.36 #Concentration of 2ʹ-dG in the liver μmol/kg liver
var_m$T_0.5      <-  38.5   #Half-life of DNA adduct in the liver in hours


#--Chemical parameters--#
var_m$Ka <- 5.0  #Absorption rate constant for uptake in the Small intestine in per H

#----Liver----#

#-first rate order constants-#
var_m$k_L_OH  <- 4.2e-02   #Scaled first rate order constant for the enzymatic oxidation of cinnamyl alcohol in the liver in umol/h

#--Michaelis menten constants--#
var_m$Km_L_CA     <-  8.5  #Km for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the liver in μM
var_m$Km_L_AO     <-  330  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the liver in μM
var_m$Km_L_GST    <-  100  #Km for enzymatic conjugation of cinnamaldehyde with GST in the liver in μM  
var_m$Km_L_GST_G  <-  1.7*10^3  #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the small intestine μM

#--Vmax values--#
var_m$Vsmax_L_CA    <-  9.7  #Scaled Vmax for enzymatic oxidation of cinnamaldehyde in the liver in μmol/h 
var_m$Vsmax_L_AO    <-  73   #Scaled Vmax for enzymatic reduction of cinnamaldehyde in the liver in μmol/h
var_m$Vsmax_L_GST   <-  37   #Scaled Vmax for enzymatic conjugation of cinnamaldehyde with GSH in the liver in μmol/h

#----Small intestines----#
#--Michaelis menten constants--#
var_m$Km_SI_CA    <- 70  #Km for enzymatic oxidation of cinnamaldehyde into cinnamic acid in the Small Intestine in μM
var_m$Km_SI_AO    <- 90  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the Small Intestine in μM
var_m$Km_SI_OH    <- 290 #Km for enzymatic oxidation of cinnamly alcOHol into cinnamaldehyde in the Small Intestine in μM
var_m$Km_SI_GST   <- 600 #Km for enzymatic conjugation of cinnamaldehye with GST in the Small Intestine in μM RAT value
var_m$Km_SI_GST_G <- 0 #Km toward cinnamaldehyde for enzymatic conjugation of cinnamaldehyde in the small intestine μM

#-Vmax values-#
var_m$Vsmax_SI_CA    <- 21 #Scaled Vmax for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the Small Intestine in μmol/h 
var_m$Vsmax_SI_AO    <- 30 #Scaled Vmax for enzymatic reduction of cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
var_m$Vsmax_SI_OH    <- 5.0 #Scaled Vmax for enzymatic Oxidation of cinnamyl alcohol into cinnamaldehyde in the Small Intestine in μmol/h 
var_m$Vsmax_SI_GST   <- 63 #Scaled Vmax for enzymatic Conjugation of cinnamaldehyde with GSH in the in the small intestine in μmol/h RAT value

#---Dose male---#
var_m$Oral_Dose       <- (Dose_in_mg_bw * var_m$BW)/ MW  * 1e+3     #The administered dose in umol 
var_m$Inhalation_Dose <- (Inhalation_dose_in_mg_bw * var_m$BW)/ MW  * 1e+3 #The inhaled dose in μmol



##------------Population specific parameters (female)--------------------##                                     
var_f$Age              <- Age                                                  #Age (years)
var_f$Height_start     <- 161.66 + 0.1319 * var_f$Age - 0.0027*var_f$Age^2     #Body height baseline (cm)
var_f$Height_cv        <- rnorm(N,0,0.039)                                     #Variation in body height
var_f$Height           <- var_f$Height_start * exp(var_f$Height_cv)            #Body height (cm)
var_f$BW_start         <- exp(2.7383+0.0091 * var_f$Height)                    #Body weight baseline (kg)
var_f$BW_cv            <- rnorm(N,0,0.188)                                     #Variation in body weight
var_f$BW               <- var_f$BW_start * exp(var_f$BW_cv)                    #Body weight (kg)
var_f$BSA              <- 0.007184 * var_f$Height^0.725 * var_f$BW^0.425       #Body surface area (m2)

#-Tissues volumes in % body weight-#

var_f$V_L       <- (1072.8 * (var_f$BSA)-345.7) / 1000                              #Volume liver tissue (l)
var_f$V_F       <- (1.61*var_f$BW)/(var_f$Height/100)-38.3                          #Volume adipose tissue (L)
var_f$V_F_min   <- 0.05 * var_f$BW                                                  #Minimum of adipose tissue should be at least 5% of body weight
var_f$V_F       <- ifelse(var_f$V_F < var_f$V_F_min, var_f$V_F_min, var_f$V_F)      #To ensure that adipose tissue is at least 5% of body weight
var_f$V_B       <-(((35.5 * var_f$Height + 2.27 * var_f$BW - 3382)/ 0.6178 )/ 1000) #Volume blood (L)
var_f$V_A       <-var_f$V_B / 3                                                     #Volume arterial blood (L)
var_f$V_V       <-var_f$V_B * (2/3)                                                 #Volume venous blood (L) 
var_f$V_SI      <-0.021 * (var_f$BW - var_f$V_F * 0.92) / 1.05                           #Volume gut tissue (L)
var_f$V_Pu      <-
var_f$V_RP      <-(2.331 * 10^-3 * var_f$Age + 0.1253 * var_f$BW^0.8477 + var_f$Height^0.3821 - 4.725) - var_f$V_SI - var_f$V_L   #Volume richly perfused tissue (L)
var_f$V_SP      <-var_f$BW - var_f$V_B - var_f$V_RP -var_f$V_SI - var_f$V_L - var_f$V_F  #Volume slowly perfused tissue (L)

#-Cardiac parameters-#

var_f$Q_C           <-var_f$BSA * 60 * (3 - 0.01 * (var_f$Age - 20))           #Cardiac output (L/h)
var_f$Q_Pu          <-var_f$Q_C                                                #Blood flow to the Lung (L/h) 
var_f$Q_SI          <-var_f$Q_C * 0.17                                         #Blood flow to the Small intestine (L/h)
var_f$Q_F           <-var_f$Q_C * 0.085                                        #Blood flow to adipose tissue (L/h)
var_f$Q_L           <-var_f$Q_C * 0.065                                        #Blood flow to liver (L/h)
var_f$Q_RP          <-0.626 * var_f$Q_C - var_f$Q_SI - var_f$Q_L               #Blood flow to richly perfused tissue (L/h)
var_f$Q_SP          <-0.374 * var_f$Q_C - var_f$Q_F    

var_f$P_V           <-rnorm(N,mean=390,sd=3)  #Pulmonary ventilation

#----GSH parameters female----#
#--GSH synthesis in μmol/kg tissue/h--#

var_f$G_SYN_L     <- 1122  #Liver 
var_f$G_SYN_SI    <- 27    #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
var_f$k_L_GLOS    <- 0.142 #Liver
var_f$k_SI_GLOS   <- 0.044 #Small intestine

#--Initial GSH concentration--#
var_f$init_GSH_L  <- 5639 * var_f$V_L  #initial GSH concentration in the liver in umol/kg
var_f$init_GSH_SI <- 1250 * var_f$V_SI  #initial GSH concentration in the small intestine in umol/kg

var_f$k_GSH <- 6.6 * 10^(-4) #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h
var_f$k_DNA <- 1.6 * 10^(-8) #The second-order rate constant of the reaction between cinnamaldehyde and 2ʹ-dG in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
var_f$C_PRO_L     <- 5319  #Liver
var_f$C_PRO_SI    <- 245   #Small intestine

#----DNA parameters----#
var_f$C_L_dG     <-  1.36 #Concentration of 2ʹ-dG in the liver μmol/kg liver
var_f$T_0.5      <-  38.5   #Half-life of DNA adduct in the liver in hours


#--Chemical parameters--#
var_f$Ka <- 5.0  #Absorption rate constant for uptake in the Small intestine in per H

#----Liver----#

#-first rate order constants-#
var_f$k_L_OH  <- 4.2e-02   #Scaled first rate order constant for the enzymatic oxidation of cinnamyl alcohol in the liver in umol/h

#--Michaelis menten constants--#
var_f$Km_L_CA     <-  8.5  #Km for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the liver in μM
var_f$Km_L_AO     <-  330  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the liver in μM
var_f$Km_L_GST    <-  100 #Km for enzymatic conjugation of cinnamaldehyde with GST in the liver in μM  
var_f$Km_L_GST_G  <-  1.7*10^3 #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the liver (μM)

#--Vmax values--#
var_f$Vsmax_L_CA    <-  9.7  #Scaled Vmax for enzymatic oxidation of cinnamaldehyde in the liver in μmol/h 
var_f$Vsmax_L_AO    <-  73   #Scaled Vmax for enzymatic reduction of cinnamaldehyde in the liver in μmol/h
var_f$Vsmax_L_GST   <-  37   #Scaled Vmax for enzymatic conjugation of cinnamaldehyde with GSH in the liver in μmol/h

#----Small intestines----#
#--Michaelis menten constants--#
var_f$Km_SI_CA    <- 70  #Km for enzymatic oxidation of cinnamaldehyde into cinnamic acid in the Small Intestine in μM
var_f$Km_SI_AO    <- 90  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the Small Intestine in μM
var_f$Km_SI_OH    <- 290 #Km for enzymatic oxidation of cinnamly alcOHol into cinnamaldehyde in the Small Intestine in μM
var_f$Km_SI_GST   <- 600 #Km for enzymatic conjugation of cinnamaldehye with GST in the Small Intestine in μM (RAT value)
var_f$Km_SI_GST_G <- 0  #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the small intestine (μM)

#-Vmax values-#
var_f$Vsmax_SI_CA    <- 21 #Scaled Vmax for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the Small Intestine in μmol/h 
var_f$Vsmax_SI_AO    <- 30 #Scaled Vmax for enzymatic reduction of cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
var_f$Vsmax_SI_OH    <- 5.0 #Scaled Vmax for enzymatic Oxidation of cinnamyl alcohol into cinnamaldehyde in the Small Intestine in μmol/h 
var_f$Vsmax_SI_GST   <- 63 #Scaled Vmax for enzymatic Conjugation of cinnamaldehyde with GSH in the in the small intestine in μmol/h (RAT value)

#---Dose female---#
var_f$Oral_Dose <- (Dose_in_mg_bw * var_f$BW)/ MW  * 1e+3     #The administered dose in umol 
var_f$Inhalation_Dose <- (Inhalation_dose_in_mg_bw * var_f$BW)/ MW  * 1e+3 #The inhaled d


#Combine datasets Male and Female for PBPK model
phys <- rbind(var_m,var_f)

#ONLY MALE
#phys <- var_m

#ONLY FEMALE
#  phys <- var_mf

P_F<-phys$P_F
P_L<-phys$P_L
P_SI<-phys$P_SI
P_B<-phys$P_B
P_RP<-phys$P_RP
P_SP<-phys$P_SP
P_OH_F<-phys$P_OH_F
P_OH_L<-phys$P_OH_L
P_OH_SI<-phys$P_OH_SI
P_OH_RP<-phys$P_OH_RP
P_OH_SP<-phys$P_OH_SP
P_OH_Pu<-phys$P_OH_Pu
Age <- phys$Age
Height_start <- phys$Height_start
Height_cv <- phys$Height_cv
Height <- phys$Height
BW_start <- phys$BW_start
BW_cv <- phys$BW_cv
BW<-phys$BW
BSA <- phys$BSA
V_L<-phys$V_L
V_F<-phys$V_F
V_F_min<-phys$V_F_min
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
Q_Pu<-phys$P_Pu
Q_RP<-phys$Q_RP
Q_SP<-phys$Q_SP
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
Oral_dose<- phys$Oral_Dose
Inhalation_Dose<-phys$Inhalation_Dose

#Collection of all parameters so they can be entered in the function
parameters=cbind(Volume_exposure_chamber,
                 P_F,
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
                 BW,
                 V_F,
                 V_L,
                 V_SI,
                 V_A,
                 V_V,
                 V_RP,
                 V_SP,
                 V_Pu,
                 Q_C,
                 Q_F,
                 Q_L,
                 Q_SI,
                 Q_Pu,
                 Q_RP,
                 Q_SP,
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
                 Oral_Dose,
                 Inhalation_Dose,
                 Volume_exposure_chamber,
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
                 Vsmax_SI_GST)

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

#inhalation exposure  exposure
ex <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(dose = (Dose_in_mg_bw * phys$BW)/ MW  * 1e+3, dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(dose = (Inhalation_dose_in_mg_bw * phys$BW)/ MW  * 1e+3, dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(seq(from = time.0, to = time.end, by = time.frame)) 