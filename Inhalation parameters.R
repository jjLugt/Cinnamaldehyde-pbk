#author: Joris Jean van der Lugt
#date: 05-08-2022
#Inhalation parameters 
library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(truncnorm)
library(reshape2)
library(sensitivity)
library(PKNCA)

#Simulations
set.seed(15204)                       #to ensure a reproducible output if random input is used
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
interval                   <-0.1      #interval between dosing in H
time.0                     <-0        #time start dosing
time.end                   <-24        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
Oral_Dose_in_mg_bw         <-0.7    #Dose in mg/kg-bw
Inhalation_Dose_in_mg_bw   <-0      #The inhaled dose in mg/kg-bw
iv_dose_in_mg_bw           <-0       #IV administered dose in mg/kg/bw
MW                         <-132.16   #The molecular weight of Cinnamaldehyde
BW                         <-70      #Body weight in Kg
Oral_Dose                  <-(Oral_Dose_in_mg_bw * BW)/ MW  * 1e+3       #The administered dose in μmol
Inhalation_Dose            <-(Inhalation_Dose_in_mg_bw * BW)/ MW  * 1e+3 #The inhaled dose in μmol
iv_dose                    <-(iv_dose_in_mg_bw * BW)/ MW  * 1e+3    
Volume_exposure_chamber    <-10       #volume exposure chamber in L


#--Physio-chemical parameters--#
#Values calculated in QSAR calculations file
#-Cinnamaldehyde-#

P_F      <-  47.75  #Fat/Blood partition coefficient
P_L      <-  1.83   #liver/Blood partition coefficient
P_SI     <-  1.81   #Small intestine/Blood partition coefficients
P_RP     <-  1.81   #Richly perfused tissues/Blood partition coefficients
P_SP     <-  1.50   #Slowly perfused tissues/Blood partition coefficients
P_B      <-  0.29  #Blood/Air Partition Coefficient 
P_Pu     <-  1.81   #lung/Blood partition coefficient

#-Cinnamyl Alcohol-#
P_OH_F    <-  49.26 #Fat/Blood partition coefficient
P_OH_L    <-  1.18  #liver/Blood partition coefficient
P_OH_SI   <-  1.18  #Small intestine/Blood partition coefficients
P_OH_RP   <-  1.18  #Richly perfused tissues/Blood partition coefficients
P_OH_SP   <-  1.53  #Slowly perfused tissues/Blood partition coefficients
P_OH_Pu   <-  1.18  #Lung/Blood partition coefficients

#--Physiological Parameters--#

fV_fat   <- 0.2142  #Unit less fraction of total tissue volume that consist of fat
fV_liver <- 0.0257  #Unit less fraction of total tissue volume that consist of  liver
fV_si    <- 0.0091  #Unit less fraction of total tissue volume that consist of Small intestine
fV_a     <- 0.02  #Unit less fraction of total tissue volume that consist of arterial blood
fV_v     <- 0.059  #Unit less fraction of total tissue volume that consist of venous blood
fV_rp    <- 0.048  #Unit less fraction of total tissue volume that consist of Richly perfused tissue
fV_sp    <- 0.6164  #Unit less fraction of total tissue volume that consist of Slowly perfused tissue
fV_pu    <- 0.076  #Unit less fraction of total tissue volume that consist of

V_adjust <- fV_fat + fV_liver + fV_si + fV_a + fV_v + fV_rp + fV_sp +fV_pu

#-Tissues volumes in L #
V_F      <- (fV_fat / V_adjust)* BW    #Fat
V_L      <- (fV_liver / V_adjust)* BW  #Liver
V_SI     <- (fV_si / V_adjust)* BW     #Small intestine
V_A      <- (fV_a / V_adjust)* BW  #Arterial Blood
V_V      <- (fV_v / V_adjust)* BW  #Venous Blood
V_RP     <- (fV_rp / V_adjust)* BW  #Richly perfused (RP)
V_SP     <- (fV_sp / V_adjust)* BW   #Slowly perfused (SP)
V_Pu     <- (fV_pu / V_adjust)* BW    #Lung

#-Cardiac parameters-#

Q_C       <-   390   #Cardiac output in L/h

fQC_fat   <- 0.052  #Unit less fraction of total Cardiac output that goes to fat
fQC_liver <- 0.141  #Unit less fraction of total Cardiac output that goes to the liver
fQC_si    <- 0.086  #Unit less fraction of total Cardiac output that goes to the Small intestine
fQC_rp    <- 0.473  #Unit less fraction of total Cardiac output that goes to the Richly perfused tissue
fQC_sp    <- 0.248  #Unit less fraction of total Cardiac output that goes to the Slowly perfused tissue

Q_adjust <- fQC_fat + fQC_liver + fQC_si + fQC_rp + fQC_sp

#-Blood flow to tissues in l#
Q_F      <- (fQC_fat / Q_adjust)* Q_C    #Fat
Q_L      <- (fQC_liver / Q_adjust)* Q_C  #Liver
Q_SI     <- (fQC_si / Q_adjust)* Q_C     #Small intestine
Q_RP     <- (fQC_rp / Q_adjust)* Q_C  #Richly perfused (RP)
Q_SP     <- (fQC_sp / Q_adjust)* Q_C   #Slowly perfused (SP)
Q_Pu     <- Q_C    #Lung

#inhalation parameters
P_V     <-   540             #Alveolar ventilation L/h

#----GSH parameters----#
#--GSH synthesis in μmol/kg tissue/h--#

G_SYN_L     <- 860 * V_L * 0.9       #Liver 
G_SYN_SI    <- 78   * V_L * 0.9        #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
k_L_GLOS    <- 0.142        #Liver
k_SI_GLOS   <- 0.044        #Small intestine

#--Initial GSH concentration--#
init_GSH_L  <- 6120 * V_L   #initial amount of GSH in the liver in μmol
init_GSH_SI <- 1780 * V_SI  #initial amount of GSH in the small intestine in μmol

k_GSH <- 6.6 * 10^(-4)      #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
C_PRO_L     <- 5319   * V_L      #Liver
C_PRO_SI    <- 245    * V_SI     #Small intestine

#--Chemical parameters--#
Ka <- 0.62                   #Absorption rate constant for uptake in the Small intestine in per H

#----Liver----#
S9_scaling_L <- 143 * (V_L * 1000) #scaling factor for S9 fraction per g tissue

#----Liver----#

#-first rate order constants-#
k_L_OH  <- 4.2*10^-2 * 60 / 1000 * S9_scaling_L        #Scaled first rate order constant for the enzymatic oxidation of cinnamyl alcohol in the liver in umol/h

#--Michaelis menten constants--#
Km_L_CA     <-  8.5         #Km for enzymatic oxidation of Cinnamaldehyde into Cinnamic acid in the liver in μM
Km_L_AO     <-  330         #Km for enzymatic reduction of Cinnamaldehyde into cinnamyl alcOHol in the liver in μM
Km_L_GST    <-  1.7*10^3    #Km for enzymatic conjugation of Cinnamaldehyde with GST in the liver in μM  
Km_L_GST_G  <-  100         #Km toward GSH for enzymatic conjugation of Cinnamaldehyde in the liver (μM)

#--Vmax values--#
Vsmax_L_CA    <-  9.7 * 60 / 1000 * S9_scaling_L   #Scaled Vmax for enzymatic oxidation of Cinnamaldehyde in the liver in μmol/h 
Vsmax_L_AO    <-  73  * 60 / 1000 * S9_scaling_L      #Scaled Vmax for enzymatic reduction of Cinnamaldehyde in the liver in μmol/h
Vsmax_L_GST   <-  37  * 60 / 1000 * S9_scaling_L       #Scaled Vmax for enzymatic conjugation of Cinnamaldehyde with GSH in the liver in μmol/h

#----Small intestines----#
S9_scaling_SI <-  11.4 * (V_SI *1000) #scaling factor fraction S9 protein per g tissue

#--Michaelis Menten constants--#
Km_SI_CA    <- 70          #Km for enzymatic oxidation of Cinnamaldehyde into Cinnamic acid in the Small Intestine in μM
Km_SI_AO    <- 90          #Km for enzymatic reduction of Cinnamaldehyde into Cinnamyl alcOHol in the Small Intestine in μM
Km_SI_OH    <- 290         #Km for enzymatic oxidation of Cinnamly alcOHol into Cinnamaldehyde in the Small Intestine in μM

#-Vmax values-#
Vsmax_SI_CA    <- 21   * 60/1000 * S9_scaling_SI   #Scaled Vmax for enzymatic oxidation of Cinnamaldehyde into Cinnamic acid in the Small Intestine in μmol/h 
Vsmax_SI_AO    <- 30   * 60/1000 * S9_scaling_SI   #Scaled Vmax for enzymatic reduction of Cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
Vsmax_SI_OH    <- 5.0  * 60/1000 * S9_scaling_SI   #Scaled Vmax for enzymatic Oxidation of Cinnamyl alcohol into Cinnamaldehyde in the Small Intestine in μmol/h 

#Collection of all parameters so they can be entered in the function
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
                 C_PRO_L,
                 C_PRO_SI,
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
                 Vsmax_SI_CA,
                 Vsmax_SI_AO,
                 Vsmax_SI_OH,
                 Volume_exposure_chamber)


#defining the begin situation of the model Inhalation variation 
inits <- c("A_GI"         =0,
           "A_Exhalation" =0,
           "A_inhalation_Dose" =0,
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
           "A_OH_M_L_C_A" =0,
           "A_OH_L"       =0,
           "A_L"          =0,
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);




#inhalation exposure  exposure
ex <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(dose = Oral_Dose, dur=0.01, cmt="A_GI", nbr.doses=1)%>%
  et(dose = Inhalation_Dose, dur=0.001, cmt="A_inhalation_Dose", nbr.doses=nbr.doses)%>%
  et(dose= init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(dose= init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)%>%
  et(dose=iv_dose, dur=0.005,cmt="A_V",nbr.doses=1)%>%
  et(seq(from = time.0, to = time.end, by = time.frame))

