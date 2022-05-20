#author: Joris Jean van der Lugt
#date: 20-05-2021
#Inhalation parameters 



library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(shiny)
library(truncnorm)
library(reshape2)

#Simulations
set.seed(15204)         #to ensure a reproducible output
amount.units <-"umol"
time.units   <-"h"
nbr.doses    <-1        #number of doses
time.0       <-0        #time start dosing
time.end     <-8        #time end of simulation
time.frame   <-0.01     #time steps of simulation
Dose_in_mg   <-0     #Dose in mg/kg-bw
MW           <-132.16   #The molecular weight of Cinnamaldehyde
inhalation_dose_in_mg <- 100  #The inhaled dose in mg/kg
DOSE         <-(Dose_in_mg * 70)/ MW  * 1e+3     #The administered dose in umol 



RM_L_DA <- 0 
RM_Lc_GSH  <- 0 
RM_SI_AG_GST <- 0
RM_SI_AG_CHEM <- 0
RM_SIc_GSH <- 0 


#--Physico-chemical parameters--#
#-Cinnamaldehyde-#

P_F      <-  39.3 #Fat/Blood partition coefficient
P_L      <-  2.04 #liver/Blood partition coefficient
P_SI     <-  2.04 #Small intestine/Blood partition coefficients
P_RP     <-  2.04 #Richly perfused tissues/Blood partition coefficients
P_SP     <-  1.57 #Slowely perfused tissues/Blood partition coefficients
P_B      <- 274.84 # Blood/Air Partition Coefficient (unitless)
P_PU     <-  2.04  #lung/Blood partition coefficient

#-Cinnamyl AlcOHol-#
P_OH_F    <-  40.5 #Fat/Blood partition coefficient
P_OH_L    <-  2.09 #liver/Blood partition coefficient
P_OH_SI   <-  2.09 #Small intestine/Blood partition coefficients
P_OH_RP   <-  2.09 #Richly perfused tissues/Blood partition coefficients
P_OH_SP   <-  1.60 #Slowly perfused tissues/Blood partition coefficients
P_OH_PU   <-  2.09 #Lung/Blood partition coefficients
#--Pyshiological Parameters--#
BW      <- 70    #Body weight in Kg

#-Tissues volumes in % body weight-#

V_F      <- 21.4  #Fat
V_L      <- 2.6   #Liver
V_SI     <- 0.9   #Small intestine
V_A      <- 2.0   #Arterial Blood
V_V      <- 5.9   #Venous Blood
V_RP     <- 3.3   #Richly perfused 
V_SP     <- 51.7  #Slowly perfused 
V_PU     <- 0.8   #Lung volume 

#-Cardiac parameters-#

Q_C      <- 310    #Cardiac output in L/h

#-Blood flow to tissues in % cardiac output-#

Q_F      <- 5.2    #Fat
Q_L      <- 14.1   #Liver
Q_SI     <- 8.6    #Small intestine
Q_RP     <- 47.3   #Richly perfused (RP)
Q_SP     <- 24.8   #Slowly perfused (SP)
Q_PU     <- 100    #lung blood flow 


#inhalation parameters
Q_P = 300  # aveolar ventilation l/h

#----GSH parameters----#
#--GSH synthesis in umol/kg tissue/h--#

G_SYN_L     <- 1122  #Liver 
G_SYN_SI    <- 27    #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
k_L_GLOS    <- 0.142 #Liver
k_SI_GLOS   <- 0.044 #Small intestine

#--Initial GSH concentration--#
init_GSH_L  <- 5639 * V_L   #initial amount of GSH in the liver in umol/kg
init_GSH_SI <- 1250 * V_SI  #initial amount of GSH in the small intestine in umol/kg

k_GSH <- 6.6 * 10^(-4) #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h
k_DNA <- 1.6 * 10^(-8) #The second-order rate constant of the reaction between cinnamaldehyde and 2ʹ-dG in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
C_PRO_L     <- 5319  #Liver
C_PRO_SI    <- 245   #Small intestine

#----DNA parameters----#
C_L_dG     <-  1.36 #Concentration of 2ʹ-dG in the liver μmol/kg liver
T_0.5      <-  38.5   #Half-life of DNA adduct in the liver in hours


#--Chemical parameters--#
Ka <- 5.0  #Absorption rate constant for uptake in the Small intestine in per H

#----Liver----#

#-first rate order constants-#
k_L_OH  <- 4.2e-02   #Scaled first rate order constant for the enzymatic oxidation of cinnamyl alcohol in the liver in umol/h

#--Michaelis menten constants--#
Km_L_CA     <-  8.5  #Km for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the liver in μM
Km_L_AO     <-  330  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the liver in μM
Km_L_GST    <-  100 #Km for enzymatic conjugation of cinnamaldehyde with GST in the liver in μM  
Km_L_GST_G  <-  1.7*10^3 #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the liver (μM)

#--Vmax values--#
Vsmax_L_CA    <-  9.7  #Scaled Vmax for enzymatic oxidation of cinnamaldehyde in the liver in μmol/h 
Vsmax_L_AO    <-  73   #Scaled Vmax for enzymatic reduction of cinnamaldehyde in the liver in μmol/h
Vsmax_L_GST   <-  37   #Scaled Vmax for enzymatic conjugation of cinnamaldehyde with GSH in the liver in μmol/h

#----Small intestines----#
#--Michaelis menten constants--#
Km_SI_CA    <- 70  #Km for enzymatic oxidation of cinnamaldehyde into cinnamic acid in the Small Intestine in μM
Km_SI_AO    <- 90  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the Small Intestine in μM
Km_SI_OH    <- 290 #Km for enzymatic oxidation of cinnamly alcOHol into cinnamaldehyde in the Small Intestine in μM
Km_SI_GST   <- 600 #Km for enzymatic conjugation of cinnamaldehye with GST in the Small Intestine in μM (RAT value)
Km_SI_GST_G <- 0  #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the small intestine (μM)

#-Vmax values-#
Vsmax_SI_CA    <- 21 #Scaled Vmax for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the Small Intestine in μmol/h 
Vsmax_SI_AO    <- 30 #Scaled Vmax for enzymatic reduction of cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
Vsmax_SI_OH    <- 5.0 #Scaled Vmax for enzymatic Oxidation of cinnamyl alcohol into cinnamaldehyde in the Small Intestine in μmol/h 
Vsmax_SI_GST   <- 63 #Scaled Vmax for enzymatic Conjugation of cinnamaldehyde with GSH in the in the small intestine in μmol/h (RAT value)

#Collection of all parameters so they can be entered in the function
parameters=cbind(RM_L_DA=RM_L_DA,  
                 RM_Lc_GSH=RM_Lc_GSH, 
                 RM_SI_AG_GST=RM_SI_AG_GST,
                 RM_SI_AG_CHEM=RM_SI_AG_CHEM,
                 RM_SIc_GSH=RM_SIc_GSH,
                 P_F=P_F,
                 P_L=P_L,
                 P_SI=P_SI,
                 P_RP=P_RP,
                 P_SP=P_SP,
                 P_B=P_B,
                 P_PU=P_PU,
                 P_OH_F=P_OH_F,
                 P_OH_L=P_OH_L,
                 P_OH_SI=P_OH_SI,
                 P_OH_RP=P_OH_RP,
                 P_OH_SP=P_OH_SP,
                 P_OH_PU=P_OH_PU,
                 BW=BW,
                 V_F=V_F,
                 V_L=V_L,
                 V_SI=V_SI,
                 V_A=V_A,
                 V_V=V_V,
                 V_RP=V_RP,
                 V_SP=V_SP,
                 V_PU=V_PU,
                 Q_C=Q_C,
                 Q_F=Q_F,
                 Q_L=Q_L,
                 Q_SI=Q_SI,
                 Q_RP=Q_RP,
                 Q_SP=Q_SP,
                 Q_PU=Q_PU,
                 Q_P=Q_P,
                 G_SYN_L=G_SYN_L,
                 G_SYN_SI=G_SYN_SI,
                 k_L_GLOS=k_L_GLOS,
                 k_SI_GLOS=k_SI_GLOS,
                 init_GSH_L=init_GSH_L,
                 init_GSH_SI=init_GSH_SI,
                 k_GSH=k_GSH,
                 k_DNA=k_DNA,
                 C_PRO_L=C_PRO_L,
                 C_PRO_SI=C_PRO_SI,
                 C_L_dG=C_L_dG,
                 T_0.5=T_0.5,
                 DOSE=DOSE,
                 inhalation_dose=inhalation_dose,
                 Ka=Ka,
                 k_L_OH = k_L_OH,
                 Km_L_CA=Km_L_CA,
                 Km_L_AO=Km_L_AO,
                 Km_L_GST=Km_L_GST,
                 Km_L_GST_G=Km_L_GST_G,
                 Vsmax_L_CA=Vsmax_L_CA,
                 Vsmax_L_AO=Vsmax_L_AO,
                 Vsmax_L_GST=Vsmax_L_GST,
                 Km_SI_CA=Km_SI_CA,
                 Km_SI_AO=Km_SI_AO,
                 Km_SI_OH=Km_SI_OH,
                 Km_SI_GST=Km_SI_GST,
                 Km_SI_GST_G=Km_SI_GST_G,
                 Vsmax_SI_CA=Vsmax_SI_CA,
                 Vsmax_SI_AO=Vsmax_SI_AO,
                 Vsmax_SI_OH=Vsmax_SI_OH,
                 Vsmax_SI_GST=Vsmax_SI_GST)


#defining the begin situation of the model Inhalation varation 
inits <- c("A_GI"         = 0 ,
           "A_I"          = 0,
           "A_X"          = 0,
           "A_OH_PU"     =0,
           "A_V"          = 0,
           "A_OH_V"       = 0,
           "A_F"          = 0,
           "A_OH_F"       = 0,
           "AM_L_CA"      = 0,
           "AM_L_AO"      = 0,
           "AM_L_AG_GST"  = 0,
           "AM_L_AG_CHEM" = 0,
           "AM_L_AP"      = 0,
           "AM_L_DA_FORM" = 0,
           "AM_L_DA"      = 0,
           "A_OH_M_L_C_A" = 0,
           "A_OH_L"       = 0,
           "A_L"          = 0,
           "AM_Lc_GSH"    =init_GSH_L, 
           "AM_SI_CA"     = 0,
           "AM_SI_AO"     = 0,
           "AM_SI_AG_GST" = 0,
           "AM_SI_AG_CHEM"= 0,
           "AM_SI_AP"     = 0,
           "A_OH_M_SI_C_A"= 0,
           "A_OH_SI"      = 0,
           "A_SI"         = 0,
           "AM_SIc_GSH"   =init_GSH_SI,
           "A_RP"         = 0,
           "A_OH_RP"      = 0,
           "A_SP"          = 0,
           "A_OH_SP"      = 0
);




#inhalation exposure  exposure
inhalation_dose <- (inhalation_dose_in_mg * 70)/ MW * 1e+3  #The inhaled dose in umol
ex <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(dose = DOSE, dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(dose = inhalation_dose, dur = 0.01, cmt="A_I", nbr.doses=nbr.doses)%>%
  et(seq(from = time.0, to = time.end, by = time.frame)) 