#author: Joris Jean van der Lugt
#date: 20-05-2021
#Single person rxode model parameters 
#Simulations

library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(shiny)
library(truncnorm)
library(reshape2)

set.seed(15204)         #to ensure a reproducible output
amount.units <-"μmol"
time.units   <-"h"
nbr.doses    <-1        #number of doses
time.0       <-0        #time start dosing
time.end     <-8        #time end of simulation
time.frame   <-0.01     #time steps of simulation
Dose_in_mg   <-100   #Dose in mg/kg-bw
MW           <-132.16   #The molecular weight of Cinnamaldehyde
DOSE         <-(Dose_in_mg * 70)/ MW  * 1e+3     #The administered dose in μmol 


#--Physico-chemical parameters--#
#-Cinnamaldehyde-#

P_F      <-  39.3 #Fat/Blood partition coefficient
P_L      <-  2.04 #Fat/Blood partition coefficient
P_SI     <-  2.04 #Small intestine/Blood partition coefficients
P_RP     <-  2.04 #Richly perfused tissues/Blood partition coefficients
P_SP     <-  1.57 #Slowely perfused tissues/Blood partition coefficients

#-Cinnamyl AlcOHol-#
P_OH_F    <-  40.5 #Fat/Blood partition coefficient
P_OH_L    <-  2.09 #Fat/Blood partition coefficient
P_OH_SI   <-  2.09 #Small intestine/Blood partition coefficients
P_OH_RP   <-  2.09 #Richly perfused tissues/Blood partition coefficients
P_OH_SP   <-  1.60 #Slowly perfused tissues/Blood partition coefficients

#--Pyshiological Parameters--#
BW      <- 70    #Body weight in Kg

#-Tissues volumes in % body weight-#

V_F      <- 0.214 * BW  #Fat
V_L      <- 0.026 * BW   #Liver
V_SI     <- 0.009 * BW  #Small intestine
V_A      <- 0.02  * BW #Arterial Blood
V_V      <- 0.059 * BW  #Venous Blood
V_RP     <- 0.052 * BW  #Richly perfused 
V_SP     <- 0.6 * BW #Slowly perfused 

#-Cardiac parameters-#

Q_C      <- 310    #Cardiac output in L/h

#-Blood flow to tissues in % cardiac output-#

Q_F      <- 0.052 * Q_C    #Fat
Q_L      <- 0.141 * Q_C   #Liver
Q_SI     <- 0.087 * Q_C    #Small intestine
Q_RP     <- 0.472 * Q_C   #Richly perfused (RP)
Q_SP     <- 0.248 * Q_C   #Slowly perfused (SP)

#----GSH parameters----#
#--GSH synthesis in μmol/kg tissue/h--#

G_SYN_L     <- 1122  #Liver 
G_SYN_SI    <- 270    #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
k_L_GLOS    <- 0.142 #Liver
k_SI_GLOS   <- 0.044 #Small intestine

#--Initial GSH concentration--#
init_GSH_L  <- 5639 * V_L   #initial amount of GSH in the liver in μmol/kg
init_GSH_SI <- 1250 * V_SI  #initial amount of GSH in the small intestine in μmol/kg

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
k_L_OH  <- 4.2e-02   #Scaled first rate order constant for the enzymatic oxidation of cinnamyl alcohol in the liver in μmol/h

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
parameters=cbind(P_F,
                 P_L,
                 P_SI,
                 P_RP,
                 P_SP,
                 P_OH_F,
                 P_OH_L,
                 P_OH_SI,
                 P_OH_RP,
                 P_OH_SP,
                 BW,
                 V_F,
                 V_L,
                 V_SI,
                 V_A,
                 V_V,
                 V_RP,
                 V_SP,
                 Q_C,
                 Q_F,
                 Q_L,
                 Q_SI,
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
                 DOSE,
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


#defining the begin situation of the model (in this case no chemical present in the organs)
inits <- c("A_GI"         = 0,
           "A_V"          = 0,
           "A_OH_V"       = 0,
           "A_A"          = 0,
           "A_OH_A"       = 0,
           "A_F"          = 0,
           "A_OH_F"       = 0,
           "AM_L_CA"      = 0,
           "AM_L_AP"      = 0,
           "AM_L_AO"      = 0,
           "A_OH_M_L_C_A" = 0,
           "A_OH_L"       = 0,
           "A_L"          = 0,
           "AM_Lc_GSH"    = init_GSH_L,
           "AM_L_AG_GST"  = 0,
           "AM_L_AG_CHEM" = 0,
           "AM_L_DA_FORM" = 0,
           "AM_L_DA"      = 0,
           "A_OH_SI"      = 0,
           "AM_SI_AO"     = 0,
           "A_OH_M_SI_C_A"= 0,
           "A_SI"         = 0,
           "AM_SI_CA"     = 0,
           "AM_SI_AP"     = 0,
           "AM_SIc_GSH"   = init_GSH_SI,
           "AM_SI_AG_GST" = 0,
           "AM_SI_AG_CHEM"= 0,
           "A_RP"         = 0,
           "A_OH_RP"      = 0,
           "A_SP"         = 0,
           "A_OH_SP"      = 0
);



#Step 3 exposure
ex <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(dose = DOSE, dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(seq(from = time.0, to = time.end, by = time.frame)) 
