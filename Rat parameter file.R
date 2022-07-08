#author: Joris Jean van der Lugt
#date: 17-06-2022
#Rat cinnamaldehyde pbk Model adapted from:  "Dose-dependent DNA adduct formation by cinnamaldehyde and other food-borne α,β-unsaturated aldehydes predicted by physiologically based in silico modelling"

library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(shiny)
library(truncnorm)
library(reshape2)
library(plotly)

#Simulations
set.seed(15204)         #to ensure a reproducible output
amount.units <-"umol"
time.units   <-"h"
nbr.doses    <-1        #number of doses
time.0       <-0        #time start dosing
time.end     <-8        #time end of simulation
time.frame   <-0.1     #time steps of simulation
#N           <-1000     #Number of males
#NF          <-1000     #Number of females
Dose_in_mg   <-50     #Dose in mg/kg-bw
MW           <-132.16   #The molecular weight of Cinnamaldehyde
DOSE         <-(Dose_in_mg * 0.25)/ MW  * 1e+3     #The administered dose in umol 


#--Physico-chemical parameters--#
#-Cinnamaldehyde-#

P_F      <-  14.2 #Fat/Blood partition coefficient
P_L      <-  1.21 #Fat/Blood partition coefficient
P_SI     <-  1.21 #Small intestine/Blood partition coefficients
P_RP     <-  1.21 #Richly perfused tissues/Blood partition coefficients
P_SP     <-  0.57 #Slowely perfused tissues/Blood partition coefficients

#-Cinnamyl AlcOHol-#
P_OH_F    <-  14.6 #Fat/Blood partition coefficient
P_OH_L    <-  1.22 #Fat/Blood partition coefficient
P_OH_SI   <-  1.22 #Small intestine/Blood partition coefficients
P_OH_RP   <-  1.22 #Richly perfused tissues/Blood partition coefficients
P_OH_SP   <-  0.57 #Slowly perfused tissues/Blood partition coefficients

#--Pyshiological Parameters--#
BW      <- 0.25    #Body weight in Kg

#-Tissues volumes in % body weight-#

V_F      <- 7  #Fat
V_L      <- 3.4   #Liver
V_SI     <- 1.4   #Small intestine
V_A      <- 1.9   #Arterial Blood
V_V      <- 5.6   #Venous Blood
V_RP     <- 4.2   #Richly perfused 
V_SP     <- 67.6  #Slowly perfused 

#-Cardiac parameters-#

Q_C      <- 5.4    #Cardiac output in L/h

#-Blood flow to tissues in % cardiac output-#

Q_F      <- 0.07    #Fat
Q_L      <- 0.13   #Liver
Q_SI     <- 0.12    #Small intestine
Q_RP     <- 0.64   #Richly perfused (RP)
Q_SP     <- 0.17   #Slowly perfused (SP)

#----GSH parameters----#
#--GSH synthesis in umol/kg tissue/h--#

G_SYN_L     <- 6.6  #Liver 
G_SYN_SI    <- 0.25    #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
k_L_GLOS    <- 0.142 #Liver
k_SI_GLOS   <- 0.044 #Small intestine

#--Initial GSH concentration--#
init_GSH_L  <- 6120   #initial GSH concentration in the liver in umol/kg
init_GSH_SI <- 1780   #initial GSH concentration in the small intestine in umol/kg

k_GSH <- 6.6 * 10^(-4) #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h
k_DNA <- 1.6 * 10^(-8) #The second-order rate constant of the reaction between cinnamaldehyde and 2ʹ-dG in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
C_PRO_L     <- 3000  #Liver
C_PRO_SI    <- 774   #Small intestine

#----DNA parameters----#
C_L_dG     <-  1.36 #Concentration of 2ʹ-dG in the liver μmol/kg liver
T_0.5      <-  38.5   #Half-life of DNA adduct in the liver in hours

#--Chemical parameters--#
Ka <- 5.0  #Absorption rate constant for uptake in the Small intestine in per H

#----Liver----#

#First order rate constants
k_L_CA   <-  7.4*10^(-4)  #Scaled first-order rate constant for enzymatic oxidation of cinnamaldehyde in the liver (μmol/h)
k_L_GST  <-  6.2*10^(-2)  #Scaled first-order rate constant for enzymatic conjugation of cinnamaldehyde with GSH in the liver (μmol/h)

#--Michaelis menten constants--#

Km_L_CA     <-  8.5  #Km for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the liver in μM
Km_L_AO     <-  120  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the liver in μM
Km_L_GST    <-  100  #Km for enzymatic conjugation of cinnamaldehyde with GST in the liver in μM  
Km_L_GST_G  <-  100  #Km for enzymatic conjugation of cinnamaldehyde with GST in the liver in μM 
Km_L_OH     <-  1300 #Km for enzymatic oxidationf of cinnamyl alcohol to cinnamaldehyde in the liver in umol

#--Vmax values--#

Vsmax_L_CA    <- 9.7  #Scaled Vmax for enzymatic oxidation of cinnamaldehyde in the liver in μmol/h 
Vsmax_L_AO    <- 29   #Scaled Vmax for enzymatic reduction of cinnamaldehyde in the liver in μmol/h
Vsmax_L_GST   <- 37   #Scaled Vmax for enzymatic conjugation of cinnamaldehyde with GSH in the liver in μmol/h
Vsmax_L_GST_G <- 100  #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the Liver (μM RAT value)
Vsmax_L_OH    <- 15   #Scaled Vmax for enzymatic oxidation of cinnamyl alcohol to cinnamaldehyde in the liver in umol/h

#----Small intestines----#
#First order rate constants
k_SI_CA  <- 3.9*10^(-3)            #Scaled first-order rate constant for enzymatic oxidation of cinnamaldehyde in the small intestine (μmol/h)

#--Michaelis menten constants SI--#
Km_SI_CA    <- 70  #Km for enzymatic oxidation of cinnamaldehyde into cinnamic acid in the Small Intestine in μM
Km_SI_AO    <- 75  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the Small Intestine in μM
Km_SI_OH    <- 290 #Km for enzymatic oxidation of cinnamly alcOHol into cinnamaldehyde in the Small Intestine in μM
Km_SI_GST   <- 600 #Km for enzymatic conjugation of cinnamaldehye with GST in the Small Intestine in μM (RAT value)
Km_SI_GST_G <- 100  #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the small intestine (μM

#-Vmax values-#
Vsmax_SI_CA    <- 21 #Scaled Vmax for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the Small Intestine in μmol/h 
Vsmax_SI_AO    <- 5.8 #Scaled Vmax for enzymatic reduction of cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
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
                 k_L_CA,
                 k_L_GST,
                 Km_L_CA,
                 Km_L_AO,
                 Km_L_GST,
                 Km_L_GST_G,
                 Km_L_OH,
                 Vsmax_L_CA,
                 Vsmax_L_AO,
                 Vsmax_L_GST,
                 Vsmax_L_OH,
                 k_SI_CA,
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
           "AA_A"         = 0,
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
           "A_SP"         = 0,
           "A_OH_SP"      = 0
);



#Step 3 exposure
ex <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(dose = DOSE, dur=0.1, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(seq(from = time.0, to = time.end, by = time.frame)) 

