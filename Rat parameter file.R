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
set.seed(15204)                       #to ensure a reproducible output
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-4        #time end of simulation
time.frame                 <-0.1      #time steps of simulation
Oral_dose_in_mg_bw         <-0     #Dose in mg/kg-bw
Inhalation_dose_in_mg_bw   <-0        #The inhaled dose in mg/kg-bw
iv_dose_in_mg_bw           <-20       #IV administered dose in mg/kg/bw
MW                         <-132.16   #The molecular weight of Cinnamaldehyde
BW                         <- 0.25    #Body weight in Kg
Oral_Dose                  <-(Oral_dose_in_mg_bw * BW)/ MW  * 1e+3       #The administered dose in μmol
Inhalation_Dose            <-(Inhalation_dose_in_mg_bw * BW)/ MW  * 1e+3 #The inhaled dose in μmol
iv_dose                    <-(iv_dose_in_mg_bw * BW)/ MW  * 1e+3         #the administered IV dose 
Volume_exposure_chamber    <-10       #volume exposure chamber in L


#--Physio-chemical parameters--#
#-Cinnamaldehyde-#

P_F      <-  17.42#Fat/Blood partition coefficient
P_L      <-  1.18 #Fat/Blood partition coefficient 
P_SI     <-  1.18 #Small intestine/Blood partition coefficients
P_RP     <-  1.18 #Richly perfused tissues/Blood partition coefficients
P_SP     <-  0.39 #Slowly perfused tissues/Blood partition coefficients
P_B      <-  0.29 #Blood air partition coefficients
P_Pu     <-  0.49 #Lung/Blood partition coefficients

#-Cinnamyl Alcohol-#
P_OH_F    <- 17.65 #Fat/Blood partition coefficient
P_OH_L    <-  1.18 #Fat/Blood partition coefficient
P_OH_SI   <-  1.18 #Small intestine/Blood partition coefficients
P_OH_RP   <-  1.18 #Richly perfused tissues/Blood partition coefficients
P_OH_SP   <-  0.49 #Slowly perfused tissues/Blood partition coefficients
P_OH_Pu   <-  1.18 #Lung/Blood partition coefficients


#--Physiological Parameters--#

#-Tissues volumes in L-#

V_F      <- 0.07  * BW  #Fat
V_L      <- 0.034 * BW  #Liver
V_SI     <- 0.014 * BW  #Small intestine
V_A      <- 0.019 * BW  #Arterial Blood
V_V      <- 0.059 * BW  #Venous Blood
V_RP     <- 0.037 * BW  #Richly perfused 
V_SP     <- 0.676 * BW  #Slowly perfused
V_Pu     <- 0.005 * BW  #Lung

#-Cardiac parameters-#

Q_C      <- 5.4    #Cardiac output in L/h

#-Blood flow to tissues in % cardiac output-#

Q_F      <- 0.07 * Q_C   #Fat
Q_L      <- 0.13 * Q_C #Liver
Q_SI     <- 0.12 * Q_C  #Small intestine
Q_RP     <- 0.64 * Q_C  #Richly perfused
Q_SP     <- 0.17 * Q_C  #Slowly perfused
Q_Pu     <- Q_C    #Blood flow through the lung

P_V      <- 0.75 #Pulmonary ventilation in L/h 

#----GSH parameters----#
#--GSH synthesis in umol/per organ/h--#

G_SYN_L     <- 869 * V_L * 0.9  #Liver 
G_SYN_SI    <- 78  * V_SI * 0.9  #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
k_L_GLOS    <- 0.142 #Liver
k_SI_GLOS   <- 0.044 #Small intestine

#--Initial GSH concentration--#
init_GSH_L  <- 6120 * V_L  #initial GSH concentration in the liver in umol
init_GSH_SI <- 1780 * V_SI #initial GSH concentration in the small intestine in umol

k_GSH <- 6.6 * 10^(-4) #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
C_PRO_L     <- 5319  * V_L  #Liver
C_PRO_SI    <- 245   * V_SI #Small intestine

#--Chemical parameters--#
Ka <- 1.37  #Absorption rate constant for uptake in the Small intestine in per H

#----Liver----#
S9_scaling_L <- 143 * (V_L * 1000) #scaling factor for S9 fraction per g tissue


#First order rate constants
k_L_CA   <-  7.4*10^(-3)* 60 / 1000 * S9_scaling_L   #Scaled first-order rate constant for enzymatic oxidation of cinnamaldehyde in the liver (L/h)
k_L_GST  <-  6.2*10^(-2)* 60 / 1000 * S9_scaling_L  #Scaled first-order rate constant for enzymatic conjugation of cinnamaldehyde with GSH in the liver (L/h)

#--Michaelis menten constants--#

Km_L_AO     <-  120  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the liver in μM
Km_L_GST_G  <-  100  #Km for enzymatic conjugation of cinnamaldehyde with GST in the liver in μM 
Km_L_OH     <-  1300 #Km for enzymatic oxidation of cinnamyl alcohol to cinnamaldehyde in the liver in μM

#--Vmax values--#
Vsmax_L_AO    <- 29  * 60 / 1000 * S9_scaling_L      #Scaled Vmax for enzymatic reduction of cinnamaldehyde in the liver in μmol/h
Vsmax_L_GST_G <- 100  *60 / 1000 * S9_scaling_L      #Scaled Vmax toward GSH for enzymatic conjugation of cinnamaldehyde in the Liver (μM RAT value)
Vsmax_L_OH    <- 15   *60 / 1000 * S9_scaling_L      #Scaled Vmax for enzymatic oxidation of cinnamyl alcohol to cinnamaldehyde in the liver in μmol/h

#----Small intestines----#
S9_scaling_SI <- 11.4 * V_SI *1000      #scaling factor fraction S9 protein per g tissue

#First order rate constants
k_SI_CA  <- 3.9*10^(-3) * 60/1000 * S9_scaling_SI #Scaled first-order rate constant for enzymatic oxidation of cinnamaldehyde in the small intestine (μmol/h)

#--Michaelis menten constants SI--#
Km_SI_CA    <- 0   #Km for enzymatic oxidation of cinnamaldehyde into cinnamic acid in the Small Intestine in μM
Km_SI_AO    <- 75  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the Small Intestine in μM
Km_SI_GST   <- 600 #Km for enzymatic conjugation of cinnamaldehye with GST in the Small Intestine in μM (RAT value)
Km_SI_GST_G <- 100 #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the small intestine (μM)

#-Vmax values-#
Vsmax_SI_AO    <- 5.8 * 60 / 1000 * S9_scaling_SI  #Scaled Vmax for enzymatic reduction of cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
Vsmax_SI_GST   <- 63  * 60 / 1000 * S9_scaling_SI  #Scaled Vmax for enzymatic Conjugation of cinnamaldehyde with GSH in the in the small intestine in μmol/h (RAT value)

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
                 k_L_GST,
                 Km_L_OH,
                 Km_L_AO,
                 Km_L_GST_G,
                 k_L_CA,
                 Vsmax_L_AO,
                 Vsmax_L_OH,
                 k_SI_CA,
                 Km_SI_CA,
                 Km_SI_AO,
                 Km_SI_GST,
                 Km_SI_GST_G,
                 Vsmax_SI_AO,
                 Vsmax_SI_GST,
                 S9_scaling_L,
                 S9_scaling_SI)


#defining the begin situation of the model Inhalation variation 
inits <- c("A_GI"         =0,
           "A_V"          =0,
           "A_OH_V"       =0,
           "A_A"          =0,
           "A_OH_A"       =0,
           "A_inhalation_dose" =0,
           "A_Exhalation" =0,
           "A_Pu"         =0,
           "A_OH_Pu"      =0,
           "A_F"          =0,
           "A_OH_F"       =0,
           "A_L"          =0,
           "A_OH_L"       =0,
           "AM_L_GST"     =0,
           "AM_L_CHEM"    =0,
           "AM_L_CA"      =0,
           "AM_L_AP"      =0,
           "AM_Lc_GSH"    =0,
           "AM_L_AO"      =0,
           "AM_OH_L_C_A"  =0,
           "A_SI"         =0,
           "A_OH_SI"      =0,
           "AM_SI_GST"    =0,
           "AM_SI_CHEM"   =0,
           "AM_SI_CA"     =0,
           "AM_SI_AP"     =0,
           "AM_SIc_GSH"   =0,
           "AM_SI_AO"     =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);




#inhalation exposure  exposure
ex <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(dose = Oral_Dose, dur=0.001, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(dose=iv_dose, dur=0.005,cmt="A_V",nbr.doses=1)%>%
  et(dose= init_GSH_L, cmt="AM_Lc_GSH", nbr.doses=1)%>%
  et(dose= init_GSH_SI, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(dose = Inhalation_Dose, dur=0.001, cmt="A_inhalation_Dose", nbr.doses=nbr.doses)%>%
  et(seq(from = time.0, to = time.end, by = time.frame)) 
