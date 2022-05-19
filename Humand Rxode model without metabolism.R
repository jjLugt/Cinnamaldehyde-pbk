#author: Joris Jean van der Lugt
#date: 28-10-2021
#Human cinnamaldehyde pbk Model adapted from:  "Dose-dependent DNA adduct formation by cinnamaldehyde and other food-borne α,β-unsaturated aldehydes predicted by physiologically based in silico modelling"

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
Dose_in_mg   <-100    #Dose in mg/kg-bw
MW           <-132.16   #The molecular weight of Cinnamaldehyde
DOSE         <-(Dose_in_mg * 70)/ MW  * 1e+3     #The administered dose in umol 


Rin <-0
R_F<-0
R_L<-0
R_SI<-0
R_V<-0
R_SP<-0
R_RP<-0

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

V_F      <- 21.4  #Fat
V_L      <- 2.6   #Liver
V_SI     <- 0.9   #Small intestine
V_A      <- 2.0   #Arterial Blood
V_V      <- 5.9   #Venous Blood
V_RP     <- 4.1   #Richly perfused 
V_SP     <- 51.7  #Slowly perfused 

#-Cardiac parameters-#

Q_C      <- 310    #Cardiac output in L/h

#-Blood flow to tissues in % cardiac output-#

Q_F      <- 5.2    #Fat
Q_L      <- 14.1   #Liver
Q_SI     <- 8.6    #Small intestine
Q_RP     <- 47.3   #Richly perfused (RP)
Q_SP     <- 24.8   #Slowly perfused (SP)

#----GSH parameters----#
#--GSH synthesis in umol/kg tissue/h--#

G_SYN_L     <- 1122  #Liver 
G_SYN_SI    <- 27    #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
k_L_GLOS    <- 0.142 #Liver
k_SI_GLOS   <- 0.044 #Small intestine

#--Initial GSH concentration--#
init_GSH_L  <- 5639 * V_L   #initial amount of GSH in the liver in umol
init_GSH_SI <- 1250 * V_SI  #initial amount of GSH in the small intestine in umol

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
parameters=cbind(Rin <-Rin,
                 R_F<-R_F,
                 R_L<-R_L,
                 R_SI<-R_SI,
                 R_V<-R_V,
                 R_SP<-R_SP,
                 R_RP<-R_SP,
                 P_F=P_F,
                 P_L=P_L,
                 P_SI=P_SI,
                 P_RP=P_RP,
                 P_SP=P_SP,
                 P_OH_F=P_OH_F,
                 P_OH_L=P_OH_L,
                 P_OH_SI=P_OH_SI,
                 P_OH_RP=P_OH_RP,
                 P_OH_SP=P_OH_SP,
                 BW=BW,
                 V_F=V_F,
                 V_L=V_L,
                 V_SI=V_SI,
                 V_A=V_A,
                 V_V=V_V,
                 V_RP=V_RP,
                 V_SP=V_SP,
                 Q_C=Q_C,
                 Q_F=Q_F,
                 Q_L=Q_L,
                 Q_SI=Q_SI,
                 Q_RP=Q_RP,
                 Q_SP=Q_SP,
                 DOSE=DOSE,
                 Ka=Ka
                 )


#defining the begin situation of the model (in this case no chemical present in the organs)
inits <- c("A_GI"         = 0,
           "A_V"          = 0,
           "A_A"          =0,
           "A_F"          = 0,
           "A_L"          = 0,
           "A_SI"         = 0,
           "A_RP"         =0,
           "A_SP"         =0
);



#Step 3 exposure
ex <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(dose = DOSE, dur=0.001, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(seq(from = time.0, to = time.end, by = time.frame)) 


PBK_Cinnamaldehyde <- RxODE({
  
  #--Defining the compartments of the model--#
  
  #rate of change in cinnamaldehyde concentration in the GI cavity in umol-#
  Rin            <- -Ka * A_GI;
  
  #-Blood concentrations-#
  C_V            <- A_V       / V_V ;           #Concentration of Cinnamaldehyde in Venous blood in umol/l
  C_A            <- A_A     / V_A                          #Concentration of cinnamaldehyde in Arterial blood   
  
  
  #-Concentration in fat-#
  C_F            <- A_F       / V_F;                    #Concentration in Fat in umol/kg
  C_V_F          <- C_F       / P_F;                    #Concentration of cinnamaldehyde in venous blood leaving Fat in umol/l
  R_F            <- Q_F * (C_A - C_V_F); #rate of change in cinnamaldehyde concentration in the Fat in umol
  
  
  #-Concentration in richly perfused tissue-#
  C_RP           <- A_RP      / V_RP;                   #Concentration of Cinnamaldehyde in RP tissue in umol/kg
  C_V_RP         <- C_RP      / P_RP;                   #concentration of Cinnamaldehyde in RP venous blood leaving the tissue in umol/L
  R_RP           <- Q_RP * (C_A - C_V_RP); #rate of change in cinnamaldehyde concentration in the RP in umol-#
  
  #-Concentration in slowly perfused tissue-#
  C_SP           <- A_SP     / V_SP;                   #Concentration in SP tissue in umol/kg
  C_V_SP         <- C_SP      / P_SP;                   #Concentration of Cinnamaldehyde in venous blood leaving the SP in umol/l
  R_SP           <- Q_SP * (C_A - C_V_SP); #rate of change in cinnamaldehyde concentration in the SP tissue in umol#
  
  #-Concentration in the Small intestine-#
  C_SI           <- A_SI      / V_SI;                   #Concentration Cinnamaldehyde in the Small intestine in umol/kg
  C_V_SI         <- C_SI      / P_SI;                   #Concentration of cinnamaldehyde in venous blood leaving the Small intestine in umol/l
  R_SI           <- Q_SI * (C_A - C_V_SI) +Ka *A_GI #rate of change in cinnamaldehyde concentration in the SI in umol-#
  
  #-Concentration in Liver-#
  C_L            <- A_L       / V_L;                    #Concentration Cinnamaldehyde in the Liver in umol/kg
  C_V_L          <- C_L       / P_L;                    #Concentration of cinnamaldehyde in venous blood leaving the Liver in umol/l
  R_L            <- Q_L * C_A + Q_SI * C_V_SI - (Q_L + Q_SI) * C_V_L; #rate of change in cinnamaldehyde concentration in the liver in umol-#
  
  R_V            <- Q_F * C_V_F + (Q_L + Q_SI) * C_V_L + Q_RP * C_V_RP + Q_SP * C_V_SP - Q_C * C_V; #rate of change in cinnamaldehyde concentration in the venous blood in umol
  R_A            <- Q_C * C_V - (Q_F * C_A + Q_L * C_A + Q_SI * C_A + Q_RP * C_A + Q_SP * C_A);
  #--Differential equations--#
  
  #-amount of Cinnamaldehyde in GI cavity in umol
  d/dt(A_GI)     <- Rin;
  
  #-Amount of Cinnamaldehyde in Venous blood in umol 
  
  d/dt(A_V)      <- R_V;
  d/dt(A_A)      <-R_A;
  
  #-Fat-#
  #Amount of Cinnamaldehyde in the Fat in umol 
  d/dt(A_F)      <- R_F;                         

  
  #-Amount of cinnamaldehyde in the liver in umol-#
  d/dt(A_L)      <- R_L; 
  
  #-Amount of Cinnamaldehyde in the Small intestine-#
  d/dt(A_SI)     <- R_SI;
  
  #Richly perfused Tissue#
  #Amount of Cinnamaldehyde in the RP tissue in umol 
  d/dt(A_RP)     <- R_RP;           
  
  
  #-Slowly perfused Tissue-#
  d/dt(A_SP)     <- R_SP;         #Amount of Cinnamaldehyde in the SP tissue in umol 
  
  
})

solve.pbk_nonpop <- solve(PBK_Cinnamaldehyde, parameters, events = ex, inits, cores=4) #Solve the PBPK model


solve.pbk_nonpop <- solve.pbk_nonpop/70 * MW /1e+3
mass_df <-solve.pbk_nonpop[,c(14:20)]
mass_at_t <- data.frame(mass=as.numeric())


for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(solve.pbk_nonpop[i,])
}



plot(mass_at_t[,1])