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

#--Chemical parameters--#
Ka <- 5.0  #Absorption rate constant for uptake in the Small intestine in per H



#Collection of all parameters so they can be entered in the function
parameters=cbind(P_F,
                 P_L,
                 P_SI,
                 P_RP,
                 P_SP,
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
                 DOSE,
                 Ka
                 )


#defining the begin situation of the model (in this case no chemical present in the organs)
inits <- c("A_GI"         =0,
           "A_V"          =0,
           "A_A"          =0,
           "A_F"          =0,
           "A_L"          =0,
           "A_SI"         =0,
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
  C_V            <- A_V     / V_V ;           #Concentration of Cinnamaldehyde in Venous blood in umol/l
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
  R_SI           <- Q_SI * (C_A - C_V_SI) -Rin #rate of change in cinnamaldehyde concentration in the SI in umol-#
  
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
  
  d/dt(A_V)      <-R_V;
  
  #-Amount of Cinnamaldehyde in Artblood in umol 
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


solve.pbk_nonpop <- solve.pbk_nonpop/BW * MW /1e+3
mass_df <-solve.pbk_nonpop[,c(22:29)]
mass_at_t <- data.frame(mass=as.numeric())


for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])