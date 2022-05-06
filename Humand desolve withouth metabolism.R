#author: Joris Jean van der Lugt
#date: 28-10-2021
#Human cinnamaldehyde pbk Model adapted from:  "Dose-dependent DNA adduct formation by cinnamaldehyde and other food-borne α,β-unsaturated aldehydes predicted by physiologically based in silico modelling"

library(deSolve)
library(tidyverse)
library(readxl)
library(readr)
library(shiny)

#---- parameter definitions ----#
RM_L_DA <- 0 
RM_Lc_GSH  <- 0 
RM_SI_AG_GST <- 0
RM_SI_AG_CHEM <- 0
RM_SIc_GSH <- 0 
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
init_GSH_L  <- 5639 * V_L  #initial GSH concentration in the liver in umol/kg
init_GSH_SI <- 1250 *V_SI  #initial GSH concentration in the small intestine in umol/kg

k_GSH <- 6.6 * 10^(-4) #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h
k_DNA <- 1.6 * 10^(-8) #The second-order rate constant of the reaction between cinnamaldehyde and 2ʹ-dG in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
C_PRO_L     <- 3000  #Liver
C_PRO_SI    <- 774   #Small intestine

#----DNA parameters----#
C_L_dG     <-  1.36 #Concentration of 2ʹ-dG in the liver μmol/kg liver
T_0.5      <-  38.5   #Half-life of DNA adduct in the liver in hours


#---Exposure parameters---#
Dose_in_mg <-250   #Dose in mg/kg-bw
MW         <-132.16 #The molecular weight of Cinnamaldehyde
DOSE       <-(Dose_in_mg * 70)/ MW  * 1e+6     #The administered dose in umol 



#exposure Time definition 
times=seq(0,8,by=0.01) #from 0 to 8 hours in increments of 0.01 h

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
parameters=c(P_F, P_L, P_SI, P_RP, P_SP, P_OH_F, P_OH_L, P_OH_SI, P_OH_RP, P_OH_SP, BW, V_F, V_L, V_SI, V_A, V_V, V_SI, V_RP, V_SP, Q_C, Q_F, Q_L, Q_SI, Q_RP, Q_SP, G_SYN_L, G_SYN_SI, 
             k_L_GLOS, k_SI_GLOS, init_GSH_L, init_GSH_SI, k_GSH, k_DNA, C_PRO_L, C_PRO_SI, C_L_dG, T_0.5, DOSE, Ka, Km_L_CA, Km_L_CA, Km_L_GST, Km_L_GST_G, Vsmax_L_CA, Vsmax_L_AO,
             Vsmax_L_GST,  Km_SI_CA, Km_SI_AO, Km_SI_OH, Km_SI_GST, Km_SI_GST_G, Vsmax_SI_CA, Vsmax_SI_AO, Vsmax_SI_GST)


#----Compartment definitions,and definition of the model----#
PBK_Cinnamaldehyde=function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #--Defining the compartments of the model--#
    #-Concentration in fat-#
    C_F            <- A_F       / V_F                    #Concentration in Fat in umol/kg
    C_V_F          <- C_F       / P_F                    #Concentration of cinnamaldehyde in venous blood leaving Fat in umol/l
    
    #-Concentration in Liver-#
    C_L            <- A_L       / V_L                    #Concentration Cinnamaldehyde in the Liver in umol/kg
    C_V_L          <- C_L       / P_L                    #Concentration of cinnamaldehyde in venous blood leaving the Liver in umol/l
  
     #-Concentration in the Small intestine-#
    C_SI           <- A_SI      / V_SI                   #Concentration Cinnamaldehyde in the Small intestine in umol/kg
    C_V_SI         <- C_SI      / P_SI                   #Concentration of cinnamaldehyde in venous blood leaving the Small intestine in umol/l
    
   
    #-Concentration in richly perfused tissue-#
    C_RP           <- A_RP      / V_RP                   #Concentration of Cinnamaldehyde in RP tissue in umol/kg
    C_V_RP         <- C_RP      / P_RP                   #concentration of Cinnamaldehyde in RP venous blood leaving the tissue in umol/L

    
    #-Concentration in slowly perfused tissue-#
    C_SP            <- A_SP     / V_SP                   #Concentration in SP tissue in umol/kg
    C_V_SP         <- C_SP      / P_SP                   #Concentration of Cinnamaldehyde in venous blood leaving the SP in umol/l

    
    #-Blood concentrations-#
    C_V            <- A_V       / ( V_A + V_V)           #Concentration of Cinnamaldehyde in Venous blood in umol/l
    C_A            <- C_V
    #--Differential equations--#
    
    #-amount of Cinnamaldehyde in GI cavity in umol
    dA_GI          <- -Ka * A_GI  
    
    #-Amount of Cinnamaldehyde in Venous blood in umol 
    dA_V           <- Q_F * C_V_F + (Q_L + Q_SI) * C_V_L + Q_RP * C_V_RP + Q_SP * C_V_SP - Q_C * C_V
  
    
    #-Fat-#
    dA_F           <- Q_F * (C_A - C_V_F)                            #Amount of Cinnamaldehyde in the Fat in umol 
  
    #-Liver-#
    #-Differential equations for metabolism in the liver-#
    #-Amount of cinnamaldehyde in the liver in umol-#
    dA_L           <-  Q_L * C_A + Q_SI * C_V_SI - (Q_L + Q_SI) * C_V_L# - (RM_L_CA #+ RM_L_AO + RM_L_AG_GST + RM_L_AG_CHEM + RM_L_AP + RM_L_DA_FORM)# + R_OH_M_L_C_A            #amount in mg/h in time in liver
    
    #-Amount of Cinnamaldehyde in the Small intestine-#
    dA_SI          <- Q_SI * (C_A - C_V_SI) + Ka * A_GI #- (RM_SI_CA #+ RM_SI_AO + RM_SI_AG_GST + RM_SI_AP+ RM_L_AG_CHEM)# + R_OH_M_SI_C_A
    
    #Richly perfused Tissue_#
    dA_RP    <- Q_RP * (C_A - C_V_RP)             #Amount of Cinnamaldehyde in the RP tissue in umol 
    
    #-Slowly perfused Tissue-#
    dA_SP     <- Q_SP * (C_A - C_V_SP)             #Amount of Cinnamaldehyde in the SP tissue in umol 
    
    list(c(dA_GI,dA_V,dA_F,dA_L,dA_SI,dA_RP,dA_SP))

  })
} 

#defining the begin situation of the model (in this case no prior Cinnamaldehyde present in the organs)
state <- c("A_GI"         = DOSE ,
           "A_V"          = 0,
           "A_F"          = 0,
           "A_L"          = 0,
           "A_SI"         = 0,
           "A_RP"         = 0,
           "A_SP"          = 0
)
#Running the model---------------------------------------------------------------------------------------------------------------------------------- 
PBK_Cinnamaldehyde_results = ode(y= state,times=times,func=PBK_Cinnamaldehyde, parms = parameters, hini = 0.1, method = "euler")

df_pbk_results <- as_tibble(PBK_Cinnamaldehyde_results)
df_Cinnamaldehyde_results <- as_tibble(PBK_Cinnamaldehyde_results[,c(2:17, 19:26, 28:31)])
df_Cinnamyl_alcohol_results <- as_tibble(PBK_Cinnamaldehyde_results[,c(1,4,6,14,15,23,24,28,30)])
df_Carboxylic_acid_results <- as_tibble(PBK_Cinnamaldehyde_results[,c(1,7,18)])


mass_at_t <- rowSums(df_pbk_results[800,])
mass_at_t/1e+6 - (DOSE / 1e+6)

#plotting the results of the differential equations----------------------------------------------------------------------------------------------------------------

#plotting the concentration in the chamber to evaluate it
#---Blood plot--#
pL_GSH = ggplot(df_pbk_results, aes(time, AM_Lc_GSH)) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of GSH in the liver")
pL_GSH + scale_y_log10()

pA_L = ggplot(df_pbk_results, aes(time, A_L )) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of Cinnamaldehyde in the liver")
pA_L 

pA_SP = ggplot(df_pbk_results, aes(time, A_SP )) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of Cinnamaldehyde in Slowely perfused tissue")
pA_SP

pA_V = ggplot(df_pbk_results, aes(time, A_V )) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of Cinnamaldehyde in Venous blood")
pA_V

pAA_A = ggplot(df_pbk_results, aes(time, AA_A )) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of Cinnamaldehyde in Arterial blood")
pAA_A