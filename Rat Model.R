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


#---Exposure parameters---#
Dose_in_mg <-250   #Dose in mg/kg-bw
MW         <-132.16 #The molecular weight of Cinnamaldehyde
DOSE       <-(Dose_in_mg * 0.25)/ MW  * 1e+6     #The administered dose in umol 


#exposure Time definition 
times=seq(0,8,by=0.01) #from 0 to 8 hours in increments of 0.01 h

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
Km_L_GST_G  <-  100  #??????
Km_L_OH     <-  1300 #Km for enzymatic oxidationf of cinnamyl alcohol to cinnamaldehyde in the liver in umol

#--Vmax values--#

Vsmax_L_CA    <- 9.7  #Scaled Vmax for enzymatic oxidation of cinnamaldehyde in the liver in μmol/h 
Vsmax_L_AO    <- 29   #Scaled Vmax for enzymatic reduction of cinnamaldehyde in the liver in μmol/h
Vsmax_L_GST   <- 37   #Scaled Vmax for enzymatic conjugation of cinnamaldehyde with GSH in the liver in μmol/h
Vsmax_L_GST_G <- 100  #???????????
Vsmax_L_OH    <- 15   #Scaled Vmax for enzymatic oxidation of cinnamyl alcohol to cinnamaldehyde in the liver in umol/h

#----Small intestines----#
#First order rate constants
k_SI_CA  <- 3.9*10^(-3)            #Scaled first-order rate constant for enzymatic oxidation of cinnamaldehyde in the small intestine (μmol/h)

#--Michaelis menten constants SI--#
Km_SI_CA    <- 70  #Km for enzymatic oxidation of cinnamaldehyde into cinnamic acid in the Small Intestine in μM
Km_SI_AO    <- 75  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the Small Intestine in μM
Km_SI_OH    <- 290 #Km for enzymatic oxidation of cinnamly alcOHol into cinnamaldehyde in the Small Intestine in μM
Km_SI_GST   <- 600 #Km for enzymatic conjugation of cinnamaldehye with GST in the Small Intestine in μM (RAT value)
Km_SI_GST_G <- 100  #?????????

#-Vmax values-#
Vsmax_SI_CA    <- 21 #Scaled Vmax for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the Small Intestine in μmol/h 
Vsmax_SI_AO    <- 5.8 #Scaled Vmax for enzymatic reduction of cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
Vsmax_SI_OH    <- 5.0 #Scaled Vmax for enzymatic Oxidation of cinnamyl alcohol into cinnamaldehyde in the Small Intestine in μmol/h 
Vsmax_SI_GST   <- 63 #Scaled Vmax for enzymatic Conjugation of cinnamaldehyde with GSH in the in the small intestine in μmol/h (RAT value)

#Collection of all parameters so they can be entered in the function
parameters=c(P_F, P_L, P_SI, P_RP, P_SP, P_OH_F, P_OH_L, P_OH_SI, P_OH_RP, P_OH_SP, BW, V_F, V_L, V_SI, V_A, V_V, V_SI, V_RP, V_SP, Q_C, Q_F, Q_L, Q_SI, Q_RP, Q_SP, G_SYN_L, G_SYN_SI, 
             k_L_GLOS, k_SI_GLOS, init_GSH_L, init_GSH_SI, k_GSH, k_DNA, C_PRO_L, C_PRO_SI, C_L_dG, T_0.5, DOSE, Ka, k_L_CA, k_L_GST, Km_L_CA, Km_L_CA, Km_L_GST, Km_L_GST_G, Vsmax_L_CA, Vsmax_L_AO,
             Vsmax_L_GST, Vsmax_L_GST_G, Vsmax_L_OH, k_SI_CA, Km_SI_CA, Km_SI_AO, Km_SI_OH, Km_SI_GST, Km_SI_GST_G, Vsmax_SI_CA, Vsmax_SI_AO, Vsmax_SI_GST)


#----Compartment definitions,and definition of the model----#
PBK_Cinnamaldehyde=function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #--Defining the compartments of the model--#
    
    #-Concentration in fat-#
    C_F            <- A_F       / V_F                    #Concentration in Fat in umol/kg
    C_V_F          <- C_F       / P_F                    #Concentration of cinnamaldehyde in venous blood leaving Fat in umol/l
   C_OH_F         <- A_OH_F    / V_F                    #Concentration of Cinnamyl alcOHol in Fat in umol/kg
   C_OH_V_F       <- C_OH_F    / P_OH_F                 #Concentration of Cinnamyl alcOHol in venous blood leaving Fat in umol/l
    
    #-Concentration in Liver-#
    C_L            <- A_L       / V_L                    #Concentration Cinnamaldehyde in the Liver in umol/kg
    C_V_L          <- C_L       / P_L                    #Concentration of cinnamaldehyde in venous blood leaving the Liver in umol/l
    C_OH_L         <- A_OH_L    / V_L                    #Concentration of Cinnamyl alcOHol in the Liver in umol/kg
    C_OH_V_L       <- C_OH_L   / P_OH_L                 #Concentration of Cinnamyl alcOHol in venous blood leaving the Liver in umol/l
    
    RM_L_CA       <- k_L_CA * C_V_L         #Amount of Cinnamaldehyde oxidized to carboxylic acid in the liver in umol
    RM_L_AO       <- Vsmax_L_AO * C_V_L / (Km_L_AO + C_V_L )        #Amount of Cinnamaldehyde reduced to cinnamyl alcOHol in the liver in umol
    
    #-GSH in the liver-#
    C_Lc_GSH       <- AM_Lc_GSH / V_L                    #Concentration of GSH in the liver cytosol in umol/l  
    
    RM_L_AG_GST   <- k_L_GST * C_V_L * C_Lc_GSH/(Km_L_GST_G + C_Lc_GSH) #Amount of cinnamaldehyde metabolized with GSH in the liver to conjugate GST   #KM_L_GST_G is stil unknown
    RM_L_AG_CHEM  <- k_GSH * C_V_L * C_Lc_GSH * V_L                 #Amount of Cinnamaldehyde chemically bound in liver to GSH in umol
    RM_L_AP       <- k_GSH * C_V_L * C_PRO_L * V_L                  #Amount of Cinnamaldehyde proteins adducts in the liver in umol
    RM_L_DA_FORM  <- k_DNA * C_V_L * C_L_dG * V_L                   #Formation of DNA adduct in the liver 
    RM_L_DA       <- RM_L_DA_FORM - RM_L_DA * (log(2)/T_0.5)        #Amount of DNA adduct in the liver
    R_OH_M_L_C_A  <- Vsmax_L_OH *C_V_L/(Km_L_OH + C_V_L)                                #Amount of Cinnamyl alcOHol oxidized to cinnamaldehyde in the liver in umol
    RM_Lc_GSH     <- G_SYN_L * V_L * 0.9 - (RM_L_AG_GST + RM_L_AG_CHEM + k_L_GLOS * RM_Lc_GSH)  #Amount of GSH in the liver cytosol
  
     #-Concentration in the Small intestine-#
    C_SI           <- A_SI      / V_SI                   #Concentration Cinnamaldehyde in the Small intestine in umol/kg
    C_V_SI         <- C_SI      / P_SI                   #Concentration of cinnamaldehyde in venous blood leaving the Small intestine in umol/l
    C_OH_SI        <- A_OH_SI   / V_SI                   #Concentration of Cinnamyl alcOHol in the Small intestine in umol/kg
    C_OH_V_SI      <- C_OH_SI   / P_OH_SI                #Concentration of Cinnamyl alcOHol in venous blood leaving the Small intestine  in umol/l
    
    #Cinnamaldehyde metabolism in the SI
    RM_SI_CA      <- k_SI_CA * C_V_SI    #Ammount of Cinnamaldehyde enzymatically oxidized cabroxylic acid in the small intestine in umol
    RM_SI_AO      <- Vsmax_SI_AO * C_V_SI / (Km_SI_AO + C_V_SI)    #Ammount of Cinnamaldehyde reduced to cinnamyl alcOHol in the small intestine in umol
    
    #-Concentration of GSH in the Small Intestine-#
    RM_SIc_GSH    <- G_SYN_SI * V_SI * 0.9 - (RM_SI_AG_GST + RM_SI_AG_CHEM + k_SI_GLOS * RM_SIc_GSH) 
    C_SIc_GSH      <- RM_SIc_GSH/ V_SI                   #Concentration of GSH in the Small Intestine in umol/l
    
    #Cinnamaldehyde metabolism in the SI
    RM_SI_AG_GST  <- Vsmax_SI_GST * C_V_SI * C_SIc_GSH / (Km_SI_GST_G * C_V_SI + Km_SI_GST * C_SIc_GSH + C_SIc_GSH * C_V_SI)  #-amount of cinnamaldehyde metabolized in the small intestine to GSH conjugate by GST in umol
    R_OH_M_SI_C_A <- Vsmax_SI_OH * C_OH_V_SI/(Km_SI_OH + C_OH_V_SI)#Amount of Cinnamyl alcOHol enzymatically oxidized to cinnamaldehyde in the small intestine in umol  umol 
    RM_SI_AG_CHEM <- k_GSH * C_V_SI * C_SIc_GSH * V_SI             #Amount of Cinnamaldehyde bound in the small intestine to GSH in umol
    RM_SI_AP      <- k_GSH * C_V_SI * C_PRO_SI * V_SI              #Amount of Cinnamaldehyde protein adducts in the small intestine in umol
   
    #-Concentration in richly perfused tissue-#
    C_RP           <- A_RP      / V_RP                   #Concentration of Cinnamaldehyde in RP tissue in umol/kg
    C_V_RP         <- C_RP      / P_RP                   #concentration of Cinnamaldehyde in RP venous blood leaving the tissue in umol/L
    C_OH_RP        <- A_OH_RP   / V_RP                   #Concentration of Cinnamyl alcOHol in the RP tissue in umol/kg
    C_OH_V_RP      <- C_OH_RP   / P_OH_RP                #Concentration of Cinnamyl alcOHol in venous blood leaving the RP tissue in umol/l 
    
    #-Concentration in slowly perfused tissue-#
    C_SP            <- A_SP     / V_SP                   #Concentration in SP tissue in umol/kg
    C_V_SP         <- C_SP      / P_SP                   #Concentration of Cinnamaldehyde in venous blood leaving the SP in umol/l
    C_OH_SP        <- A_OH_SP   / V_SP                   #Concentration of Cinnamyl alcOHol in the SP tissue in umol/kg
    C_OH_V_SP      <- C_OH_SP   / P_OH_SP                #Concentration of Cinnamyl alcOHol in venous blood leaving the SP tissue in umol/l 
    
    #-Blood concentrations-#
    C_V            <- A_V       / ( V_A + V_V)           #Concentration of Cinnamaldehyde in Venous blood in umol/l
    C_A            <- C_V 
    C_OH_V         <- A_OH_V    / ( V_A + V_V)           #Concentration of Cinnamyl alcOHol in venous Blood in umol/l
    C_OH_A         <- C_OH_V 
    #--Differential equations--#
    
    #-amount of Cinnamaldehyde in GI cavity in umol
    dA_GI          <- -Ka * A_GI  
    
    #-Amount of Cinnamaldehyde in Venous blood in umol 
    dA_V           <- Q_F * C_V_F + (Q_L + Q_SI) * C_V_L + Q_RP * C_V_RP + Q_SP * C_V_SP - Q_C*C_V
    dAA_A          <- C_A
    
    #-Amount of Cinnamyl alcOHol in venous blood in umol 
    dA_OH_V        <- Q_F * C_OH_V_F + (Q_L + Q_SI) * C_OH_V_L + Q_RP * C_OH_V_RP + Q_SP * C_OH_V_SP - Q_C * C_OH_V 
    
    #-Fat-#
    dA_F           <- Q_F * (C_A - C_V_F)                            #Amount of Cinnamaldehyde in the Fat in umol 
    dA_OH_F        <- Q_F * (C_OH_A - C_OH_V_F)                      #Amount of Cinnamyl alcOHol in Fat in umol
    
    #-Liver-#
    #-Differential equations for metabolism in the liver-#
    dAM_L_CA       <- RM_L_CA        #Amount of Cinnamaldehyde oxidized to carboxylic acid in the liver in umol
    
    dAM_L_AO       <- RM_L_AO        #Amount of Cinnamaldehyde reduced to cinnamyl alcOHol in the liver in umol
    
    dAM_L_AG_GST   <- RM_L_AG_GST    #Amount of cinnamaldehyde metabolized with GSH in the liver to conjugate GST   #KM_L_GST_G is stil unknown
    
    dAM_L_AG_CHEM  <- RM_SI_AG_CHEM               #Amount of Cinnamaldehyde chemically bound in liver to GSH in umol
    
    dAM_L_AP       <- RM_L_AP                 #Amount of Cinnamaldehyde proteins adducts in the liver in umol
    
    dAM_L_DA_FORM  <- RM_L_DA_FORM                 #Formation of DNA adduct in the liver 
    
    dAM_L_DA       <- RM_L_DA        #Amount of DNA adduct in the liver
    
    dA_OH_M_L_C_A  <- R_OH_M_L_C_A          #Amount of Cinnamyl alcOHol oxidized to cinnamaldehyde in the liver in umol
    
    dA_OH_L        <- Q_L * C_OH_A + Q_SI *C_OH_V_SI - (Q_L+Q_SI) * C_OH_V_L + RM_L_AO - R_OH_M_L_C_A # Amount of Cinnamyl alcOHol in the liver in umol 
    
    #-Amount of cinnamaldehyde in the liver in umol-#
    dA_L           <-  Q_L * C_A + Q_SI + C_V_SI - (Q_L + Q_SI) * C_V_L - (RM_L_CA + RM_L_AO + RM_L_AG_GST + RM_L_AG_CHEM + RM_L_AP + RM_L_DA_FORM + R_OH_M_L_C_A)            #amount in mg/h in time in liver
    
    #--GSH in the Liver cytosol--#
    
    dAM_Lc_GSH     <- RM_Lc_GSH  #Amount of GSH in the liver cytosol
    
    
    #--Small intestine--#
    
    dAM_SI_CA      <- RM_SI_CA
    
    dAM_SI_AO      <- RM_SI_AO 
    
    dAM_SI_AG_GST  <- RM_SI_AG_GST         #amount of cinnamaldehyde metabolized in the small intestine to GSH conjugate by GST in umol
    
    dAM_SI_AG_CHEM <- RM_SI_AG_CHEM        #Amount of Cinnamaldehyde bound in the small intestine to GSH in umol
    
    dAM_SI_AP      <- RM_SI_AP             #Amount of Cinnamaldehyde protein adducts in the small intestine in umol
    
    dA_OH_M_SI_C_A <- R_OH_M_SI_C_A #Amount of Cinnamyl alcOHol enzymatically oxidized to cinnamaldehyde in the small intestine in umol  umol 
    
    dA_OH_SI       <- Q_SI * (C_OH_A - C_OH_V_SI) + RM_SI_AO - R_OH_M_SI_C_A #Amount of Cinnamyl alcOHol in the Small intestine in umol 
    
    #-Amount of Cinnamaldehyde in the Small intestine-#
    dA_SI          <- Q_SI * (C_A - C_V_SI) + Ka *A_GI - (RM_SI_CA + RM_SI_AO + RM_SI_AG_GST + RM_SI_AP+ RM_L_AG_CHEM + R_OH_M_SI_C_A)
    
    #-Amount of GSH in the Small intestie cytosol-#
    dAM_SIc_GSH    <- RM_SIc_GSH
    
    #Richly perfused Tissue_#
    dA_RP    <- Q_RP * (C_A - C_V_RP)             #Amount of Cinnamaldehyde in the RP tissue in umol 
    dA_OH_RP <- Q_RP * (C_OH_A - C_OH_V_RP)       #Amount of Cinnamyl alcOHol in the RP tissue in umol
    
    #-Slowly perfused Tissue-#
    dA_SP     <- Q_SP * (C_A - C_V_SP)             #Amount of Cinnamaldehyde in the SP tissue in umol 
    dA_OH_SP  <- Q_SP * (C_OH_A - C_OH_V_SP)       #Amount of Cinnamyl alcOHol in the SP tissue in umol
    
    list(c(dA_GI,dA_V,dAA_A,dA_OH_V,dA_F,dA_OH_F,dAM_L_CA,dAM_L_AO,dAM_L_AG_GST,dAM_L_AG_CHEM,dAM_L_AP,dAM_L_DA_FORM,dAM_L_DA, 
           dA_OH_M_L_C_A,dA_OH_L,dA_L,dAM_Lc_GSH,dAM_SI_CA,dAM_SI_AO,dAM_SI_AG_GST,dAM_SI_AG_CHEM,dAM_SI_AP,dA_OH_M_SI_C_A,
           dA_OH_SI,dA_SI,dAM_SIc_GSH, dA_RP,dA_OH_RP,dA_SP,dA_OH_SP))

  })
} 

#defining the begin situation of the model (in this case no chemical present in the organs)
state <- c("A_GI"         = DOSE ,
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
           "A_OH_M_L_C_A"  = 0,
           "A_OH_L"       = 0,
           "A_L"          = 0,
           "AM_Lc_GSH"    =init_GSH_L, 
           "AM_SI_CA"     = 0,
           "AM_SI_AO"     = 0,
           "AM_SI_AG_GST" = 0,
           "AM_SI_AG_CHEM"= 0,
           "AM_SI_AP"     = 0,
           "A_OH_M_SI_C_A" = 0,
           "A_OH_SI"      = 0,
           "A_SI"         = 0,
           "AM_SIc_GSH"   =init_GSH_SI,
           "A_RP"         = 0,
           "A_OH_RP"      = 0,
           "A_SP"         = 0,
           "A_OH_SP"      = 0
)
#Running the model---------------------------------------------------------------------------------------------------------------------------------- 
PBK_Cinnamaldehyde_results = ode(y= state,times=times,func=PBK_Cinnamaldehyde, parms = parameters, hini = 0.01, method = "euler")

df_pbk_results <- as_tibble(PBK_Cinnamaldehyde_results)
df_Cinnamaldehyde_results <- as_tibble(PBK_Cinnamaldehyde_results[,c(1:3,5,16,25,27,29)])
df_Cinnamyl_alcohol_results <- as_tibble(PBK_Cinnamaldehyde_results[,c(1,4,6,14,15,23,24,28,30)])
df_Carboxylic_acid_results <- as_tibble(PBK_Cinnamaldehyde_results[,c(1,7,18)])



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