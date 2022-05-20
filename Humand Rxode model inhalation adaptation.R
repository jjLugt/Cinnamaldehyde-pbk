#author: Joris Jean van der Lugt
#date: 20-05-2021
#Human cinnamaldehyde inhalation pbk Model adapted from:  "Dose-dependent DNA adduct formation by cinnamaldehyde and other food-borne α,β-unsaturated aldehydes predicted by physiologically based in silico modelling"
library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(shiny)
library(truncnorm)
library(reshape2)



PBK_Cinnamaldehyde <- RxODE({
  
  #--Defining the compartments of the model--#
  
  #Concentration in lung
  C_I          <-  A_I / Q_P  ;                         #inhaled concentration umol/l
  C_PU         <- (A_I - A_X)/ V_PU;                    #concentration in pulmonary tissue
  C_V_PU       <- C_PU       / P_PU;                    #Concentration of cinnamaldehyde in venous blood leaving Fat in umol/l
  C_OH_PU      <- A_OH_PU    / V_PU;                    #Concentration of Cinnamyl alcOHol in Fat in umol/kg
  C_OH_V_PU    <- C_OH_PU    / P_OH_RP;                 #Concentration of Cinnamyl alcOHol in venous blood leaving Fat in umol/l
  
  
  C_X <- C_PU / P_B ;                      #exhaled concentration mg/l
  
  
  #-Concentration in fat-#
  C_F            <- A_F       / V_F;                    #Concentration in Fat in umol/kg
  C_V_F          <- C_F       / P_F;                    #Concentration of cinnamaldehyde in venous blood leaving Fat in umol/l
  C_OH_F         <- A_OH_F    / V_F;                    #Concentration of Cinnamyl alcOHol in Fat in umol/kg
  C_OH_V_F       <- C_OH_F    / P_OH_F;                 #Concentration of Cinnamyl alcOHol in venous blood leaving Fat in umol/l
  
  #-Concentration in Liver-#
  C_L            <- A_L       / V_L;                    #Concentration Cinnamaldehyde in the Liver in umol/kg
  C_V_L          <- C_L       / P_L;                    #Concentration of cinnamaldehyde in venous blood leaving the Liver in umol/l
  C_OH_L         <- A_OH_L    / V_L;                    #Concentration of Cinnamyl alcOHol in the Liver in umol/kg
  C_OH_V_L       <- C_OH_L   / P_OH_L;                 #Concentration of Cinnamyl alcOHol in venous blood leaving the Liver in umol/l
  
  RM_L_CA       <- Vsmax_L_CA * C_V_L / (Km_L_CA + C_V_L);         #Amount of Cinnamaldehyde oxidized to carboxylic acid in the liver in umol
  RM_L_AO       <- Vsmax_L_AO * C_V_L / (Km_L_AO + C_V_L );        #Amount of Cinnamaldehyde reduced to cinnamyl alcOHol in the liver in umol
  
  #-GSH in the liver-#
  C_Lc_GSH       <- AM_Lc_GSH / V_L;                    #Concentration of GSH in the liver cytosol in umol/l  
  
  RM_L_AG_GST   <- Vsmax_L_GST * C_V_L * C_Lc_GSH/(Km_L_GST_G * C_V_L + Km_L_GST * C_Lc_GSH + C_Lc_GSH * C_V_L); #Amount of cinnamaldehyde metabolized with GSH in the liver to conjugate GST   #KM_L_GST_G is stil unknown
  RM_L_AG_CHEM  <- k_GSH * C_V_L * C_Lc_GSH * V_L;                 #Amount of Cinnamaldehyde chemically bound in liver to GSH in umol
  RM_L_AP       <- k_GSH * C_V_L * C_PRO_L * V_L;                  #Amount of Cinnamaldehyde proteins adducts in the liver in umol
  RM_L_DA_FORM  <- k_DNA * C_V_L * C_L_dG * V_L;                   #Formation of DNA adduct in the liver 
  RM_L_DA       <- RM_L_DA_FORM - RM_L_DA * (log(2)/T_0.5);        #Amount of DNA adduct in the liver
  R_OH_M_L_C_A  <- k_L_OH * C_OH_V_L;                                #Amount of Cinnamyl alcOHol oxidized to cinnamaldehyde in the liver in umol
  RM_Lc_GSH     <- G_SYN_L * V_L * 0.9 - (RM_L_AG_GST + RM_L_AG_CHEM + k_L_GLOS * RM_Lc_GSH);  #Amount of GSH in the liver cytosol
  
  #-Concentration in the Small intestine-#
  C_SI           <- A_SI      / V_SI;                   #Concentration Cinnamaldehyde in the Small intestine in umol/kg
  C_V_SI         <- C_SI      / P_SI;                   #Concentration of cinnamaldehyde in venous blood leaving the Small intestine in umol/l
  C_OH_SI        <- A_OH_SI   / V_SI;                   #Concentration of Cinnamyl alcOHol in the Small intestine in umol/kg
  C_OH_V_SI      <- C_OH_SI   / P_OH_SI;                #Concentration of Cinnamyl alcOHol in venous blood leaving the Small intestine  in umol/l
  
  #Cinnamaldehyde metabolism in the SI
  RM_SI_CA      <- Vsmax_SI_CA * C_V_SI / (Km_SI_CA + C_V_SI);    #Ammount of Cinnamaldehyde enzymatically oxidized cabroxylic acid in the small intestine in umol
  RM_SI_AO      <- Vsmax_SI_AO * C_V_SI / (Km_SI_AO + C_V_SI);    #Ammount of Cinnamaldehyde reduced to cinnamyl alcOHol in the small intestine in umol
  
  #-Concentration of GSH in the Small Intestine-#
  RM_SIc_GSH    <- G_SYN_SI * V_SI * 0.9 - (RM_SI_AG_GST + RM_SI_AG_CHEM + k_SI_GLOS * RM_SIc_GSH); 
  C_SIc_GSH      <- RM_SIc_GSH/ V_SI;                   #Concentration of GSH in the Small Intestine in umol/l
  
  #Cinnamaldehyde metabolism in the SI
  RM_SI_AG_GST  <- Vsmax_SI_GST * C_V_SI * C_SIc_GSH / (Km_SI_GST_G * C_V_SI + Km_SI_GST * C_SIc_GSH + C_SIc_GSH * C_V_SI);  #-amount of cinnamaldehyde metabolized in the small intestine to GSH conjugate by GST in umol
  R_OH_M_SI_C_A <- Vsmax_SI_OH * C_OH_V_SI/(Km_SI_OH + C_OH_V_SI);#Amount of Cinnamyl alcOHol enzymatically oxidized to cinnamaldehyde in the small intestine in umol  umol 
  RM_SI_AG_CHEM <- k_GSH * C_V_SI * C_SIc_GSH * V_SI;             #Amount of Cinnamaldehyde bound in the small intestine to GSH in umol
  RM_SI_AP      <- k_GSH * C_V_SI * C_PRO_SI * V_SI;              #Amount of Cinnamaldehyde protein adducts in the small intestine in umol
  
  #-Concentration in richly perfused tissue-#
  C_RP           <- A_RP      / V_RP;                   #Concentration of Cinnamaldehyde in RP tissue in umol/kg
  C_V_RP         <- C_RP      / P_RP;                   #concentration of Cinnamaldehyde in RP venous blood leaving the tissue in umol/L
  C_OH_RP        <- A_OH_RP   / V_RP;                   #Concentration of Cinnamyl alcOHol in the RP tissue in umol/kg
  C_OH_V_RP      <- C_OH_RP   / P_OH_RP;                #Concentration of Cinnamyl alcOHol in venous blood leaving the RP tissue in umol/l 
  
  #-Concentration in slowly perfused tissue-#
  C_SP            <- A_SP     / V_SP;                   #Concentration in SP tissue in umol/kg
  C_V_SP         <- C_SP      / P_SP;                   #Concentration of Cinnamaldehyde in venous blood leaving the SP in umol/l
  C_OH_SP        <- A_OH_SP   / V_SP;                   #Concentration of Cinnamyl alcOHol in the SP tissue in umol/kg
  C_OH_V_SP      <- C_OH_SP   / P_OH_SP;                #Concentration of Cinnamyl alcOHol in venous blood leaving the SP tissue in umol/l 
  
  #-Blood concentrations-#
  C_V            <- A_V       / ( V_A + V_V);           #Concentration of Cinnamaldehyde in Venous blood in umol/l
  C_A            <- C_V; 
  C_OH_V         <- A_OH_V    / ( V_A + V_V);           #Concentration of Cinnamyl alcOHol in venous Blood in umol/l
  C_OH_A         <- C_OH_V; 
  #--Differential equations--#
  
  #-amount of Cinnamaldehyde in GI cavity in umol
  d/dt(A_GI)          <- -Ka * A_GI;  
  
  #amount of cinnamaldehyde inhaled/exhaled  in umol
  d/dt(A_I)           <- (Q_P * C_I) - A_X ;
  d/dt(A_X)           <- Q_P * C_X ;
  d/dt(A_OH_PU)       <- Q_PU * (C_OH_A - C_OH_V_PU);
  
  
  #-Amount of Cinnamaldehyde in Venous blood in umol 
  d/dt(A_V)           <- Q_F * C_V_F + (Q_L + Q_SI) * C_V_L + Q_RP * C_V_RP + Q_SP * C_V_SP + Q_PU * C_PU - Q_C*C_V;
  #d/dt(AA_A)          <- 
  
  #-Amount of Cinnamyl alcOHol in venous blood in umol 
  d/dt(A_OH_V)        <- Q_F * C_OH_V_F + (Q_L + Q_SI) * C_OH_V_L + Q_RP * C_OH_V_RP + Q_SP * C_OH_V_SP + Q_PU * C_OH_V_PU - Q_C * C_OH_V; 
  
  #-Fat-#
  d/dt(A_F)           <- Q_F * (C_A - C_V_F);                            #Amount of Cinnamaldehyde in the Fat in umol 
  d/dt(A_OH_F)        <- Q_F * (C_OH_A - C_OH_V_F);                      #Amount of Cinnamyl alcOHol in Fat in umol
  
  #-Liver-#
  #-Differential equations for metabolism in the liver-#
  d/dt(AM_L_CA)       <- RM_L_CA;       #Amount of Cinnamaldehyde oxidized to carboxylic acid in the liver in umol
  
  d/dt(AM_L_AO)       <- RM_L_AO;        #Amount of Cinnamaldehyde reduced to cinnamyl alcOHol in the liver in umol
  
  d/dt(AM_L_AG_GST)   <- RM_L_AG_GST;    #Amount of cinnamaldehyde metabolized with GSH in the liver to conjugate GST   #KM_L_GST_G is stil unknown
  
  d/dt(AM_L_AG_CHEM)  <- RM_SI_AG_CHEM;               #Amount of Cinnamaldehyde chemically bound in liver to GSH in umol
  
  d/dt(AM_L_AP)       <- RM_L_AP;                 #Amount of Cinnamaldehyde proteins adducts in the liver in umol
  
  d/dt(AM_L_DA_FORM)  <- RM_L_DA_FORM;                 #Formation of DNA adduct in the liver 
  
  d/dt(AM_L_DA)       <- RM_L_DA;        #Amount of DNA adduct in the liver
  
  d/dt(A_OH_M_L_C_A)  <- R_OH_M_L_C_A;          #Amount of Cinnamyl alcOHol oxidized to cinnamaldehyde in the liver in umol
  
  d/dt(A_OH_L)        <- Q_L * C_OH_A + Q_SI *C_OH_V_SI - (Q_L+Q_SI) * C_OH_V_L + RM_L_AO - R_OH_M_L_C_A; # Amount of Cinnamyl alcOHol in the liver in umol 
  
  #-Amount of cinnamaldehyde in the liver in umol-#
  d/dt(A_L)           <-  Q_L * C_A + Q_SI * C_V_SI - (Q_L + Q_SI) * C_V_L - (RM_L_CA + RM_L_AO + RM_L_AG_GST + RM_L_AG_CHEM + RM_L_AP + RM_L_DA_FORM + R_OH_M_L_C_A);            #amount in mg/h in time in liver
  
  #--GSH in the Liver cytosol--#
  
  d/dt(AM_Lc_GSH)     <- RM_Lc_GSH;  #Amount of GSH in the liver cytosol
  
  
  #--Small intestine--#
  
  d/dt(AM_SI_CA)      <- RM_SI_CA;
  
  d/dt(AM_SI_AO)      <- RM_SI_AO; 
  
  d/dt(AM_SI_AG_GST)  <- RM_SI_AG_GST;         #amount of cinnamaldehyde metabolized in the small intestine to GSH conjugate by GST in umol
  
  d/dt(AM_SI_AG_CHEM) <- RM_SI_AG_CHEM;        #Amount of Cinnamaldehyde bound in the small intestine to GSH in umol
  
  d/dt(AM_SI_AP)      <- RM_SI_AP;             #Amount of Cinnamaldehyde protein adducts in the small intestine in umol
  
  d/dt(A_OH_M_SI_C_A) <- R_OH_M_SI_C_A; #Amount of Cinnamyl alcOHol enzymatically oxidized to cinnamaldehyde in the small intestine in umol  umol 
  
  d/dt(A_OH_SI)       <- Q_SI * (C_OH_A - C_OH_V_SI) + RM_SI_AO - R_OH_M_SI_C_A; #Amount of Cinnamyl alcOHol in the Small intestine in umol 
  
  #-Amount of Cinnamaldehyde in the Small intestine-#
  d/dt(A_SI)          <- Q_SI * (C_A - C_V_SI) + Ka *A_GI - (RM_SI_CA + RM_SI_AO + RM_SI_AG_GST + RM_SI_AP+ RM_L_AG_CHEM + R_OH_M_SI_C_A);
  
  #-Amount of GSH in the Small intestie cytosol-#
  d/dt(AM_SIc_GSH)    <- RM_SIc_GSH;
  
  #Richly perfused Tissue_#
  d/dt(A_RP)    <- Q_RP * (C_A - C_V_RP);             #Amount of Cinnamaldehyde in the RP tissue in umol 
  d/dt(A_OH_RP) <- Q_RP * (C_OH_A - C_OH_V_RP);       #Amount of Cinnamyl alcOHol in the RP tissue in umol
  
  #-Slowly perfused Tissue-#
  d/dt(A_SP)     <- Q_SP * (C_A - C_V_SP);             #Amount of Cinnamaldehyde in the SP tissue in umol 
  d/dt(A_OH_SP)  <- Q_SP * (C_OH_A - C_OH_V_SP);       #Amount of Cinnamyl alcOHol in the SP tissue in umol
})

print(PBK_Cinnamaldehyde)
solve.pbk_nonpop <- solve(PBK_Cinnamaldehyde, parameters, events = ex, inits, cores=4) #Solve the PBPK model

