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


PBK_Cinnamaldehyde <- RxODE({
    
    #--Defining the compartments of the model--#
   
    Rin            <- -Ka * A_GI;                         #Rate of change in cinnamaldehyde concentration in the GI cavity in umol/h
    
    #-------------------Blood------------------------#
    #Cinnamaldehyde#
    C_V            <- A_V       / V_V ;                   #Concentration of Cinnamaldehyde in Venous blood in umol/l
    C_A            <- A_A       / V_A                     #Concentration of cinnamaldehyde in Arterial blood   
    
    #Cinnamyl alcohol#
    C_OH_V         <- A_OH_V    /  V_V;                   #Concentration of Cinnamyl alcOHol in venous Blood in umol/l
    C_OH_A         <- A_OH_A    /  V_A;                   #Concentration of Cinnamyl alcOHol in Arterial Blood in umol/l 
    
    #-----------FAT---------------#
    #Cinnamaldehyde#
    C_F            <- A_F       / V_F;                    #Concentration in Fat in umol/kg
    C_V_F          <- C_F       / P_F;                    #Concentration of cinnamaldehyde in venous blood leaving Fat in umol/l
    R_F            <- Q_F * (C_A - C_V_F);                #rate of change in cinnamaldehyde concentration in the Fat in umol/h
    #Cinnamyl alcohol# 
    C_OH_F         <- A_OH_F    / V_F;                    #Concentration of Cinnamyl alcOHol in Fat in umol/kg
    C_OH_V_F       <- C_OH_F    / P_OH_F;                 #Concentration of Cinnamyl alcOHol in venous blood leaving Fat in umol/l
    R_OH_F         <- Q_F * (C_OH_A - C_OH_V_F);          #Rate of change in Cinnamyl alcOHol concentration in the Fat in umol/h
    
    #-----------------------Richly perfused tissue-----------------#
    #Cinnamaldehyde#
    C_RP           <- A_RP      / V_RP;                   #Concentration of Cinnamaldehyde in RP tissue in umol/kg
    C_V_RP         <- C_RP      / P_RP;                   #concentration of Cinnamaldehyde in RP venous blood leaving the tissue in umol/L
    R_RP           <- Q_RP * (C_A - C_V_RP);              #rate of change in cinnamaldehyde concentration in the RP in umol/h
    #Cinnamyl alcohol#
    C_OH_RP        <- A_OH_RP   / V_RP;                   #Concentration of Cinnamyl alcOHol in the RP tissue in umol/kg
    C_OH_V_RP      <- C_OH_RP   / P_OH_RP;                #Concentration of Cinnamyl alcOHol in venous blood leaving the RP tissue in umol/l 
    R_OH_RP        <- Q_RP * (C_OH_A - C_OH_V_RP);        #Rate of change in the Cinnamyl alcohol concentration in richly perfused tissue in umol/h
    
    #-------------------Slowly perfused tissue--------------------#  
    #Cinnamaldehyde#
    C_SP           <- A_SP     / V_SP;                                                         #Concentration in SP tissue in umol/kg
    C_V_SP         <- C_SP     / P_SP;                                                         #Concentration of Cinnamaldehyde in venous blood leaving the SP in umol/l
    R_SP           <- Q_SP * (C_A - C_V_SP);                                                   #rate of change in cinnamaldehyde concentration in the SP tissue in umol/h
    #Cinnamyl alcohol
    C_OH_SP        <- A_OH_SP   / V_SP;                   #Concentration of Cinnamyl alcOHol in the SP tissue in umol/kg
    C_OH_V_SP      <- C_OH_SP   / P_OH_SP;                #Concentration of Cinnamyl alcOHol in venous blood leaving the SP tissue in umol/l 
    R_OH_SP        <- Q_SP * (C_OH_A - C_OH_V_SP);        #Rate of change in the Cinnamyl alcohol concentration in slowly perfused tissue in umol/h
    
    #-----------------------Small intestine-----------------------#
    #Cinnamaldehyde#
    C_SI           <- A_SI      / V_SI;                             #Concentration Cinnamaldehyde in the Small intestine in umol/kg
    C_V_SI         <- C_SI      / P_SI;                              #Concentration of cinnamaldehyde in venous blood leaving the Small intestine in umol/l
    RM_SI_CA       <- k_SI_CA * C_V_SI;                              #Rate of Cinnamaldehyde enzymatically oxidized cabroxylic acid in the small intestine in umol/h
    RM_SI_AP       <- k_GSH * C_V_SI * C_PRO_SI * V_SI;              #Rate of Cinnamaldehyde protein adducts in the small intestine in umol/h
    
    #GSH#
    C_SIc_GSH      <- AM_SIc_GSH / V_SI;                                                       #Concentration of GSH in the Small Intestine in umol/l
    RM_SI_AG_GST   <- Vsmax_SI_GST * C_V_SI * C_SIc_GSH / (Km_SI_GST_G * C_V_SI + Km_SI_GST * C_SIc_GSH + C_SIc_GSH * C_V_SI);  #-amount of cinnamaldehyde metabolized in the small intestine to GSH conjugate by GST in umol
    RM_SI_AG_CHEM  <- k_GSH * C_V_SI * C_SIc_GSH * V_SI;                                       #Rate of Cinnamaldehyde binding in the small intestine to GSH in umol/h
    RM_SIc_GSH     <- G_SYN_SI * V_SI * 0.9 - (RM_SI_AG_GST + RM_SI_AG_CHEM + k_SI_GLOS * AM_SIc_GSH);         #Rate of  GSH concentration in the Smal intesinte cytosol umol/h
    
    #Cinnamyl alchol#
    C_OH_SI        <- A_OH_SI   / V_SI;                                      #Concentration of Cinnamyl alcOHol in the Small intestine in umol/l
    C_OH_V_SI      <- C_OH_SI   / P_OH_SI;                                   #Concentration of Cinnamyl alcOHol in venous blood leaving the Small intestine  in umol/l
    RM_SI_AO       <- Vsmax_SI_AO * C_V_SI / (Km_SI_AO + C_V_SI);             #Ammount of Cinnamaldehyde reduced to cinnamyl alcOHol in the small intestine in umol
    R_OH_M_SI_C_A  <- Vsmax_SI_OH * C_OH_V_SI/(Km_SI_OH + C_OH_V_SI);         #Rate of Cinnamyl alcOHol enzymatically oxidized to cinnamaldehyde in the small intestine in umol 
    R_OH_SI        <- Q_SI * (C_OH_A - C_OH_V_SI) + RM_SI_AO - R_OH_M_SI_C_A;  #Rate of Cinnamyl alcohol concentration change in the small intestine in umol/h
    
    #Over all output small intestine#
    R_SI           <- Q_SI * (C_A - C_V_SI) -Rin - (RM_SI_CA + RM_SI_AP + RM_SI_AG_GST + RM_SI_AG_CHEM + RM_SI_AO ) + R_OH_M_SI_C_A;        #Rate of change in cinnamaldehyde concentration in the SI in umol/h
    
    #---------------Liver-----------------------------------#
    #-Cinnamaldehyde-#
    C_L            <- A_L       / V_L;                                                         #Concentration Cinnamaldehyde in the Liver in umol/kg
    C_V_L          <- C_L       / P_L;                                                         #Concentration of cinnamaldehyde in venous blood leaving the Liver in umol/l
    RM_L_CA        <- k_L_CA * C_V_L;                                 #Rate of Cinnamaldehyde oxidation to carboxylic acid in the liver in umol/h
    RM_L_AP        <- k_GSH * C_V_L * C_PRO_L * V_L;                                           #Rate of Cinnamaldehyde proteins adducts formation in the liver in umol/h
    
    #GSH#
    C_Lc_GSH       <- AM_Lc_GSH / V_L;                                                          #Concentration of GSH in the liver cytosol in umol/l  
    RM_L_AG_GST    <- k_L_GST * C_V_L * C_Lc_GSH/(Km_L_GST_G + C_Lc_GSH); #Amount of cinnamaldehyde metabolized with GSH in the liver to conjugate GST   #KM_L_GST_G is stil unknown
    RM_L_AG_CHEM   <- k_GSH * C_V_L * C_Lc_GSH * V_L;                                           #Amount of Cinnamaldehyde chemically bound in liver to GSH in umol
    RM_Lc_GSH      <- G_SYN_L * V_L * 0.9 - (RM_L_AG_GST + RM_L_AG_CHEM + k_L_GLOS * AM_Lc_GSH);#Rate of change in GSH concentration in the liver cytosol umol/h
    RM_L_DA_FORM   <- k_DNA * C_V_L * C_L_dG * V_L;                                             #Rate of DNA adduct formation in the liver 
    RM_L_DA        <- RM_L_DA_FORM - AM_L_DA * (log(2)/T_0.5);                                  #Rate of DNA adduct removal in the liver
    
    #Cinnamyl alcohol
    C_OH_L         <- A_OH_L    / V_L;                    #Concentration of Cinnamyl alcOHol in the Liver in umol/kg
    C_OH_V_L       <- C_OH_L   / P_OH_L;                  #Concentration of Cinnamyl alcOHol in venous blood leaving the Liver in umol/l
    RM_L_AO        <- Vsmax_L_AO * C_V_L / (Km_L_AO + C_V_L );        #Rate of Cinnamaldehyde reduced to cinnamyl alcOHol in the liver in umol/h
    R_OH_M_L_C_A   <- Vsmax_L_OH *C_V_L/(Km_L_OH + C_V_L);                               #Rate of Cinnamyl alcOHol oxidized to cinnamaldehyde in the liver in umo/h
    R_OH_L         <- Q_L * C_OH_A + Q_SI * C_OH_V_SI + RM_L_AO - (Q_L + Q_SI) * C_OH_V_L - R_OH_M_L_C_A;#Rate of in  Cinnamyl alcOHol concentration in the liver in umol
    
    #Over al output Liver#
    R_L            <- Q_L * C_A + Q_SI * C_V_SI - (Q_L + Q_SI) * C_V_L - (RM_L_CA + RM_L_AP + RM_L_AG_GST + RM_L_AG_CHEM + RM_L_DA_FORM + RM_L_AO) + R_OH_M_L_C_A ; #Rate of change in Cinnamaldehyde concentration in the liver in umol/h
   
    #----------------Blood------------------------#
    #Cinnamaldehyde#
    R_V            <- Q_F * C_V_F + (Q_L + Q_SI) * C_V_L + Q_RP * C_V_RP + Q_SP * C_V_SP - Q_C * C_V; #Rate of change in Cinnamaldehyde concentration in the venous blood in umol
    R_A            <- Q_C * C_V - (Q_F * C_A + Q_L * C_A + Q_SI * C_A + Q_RP * C_A + Q_SP * C_A);     #Rate of change in Cinnamaldehyde concentratin in the arterial blood in umol 
    #Cinnamyl alcohol
    R_OH_V         <- Q_F * C_OH_V_F + (Q_L + Q_SI) * C_OH_V_L + Q_RP * C_OH_V_RP + Q_SP * C_OH_V_SP - Q_C * C_OH_V; 
    R_OH_A         <- Q_C * C_OH_V - (Q_F * C_OH_A + Q_L * C_OH_A + Q_SI * C_OH_A + Q_RP * C_OH_A + Q_SP * C_OH_A);
    
    
    
    #--Differential equations--#
    
    #-amount of Cinnamaldehyde in GI cavity in umol
    d/dt(A_GI)          <- Rin;           #amount of Cinnamaldehyde in GI cavity in umol
    
    #--------Blood----------#
    #-Venous-#
    d/dt(A_V)      <- R_V;            #Amount of Cinnamaldehyde in Venous blood in umol 
    d/dt(A_OH_V)   <- R_OH_V;        #Amount of Cinnamyl alcohol in Venous blood in umol 
    #-Arterial-#
    d/dt(A_A)      <- R_A;            #Amount of Cinnamaldehyde in Arterial blood in umol 
    d/dt(A_OH_A)   <- R_OH_A;        #Amount of Cinnamyl alcohol in Arterial blood in umol
    
    #-Fat calculations-#
    d/dt(A_F)      <- R_F;           #Amount of Cinnamaldehyde in the Fat in umol                      
    d/dt(A_OH_F)   <- R_OH_F;        #Amount of Cinnamyl alcOHol in Fat in umol
    
    #--------------------Liver calculations-------------#
    d/dt(A_L)           <- R_L;           #Amount of cinnamaldehyde in the liver in umol
    d/dt(AM_L_CA)       <- RM_L_CA;       #Amount of Cinnamaldehyde oxidized to carboxylic acid in the liver in umol
    d/dt(AM_L_AP)       <- RM_L_AP;       #Amount of Cinnamaldehyde proteins adducts in the liver in umol
    
    #--GSH in the Liver cytosol--#
    
    d/dt(AM_Lc_GSH)     <- RM_Lc_GSH;     #Amount of GSH in the liver cytosol
    d/dt(AM_L_AG_GST)   <- RM_L_AG_GST;   #Amount of cinnamaldehyde metabolized with GSH in the liver to conjugate GST   
    d/dt(AM_L_AG_CHEM)  <- RM_L_AG_CHEM;  #Amount of Cinnamaldehyde chemically bound in liver to GSH in umol
    
    #--Cinnamyl alcohol--#
    d/dt(A_OH_M_L_C_A)  <- R_OH_M_L_C_A;   #Amount of Cinnamyl alcOHol oxidized to cinnamaldehyde in the liver in umol
    d/dt(AM_L_AO)       <- RM_L_AO;        #Amount of Cinnamaldehyde reduced to cinnamyl alcOHol in the liver in umol
    d/dt(A_OH_L)        <- R_OH_L;         #Amount of Cinnamyl alcOHol in the liver in umol 
    
    #---DNA adduct formation--#
    d/dt(AM_L_DA_FORM)  <- RM_L_DA_FORM;  #Amount of DNA adduct in the liver 
    d/dt(AM_L_DA)       <- RM_L_DA;       #Amount of DNA adduct in the liver
    
    #-----------Si calculations--------------#
    d/dt(A_SI)          <- R_SI;          #Amount of Cinnamaldehyde in the Small intestine-#
    d/dt(AM_SI_CA)      <- RM_SI_CA;      #Amount of Cinnamaldehyde enzymatically oxidized cabroxylic acid in the small intestine in umol
    d/dt(AM_SI_AP)      <- RM_SI_AP;      #Amount of Cinnamaldehyde protein adducts in the small intestine in umol
    
    #-GSH in the Small intestine cytosol-#
    d/dt(AM_SI_AG_GST)  <- RM_SI_AG_GST;  #Amount of cinnamaldehyde metabolized in the small intestine to GSH conjugate by GST in umol
    d/dt(AM_SI_AG_CHEM) <- RM_SI_AG_CHEM; #Amount of Cinnamaldehyde bound in the small intestine to GSH in umol
    d/dt(AM_SIc_GSH)    <- RM_SIc_GSH;
    
    #--Cinnamyl alcohol--#
    d/dt(AM_SI_AO)      <- RM_SI_AO; 
    d/dt(A_OH_M_SI_C_A) <- R_OH_M_SI_C_A; #Amount of Cinnamyl alcOHol enzymatically oxidized to cinnamaldehyde in the small intestine in umol  umol 
    d/dt(A_OH_SI)       <- R_OH_SI;  #Amount of Cinnamyl alcOHol in the Small intestine in umol 
    
    #--------Richly perfused Tissue---------#
    #cinnamaldehyde#
    d/dt(A_RP)          <- R_RP;           #Amount of Cinnamaldehyde in the RP tissue in umol
    #Cinnamyl alchol#
    d/dt(A_OH_RP)       <-R_OH_RP;        #Amount of Cinnamyl alcOHol in the RP tissue in umol
    
    #---------------Slowly perfused Tissue----------------#
    #Cinnamaldehyde
    d/dt(A_SP)          <- R_SP;           #Amount of Cinnamaldehyde in the SP tissue in umol 
    #Cinnamyl alchol
    d/dt(A_OH_SP)       <-R_OH_SP;        #Amount of Cinnamyl alcOHol in the SP tissue in umol
    
    
  })

print(PBK_Cinnamaldehyde)
solve.pbk <- solve(PBK_Cinnamaldehyde, parameters, events = ex, inits, cores=4) #Solve the PBPK model




