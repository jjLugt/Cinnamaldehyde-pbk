#author: Joris Jean van der Lugt
#date: 28-10-2021
#Human cinnamaldehyde pbk Model adapted from:  "Dose-dependent DNA adduct formation by cinnamaldehyde and other food-borne α,β-unsaturated aldehydes predicted by physiologically based in silico modelling"
PBK_Cinnamaldehyde <- RxODE({
  
  #--Defining the compartments of the model--#
  
  Rin            <- -Ka * A_GI;                   #Rate of change in Cinnamaldehyde concentration in the GI cavity in μmol/h
  
  #-------------------Blood------------------------#
  #Cinnamaldehyde#
  C_V            <- A_V       / V_V ;             #Concentration of Cinnamaldehyde in Venous blood in μmol/l
  C_A            <- A_A       / V_A;              #Concentration of Cinnamaldehyde in Arterial blood   
  
  #Cinnamyl alcohol#
  C_OH_V         <- A_OH_V    /  V_V;             #Concentration of Cinnamyl alcOHol in venous Blood in μmol/l
  C_OH_A         <- A_OH_A    /  V_A;             #Concentration of Cinnamyl alcOHol in Arterial Blood in μmol/l 
  
  #----Inhalation------#
  #Cinnamaldehyde
  
  R_Inhalation  <- -P_V *(A_inhalation_Dose/Volume_exposure_chamber);   #Rate of Cinnamaldehyde inhalation in μmol/h
  
  C_P_Art       <-(Q_Pu * C_V + -R_Inhalation )/(Q_Pu + P_V/P_B);        #Concentration of Cinnamaldehyde  in Arterial blood leaving the lung in  μmol/l
  
  C_Alvair      <- C_P_Art/P_B;                                          #Concentration of Cinnamaldehyde exhaled in μmol/l 
  
  R_Exhalation  <- P_V * C_Alvair;                                        #Rate of Cinnamaldehyde exhalation in μmol/h
  
  
  #----Pulmonary tissue 
  C_Pu           <- A_Pu / V_Pu;                         #Concentration in pulmonary tissue in μmol/l
  C_A_Pu         <- C_Pu / P_Pu;                         #Concentration of Cinnamaldehyde in Arterial blood leaving the lung in μmol/l
  R_Pu           <- Q_Pu * (C_P_Art - C_A_Pu);           #Rate of change in the concentration of Cinnamaldehyde  in μmol/h
  
  #Cinnamyl alcohol
  C_OH_Pu        <- A_OH_Pu    / V_Pu;                  #Concentration of Cinnamyl alcOHol in lung in μmol/kg
  C_OH_V_Pu      <- C_OH_Pu   / P_OH_Pu;                #Concentration of Cinnamyl alcOHol in venous blood leaving the lung in μmol/l
  R_OH_Pu        <- Q_Pu * (C_OH_V - C_OH_V_Pu);        #Rate of change in Cinnamyl alcOHol concentration in the lung in μmol/h
  
  
  #-----------FAT---------------#
  #Cinnamaldehyde#
  C_F            <- A_F       / V_F;              #Concentration in Fat in μmol/kg
  C_V_F          <- C_F       / P_F;              #Concentration of cinnamaldehyde in venous blood leaving Fat in μmol/l
  R_F            <- Q_F * (C_A - C_V_F);          #rate of change in cinnamaldehyde concentration in the Fat in μmol/h
  
  #Cinnamyl alcohol# 
  C_OH_F         <- A_OH_F    / V_F;              #Concentration of Cinnamyl alcOHol in Fat in μmol/kg
  C_OH_V_F       <- C_OH_F    / P_OH_F;           #Concentration of Cinnamyl alcOHol in venous blood leaving Fat in μmol/l
  R_OH_F         <- Q_F * (C_OH_A - C_OH_V_F);    #Rate of change in Cinnamyl alcOHol concentration in the Fat in μmol/h
  
  
  #-----------------------Richly perfused tissue-----------------#
  #Cinnamaldehyde#
  C_RP           <- A_RP      / V_RP;             #Concentration of Cinnamaldehyde in RP tissue in μmol/kg
  C_V_RP         <- C_RP      / P_RP;             #concentration of Cinnamaldehyde in RP venous blood leaving the tissue in μmol/L
  R_RP           <- Q_RP * (C_A - C_V_RP);        #rate of change in cinnamaldehyde concentration in the RP in μmol/h
  
  #Cinnamyl alcohol#
  C_OH_RP        <- A_OH_RP   / V_RP;             #Concentration of Cinnamyl alcOHol in the RP tissue in μmol/kg
  C_OH_V_RP      <- C_OH_RP   / P_OH_RP;          #Concentration of Cinnamyl alcOHol in venous blood leaving the RP tissue in μmol/l 
  R_OH_RP        <- Q_RP * (C_OH_A - C_OH_V_RP);  #Rate of change in the Cinnamyl alcohol concentration in richly perfused tissue in μmol/h
  
  #-------------------Slowly perfused tissue--------------------#  
  #Cinnamaldehyde#
  C_SP           <- A_SP     / V_SP;              #Concentration in SP tissue in μmol/kg
  C_V_SP         <- C_SP     / P_SP;              #Concentration of Cinnamaldehyde in venous blood leaving the SP in μmol/l
  R_SP           <- Q_SP * (C_A - C_V_SP);        #rate of change in cinnamaldehyde concentration in the SP tissue in μmol/h
 
  #Cinnamyl alcohol
  C_OH_SP        <- A_OH_SP   / V_SP;             #Concentration of Cinnamyl alcOHol in the SP tissue in μmol/kg
  C_OH_V_SP      <- C_OH_SP   / P_OH_SP;          #Concentration of Cinnamyl alcOHol in venous blood leaving the SP tissue in μmol/l 
  R_OH_SP        <- Q_SP * (C_OH_A - C_OH_V_SP);  #Rate of change in the Cinnamyl alcohol concentration in slowly perfused tissue in μmol/h
  
  #-----------------------Small intestine-----------------------#
  #Cinnamaldehyde#
  C_SI           <- A_SI      / V_SI;             #Concentration Cinnamaldehyde in the Small intestine in μmol/kg
  C_V_SI         <- C_SI      / P_SI;             #Concentration of cinnamaldehyde in venous blood leaving the Small intestine in μmol/l
  C_SIc_GSH      <- AM_SIc_GSH / V_SI;            #Concentration of GSH in the Small Intestine in μmol/l
  
  #Cinnamyl alcohol#
  C_OH_SI        <- A_OH_SI   / V_SI;             #Concentration of Cinnamyl alcOHol in the Small intestine in μmol/l
  C_OH_V_SI      <- C_OH_SI   / P_OH_SI;          #Concentration of Cinnamyl alcOHol in venous blood leaving the Small intestine  in μmol/l
 
  #SI metabolism
  RM_SI_GST      <- Vsmax_SI_GST * C_V_SI * C_SIc_GSH / (Km_SI_GST_G * C_V_SI + Km_SI_GST * C_SIc_GSH + C_SIc_GSH * C_V_SI);  #Amount of cinnamaldehyde metabolized in the small intestine to GSH conjugate by GST in μmol
  RM_SI_CHEM     <- k_GSH * C_V_SI * C_SIc_GSH * V_SI;             #Rate of Cinnamaldehyde binding in the small intestine to GSH in μmol/h
  RM_SI_CA       <- k_SI_CA * C_V_SI;                              #Rate of Cinnamaldehyde enzymatically oxidized to carboxylic acid in the small intestine in μmol/h
  RM_SI_AP       <- k_GSH * C_V_SI * C_PRO_SI;                     #Rate of Cinnamaldehyde protein adducts in the small intestine in μmol/h
  RM_SIc_GSH     <- G_SYN_SI -(RM_SI_GST + RM_SI_CHEM + k_SI_GLOS * AM_SIc_GSH);      #Rate of  GSH concentration in the small intestine cytosol μmol/h
  RM_SI_AO       <- (Vsmax_SI_AO * C_V_SI) / (Km_SI_AO + C_V_SI);                             #Rate of Cinnamaldehyde reduction to cinnamyl alcOHol in the small intestine in μmol
  
  #Output small intestine#
  R_OH_SI        <- Q_SI * (C_OH_A - C_OH_V_SI) + RM_SI_AO        #Rate of Cinnamyl alcohol concentration change in the small intestine in μmol/h
  
  R_SI           <- Q_SI * (C_A - C_V_SI) - Rin - RM_SI_CA - RM_SI_AP - RM_SI_GST - RM_SI_CHEM - RM_SI_AO # #Rate of change in cinnamaldehyde concentration in the SI in μmol/h
 
  
   #---------------Liver-----------------------------------#
  #-Cinnamaldehyde-#
  C_L            <- A_L       / V_L;                    #Concentration Cinnamaldehyde in the Liver in μmol/kg
  C_V_L          <- C_L       / P_L;                    #Concentration of cinnamaldehyde in venous blood leaving the Liver in μmol/l
  C_Lc_GSH       <- AM_Lc_GSH / V_L;                    #Concentration of GSH in the liver cytosol in μmol/kg  
  
  #Cinnamyl alcohol
  C_OH_L         <- A_OH_L    / V_L;                    #Concentration of Cinnamyl alcOHol in the Liver in μmol/kg
  C_OH_V_L       <- C_OH_L   / P_OH_L;                  #Concentration of Cinnamyl alcOHol in venous blood leaving the Liver in μmol/l
  
  
  #Liver metabolism
  RM_L_GST       <- k_L_GST * C_V_L * C_Lc_GSH/(Km_L_GST_G + C_Lc_GSH); #Amount of cinnamaldehyde metabolized with GSH in the liver to conjugate GST   
  RM_L_CHEM      <- k_GSH * C_V_L * C_Lc_GSH * V_L;     #Rate of Cinnamaldehyde chemically binding to GSH in the liver in μmol/h
  RM_L_CA        <- k_L_CA * C_V_L;                     #Rate of Cinnamaldehyde oxidation to carboxylic acid in the liver in μmol/h
  RM_L_AP        <- k_GSH * C_V_L * C_PRO_L;            #Rate of Cinnamaldehyde proteins adducts formation in the liver in μmol/h
  RM_Lc_GSH      <- G_SYN_L - (RM_L_GST + RM_L_CHEM + k_L_GLOS * AM_Lc_GSH);      #Rate of change in GSH concentration in the liver cytosol μmol/h
  RM_L_AO        <- (Vsmax_L_AO * C_V_L) / (Km_L_AO + C_V_L );                #Rate of Cinnamaldehyde reduced to cinnamyl alcOHol in the liver in μmol/h
  RM_OH_L_C_A    <- (Vsmax_L_OH * C_OH_V_L) /(Km_L_OH + C_OH_V_L);            #Rate of Cinnamyl alcOHol oxidized to cinnamaldehyde in the liver in μmol/h
  
  
  
  #Over all output Liver#
  R_OH_L         <- Q_L * C_OH_A + Q_SI * C_OH_V_SI  - (Q_L + Q_SI) * C_OH_V_L + RM_L_AO - RM_OH_L_C_A;  #Rate of in  Cinnamyl alcOHol concentration in the liver in μmol
  
  R_L            <- Q_L * C_A + Q_SI * C_V_SI + RM_OH_L_C_A - (Q_L + Q_SI) * C_V_L - RM_L_CA - RM_L_AP - RM_L_GST - RM_L_CHEM - RM_L_AO #Rate of change in Cinnamaldehyde concentration in the liver in μmol/h
  
  #----------------Blood------------------------#
  #Cinnamyl alcohol
  R_OH_V         <- Q_F * C_OH_V_F + (Q_L + Q_SI) * C_OH_V_L + Q_RP * C_OH_V_RP + Q_SP * C_OH_V_SP - Q_C * C_OH_V; 
  
  R_OH_A         <- Q_C * C_OH_V_Pu - (Q_F * C_OH_A + Q_L * C_OH_A + Q_SI * C_OH_A + Q_RP * C_OH_A + Q_SP * C_OH_A);
  
  #Cinnamaldehyde#
  R_V            <- Q_F * C_V_F + (Q_L + Q_SI) * C_V_L + Q_RP * C_V_RP + Q_SP * C_V_SP - Q_C * C_V;    #Rate of change in Cinnamaldehyde concentration in the venous blood in μmol
  
  R_A            <- Q_C * C_A_Pu - (Q_F * C_A + Q_L * C_A + Q_SI * C_A + Q_RP * C_A + Q_SP * C_A );       #Rate of change in Cinnamaldehyde concentration in the arterial blood in μmol 
  
 
  #----------------------------------------------Differential equations-------------------------------------------------------------------------------#
  
  #-GI tract-#
  d/dt(A_GI)       <- Rin;      #Amount of Cinnamaldehyde in GI cavity in μmol
  
  #--------Blood----------#
  #-Venous-#
  d/dt(A_V)        <- R_V;      #Amount of Cinnamaldehyde in Venous blood in μmol 
  #Cinnamyl alcohol
  d/dt(A_OH_V)     <- R_OH_V;   #Amount of Cinnamyl alcohol in Venous blood in μmol 
  #-Arterial-#
  d/dt(A_A)        <- R_A;      #Amount of Cinnamaldehyde in Arterial blood in μmol 
  #Cinnamyl alcohol
  d/dt(A_OH_A)     <- R_OH_A;   #Amount of Cinnamyl alcohol in Arterial blood in μmol
  
  #-----Lung---#
  d/dt(A_inhalation_Dose) <- R_Inhalation       #Amount of Cinnamaldehyde in the exposure chamber in μmol
  d/dt(A_Exhalation)      <- R_Exhalation       #Amount of Cinnamaldehyde exhaled in μmol
  d/dt(A_Pu)              <- R_Pu;              #Amount of Cinnamaldehyde in the lung in μmol
  d/dt(A_OH_Pu)           <- R_OH_Pu;           #Amount of Cinnamyl alcohol in the lung in μmol
  
  #-Fat calculations-#
  d/dt(A_F)        <- R_F;      #Amount of Cinnamaldehyde in the Fat in μmol                      
  #Cinnamyl alcohol
  d/dt(A_OH_F)     <- R_OH_F;   #Amount of Cinnamyl alcOHol in Fat in μmol
  
  #--------Richly perfused Tissue---------#
  #Cinnamaldehyde#
  d/dt(A_RP)        <- R_RP;    #Amount of Cinnamaldehyde in the RP tissue in μmol
  #Cinnamyl alcohol#
  d/dt(A_OH_RP)     <- R_OH_RP; #Amount of Cinnamyl alcOHol in the RP tissue in μmol
  
  #---------------Slowly perfused Tissue----------------#
  #Cinnamaldehyde
  d/dt(A_SP)        <- R_SP;    #Amount of Cinnamaldehyde in the SP tissue in μmol 
  #Cinnamyl alcohol
  d/dt(A_OH_SP)     <- R_OH_SP; #Amount of Cinnamyl alcohol in the SP tissue in μmol
  
  #-----------Si calculations--------------#
  d/dt(AM_SI_GST) <- RM_SI_GST; #Amount of Cinnamaldehyde metabolized in the small intestine to GSH conjugate by GST in μmol
  d/dt(AM_SI_CHEM)<- RM_SI_CHEM;#Amount of Cinnamaldehyde bound in the small intestine to GSH in μmol
  d/dt(AM_SI_CA)  <- RM_SI_CA;  #Amount of Cinnamaldehyde enzymatically oxidized Cinnamic acid in the small intestine in μmol
  d/dt(AM_SI_AP)  <- RM_SI_AP;  #Amount of Cinnamaldehyde protein adducts in the small intestine in μmol
  d/dt(AM_SI_AO)  <- RM_SI_AO;  #Amount of Cinnamaldehyde reduced to cinnamyl alcohol in the SI in μmol
  
  d/dt(A_SI)      <- R_SI;      #Amount of Cinnamaldehyde in the small intestine-#

  #GSH
  d/dt(AM_SIc_GSH)<- RM_SIc_GSH;#Amount of GSH in the small intestine cytosol
  d/dt(A_OH_SI)   <- R_OH_SI;   #Amount of Cinnamyl alcohol in the small intestine in μmol 
  
  #--------------------Liver calculations-------------#
  d/dt(AM_L_GST) <- RM_L_GST;   #Amount of Cinnamaldehyde metabolized with GSH in the liver to conjugate GST 
  d/dt(AM_L_CHEM)<- RM_L_CHEM;  #Amount of Cinnamaldehyde chemically bound in liver to GSH in μmol
  d/dt(AM_L_CA)  <- RM_L_CA;    #Amount of Cinnamaldehyde oxidized to cinnamic acid in the liver in μmol
  d/dt(AM_L_AP)  <- RM_L_AP;    #Amount of Cinnamaldehyde proteins adduct in the liver in μmol
  d/dt(AM_L_AO)  <- RM_L_AO;    #Amount of Cinnamaldehyde reduced to cinnamyl alcOHol in the liver in μmol
  
  d/dt(A_L)      <- R_L;        #Amount of Cinnamaldehyde in the liver in μmol
 
  #GSH
  d/dt(AM_Lc_GSH) <- RM_Lc_GSH; #Amount of GSH in the liver cytosol
  d/dt(A_OH_L)    <- R_OH_L;    #Amount of Cinnamyl alcOHol in the liver in μmol
  d/dt(AM_OH_L_C_A)  <- RM_OH_L_C_A;  #Amount of Cinnamyl alcOHol oxidized to Cinnamaldehyde in the liver in μmol
  
})

solve.pbk_rat <- solve(PBK_Cinnamaldehyde, parameters, events = ex, inits) #Solve the PBPK model
