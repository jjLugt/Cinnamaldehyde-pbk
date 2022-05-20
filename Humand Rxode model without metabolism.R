#author: Joris Jean van der Lugt
#date: 20-05-2021
#Human cinnamaldehyde pbk Model adapted from:  "Dose-dependent DNA adduct formation by cinnamaldehyde and other food-borne α,β-unsaturated aldehydes predicted by physiologically based in silico modelling"
#base model without metabolism
PBK_Cinnamaldehyde <- RxODE({
  
  #--Defining the compartments of the model--#
  
  Rin            <- -Ka * A_GI;                         #rate of change in cinnamaldehyde concentration in the GI cavity in umol/h
  
  #-Blood concentrations-#
  C_V            <- A_V     / V_V ;                     #Concentration of Cinnamaldehyde in Venous blood in umol/l
  C_A            <- A_A     / V_A                       #Concentration of cinnamaldehyde in Arterial blood   
  
  
  #-Concentration in fat-#
  C_F            <- A_F       / V_F;                    #Concentration in Fat in umol/kg
  C_V_F          <- C_F       / P_F;                    #Concentration of cinnamaldehyde in venous blood leaving Fat in umol/l
  R_F            <- Q_F * (C_A - C_V_F);                #rate of change in cinnamaldehyde concentration in the Fat in umol/h
  
  
  #-Concentration in richly perfused tissue-#
  C_RP           <- A_RP      / V_RP;                   #Concentration of Cinnamaldehyde in RP tissue in umol/kg
  C_V_RP         <- C_RP      / P_RP;                   #concentration of Cinnamaldehyde in RP venous blood leaving the tissue in umol/L
  R_RP           <- Q_RP * (C_A - C_V_RP);              #rate of change in cinnamaldehyde concentration in the RP in umol/h
  
  #-Concentration in slowly perfused tissue-#
  C_SP           <- A_SP     / V_SP;                    #Concentration in SP tissue in umol/kg
  C_V_SP         <- C_SP      / P_SP;                   #Concentration of Cinnamaldehyde in venous blood leaving the SP in umol/l
  R_SP           <- Q_SP * (C_A - C_V_SP);              #rate of change in cinnamaldehyde concentration in the SP tissue in umol/h
  
  #-Concentration in the Small intestine-#
  C_SI           <- A_SI      / V_SI;                   #Concentration Cinnamaldehyde in the Small intestine in umol/kg
  C_V_SI         <- C_SI      / P_SI;                   #Concentration of cinnamaldehyde in venous blood leaving the Small intestine in umol/l
  RM_SI_CA       <- Vsmax_SI_CA * C_V_SI / (Km_SI_CA + C_V_SI);                  #Rate of Cinnamaldehyde enzymatically oxidized cabroxylic acid in the small intestine in umol/h
  RM_SI_AP       <- k_GSH * C_V_SI * C_PRO_SI * V_SI;                            #Amount of Cinnamaldehyde protein adducts in the small intestine in umol
  R_SI           <- Q_SI * (C_A - C_V_SI) -Rin - (RM_SI_CA + RM_SI_AP)           #Rate of change in cinnamaldehyde concentration in the SI in umol
  
  #-Concentration in Liver-#
  C_L            <- A_L       / V_L;                    #Concentration Cinnamaldehyde in the Liver in umol/kg
  C_V_L          <- C_L       / P_L;                    #Concentration of cinnamaldehyde in venous blood leaving the Liver in umol/l
  RM_L_CA        <- Vsmax_L_CA * C_V_L / (Km_L_CA + C_V_L);                                  #Rate of Cinnamaldehyde oxidation to carboxylic acid in the liver in umol/h
  RM_L_AP        <- k_GSH * C_V_L * C_PRO_L * V_L;                                           #Rate of Cinnamaldehyde proteins adducts formation in the liver in umol/h
  R_L            <- Q_L * C_A + Q_SI * C_V_SI - (Q_L + Q_SI) * C_V_L - (RM_L_CA + RM_L_AP) ; #Rate of change in Cinnamaldehyde concentration in the liver in umol/h
  
  #Rate of Cinnamaldehyde in Blood
  R_V            <- Q_F * C_V_F + (Q_L + Q_SI) * C_V_L + Q_RP * C_V_RP + Q_SP * C_V_SP - Q_C * C_V; #Rate of change in Cinnamaldehyde concentration in the venous blood in umol
  R_A            <- Q_C * C_V - (Q_F * C_A + Q_L * C_A + Q_SI * C_A + Q_RP * C_A + Q_SP * C_A);     #Rate of change in Cinnamaldehyde concentratin in the arterial blood in umol 
  
  
  
  
  
  #--Differential equations--#
  
  #-GI tract-#
  d/dt(A_GI)     <- Rin;           #amount of Cinnamaldehyde in GI cavity in umol
  
  
  #-Venous blood calculations-#
  d/dt(A_V)      <-R_V;            #Amount of Cinnamaldehyde in Venous blood in umol 
  
  #-Arterial blood calculations-#
  d/dt(A_A)      <-R_A;            #Amount of Cinnamaldehyde in Arterial blood in umol 
  
  #-Fat calculations-#
  d/dt(A_F)      <- R_F;           #Amount of Cinnamaldehyde in the Fat in umol                      

  
  #-Liver calculations-#
  d/dt(A_L)      <- R_L;           #Amount of cinnamaldehyde in the liver in umol
  d/dt(AM_L_CA)  <- RM_L_CA;       #Amount of Cinnamaldehyde oxidized to carboxylic acid in the liver in umol
  d/dt(AM_L_AP)  <- RM_L_AP;       #Amount of Cinnamaldehyde proteins adducts in the liver in umol
 
  #-Si calculations-#
  d/dt(A_SI)     <- R_SI;          #Amount of Cinnamaldehyde in the Small intestine-#
  d/dt(AM_SI_CA) <- RM_SI_CA;      #Amount of Cinnamaldehyde enzymatically oxidized cabroxylic acid in the small intestine in umol
  d/dt(AM_SI_AP) <- RM_SI_AP;      #Amount of Cinnamaldehyde protein adducts in the small intestine in umol
  
  #-Richly perfused Tissue calculations-#
  d/dt(A_RP)     <- R_RP;          #Amount of Cinnamaldehyde in the RP tissue in umol
  
  
  #-Slowly perfused Tissue calculaitons-#
  d/dt(A_SP)     <- R_SP;          #Amount of Cinnamaldehyde in the SP tissue in umol 
  
  
})

solve.pbk_nonpop <- solve(PBK_Cinnamaldehyde, parameters, events = ex, inits, cores=4) #Solve the PBPK model


