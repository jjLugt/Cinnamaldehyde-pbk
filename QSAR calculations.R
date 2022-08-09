#author: Joris Jean van der Lugt
#date: 05-08-2022
#De Jongh et al 2011 Qsar calculation 

#EPIsuite Log Kow QSAR
Log_Kow_Benzaldehyde    <- 1.71
Log_Kow_Cinnamaldehyde  <- 1.82
Log_Kow_Cinnamylalcohol <- 1.84

#EPIsuite Log Koa QSAR
Log_Koa_Benzaldehyde    <- 4.79
Log_Koa_Cinnamaldehyde  <- 6

#EPIsuite Henry law constant
HL_Benzaldehyde   <- 2.67E-05 #atm-m3/mole
HL_Cinnamaldehyde <- 3.54E-06 #atm-m3/mole

#Water and lipid fractions of human tissues source: Table 1 De Jongh et al 1997
Blood_Water_Fraction   <- 0.83
Blood_Lipid_Fraction   <- 0.0056
Adipose_Water_Fraction <- 0.200
Adipose_Lipid_Fraction <- 0.800
Liver_Water_Fraction   <- 0.711
Liver_Lipid_Fraction   <- 0.049
Muscle_Water_Fraction  <- 0.792
Muscle_Lipid_Fraction  <- 0.031
Brain_Water_Fraction   <- 0.775
Brain_Lipid_Fraction   <- 0.133

#RAT
Blood_Water_Fraction_R   <- 0.820
Blood_Lipid_Fraction_R   <- 0.0036
Adipose_Lipid_Fraction_R <- 0.200
Adipose_Lipid_Fraction_R <- 0.800
Liver_Water_Fraction_R   <- 0.740
Liver_Lipid_Fraction_R   <- 0.065
Muscle_Water_Fraction_R  <- 0.720
Muscle_Lipid_Fraction_R  <- 0.043

#A and B values for the Tissue partition coefficients for Humans Table 2 de Jongh et al 1997
#A
Fat_Blood_A     <- 1.03 #+- 0.03
Liver_Blood_A   <- 0.81 #+- 0.02
Muscle_Blood_A  <- 0.81 #+- 0.05
Kidney_Blood_A  <- 0.57 #+- 0.05
Brain_Blood_A   <- 0.48 #+- 0.03

Fat_Blood_B     <- -0.38 #+- 0.25
Liver_Blood_B   <- -0.35 #+- 0.10
Muscle_Blood_B  <- -0.22 #+- 0.02
Kidney_Blood_B  <- -0.19 #+- 0.01
Brain_Blood_B   <- -0.21 #+- 0.02

#Rat
#A
Fat_Blood_A_R     <- 0.70 #+- 0.02
Liver_Blood_A_R   <- 0.44 #+- 0.03
Muscle_Blood_A_R <- 0.29 #+- 0.06
#B
Fat_Blood_B_R     <- -0.02 #+- 0.73
Liver_Blood_B_R   <- -0.19 #+- 0.11
Muscle_Blood_B_R  <- -0.55 #+- 0.06

#--------------lipid/water partition coefficients

#Fat/Blood
Plw_Fat_Blood_Benzaldehyde <- Log_Kow_Benzaldehyde ^ Fat_Blood_A
Plw_Fat_Blood_Cinnamaldehyde <- Log_Kow_Cinnamaldehyde ^ Fat_Blood_A
Plw_Fat_Blood_Cinnamylalcohol <- Log_Kow_Cinnamylalcohol ^ Fat_Blood_A

#Liver/Blood
Plw_Liver_Blood_Benzaldehyde <- Log_Kow_Benzaldehyde ^ Liver_Blood_A
Plw_Liver_Blood_Cinnamaldehyde <- Log_Kow_Cinnamaldehyde ^ Liver_Blood_A
Plw_Liver_Blood_Cinnamylalcohol <- Log_Kow_Cinnamylalcohol ^ Liver_Blood_A

#Muscle/Slowly perfused tissue
Plw_Muscle_Blood_Benzaldehyde <- Log_Kow_Benzaldehyde ^ Muscle_Blood_A
Plw_Muscle_Blood_Cinnamaldehyde <- Log_Kow_Cinnamaldehyde ^ Muscle_Blood_A
Plw_Muscle_Blood_Cinnamylalcohol <- Log_Kow_Cinnamylalcohol ^ Muscle_Blood_A

#Muscle/Slowly perfused tissue
Plw_Brain_Blood_Benzaldehyde <- Log_Kow_Benzaldehyde ^ Brain_Blood_A
Plw_Brain_Blood_Cinnamaldehyde <- Log_Kow_Cinnamaldehyde ^ Brain_Blood_A
Plw_Brain_Blood_Cinnamylalcohol <- Log_Kow_Cinnamylalcohol ^ Brain_Blood_A


#----------------Partition coefficient 
#Fat/Blood 
P_F_Benzaldehyde   <-  (Adipose_Lipid_Fraction * Plw_Fat_Blood_Benzaldehyde + Adipose_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Fat_Blood_Benzaldehyde + Blood_Water_Fraction)+ Fat_Blood_B
P_F_Cinnamaldehyde <-  (Adipose_Lipid_Fraction * Plw_Fat_Blood_Cinnamaldehyde + Adipose_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Fat_Blood_Cinnamaldehyde + Blood_Water_Fraction) + Fat_Blood_B
P_F_Cinnamylalcohol <-  (Adipose_Lipid_Fraction * Plw_Fat_Blood_Cinnamylalcohol + Adipose_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Fat_Blood_Cinnamylalcohol + Blood_Water_Fraction) + Fat_Blood_B

#Liver/Blood 
P_L_Benzaldehyde   <-  (Liver_Lipid_Fraction * Plw_Liver_Blood_Benzaldehyde + Liver_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Liver_Blood_Benzaldehyde + Blood_Water_Fraction)+ Liver_Blood_B
P_L_Cinnamaldehyde <-  (Liver_Lipid_Fraction * Plw_Liver_Blood_Cinnamaldehyde + Liver_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Liver_Blood_Cinnamaldehyde + Blood_Water_Fraction) + Liver_Blood_B
P_L_Cinnamylalcohol <-  (Liver_Lipid_Fraction * Plw_Liver_Blood_Cinnamylalcohol + Liver_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Liver_Blood_Cinnamylalcohol + Blood_Water_Fraction) + Liver_Blood_B

#Muscle/Blood
P_SP_Benzaldehyde   <-  (Muscle_Lipid_Fraction * Plw_Muscle_Blood_Benzaldehyde + Muscle_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Muscle_Blood_Benzaldehyde + Blood_Water_Fraction)+ Muscle_Blood_B
P_SP_Cinnamaldehyde <-  (Muscle_Lipid_Fraction * Plw_Muscle_Blood_Cinnamaldehyde + Muscle_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Muscle_Blood_Cinnamaldehyde + Blood_Water_Fraction) + Muscle_Blood_B
P_SP_Cinnamylalcohol <-  (Muscle_Lipid_Fraction * Plw_Muscle_Blood_Cinnamylalcohol + Muscle_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Muscle_Blood_Cinnamylalcohol + Blood_Water_Fraction) + Muscle_Blood_B

#Brain/Blood
P_Br_Benzaldehyde   <-  (Brain_Lipid_Fraction * Plw_Brain_Blood_Benzaldehyde + Brain_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Brain_Blood_Benzaldehyde + Blood_Water_Fraction)+ Brain_Blood_B
P_Br_Cinnamaldehyde <-  (Brain_Lipid_Fraction * Plw_Brain_Blood_Cinnamaldehyde + Brain_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Brain_Blood_Cinnamaldehyde + Blood_Water_Fraction) + Brain_Blood_B
P_Br_Cinnamylalcohol <-  (Brain_Lipid_Fraction * Plw_Brain_Blood_Cinnamylalcohol + Brain_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Brain_Blood_Cinnamylalcohol + Blood_Water_Fraction) + Brain_Blood_B

#Blood:air partition coefficient calculations based on ten Berge et al 2011
P_B_Cinnamaldehyde <- 0.4445/HL_Cinnamaldehyde + 0.005189 * Log_Koa_Cinnamaldehyde
P_B_Benzaldehdye <- 0.4445/HL_Benzaldehyde + 0.005189 * Log_Koa_Benzaldehyde



