#author: Joris Jean van der Lugt
#date: 05-08-2022
#De Jongh et al 2011 Qsar calculation 

#EPIsuite Log Kow QSAR
Log_Kow_Benzaldehyde    <- 1.71
Kow_Benzaldehyde       <- 1*10^1.71
Log_Kow_Cinnamaldehyde  <- 1.82
Kow_Cinnamaldehyde      <- 1*10^1.82
Log_Kow_Cinnamylalcohol <- 1.84
Kow_Cinnamylalcohol     <- 1*10^1.84

#EPIsuite Log Koa QSAR
Log_Koa_Benzaldehyde    <- 4.79
Koa_Benzaldehyde        <- 1*10^4.79
Log_Koa_Cinnamaldehyde  <- 1.119928             #calculated using Log_Kow_Cinnamaldehyde - log(Cinnamaldehyde_Henry_DL)
Koa_Cinnamaldehyde      <- 1*10^Log_Koa_Cinnamaldehyde

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
Adipose_Water_Fraction_R <- 0.200
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
#Human
Plw_Fat_Blood_Benzaldehyde <- Kow_Benzaldehyde ^ Fat_Blood_A
Plw_Fat_Blood_Cinnamaldehyde <- Kow_Cinnamaldehyde ^ Fat_Blood_A
Plw_Fat_Blood_Cinnamylalcohol <- Kow_Cinnamylalcohol ^ Fat_Blood_A
#Rat
Plw_Fat_Blood_Benzaldehyde_R <- Kow_Benzaldehyde^Fat_Blood_A_R
Plw_Fat_Blood_Cinnamaldehyde_R <-Kow_Cinnamaldehyde^Fat_Blood_A_R
Plw_Fat_Blood_Cinnamylalcohol_R <- Kow_Cinnamylalcohol^Fat_Blood_A_R

#Liver/Blood
#Human
Plw_Liver_Blood_Benzaldehyde <- Kow_Benzaldehyde ^ Liver_Blood_A
Plw_Liver_Blood_Cinnamaldehyde <- Kow_Cinnamaldehyde ^ Liver_Blood_A
Plw_Liver_Blood_Cinnamylalcohol <- Kow_Cinnamylalcohol ^ Liver_Blood_A
#Rat
Plw_Liver_Blood_Benzaldehyde_R <-Kow_Benzaldehyde ^ Liver_Blood_A_R
Plw_Liver_Blood_Cinnamaldehyde_R <- Kow_Cinnamaldehyde ^ Liver_Blood_A_R
Plw_Liver_Blood_Cinnamylalcohol_R <- Kow_Cinnamylalcohol ^ Liver_Blood_A_R

#Muscle/Slowly perfused tissue
#Human
Plw_Muscle_Blood_Benzaldehyde <- Kow_Benzaldehyde ^ Muscle_Blood_A
Plw_Muscle_Blood_Cinnamaldehyde <- Kow_Cinnamaldehyde ^ Muscle_Blood_A
Plw_Muscle_Blood_Cinnamylalcohol <-Kow_Cinnamylalcohol ^ Muscle_Blood_A
#Rat
Plw_Muscle_Blood_Benzaldehyde_R <- Kow_Benzaldehyde ^ Muscle_Blood_A_R
Plw_Muscle_Blood_Cinnamaldehyde_R <- Kow_Cinnamaldehyde ^ Muscle_Blood_A_R
Plw_Muscle_Blood_Cinnamylalcohol_R <- Kow_Cinnamylalcohol ^ Muscle_Blood_A_R

#Brain/Richly perfused tissue
#Human
Plw_Brain_Blood_Benzaldehyde <- Kow_Benzaldehyde ^ Brain_Blood_A
Plw_Brain_Blood_Cinnamaldehyde <- Kow_Cinnamaldehyde ^ Brain_Blood_A
Plw_Brain_Blood_Cinnamylalcohol <- Kow_Cinnamylalcohol ^ Brain_Blood_A
#Rat
Plw_Brain_Blood_Benzaldehyde_R <- Kow_Benzaldehyde ^ Brain_Blood_A_R
Plw_Brain_Blood_Cinnamaldehyde_R <- Kow_Cinnamaldehyde ^ Brain_Blood_A_R
Plw_Brain_Blood_Cinnamylalcohol_R <- Kow_Cinnamylalcohol ^ Brain_Blood_A_R

#----------------Partition coefficient 
#Fat/Blood

#Human
P_F_Benzaldehyde   <-  (Adipose_Lipid_Fraction * Plw_Fat_Blood_Benzaldehyde + Adipose_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Fat_Blood_Benzaldehyde + Blood_Water_Fraction)+ Fat_Blood_B
P_F_Cinnamaldehyde <-  (Adipose_Lipid_Fraction * Plw_Fat_Blood_Cinnamaldehyde + Adipose_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Fat_Blood_Cinnamaldehyde + Blood_Water_Fraction) + Fat_Blood_B
P_F_Cinnamylalcohol <-  (Adipose_Lipid_Fraction * Plw_Fat_Blood_Cinnamylalcohol + Adipose_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Fat_Blood_Cinnamylalcohol + Blood_Water_Fraction) + Fat_Blood_B

#Rat
P_F_Benzaldehyde_R   <-  (Adipose_Lipid_Fraction_R * Plw_Fat_Blood_Benzaldehyde_R + Adipose_Water_Fraction_R)/(Blood_Lipid_Fraction_R * Plw_Fat_Blood_Benzaldehyde_R + Blood_Water_Fraction_R)+ Fat_Blood_B_R
P_F_Cinnamaldehyde_R <-  (Adipose_Lipid_Fraction_R * Plw_Fat_Blood_Cinnamaldehyde_R + Adipose_Water_Fraction_R)/(Blood_Lipid_Fraction_R * Plw_Fat_Blood_Cinnamaldehyde_R + Blood_Water_Fraction_R) + Fat_Blood_B_R
P_F_Cinnamylalcohol_R <- (Adipose_Lipid_Fraction_R * Plw_Fat_Blood_Cinnamylalcohol_R + Adipose_Water_Fraction_R)/(Blood_Lipid_Fraction_R * Plw_Fat_Blood_Cinnamylalcohol_R + Blood_Water_Fraction_R) + Fat_Blood_B_R


#----------------Liver/Blood 
#Human
P_L_Benzaldehyde   <-  (Liver_Lipid_Fraction * Plw_Liver_Blood_Benzaldehyde + Liver_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Liver_Blood_Benzaldehyde + Blood_Water_Fraction)+ Liver_Blood_B
P_L_Cinnamaldehyde <-  (Liver_Lipid_Fraction * Plw_Liver_Blood_Cinnamaldehyde + Liver_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Liver_Blood_Cinnamaldehyde + Blood_Water_Fraction) + Liver_Blood_B
P_L_Cinnamylalcohol <-  (Liver_Lipid_Fraction * Plw_Liver_Blood_Cinnamylalcohol + Liver_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Liver_Blood_Cinnamylalcohol + Blood_Water_Fraction) + Liver_Blood_B

#Rat
P_L_Benzaldehyde_R   <- (Liver_Lipid_Fraction_R * Plw_Liver_Blood_Benzaldehyde_R + Liver_Water_Fraction_R)/(Blood_Lipid_Fraction_R * Plw_Liver_Blood_Benzaldehyde_R + Blood_Water_Fraction_R)+ Liver_Blood_B_R
P_L_Cinnamaldehyde_R <- (Liver_Lipid_Fraction_R * Plw_Liver_Blood_Cinnamaldehyde_R + Liver_Water_Fraction_R)/(Blood_Lipid_Fraction_R * Plw_Liver_Blood_Cinnamaldehyde_R + Blood_Water_Fraction_R) + Liver_Blood_B_R
P_L_Cinnamylalcohol_R <- (Liver_Lipid_Fraction_R * Plw_Liver_Blood_Cinnamylalcohol_R + Liver_Water_Fraction_R)/(Blood_Lipid_Fraction_R * Plw_Liver_Blood_Cinnamylalcohol_R + Blood_Water_Fraction_R) + Liver_Blood_B_R

#----------------Muscle/Blood
#Human
P_SP_Benzaldehyde   <-  (Muscle_Lipid_Fraction * Plw_Muscle_Blood_Benzaldehyde + Muscle_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Muscle_Blood_Benzaldehyde + Blood_Water_Fraction)+ Muscle_Blood_B
P_SP_Cinnamaldehyde <-  (Muscle_Lipid_Fraction * Plw_Muscle_Blood_Cinnamaldehyde + Muscle_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Muscle_Blood_Cinnamaldehyde + Blood_Water_Fraction) + Muscle_Blood_B
P_SP_Cinnamylalcohol <-  (Muscle_Lipid_Fraction * Plw_Muscle_Blood_Cinnamylalcohol + Muscle_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Muscle_Blood_Cinnamylalcohol + Blood_Water_Fraction) + Muscle_Blood_B

#Rat
P_SP_Benzaldehyde_R   <-  (Muscle_Lipid_Fraction_R * Plw_Muscle_Blood_Benzaldehyde_R + Muscle_Water_Fraction_R)/(Blood_Lipid_Fraction_R * Plw_Muscle_Blood_Benzaldehyde_R + Blood_Water_Fraction_R)+ Muscle_Blood_B_R
P_SP_Cinnamaldehyde_R <-  (Muscle_Lipid_Fraction_R * Plw_Muscle_Blood_Cinnamaldehyde_R + Muscle_Water_Fraction_R)/(Blood_Lipid_Fraction_R * Plw_Muscle_Blood_Cinnamaldehyde_R + Blood_Water_Fraction_R) + Muscle_Blood_B_R
P_SP_Cinnamylalcohol_R <- (Muscle_Lipid_Fraction_R * Plw_Muscle_Blood_Cinnamylalcohol_R + Muscle_Water_Fraction_R)/(Blood_Lipid_Fraction_R * Plw_Muscle_Blood_Cinnamylalcohol_R + Blood_Water_Fraction_R) + Muscle_Blood_B_R

#-------------------Brain/Blood
#Human
P_Br_Benzaldehyde   <-  (Brain_Lipid_Fraction * Plw_Brain_Blood_Benzaldehyde + Brain_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Brain_Blood_Benzaldehyde + Blood_Water_Fraction)+ Brain_Blood_B
P_Br_Cinnamaldehyde <-  (Brain_Lipid_Fraction * Plw_Brain_Blood_Cinnamaldehyde + Brain_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Brain_Blood_Cinnamaldehyde + Blood_Water_Fraction) + Brain_Blood_B
P_Br_Cinnamylalcohol <- (Brain_Lipid_Fraction * Plw_Brain_Blood_Cinnamylalcohol + Brain_Water_Fraction)/(Blood_Lipid_Fraction * Plw_Brain_Blood_Cinnamylalcohol + Blood_Water_Fraction) + Brain_Blood_B


#Blood:air partition coefficient calculations based on ten Berge et al 2011
#Human

#Henry coefficient calculations
Cinnamaldehyde_vapour_pressure  <- 0.0337 #mm Hg
Cinnamaldehyde_Molecular_weight <- 132.16 #unitless 
Cinnamaldehyde_water_solubility <- 2.15E3 #mg\l
Cinnamaldehyde_gas_constant     <- 3.45E-6 #atm-m3/mole
temperature                      <- 298.15 #kelvin


Cinnamaldehyde_Henry_DL<- (Cinnamaldehyde_vapour_pressure * Cinnamaldehyde_Molecular_weight)/(Cinnamaldehyde_water_solubility * Cinnamaldehyde_gas_constant * temperature)



P_B_Cinnamaldehyde <- 0.4445/Cinnamaldehyde_Henry_DL + 0.005189 * Koa_Cinnamaldehyde   
P_B_Benzaldehdye <- 0.4445/HL_Benzaldehyde + 0.005189 * Koa_Benzaldehyde

#Rat
P_B_Cinnamaldehyde_R <- 0.4445/Cinnamaldehyde_Henry_DL + 0.005189 * Koa_Cinnamaldehyde
P_B_Benzaldehdye_R <- 0.4445/HL_Benzaldehyde + 0.005189 * Koa_Benzaldehyde



#Calculation of Ka and Fa based on the QSAR in Ans Punt
#Just change h polar surface area
#rat
R<-0.125 #radius small intestine 0.252/2
Tsi<-1.63 #small intestine transit time 
kt<-1/(Tsi/7)

#Polar surface area
T_PSA=17.07

Papp.cm.power.minus.six.per.s = 10^(-4.36 - 0.01*T_PSA)*1000000
Peff.cm.power.minus.four.per.s.human = (10^(0.4926*log10(Papp.cm.power.minus.six.per.s)- 0.1454))
Peff.cm.power.minus.four.per.s = Peff.cm.power.minus.four.per.s.human/11.04
Peff.cm.per.hr = Peff.cm.power.minus.four.per.s/10000*3600
Ka= (Peff.cm.per.hr*2)/R
Fa = 1-(1+Ka/kt)**-7

#Calculation of Ka and Fa based on the QSAR in Ans Punt
#Just change h polar sufrace area
Human
R<-1.25 #radius small intestine=2.52·π≈0.39789cm
Tsi<-1.63 #small intestine transit time 
kt<-1/(Tsi/7)

#Polar surface area
T_PSA=17.07

Papp.cm.power.minus.six.per.s = 10^(-4.36 - 0.01*T_PSA)*1000000
Peff.cm.power.minus.four.per.s.human = (10^(0.4926*log10(Papp.cm.power.minus.six.per.s)- 0.1454))
Peff.cm.power.minus.four.per.s = Peff.cm.power.minus.four.per.s.human/11.04
Peff.cm.per.hr = Peff.cm.power.minus.four.per.s/10000*3600
Ka= (Peff.cm.per.hr*2)/R
Fa = 1-(1+Ka/kt)**-7
