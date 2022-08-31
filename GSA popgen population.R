#author: Joris Jean van der Lugt
#date: 30-08-2022
#GSA popgen population
library(tidyverse)
library(readxl)
library(readr)
library(truncnorm)
library(reshape2)
library(sensitivity)
library(PKNCA)

#Simulations
set.seed(15204)         #to ensure a reproducible output
amount.units            <-"umol"
time.units              <-"h"
nbr.doses               <-1        #number of doses
time.0                  <-0        #time start dosing
time.end                <-8        #time end of simulation
time.frame              <-0.1     #time steps of simulation
N                       <-4000     #Number of males
Volume_exposure_chamber <-10       #volume exposure chamber in L
MW                      <-132.16   #The molecular weight of Cinnamaldehyde



#Importing the popgen data set 
#GSA_Base <- read_csv("Popgen GSA N=2000.csv", 
 #                   col_types = cols(Age = col_number(), 
  #                                    `Body Mass` = col_number(), Height = col_number(), 
   #                                   `Cardiac Output` = col_number(), 
    #                                  `Lung mass` = col_number(), `Lung flow` = col_number(), 
     #                                 `Liver mass` = col_number(), `Liver flow` = col_number(), 
      #                                `Small intestine mass` = col_number(), 
       #                               `Small intestine flow` = col_number(), 
        #                              `Adipose mass` = col_number(), `Adipose flow` = col_number(), 
         #                             `Liver Total flow` = col_number(), 
          #                            `Slowly Perfused mass` = col_number(), 
           #                           `Slowly Perfused flow` = col_number(), 
            #                          `Richly Perfused mass` = col_number(), 
             #                         `Richly Perfused flow` = col_number(), 
              #                        `Lung Bronchial flow` = col_number()))

#Importing the popgen data set 
GSA_Base <- GSA_popgen_n_4000 <- read_csv("GSA popgen n=4000.csv", 
                     col_types = cols(Age = col_number(), 
                                      `Body Mass` = col_number(), Height = col_number(), 
                                      `Cardiac Output` = col_number(), 
                                      `Lung mass` = col_number(), `Lung flow` = col_number(), 
                                      `Liver mass` = col_number(), `Liver flow` = col_number(), 
                                      `Small intestine mass` = col_number(), 
                                       `Small intestine flow` = col_number(), 
                                       `Adipose mass` = col_number(), `Adipose flow` = col_number(), 
                                      `Liver Total flow` = col_number(), 
                                      `Slowly Perfused mass` = col_number(), 
                                      `Slowly Perfused flow` = col_number(), 
                                      `Richly Perfused mass` = col_number(), 
                                      `Richly Perfused flow` = col_number(), 
                                     `Lung Bronchial flow` = col_number()))



colnames <-c("Age","Height","BW","V_L","V_F",
             "V_B","V_A","V_V","V_SI","V_Pu","V_RP","V_SP","Q_C","Q_F","Q_L","Q_SI","Q_RP","Q_SP")

GSA_pop_number <- length(colnames)
#create data frames for population males
GSA_pop <- matrix(NA, nrow = N, ncol =GSA_pop_number)
colnames(GSA_pop) <- colnames
GSA_pop <- as.data.frame(GSA_pop)


#--Physico-chemical parameters--#
#-Cinnamaldehyde-#
GSA_pop$P_F      <-  1.62    #Fat/Blood partition coefficient
GSA_pop$P_L      <-  0.59    #Fat/Blood partition coefficient
GSA_pop$P_SI     <-  0.59    #Small intestine/Blood partition coefficients
GSA_pop$P_RP     <-  0.59    #Richly perfused tissues/Blood partition coefficients
GSA_pop$P_B      <-  1.25E5  #Blood/Air Partition Coefficient 
GSA_pop$P_SP     <-  0.78    #Slowly perfused tissues/Blood partition coefficients
GSA_pop$P_Pu     <-  0.59    #Lung/Blood partition coefficients

#-Cinnamyl Alcohol-#
GSA_pop$P_OH_F    <-  1.64   #Fat/Blood partition coefficient
GSA_pop$P_OH_L    <-  0.59   #Fat/Blood partition coefficient
GSA_pop$P_OH_SI   <-  0.59   #Small intestine/Blood partition coefficients
GSA_pop$P_OH_RP   <-  0.59   #Richly perfused tissues/Blood partition coefficients
GSA_pop$P_OH_SP   <-  0.78   #Slowly perfused tissues/Blood partition coefficients
GSA_pop$P_OH_Pu   <-  0.59   #Lung/Blood partition coefficients

#--Physiological Parameters--#

#Population specific parameters (Male)
GSA_pop$Age              <- GSA_Base$Age                            #Age (years)           
GSA_pop$Height           <- GSA_Base$Height                         #Body height (cm)                                     
GSA_pop$BW               <- GSA_Base$'Body Mass'                    #Body weight (kg)

#-Tissues volumes in % body weight-#

GSA_pop$V_L      <-GSA_Base$`Liver mass`              #Volume liver tissue (l)
GSA_pop$V_F      <-GSA_Base$`Adipose mass`            #Volume adipose tissue (L)
GSA_pop$V_B      <-(((13.1 * GSA_pop$Height + 18.05 * GSA_pop$BW - 480) / 0.5723) / 1000)#Volume blood (L)
GSA_pop$V_A      <-GSA_pop$V_B / 3                  #Volume arterial blood (L)
GSA_pop$V_V      <-GSA_pop$V_B * (2/3)              #Volume venous blood (L) 
GSA_pop$V_SI     <-GSA_Base$`Small intestine mass`    #Volume Small intestine (L)
GSA_pop$V_Pu     <-GSA_Base$`Lung mass`
GSA_pop$V_RP     <-GSA_Base$`Richly Perfused mass`    #Volume richly perfused tissue (L)
GSA_pop$V_SP     <-GSA_Base$`Slowly Perfused mass`    #Volume slowly perfused tissue (L)

#-Cardiac parameters-#
GSA_pop$Q_C           <- GSA_Base$`Cardiac Output`           #Cardiac output (L/h)
GSA_pop$Q_SI          <- GSA_Base$`Small intestine flow`                                         #Blood flow to the gut (L/h)
GSA_pop$Q_F           <- GSA_Base$`Adipose flow`                                         #Blood flow to adipose tissue (L/h)
GSA_pop$Q_L           <- GSA_Base$`Liver flow`                                        #Blood flow to liver(L/h)
GSA_pop$Q_Pu          <- GSA_Base$`Lung flow`
GSA_pop$Q_RP          <- GSA_Base$`Richly Perfused flow`             #Blood flow to richly perfused tissue (L/h)
GSA_pop$Q_SP          <- GSA_Base$`Slowly Perfused flow`                            #Blood flow to slowly perfused tissue (L/h)
GSA_pop$check         <- GSA_pop$Q_SI + GSA_pop$Q_F + GSA_pop$Q_L + GSA_pop$Q_RP + GSA_pop$Q_SP

#Pulmonary ventilation
GSA_pop$P_V           <-rnorm(N,mean=540,sd=3)

#----GSH parameters----#
#--GSH synthesis in umol/kg tissue/h--#

GSA_pop$G_SYN_L     <- 1122  #Liver 
GSA_pop$G_SYN_SI    <- 27    #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
GSA_pop$k_L_GLOS    <- 0.142 #Liver
GSA_pop$k_SI_GLOS   <- 0.044 #Small intestine

#--Initial GSH concentration--#
GSA_pop$init_GSH_L  <- 5639 * GSA_pop$V_L  #initial GSH concentration in the liver in umol
GSA_pop$init_GSH_SI <- 1250 * GSA_pop$V_SI  #initial GSH concentration in the small intestine in umol

GSA_pop$k_GSH <- 6.6 * 10^(-4) #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h
GSA_pop$k_DNA <- 1.6 * 10^(-8) #The second-order rate constant of the reaction between cinnamaldehyde and 2ʹ-dG in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
GSA_pop$C_PRO_L     <- 5319  #Liver
GSA_pop$C_PRO_SI    <- 245   #Small intestine

#----DNA parameters----#
GSA_pop$C_L_dG     <-  1.36 #Concentration of 2ʹ-dG in the liver μmol/kg liver
GSA_pop$T_0.5      <-  38.5   #Half-life of DNA adduct in the liver in hours


#--Chemical parameters--#
GSA_pop$Ka <- 5.0  #Absorption rate constant for uptake in the Small intestine in per H

#----Liver----#

#-first rate order constants-#
GSA_pop$k_L_OH  <- 4.2e-02   #Scaled first rate order constant for the enzymatic oxidation of cinnamyl alcohol in the liver in umol/h

#--Michaelis menten constants--#
GSA_pop$Km_L_CA     <-  8.5  #Km for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the liver in μM
GSA_pop$Km_L_AO     <-  330  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the liver in μM
GSA_pop$Km_L_GST    <-  100  #Km for enzymatic conjugation of cinnamaldehyde with GST in the liver in μM  
GSA_pop$Km_L_GST_G  <-  1.7*10^3  #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the small intestine μM

#--Vmax values--#
GSA_pop$Vsmax_L_CA    <-  9.7  #Scaled Vmax for enzymatic oxidation of cinnamaldehyde in the liver in μmol/h 
GSA_pop$Vsmax_L_AO    <-  73   #Scaled Vmax for enzymatic reduction of cinnamaldehyde in the liver in μmol/h
GSA_pop$Vsmax_L_GST   <-  37   #Scaled Vmax for enzymatic conjugation of cinnamaldehyde with GSH in the liver in μmol/h

#----Small intestines----#
#--Michaelis menten constants--#
GSA_pop$Km_SI_CA    <- 70  #Km for enzymatic oxidation of cinnamaldehyde into cinnamic acid in the Small Intestine in μM
GSA_pop$Km_SI_AO    <- 90  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the Small Intestine in μM
GSA_pop$Km_SI_OH    <- 290 #Km for enzymatic oxidation of cinnamly alcOHol into cinnamaldehyde in the Small Intestine in μM
GSA_pop$Km_SI_GST   <- 600 #Km for enzymatic conjugation of cinnamaldehye with GST in the Small Intestine in μM RAT value
GSA_pop$Km_SI_GST_G <- 100 #Km toward cinnamaldehyde for enzymatic conjugation of cinnamaldehyde in the small intestine μM

#-Vmax values-#
GSA_pop$Vsmax_SI_CA    <- 21 #Scaled Vmax for enzymatic oxidation of cinnamaldehyde into Cinnamic acid in the Small Intestine in μmol/h 
GSA_pop$Vsmax_SI_AO    <- 30 #Scaled Vmax for enzymatic reduction of cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
GSA_pop$Vsmax_SI_OH    <- 5.0 #Scaled Vmax for enzymatic Oxidation of cinnamyl alcohol into cinnamaldehyde in the Small Intestine in μmol/h 
GSA_pop$Vsmax_SI_GST   <- 63 #Scaled Vmax for enzymatic Conjugation of cinnamaldehyde with GSH in the in the small intestine in μmol/h RAT value
