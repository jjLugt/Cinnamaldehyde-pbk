#author: Joris Jean van der Lugt
#date: 20-05-2021
#Global sensitivity analysis
library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(truncnorm)
library(reshape2)
library(sensitivity)
library(PKNCA)
library(ggplot2)
library(gridExtra)

#Simulations
set.seed(15204)                       #to ensure a reproducible output if random input is used
MW                         <-132.16   #The molecular weight of Cinnamaldehyde
BW                         <-0.25    #Body weight in Kg
Volume_exposure_chamber    <-10       #volume exposure chamber in L


#--Physio-chemical parameters--#
#-Cinnamaldehyde-#

P_F      <-  17.42#Fat/Blood partition coefficient
P_L      <-  1.18 #Fat/Blood partition coefficient 
P_SI     <-  1.18 #Small intestine/Blood partition coefficients
P_RP     <-  1.18 #Richly perfused tissues/Blood partition coefficients
P_SP     <-  0.39 #Slowly perfused tissues/Blood partition coefficients
P_B      <-  0.29 #Blood air partition coefficients
P_Pu     <-  1.18 #Lung/Blood partition coefficients

#-Cinnamyl Alcohol-#
P_OH_F    <- 17.65 #Fat/Blood partition coefficient
P_OH_L    <-  1.18 #Fat/Blood partition coefficient
P_OH_SI   <-  1.18 #Small intestine/Blood partition coefficients
P_OH_RP   <-  1.18 #Richly perfused tissues/Blood partition coefficients
P_OH_SP   <-  0.49 #Slowly perfused tissues/Blood partition coefficients
P_OH_Pu   <-  1.18 #Lung/Blood partition coefficients


#-Apparent first order rate constant GSH turn over(RAT?) per h-#
k_L_GLOS    <- 0.142 #Liver
k_SI_GLOS   <- 0.044 #Small intestine

#--Physiological Parameters--#

#-Unitless fraction of body weight consisting of a specific compartment-#

fV_F      <- 0.07   #Fat
fV_L      <- 0.034 #Liver
fV_SI     <- 0.014 #Small intestine
fV_A      <- 0.019 #Arterial Blood
fV_V      <- 0.059 #Venous Blood
fV_RP     <- 0.037 #Richly perfused 
fV_SP     <- 0.676  #Slowly perfused
fV_Pu     <- 0.005  #Lung


V_adjust <- fV_F + fV_L + fV_SI + fV_A + fV_V + fV_RP + fV_SP +fV_Pu


#-Tissues volumes in L #
V_F      <- (fV_F / V_adjust)* BW    #Fat
V_L      <- (fV_L / V_adjust)* BW  #Liver
V_SI     <- (fV_SI / V_adjust)* BW     #Small intestine
V_A      <- (fV_A / V_adjust)* BW  #Arterial Blood
V_V      <- (fV_V / V_adjust)* BW  #Venous Blood
V_RP     <- (fV_RP / V_adjust)* BW  #Richly perfused (RP)
V_SP     <- (fV_SP / V_adjust)* BW   #Slowly perfused (SP)
V_Pu     <- (fV_Pu / V_adjust)* BW    #Lung

#-Cardiac parameters-#

Q_C      <- 5.4    #Cardiac output in L/h

#-Unitless fraction of cardiac output-#

fQ_F      <- 0.07  #Fat
fQ_L      <- 0.13  #Liver
fQ_SI     <- 0.12  #Small intestine
fQ_RP     <- 0.64  #Richly perfused
fQ_SP     <- 0.17  #Slowly perfused


Q_adjust <- fQ_F + fQ_L + fQ_SI + fQ_RP + fQ_SP

#-Blood flow to tissues in % cardiac output-#
Q_F      <- (fQ_F / Q_adjust)* Q_C    #Fat
Q_L      <- (fQ_L / Q_adjust)* Q_C  #Liver
Q_SI     <- (fQ_SI / Q_adjust)* Q_C     #Small intestine
Q_RP     <- (fQ_RP / Q_adjust)* Q_C  #Richly perfused (RP)
Q_SP     <- (fQ_SP / Q_adjust)* Q_C   #Slowly perfused (SP)
Q_Pu     <- Q_C   #Blood flow through the lung

P_V      <- 0.75 #Pulmonary ventilation in L/h 

#--Michaelis menten constants--#

Km_L_AO     <-  120  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the liver in μM
Km_L_GST_G  <-  100  #Km for enzymatic conjugation of cinnamaldehyde with GST in the liver in μM 
Km_L_OH     <-  1300 #Km for enzymatic oxidation of cinnamyl alcohol to cinnamaldehyde in the liver in μM

#--Michaelis menten constants SI--#
Km_SI_CA    <- 0   #Km for enzymatic oxidation of cinnamaldehyde into cinnamic acid in the Small Intestine in μM
Km_SI_AO    <- 75  #Km for enzymatic reduction of cinnamaldehyde into cinnamyl alcOHol in the Small Intestine in μM
Km_SI_GST   <- 600 #Km for enzymatic conjugation of cinnamaldehye with GST in the Small Intestine in μM (RAT value)
Km_SI_GST_G <- 100 #Km toward GSH for enzymatic conjugation of cinnamaldehyde in the small intestine (μM)

k_GSH <- 6.6 * 10^(-4) #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h

Ka <- 0.2  #Absorption rate constant for uptake in the Small intestine in per H

#Al parameters that have a body weight component can only be defined after a distribution of body weight values has been made. parameters that consist of a body weight term
#multiplied by a factor should also have a distribution. this is defined below. 
#example 
#1<- a * x(tissue volume)

#1 can only be defined after x has been defined. the factor a should also have its one distribution to see the impact of this term. so this is seperately done


G_SYN_L_Factor     <- 869
G_SYN_SI_Factor   <- 78
init_GSH_L_Factor   <- 6120
init_GSH_SI_Factor <- 1780

C_PRO_L_Factor     <- 5319
C_PRO_SI_Factor   <- 245


S9_scaling_L_Factor <- 143

k_L_CA_Factor   <-  7.4*10^(-3)
k_L_GST_Factor  <-  6.2*10^(-2)

Vsmax_L_AO_Factor    <- 29
Vsmax_L_GST_G_Factor <- 100
Vsmax_L_OH_Factor    <- 15
S9_scaling_SI_Factor <- 11.4
k_SI_CA_Factor  <- 3.9*10^(-3)
Vsmax_SI_AO_Factor    <- 5.8
Vsmax_SI_GST_Factor   <- 63

#Collection of all parameters so a distribution can be made with them later on in the file #Human
dist_para <- cbind(P_F,
                   P_L,
                   P_SI,
                   P_RP,
                   P_SP,
                   P_B,
                   P_Pu,
                   P_OH_F,
                   P_OH_L,
                   P_OH_SI,
                   P_OH_RP,
                   P_OH_SP,
                   P_OH_Pu,
                   k_L_GLOS,
                   k_SI_GLOS,
                   BW,
                   fV_F,   
                   fV_L, 
                   fV_SI,    
                   fV_A,    
                   fV_V,    
                   fV_RP,    
                   fV_SP,    
                   fV_Pu,   
                   V_adjust,
                   V_F,
                   V_L,
                   V_SI,
                   V_A,
                   V_V,
                   V_RP,
                   V_SP,
                   V_Pu,
                   Q_C,
                   fQ_F,   
                   fQ_L,
                   fQ_SI,    
                   fQ_RP,    
                   fQ_SP,    
                   Q_adjust,
                   Q_Pu,
                   Q_F,
                   Q_L,
                   Q_SI,
                   Q_RP,
                   Q_SP,
                   P_V,
                   Km_L_AO,
                   Km_L_GST_G,
                   Km_L_OH,
                   Km_SI_CA,
                   Km_SI_AO,
                   Km_SI_GST,
                   Km_SI_GST_G,
                   k_GSH,
                   Ka,
                   G_SYN_L_Factor,
                   G_SYN_SI_Factor,  
                   init_GSH_L_Factor,
                   init_GSH_SI_Factor, 
                   C_PRO_L_Factor,       
                   C_PRO_SI_Factor,
                   S9_scaling_L_Factor, 
                   k_L_CA_Factor,
                   k_L_GST_Factor,
                   Vsmax_L_AO_Factor,    
                   Vsmax_L_GST_G_Factor,
                   Vsmax_L_OH_Factor,
                   S9_scaling_SI_Factor,
                   k_SI_CA_Factor,        
                   Vsmax_SI_AO_Factor,       
                   Vsmax_SI_GST_Factor)


#Creating empty vectors for use later
colnames <- c(colnames(dist_para))
par_var <- length(colnames)

#Generating a upper and lower bound for use in the analysis 
Mean <- dist_para

Lower <- Mean - 0.01 * Mean
Upper <- Mean + 0.01 * Mean

#create data frames for population
n_sim  <- 2000             #number of iterations
X1 <- matrix(NA, nrow = n_sim, ncol = par_var)
colnames(X1) <- colnames
X1 <- as.data.frame(X1)
var <- X1

X2 <- matrix(NA, nrow = n_sim, ncol = par_var)
colnames(X2) <- colnames
X2 <- as.data.frame(X2)
var <- X2

#create a uniform distribution of variables between the upper and lower limit
for(i in 1:par_var){
  X1[,i] <- runif(n_sim, min = Lower[,i], max = Upper[,i])
  X2[,i] <- runif(n_sim, min = Lower[,i], max = Upper[,i])
}


for(i in 1:nrow(X1)){
  X1[i,73] <-X1[i,27]*X1[i,57]
  colnames(X1)[73] <- "G_SYN_L"
  X1[i,74] <-X1[i,28]*X1[i,58]
  colnames(X1)[74] <- "G_SYN_SI"
  X1[i,75] <-X1[i,27]*X1[i,59]
  colnames(X1)[75] <- "init_GSH_L"
  X1[i,76] <-X1[i,28]*X1[i,60]
  colnames(X1)[76] <- "init_GSH_SI"
  X1[i,77] <-X1[i,27]*X1[i,61]
  colnames(X1)[77] <- "C_PRO_L"
  X1[i,78] <-X1[i,28]*X1[i,62]
  colnames(X1)[78] <- "C_PRO_SI"
  X1[i,79] <-(X1[i,63]*1000)*X1[i,27]
  colnames(X1)[79] <- "S9_scaling_L"
  X1[i,80] <-(X1[i,64]*60/1000)*X1[i,79]
  colnames(X1)[80] <- "k_L_CA"
  X1[i,81] <-(X1[i,65]*60/1000)*X1[i,79]
  colnames(X1)[81] <- "k_L_GST"
  X1[i,82] <-(X1[i,66]*60/1000)*X1[i,79]
  colnames(X1)[82] <- "Vsmax_L_AO"
  X1[i,83] <-(X1[i,67]*60/1000)*X1[i,79]
  colnames(X1)[83] <- "Vsmax_L_GST"
  X1[i,84] <-(X1[i,68]*1000)*X1[i,79]
  colnames(X1)[84] <- "Vsmax_L_OH"
  X1[i,85] <-(X1[i,69]*1000)*X1[i,28]
  colnames(X1)[85] <- "S9_scaling_SI"
  X1[i,86] <-(X1[i,70]*60/1000)*X1[i,85]
  colnames(X1)[86] <- "k_SI_CA"
  X1[i,87] <-(X1[i,71]*60/1000)*X1[i,85]
  colnames(X1)[87] <- "Vsmax_SI_AO"
  X1[i,88] <-(X1[i,72]*60/1000)*X1[i,85]
  colnames(X1)[88] <- "Vsmax_SI_GST"
} 

for(i in 1:nrow(X2)){
  X2[i,73] <-X2[i,27]*X2[i,57]
  colnames(X2)[73] <- "G_SYN_L" #GSH synthesis in umol/in the liver/h
  X2[i,74] <-X2[i,28]*X2[i,58]
  colnames(X2)[74] <- "G_SYN_SI"   #GSH synthesis in umol/in the SI/h
  X2[i,75] <-X2[i,27]*X2[i,59]
  colnames(X2)[75] <- "init_GSH_L"  #initial GSH concentration in the liver in umol
  X2[i,76] <-X2[i,28]*X2[i,60]
  colnames(X2)[76] <- "init_GSH_SI" #initial GSH concentration in the small intestine in umol
  X2[i,77] <-X2[i,27]*X2[i,61]
  colnames(X2)[77] <- "C_PRO_L" #Protein reactive sites in μmol/kg LIVER
  X2[i,78] <-X2[i,28]*X2[i,62]
  colnames(X2)[78] <- "C_PRO_SI" #Protein reactive sites in μmol/kg SI
  X2[i,79] <-(X2[i,63]*1000)*X2[i,27]
  colnames(X2)[79] <- "S9_scaling_L"   #scaling factor for S9 fraction per g L
  X2[i,80] <-(X2[i,64]*60/1000)*X2[i,79]
  colnames(X2)[80] <- "k_L_CA"            #Scaled first-order rate constant for enzymatic oxidation of cinnamaldehyde in the liver (L/h)
  X2[i,81] <-(X2[i,65]*60/1000)*X2[i,79]
  colnames(X2)[81] <- "k_L_GST"             #Scaled first-order rate constant for enzymatic conjugation of cinnamaldehyde with GSH in the liver (L/h)
  X2[i,82] <-(X2[i,66]*60/1000)*X2[i,79]
  colnames(X2)[82] <- "Vsmax_L_AO"        #Scaled Vmax for enzymatic reduction of cinnamaldehyde in the liver in μmol/h
  X2[i,83] <-(X2[i,67]*60/1000)*X2[i,79]
  colnames(X2)[83] <- "Vsmax_L_GST"      #Scaled Vmax toward GSH for enzymatic conjugation of cinnamaldehyde in the Liver 
  X2[i,84] <-(X2[i,68]*1000)*X2[i,79]
  colnames(X2)[84] <- "Vsmax_L_OH"       #Scaled Vmax for enzymatic oxidation of cinnamyl alcohol to cinnamaldehyde in the liver in μmol/h
  X2[i,85] <-(X2[i,69]*1000)*X2[i,28]
  colnames(X2)[85] <- "S9_scaling_SI"   #scaling factor fraction S9 protein per g SI
  X2[i,86] <-(X2[i,70]*60/1000)*X2[i,85]
  colnames(X2)[86] <- "k_SI_CA"         #Scaled first-order rate constant for enzymatic oxidation of cinnamaldehyde in the small intestine (μmol/h)
  X2[i,87] <-(X2[i,71]*60/1000)*X2[i,85]
  colnames(X2)[87] <- "Vsmax_SI_AO"      #Scaled Vmax for enzymatic reduction of cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
  X2[i,88] <-(X2[i,72]*60/1000)*X2[i,85]
  colnames(X2)[88] <- "Vsmax_SI_GST"     #Scaled Vmax for enzymatic Conjugation of cinnamaldehyde with GSH in the in the small intestine in μmol/h (RAT value)
}


#Defining empty vectors to be used in the following calculations #Human
Flow_volumes <-c()
Tissue_volumes<-c()
Total_body_weight<-c()
Total_QC<-c()
#Based on the Changes in the volume or blood fractions recalculating blood flows and tissue volumes to ensure mass balance
for(i in 1:nrow(X1)){
  #Tissue volumes calculations
  X1[i,25]<- sum(X1[i,17:24])
  X1[i,26]<- (X1[i,17]/X1[i,25])* X1[i,16]
  X1[i,27]<- (X1[i,18]/X1[i,25])* X1[i,16]
  X1[i,28]<- (X1[i,19]/X1[i,25])* X1[i,16]
  X1[i,29]<- (X1[i,20]/X1[i,25])* X1[i,16]
  X1[i,30]<- (X1[i,21]/X1[i,25])* X1[i,16]
  X1[i,31]<- (X1[i,22]/X1[i,25])* X1[i,16]
  X1[i,32]<- (X1[i,23]/X1[i,25])* X1[i,16]
  X1[i,33]<- (X1[i,24]/X1[i,25])* X1[i,16]
  #Blood flow calculations
  X1[i,40]<- sum(X1[i,35:39])
  X1[i,41]<- X1[i,34]
  X1[i,42]<- (X1[i,35]/X1[i,40])* X1[i,34]
  X1[i,43]<- (X1[i,36]/X1[i,40])* X1[i,34]
  X1[i,44]<- (X1[i,37]/X1[i,40])* X1[i,34]
  X1[i,45]<- (X1[i,38]/X1[i,40])* X1[i,34]
  X1[i,46]<- (X1[i,39]/X1[i,40])* X1[i,34]
  Flow_volumes <-append(Flow_volumes,sum(X1[i,42:46]))
  Tissue_volumes<-append(Tissue_volumes ,sum(X1[i,26:33]))
  Total_body_weight <- append(Total_body_weight ,X1[i,16])
  Total_QC<- append(Total_QC ,X1[i,34])
  names<- c("total flow","QC","Total volume","BW")
  Mass_check <- as.data.frame(cbind(Flow_volumes,Total_QC,Tissue_volumes,Total_body_weight), col.names = names)
}
#Calculations for X2
for(i in 1:nrow(X2)){
  #Tissue volumes calculations
  X2[i,25]<- sum(X2[i,17:24])
  X2[i,26]<- (X2[i,17]/X2[i,25])* X2[i,16]
  X2[i,27]<- (X2[i,18]/X2[i,25])* X2[i,16]
  X2[i,28]<- (X2[i,19]/X2[i,25])* X2[i,16]
  X2[i,29]<- (X2[i,20]/X2[i,25])* X2[i,16]
  X2[i,30]<- (X2[i,21]/X2[i,25])* X2[i,16]
  X2[i,31]<- (X2[i,22]/X2[i,25])* X2[i,16]
  X2[i,32]<- (X2[i,23]/X2[i,25])* X2[i,16]
  X2[i,33]<- (X2[i,24]/X2[i,25])* X2[i,16]
  #Blood flow calculations
  X2[i,40]<- sum(X2[i,35:39])
  X2[i,41]<- X2[i,34]
  X2[i,42]<- (X2[i,35]/X2[i,40])* X2[i,34]
  X2[i,43]<- (X2[i,36]/X2[i,40])* X2[i,34]
  X2[i,44]<- (X2[i,37]/X2[i,40])* X2[i,34]
  X2[i,45]<- (X2[i,38]/X2[i,40])* X2[i,34]
  X2[i,46]<- (X2[i,39]/X2[i,40])* X2[i,34]
}


#Removing unnecessary variables 
X1<- X1[-c(17:25,35:40,57:72)]
X2<- X2[-c(17:25,35:40,57:72)]



#the number of bootstrap replicates
n_boot <- 1000



#Sobol design
sa <- soboljansen(model=NULL, X1, X2, nboot = n_boot, conf = 0.95, events = ex)


phys <- sa$X


#Writing the result into a file so that the environment can be cleaned to conserve memory
write.csv(phys,"D:/PBK/Cinnamaldehyde-pbk\\RAT_250mg_oral_phys_0.2ka", row.names = TRUE)



#Loading extracted simulation data. 
solve.pbk.sa <- read.csv("D:/PBK/Cinnamaldehyde-pbk\\SA_rat_250mg_oral_C_Pu_corrected", row.names=1)

#Analysing the generated data set 
solve.pbk.sa=solve.pbk.sa[which(solve.pbk.sa[,"time"]==0.2|solve.pbk.sa[,"time"]==0.5|solve.pbk.sa[,"time"]==1|solve.pbk.sa[,"time"]==1.5| 
                                  solve.pbk.sa[,"time"]==2|solve.pbk.sa[,"time"]==3|solve.pbk.sa[,"time"]==4|
                                  solve.pbk.sa[,"time"]==8),]


tab1=solve.pbk.sa[which(solve.pbk.sa[,"time"]==0.2),]
tab2=solve.pbk.sa[which(solve.pbk.sa[,"time"]==0.5),]
tab3=solve.pbk.sa[which(solve.pbk.sa[,"time"]==1),]
tab4=solve.pbk.sa[which(solve.pbk.sa[,"time"]==1.5),]
tab5=solve.pbk.sa[which(solve.pbk.sa[,"time"]==2),]
tab6=solve.pbk.sa[which(solve.pbk.sa[,"time"]==3),]
tab7=solve.pbk.sa[which(solve.pbk.sa[,"time"]==4),]
tab8=solve.pbk.sa[which(solve.pbk.sa[,"time"]==8),]


SimRes = as.data.frame(matrix(NA,118000,8))
SimRes[,1]=tab1[,2]
SimRes[,2]=tab2[,2]
SimRes[,3]=tab3[,2]
SimRes[,4]=tab4[,2]
SimRes[,5]=tab5[,2]
SimRes[,6]=tab6[,2]
SimRes[,7]=tab7[,2]
SimRes[,8]=tab8[,2]


write.csv(SimRes,"D:/PBK/Cinnamaldehyde-pbk\\SimRes_rat_250mg_oral_C_Pu_corrected", row.names = TRUE)

SimRes <- read.csv("D:/PBK/Cinnamaldehyde-pbk\\SimRes_rat_250mg_oral_C_Pu_corrected", row.names=1)

#Redefining these two variables as these are also used with dist_parm creation but not all of thet variables in dist_parm are used in the SA calculation
#so using them here would create an error.
colnames <- colnames(X1)
par_var <- length(X1)


#plots with error bars
#Sobol analysis plot blood Nrow is the number of paramters in the model
t_A<-(c(0.2,0.5,1,1.5,2,3,4,8))
Main          = TI          = TI.lower           = TI.upper =   Main.lower   = Main.upper   = matrix(NA, nrow = par_var, ncol = length(t_A))  
rownames(Main)= rownames(TI)= rownames(TI.lower) = rownames(TI.upper) =  rownames(Main.lower)=  rownames(Main.upper)= colnames

t_SA <-0.2


for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    Main[,i]       <- sa$S[,1]  #First order indices
    TI[,i]        <- sa$T[,1]  #Total indices
    TI.lower[,i] <- sa$T[,4]   #Lower CL total indices
    TI.upper[,i] <- sa$T[,5]   #Upper CL total indices
    Main.lower[,i]<- sa$S[,4]   #Lower CL total indices
    Main.upper[,i]<-sa$S[,5]    #Lower CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


Main.L <- as.matrix(Main[,1:length(t_A)])       # as.matrix
TI.L  <- as.matrix(TI[,1:length(t_A)])

Main.lower_matrix <- as.matrix(Main.lower[,1:length(t_A)])    
Main.upper_matrix <- as.matrix(Main.upper[,1:length(t_A)])   
TI.lower_matrix <- as.matrix(TI.lower[,1:length(t_A)])    
TI.upper_matrix <- as.matrix(TI.upper[,1:length(t_A)])   


Main.L.t <- apply(Main.L, 1, mean, na.rm=TRUE)       #transform it to a vector
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)
Main.lower_vector<- apply(Main.lower_matrix, 1, mean, na.rm=TRUE)
Main.upper_vector<- apply(Main.upper_matrix, 1, mean, na.rm=TRUE)
TI.lower_vector<- apply(TI.lower_matrix, 1, mean, na.rm=TRUE)
TI.upper_vector<- apply(TI.upper_matrix, 1, mean, na.rm=TRUE)

sorting <- order(TI.L.t, decreasing = F)          #sorting form low to high based on total effect
TI.L.t  <- TI.L.t[sorting]
Main.L.t <- Main.L.t[sorting]
Main.lower_vector <-Main.lower_vector[sorting]
Main.upper_vector<-Main.upper_vector[sorting]
TI.lower_vector <-TI.lower_vector[sorting]
TI.upper_vector<-TI.upper_vector[sorting]

Main.L.t <- ifelse(Main.L.t <= 0, 0, Main.L.t)
tempC    <- t(cbind(Main.L.t, TI.L.t,Main.lower_vector,Main.upper_vector,TI.lower_vector,TI.upper_vector))
tempC<- as.data.frame(tempC[,c(47:57)]) #top 10 variables with the highest total effect 



Global_sa_top_ten<-as.data.frame(melt(tempC[1,]))
Global_sa_top_ten[,3:4]<-as.data.frame(melt(tempC[2,]))
Global_sa_top_ten <- Global_sa_top_ten[ ,-c(3) ]
colnames(Global_sa_top_ten)<-c("variable","main","total")
Global_sa_top_ten<-melt(Global_sa_top_ten, id="variable")
colnames(Global_sa_top_ten)<-c("variable","indices","value")


#organising the upper and lower Ci for the First order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_main_lower<-as.data.frame(melt(tempC[3,]))
Global_sa_main_lower <- Global_sa_main_lower[ ,-c(3) ]
colnames(Global_sa_main_lower)<-c("variable","lower")
Global_sa_main_lower<-melt(Global_sa_main_lower, id="variable")
colnames(Global_sa_main_lower)<-c("variable","indices","lower")

Global_sa_main_upper<-as.data.frame(melt(tempC[4,]))
Global_sa_main_upper <- Global_sa_main_upper[ ,-c(3) ]
colnames(Global_sa_main_upper)<-c("variable","upper")
Global_sa_main_upper<-melt(Global_sa_main_upper, id="variable")
colnames(Global_sa_main_upper)<-c("variable","indices","upper")

Global_sa_main<-cbind(Global_sa_main_upper,Global_sa_main_lower)
Global_sa_main <- Global_sa_main[ ,-c(1,2,4,5) ]

#organising the upper and lower Ci for the total order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_total_lower<-as.data.frame(melt(tempC[5,]))
Global_sa_total_lower <- Global_sa_total_lower[ ,-c(3) ]
colnames(Global_sa_total_lower)<-c("variable","lower")
Global_sa_total_lower<-melt(Global_sa_total_lower, id="variable")
colnames(Global_sa_total_lower)<-c("variable","indices","lower")

Global_sa_total_upper<-as.data.frame(melt(tempC[6,]))
Global_sa_total_upper <- Global_sa_total_upper[ ,-c(3) ]
colnames(Global_sa_total_upper)<-c("variable","upper")
Global_sa_total_upper<-melt(Global_sa_total_upper, id="variable")
colnames(Global_sa_total_upper)<-c("variable","indices","upper")

Global_sa_total<-cbind(Global_sa_total_upper,Global_sa_total_lower)
Global_sa_total <- Global_sa_total[ ,-c(1,2,4,5) ]

Global_sa_CI<-rbind(Global_sa_main,Global_sa_total)

Global_sa_top_ten<-cbind(Global_sa_top_ten,Global_sa_CI)

p_0.2<-ggplot(Global_sa_top_ten,(aes(fill=indices,x=variable, y=value )))+
  geom_bar(position="dodge", stat="identity")+
  ylim(-0.1,1)+
  theme_classic()+
  scale_fill_brewer(palette="Set1")+
  coord_flip()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))+
  ggtitle('12 mins')  



t_SA <-0.5



for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    Main[,i]       <- sa$S[,1]  #First order indices
    TI[,i]        <- sa$T[,1]  #Total indices
    TI.lower[,i] <- sa$T[,4]   #Lower CL total indices
    TI.upper[,i] <- sa$T[,5]   #Upper CL total indices
    Main.lower[,i]<- sa$S[,4]   #Lower CL total indices
    Main.upper[,i]<-sa$S[,5]    #Lower CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


Main.L <- as.matrix(Main[,1:length(t_A)])       # as.matrix
TI.L  <- as.matrix(TI[,1:length(t_A)])

Main.lower_matrix <- as.matrix(Main.lower[,1:length(t_A)])    
Main.upper_matrix <- as.matrix(Main.upper[,1:length(t_A)])   
TI.lower_matrix <- as.matrix(TI.lower[,1:length(t_A)])    
TI.upper_matrix <- as.matrix(TI.upper[,1:length(t_A)])   


Main.L.t <- apply(Main.L, 1, mean, na.rm=TRUE)       #transform it to a vector
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)
Main.lower_vector<- apply(Main.lower_matrix, 1, mean, na.rm=TRUE)
Main.upper_vector<- apply(Main.upper_matrix, 1, mean, na.rm=TRUE)
TI.lower_vector<- apply(TI.lower_matrix, 1, mean, na.rm=TRUE)
TI.upper_vector<- apply(TI.upper_matrix, 1, mean, na.rm=TRUE)

sorting <- order(TI.L.t, decreasing = F)          #sorting form low to high based on total effect
TI.L.t  <- TI.L.t[sorting]
Main.L.t <- Main.L.t[sorting]
Main.lower_vector <-Main.lower_vector[sorting]
Main.upper_vector<-Main.upper_vector[sorting]
TI.lower_vector <-TI.lower_vector[sorting]
TI.upper_vector<-TI.upper_vector[sorting]

Main.L.t <- ifelse(Main.L.t <= 0, 0, Main.L.t)
tempC    <- t(cbind(Main.L.t, TI.L.t,Main.lower_vector,Main.upper_vector,TI.lower_vector,TI.upper_vector))
tempC<- as.data.frame(tempC[,c(47:57)]) #top 10 variables with the highest total effect 



Global_sa_top_ten<-as.data.frame(melt(tempC[1,]))
Global_sa_top_ten[,3:4]<-as.data.frame(melt(tempC[2,]))
Global_sa_top_ten <- Global_sa_top_ten[ ,-c(3) ]
colnames(Global_sa_top_ten)<-c("variable","main","total")
Global_sa_top_ten<-melt(Global_sa_top_ten, id="variable")
colnames(Global_sa_top_ten)<-c("variable","indices","value")


#organising the upper and lower Ci for the First order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_main_lower<-as.data.frame(melt(tempC[3,]))
Global_sa_main_lower <- Global_sa_main_lower[ ,-c(3) ]
colnames(Global_sa_main_lower)<-c("variable","lower")
Global_sa_main_lower<-melt(Global_sa_main_lower, id="variable")
colnames(Global_sa_main_lower)<-c("variable","indices","lower")

Global_sa_main_upper<-as.data.frame(melt(tempC[4,]))
Global_sa_main_upper <- Global_sa_main_upper[ ,-c(3) ]
colnames(Global_sa_main_upper)<-c("variable","upper")
Global_sa_main_upper<-melt(Global_sa_main_upper, id="variable")
colnames(Global_sa_main_upper)<-c("variable","indices","upper")

Global_sa_main<-cbind(Global_sa_main_upper,Global_sa_main_lower)
Global_sa_main <- Global_sa_main[ ,-c(1,2,4,5) ]

#organising the upper and lower Ci for the total order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_total_lower<-as.data.frame(melt(tempC[5,]))
Global_sa_total_lower <- Global_sa_total_lower[ ,-c(3) ]
colnames(Global_sa_total_lower)<-c("variable","lower")
Global_sa_total_lower<-melt(Global_sa_total_lower, id="variable")
colnames(Global_sa_total_lower)<-c("variable","indices","lower")

Global_sa_total_upper<-as.data.frame(melt(tempC[6,]))
Global_sa_total_upper <- Global_sa_total_upper[ ,-c(3) ]
colnames(Global_sa_total_upper)<-c("variable","upper")
Global_sa_total_upper<-melt(Global_sa_total_upper, id="variable")
colnames(Global_sa_total_upper)<-c("variable","indices","upper")

Global_sa_total<-cbind(Global_sa_total_upper,Global_sa_total_lower)
Global_sa_total <- Global_sa_total[ ,-c(1,2,4,5) ]

Global_sa_CI<-rbind(Global_sa_main,Global_sa_total)

Global_sa_top_ten<-cbind(Global_sa_top_ten,Global_sa_CI)

p_0.5<-ggplot(Global_sa_top_ten,(aes(fill=indices,x=variable, y=value )))+
  geom_bar(position="dodge", stat="identity")+
  ylim(-0.1,1)+
  theme_classic()+
  scale_fill_brewer(palette="Set1")+
  coord_flip()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  ggtitle('30 mins') 
t_SA <-1



for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    Main[,i]       <- sa$S[,1]  #First order indices
    TI[,i]        <- sa$T[,1]  #Total indices
    TI.lower[,i] <- sa$T[,4]   #Lower CL total indices
    TI.upper[,i] <- sa$T[,5]   #Upper CL total indices
    Main.lower[,i]<- sa$S[,4]   #Lower CL total indices
    Main.upper[,i]<-sa$S[,5]    #Lower CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


Main.L <- as.matrix(Main[,1:length(t_A)])       # as.matrix
TI.L  <- as.matrix(TI[,1:length(t_A)])

Main.lower_matrix <- as.matrix(Main.lower[,1:length(t_A)])    
Main.upper_matrix <- as.matrix(Main.upper[,1:length(t_A)])   
TI.lower_matrix <- as.matrix(TI.lower[,1:length(t_A)])    
TI.upper_matrix <- as.matrix(TI.upper[,1:length(t_A)])   


Main.L.t <- apply(Main.L, 1, mean, na.rm=TRUE)       #transform it to a vector
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)
Main.lower_vector<- apply(Main.lower_matrix, 1, mean, na.rm=TRUE)
Main.upper_vector<- apply(Main.upper_matrix, 1, mean, na.rm=TRUE)
TI.lower_vector<- apply(TI.lower_matrix, 1, mean, na.rm=TRUE)
TI.upper_vector<- apply(TI.upper_matrix, 1, mean, na.rm=TRUE)

sorting <- order(TI.L.t, decreasing = F)          #sorting form low to high based on total effect
TI.L.t  <- TI.L.t[sorting]
Main.L.t <- Main.L.t[sorting]
Main.lower_vector <-Main.lower_vector[sorting]
Main.upper_vector<-Main.upper_vector[sorting]
TI.lower_vector <-TI.lower_vector[sorting]
TI.upper_vector<-TI.upper_vector[sorting]

Main.L.t <- ifelse(Main.L.t <= 0, 0, Main.L.t)
tempC    <- t(cbind(Main.L.t, TI.L.t,Main.lower_vector,Main.upper_vector,TI.lower_vector,TI.upper_vector))
tempC<- as.data.frame(tempC[,c(47:57)]) #top 10 variables with the highest total effect 



Global_sa_top_ten<-as.data.frame(melt(tempC[1,]))
Global_sa_top_ten[,3:4]<-as.data.frame(melt(tempC[2,]))
Global_sa_top_ten <- Global_sa_top_ten[ ,-c(3) ]
colnames(Global_sa_top_ten)<-c("variable","main","total")
Global_sa_top_ten<-melt(Global_sa_top_ten, id="variable")
colnames(Global_sa_top_ten)<-c("variable","indices","value")


#organising the upper and lower Ci for the First order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_main_lower<-as.data.frame(melt(tempC[3,]))
Global_sa_main_lower <- Global_sa_main_lower[ ,-c(3) ]
colnames(Global_sa_main_lower)<-c("variable","lower")
Global_sa_main_lower<-melt(Global_sa_main_lower, id="variable")
colnames(Global_sa_main_lower)<-c("variable","indices","lower")

Global_sa_main_upper<-as.data.frame(melt(tempC[4,]))
Global_sa_main_upper <- Global_sa_main_upper[ ,-c(3) ]
colnames(Global_sa_main_upper)<-c("variable","upper")
Global_sa_main_upper<-melt(Global_sa_main_upper, id="variable")
colnames(Global_sa_main_upper)<-c("variable","indices","upper")

Global_sa_main<-cbind(Global_sa_main_upper,Global_sa_main_lower)
Global_sa_main <- Global_sa_main[ ,-c(1,2,4,5) ]

#organising the upper and lower Ci for the total order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_total_lower<-as.data.frame(melt(tempC[5,]))
Global_sa_total_lower <- Global_sa_total_lower[ ,-c(3) ]
colnames(Global_sa_total_lower)<-c("variable","lower")
Global_sa_total_lower<-melt(Global_sa_total_lower, id="variable")
colnames(Global_sa_total_lower)<-c("variable","indices","lower")

Global_sa_total_upper<-as.data.frame(melt(tempC[6,]))
Global_sa_total_upper <- Global_sa_total_upper[ ,-c(3) ]
colnames(Global_sa_total_upper)<-c("variable","upper")
Global_sa_total_upper<-melt(Global_sa_total_upper, id="variable")
colnames(Global_sa_total_upper)<-c("variable","indices","upper")

Global_sa_total<-cbind(Global_sa_total_upper,Global_sa_total_lower)
Global_sa_total <- Global_sa_total[ ,-c(1,2,4,5) ]

Global_sa_CI<-rbind(Global_sa_main,Global_sa_total)

Global_sa_top_ten<-cbind(Global_sa_top_ten,Global_sa_CI)

p_1<-ggplot(Global_sa_top_ten,(aes(fill=indices,x=variable, y=value )))+
  geom_bar(position="dodge", stat="identity")+
  ylim(-0.1,1)+
  theme_classic()+
  scale_fill_brewer(palette="Set1")+
  coord_flip()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  ggtitle('60 mins') 

t_SA <-1.5



for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    Main[,i]       <- sa$S[,1]  #First order indices
    TI[,i]        <- sa$T[,1]  #Total indices
    TI.lower[,i] <- sa$T[,4]   #Lower CL total indices
    TI.upper[,i] <- sa$T[,5]   #Upper CL total indices
    Main.lower[,i]<- sa$S[,4]   #Lower CL total indices
    Main.upper[,i]<-sa$S[,5]    #Lower CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


Main.L <- as.matrix(Main[,1:length(t_A)])       # as.matrix
TI.L  <- as.matrix(TI[,1:length(t_A)])

Main.lower_matrix <- as.matrix(Main.lower[,1:length(t_A)])    
Main.upper_matrix <- as.matrix(Main.upper[,1:length(t_A)])   
TI.lower_matrix <- as.matrix(TI.lower[,1:length(t_A)])    
TI.upper_matrix <- as.matrix(TI.upper[,1:length(t_A)])   


Main.L.t <- apply(Main.L, 1, mean, na.rm=TRUE)       #transform it to a vector
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)
Main.lower_vector<- apply(Main.lower_matrix, 1, mean, na.rm=TRUE)
Main.upper_vector<- apply(Main.upper_matrix, 1, mean, na.rm=TRUE)
TI.lower_vector<- apply(TI.lower_matrix, 1, mean, na.rm=TRUE)
TI.upper_vector<- apply(TI.upper_matrix, 1, mean, na.rm=TRUE)

sorting <- order(TI.L.t, decreasing = F)          #sorting form low to high based on total effect
TI.L.t  <- TI.L.t[sorting]
Main.L.t <- Main.L.t[sorting]
Main.lower_vector <-Main.lower_vector[sorting]
Main.upper_vector<-Main.upper_vector[sorting]
TI.lower_vector <-TI.lower_vector[sorting]
TI.upper_vector<-TI.upper_vector[sorting]

Main.L.t <- ifelse(Main.L.t <= 0, 0, Main.L.t)
tempC    <- t(cbind(Main.L.t, TI.L.t,Main.lower_vector,Main.upper_vector,TI.lower_vector,TI.upper_vector))
tempC<- as.data.frame(tempC[,c(47:57)]) #top 10 variables with the highest total effect 



Global_sa_top_ten<-as.data.frame(melt(tempC[1,]))
Global_sa_top_ten[,3:4]<-as.data.frame(melt(tempC[2,]))
Global_sa_top_ten <- Global_sa_top_ten[ ,-c(3) ]
colnames(Global_sa_top_ten)<-c("variable","main","total")
Global_sa_top_ten<-melt(Global_sa_top_ten, id="variable")
colnames(Global_sa_top_ten)<-c("variable","indices","value")


#organising the upper and lower Ci for the First order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_main_lower<-as.data.frame(melt(tempC[3,]))
Global_sa_main_lower <- Global_sa_main_lower[ ,-c(3) ]
colnames(Global_sa_main_lower)<-c("variable","lower")
Global_sa_main_lower<-melt(Global_sa_main_lower, id="variable")
colnames(Global_sa_main_lower)<-c("variable","indices","lower")

Global_sa_main_upper<-as.data.frame(melt(tempC[4,]))
Global_sa_main_upper <- Global_sa_main_upper[ ,-c(3) ]
colnames(Global_sa_main_upper)<-c("variable","upper")
Global_sa_main_upper<-melt(Global_sa_main_upper, id="variable")
colnames(Global_sa_main_upper)<-c("variable","indices","upper")

Global_sa_main<-cbind(Global_sa_main_upper,Global_sa_main_lower)
Global_sa_main <- Global_sa_main[ ,-c(1,2,4,5) ]

#organising the upper and lower Ci for the total order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_total_lower<-as.data.frame(melt(tempC[5,]))
Global_sa_total_lower <- Global_sa_total_lower[ ,-c(3) ]
colnames(Global_sa_total_lower)<-c("variable","lower")
Global_sa_total_lower<-melt(Global_sa_total_lower, id="variable")
colnames(Global_sa_total_lower)<-c("variable","indices","lower")

Global_sa_total_upper<-as.data.frame(melt(tempC[6,]))
Global_sa_total_upper <- Global_sa_total_upper[ ,-c(3) ]
colnames(Global_sa_total_upper)<-c("variable","upper")
Global_sa_total_upper<-melt(Global_sa_total_upper, id="variable")
colnames(Global_sa_total_upper)<-c("variable","indices","upper")

Global_sa_total<-cbind(Global_sa_total_upper,Global_sa_total_lower)
Global_sa_total <- Global_sa_total[ ,-c(1,2,4,5) ]

Global_sa_CI<-rbind(Global_sa_main,Global_sa_total)

Global_sa_top_ten<-cbind(Global_sa_top_ten,Global_sa_CI)

p_1.5<-ggplot(Global_sa_top_ten,(aes(fill=indices,x=variable, y=value )))+
  geom_bar(position="dodge", stat="identity")+
  ylim(-0.1,1)+
  theme_classic()+
  scale_fill_brewer(palette="Set1")+
  coord_flip()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))+
  ggtitle('90 mins') 

t_SA <-2



for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    Main[,i]       <- sa$S[,1]  #First order indices
    TI[,i]        <- sa$T[,1]  #Total indices
    TI.lower[,i] <- sa$T[,4]   #Lower CL total indices
    TI.upper[,i] <- sa$T[,5]   #Upper CL total indices
    Main.lower[,i]<- sa$S[,4]   #Lower CL total indices
    Main.upper[,i]<-sa$S[,5]    #Lower CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


Main.L <- as.matrix(Main[,1:length(t_A)])       # as.matrix
TI.L  <- as.matrix(TI[,1:length(t_A)])

Main.lower_matrix <- as.matrix(Main.lower[,1:length(t_A)])    
Main.upper_matrix <- as.matrix(Main.upper[,1:length(t_A)])   
TI.lower_matrix <- as.matrix(TI.lower[,1:length(t_A)])    
TI.upper_matrix <- as.matrix(TI.upper[,1:length(t_A)])   


Main.L.t <- apply(Main.L, 1, mean, na.rm=TRUE)       #transform it to a vector
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)
Main.lower_vector<- apply(Main.lower_matrix, 1, mean, na.rm=TRUE)
Main.upper_vector<- apply(Main.upper_matrix, 1, mean, na.rm=TRUE)
TI.lower_vector<- apply(TI.lower_matrix, 1, mean, na.rm=TRUE)
TI.upper_vector<- apply(TI.upper_matrix, 1, mean, na.rm=TRUE)

sorting <- order(TI.L.t, decreasing = F)          #sorting form low to high based on total effect
TI.L.t  <- TI.L.t[sorting]
Main.L.t <- Main.L.t[sorting]
Main.lower_vector <-Main.lower_vector[sorting]
Main.upper_vector<-Main.upper_vector[sorting]
TI.lower_vector <-TI.lower_vector[sorting]
TI.upper_vector<-TI.upper_vector[sorting]

Main.L.t <- ifelse(Main.L.t <= 0, 0, Main.L.t)
tempC    <- t(cbind(Main.L.t, TI.L.t,Main.lower_vector,Main.upper_vector,TI.lower_vector,TI.upper_vector))
tempC<- as.data.frame(tempC[,c(47:57)]) #top 10 variables with the highest total effect 



Global_sa_top_ten<-as.data.frame(melt(tempC[1,]))
Global_sa_top_ten[,3:4]<-as.data.frame(melt(tempC[2,]))
Global_sa_top_ten <- Global_sa_top_ten[ ,-c(3) ]
colnames(Global_sa_top_ten)<-c("variable","main","total")
Global_sa_top_ten<-melt(Global_sa_top_ten, id="variable")
colnames(Global_sa_top_ten)<-c("variable","indices","value")


#organising the upper and lower Ci for the First order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_main_lower<-as.data.frame(melt(tempC[3,]))
Global_sa_main_lower <- Global_sa_main_lower[ ,-c(3) ]
colnames(Global_sa_main_lower)<-c("variable","lower")
Global_sa_main_lower<-melt(Global_sa_main_lower, id="variable")
colnames(Global_sa_main_lower)<-c("variable","indices","lower")

Global_sa_main_upper<-as.data.frame(melt(tempC[4,]))
Global_sa_main_upper <- Global_sa_main_upper[ ,-c(3) ]
colnames(Global_sa_main_upper)<-c("variable","upper")
Global_sa_main_upper<-melt(Global_sa_main_upper, id="variable")
colnames(Global_sa_main_upper)<-c("variable","indices","upper")

Global_sa_main<-cbind(Global_sa_main_upper,Global_sa_main_lower)
Global_sa_main <- Global_sa_main[ ,-c(1,2,4,5) ]

#organising the upper and lower Ci for the total order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_total_lower<-as.data.frame(melt(tempC[5,]))
Global_sa_total_lower <- Global_sa_total_lower[ ,-c(3) ]
colnames(Global_sa_total_lower)<-c("variable","lower")
Global_sa_total_lower<-melt(Global_sa_total_lower, id="variable")
colnames(Global_sa_total_lower)<-c("variable","indices","lower")

Global_sa_total_upper<-as.data.frame(melt(tempC[6,]))
Global_sa_total_upper <- Global_sa_total_upper[ ,-c(3) ]
colnames(Global_sa_total_upper)<-c("variable","upper")
Global_sa_total_upper<-melt(Global_sa_total_upper, id="variable")
colnames(Global_sa_total_upper)<-c("variable","indices","upper")

Global_sa_total<-cbind(Global_sa_total_upper,Global_sa_total_lower)
Global_sa_total <- Global_sa_total[ ,-c(1,2,4,5) ]

Global_sa_CI<-rbind(Global_sa_main,Global_sa_total)

Global_sa_top_ten<-cbind(Global_sa_top_ten,Global_sa_CI)

p_2<-ggplot(Global_sa_top_ten,(aes(fill=indices,x=variable, y=value )))+
  geom_bar(position="dodge", stat="identity")+
  ylim(-0.1,1)+
  theme_classic()+
  scale_fill_brewer(palette="Set1")+
  coord_flip()+
  geom_errorbar(aes(ymin=abs(lower), ymax=abs(upper)), width=.2,
                position=position_dodge(.9))+
  ggtitle('2 hours') 

t_SA <-4



for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    Main[,i]       <- sa$S[,1]  #First order indices
    TI[,i]        <- sa$T[,1]  #Total indices
    TI.lower[,i] <- sa$T[,4]   #Lower CL total indices
    TI.upper[,i] <- sa$T[,5]   #Upper CL total indices
    Main.lower[,i]<- sa$S[,4]   #Lower CL total indices
    Main.upper[,i]<-sa$S[,5]    #Lower CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


Main.L <- as.matrix(Main[,1:length(t_A)])       # as.matrix
TI.L  <- as.matrix(TI[,1:length(t_A)])

Main.lower_matrix <- as.matrix(Main.lower[,1:length(t_A)])    
Main.upper_matrix <- as.matrix(Main.upper[,1:length(t_A)])   
TI.lower_matrix <- as.matrix(TI.lower[,1:length(t_A)])    
TI.upper_matrix <- as.matrix(TI.upper[,1:length(t_A)])   


Main.L.t <- apply(Main.L, 1, mean, na.rm=TRUE)       #transform it to a vector
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)
Main.lower_vector<- apply(Main.lower_matrix, 1, mean, na.rm=TRUE)
Main.upper_vector<- apply(Main.upper_matrix, 1, mean, na.rm=TRUE)
TI.lower_vector<- apply(TI.lower_matrix, 1, mean, na.rm=TRUE)
TI.upper_vector<- apply(TI.upper_matrix, 1, mean, na.rm=TRUE)

sorting <- order(TI.L.t, decreasing = F)          #sorting form low to high based on total effect
TI.L.t  <- TI.L.t[sorting]
Main.L.t <- Main.L.t[sorting]
Main.lower_vector <-Main.lower_vector[sorting]
Main.upper_vector<-Main.upper_vector[sorting]
TI.lower_vector <-TI.lower_vector[sorting]
TI.upper_vector<-TI.upper_vector[sorting]

Main.L.t <- ifelse(Main.L.t <= 0, 0, Main.L.t)
tempC    <- t(cbind(Main.L.t, TI.L.t,Main.lower_vector,Main.upper_vector,TI.lower_vector,TI.upper_vector))
tempC<- as.data.frame(tempC[,c(47:57)]) #top 10 variables with the highest total effect 



Global_sa_top_ten<-as.data.frame(melt(tempC[1,]))
Global_sa_top_ten[,3:4]<-as.data.frame(melt(tempC[2,]))
Global_sa_top_ten <- Global_sa_top_ten[ ,-c(3) ]
colnames(Global_sa_top_ten)<-c("variable","main","total")
Global_sa_top_ten<-melt(Global_sa_top_ten, id="variable")
colnames(Global_sa_top_ten)<-c("variable","indices","value")


#organising the upper and lower Ci for the First order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_main_lower<-as.data.frame(melt(tempC[3,]))
Global_sa_main_lower <- Global_sa_main_lower[ ,-c(3) ]
colnames(Global_sa_main_lower)<-c("variable","lower")
Global_sa_main_lower<-melt(Global_sa_main_lower, id="variable")
colnames(Global_sa_main_lower)<-c("variable","indices","lower")

Global_sa_main_upper<-as.data.frame(melt(tempC[4,]))
Global_sa_main_upper <- Global_sa_main_upper[ ,-c(3) ]
colnames(Global_sa_main_upper)<-c("variable","upper")
Global_sa_main_upper<-melt(Global_sa_main_upper, id="variable")
colnames(Global_sa_main_upper)<-c("variable","indices","upper")

Global_sa_main<-cbind(Global_sa_main_upper,Global_sa_main_lower)
Global_sa_main <- Global_sa_main[ ,-c(1,2,4,5) ]

#organising the upper and lower Ci for the total order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_total_lower<-as.data.frame(melt(tempC[5,]))
Global_sa_total_lower <- Global_sa_total_lower[ ,-c(3) ]
colnames(Global_sa_total_lower)<-c("variable","lower")
Global_sa_total_lower<-melt(Global_sa_total_lower, id="variable")
colnames(Global_sa_total_lower)<-c("variable","indices","lower")

Global_sa_total_upper<-as.data.frame(melt(tempC[6,]))
Global_sa_total_upper <- Global_sa_total_upper[ ,-c(3) ]
colnames(Global_sa_total_upper)<-c("variable","upper")
Global_sa_total_upper<-melt(Global_sa_total_upper, id="variable")
colnames(Global_sa_total_upper)<-c("variable","indices","upper")

Global_sa_total<-cbind(Global_sa_total_upper,Global_sa_total_lower)
Global_sa_total <- Global_sa_total[ ,-c(1,2,4,5) ]

Global_sa_CI<-rbind(Global_sa_main,Global_sa_total)

Global_sa_top_ten<-cbind(Global_sa_top_ten,Global_sa_CI)

p_4<-ggplot(Global_sa_top_ten,(aes(fill=indices,x=variable, y=value )))+
  geom_bar(position="dodge", stat="identity")+
  ylim(-0.1,1)+
  theme_classic()+
  scale_fill_brewer(palette="Set1")+
  coord_flip()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))+
  ggtitle('4 hours') 


t_SA <-8



for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    Main[,i]       <- sa$S[,1]  #First order indices
    TI[,i]        <- sa$T[,1]  #Total indices
    TI.lower[,i] <- sa$T[,4]   #Lower CL total indices
    TI.upper[,i] <- sa$T[,5]   #Upper CL total indices
    Main.lower[,i]<- sa$S[,4]   #Lower CL total indices
    Main.upper[,i]<-sa$S[,5]    #Lower CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


Main.L <- as.matrix(Main[,1:length(t_A)])       # as.matrix
TI.L  <- as.matrix(TI[,1:length(t_A)])

Main.lower_matrix <- as.matrix(Main.lower[,1:length(t_A)])    
Main.upper_matrix <- as.matrix(Main.upper[,1:length(t_A)])   
TI.lower_matrix <- as.matrix(TI.lower[,1:length(t_A)])    
TI.upper_matrix <- as.matrix(TI.upper[,1:length(t_A)])   


Main.L.t <- apply(Main.L, 1, mean, na.rm=TRUE)       #transform it to a vector
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)
Main.lower_vector<- apply(Main.lower_matrix, 1, mean, na.rm=TRUE)
Main.upper_vector<- apply(Main.upper_matrix, 1, mean, na.rm=TRUE)
TI.lower_vector<- apply(TI.lower_matrix, 1, mean, na.rm=TRUE)
TI.upper_vector<- apply(TI.upper_matrix, 1, mean, na.rm=TRUE)

sorting <- order(TI.L.t, decreasing = F)          #sorting form low to high based on total effect
TI.L.t  <- TI.L.t[sorting]
Main.L.t <- Main.L.t[sorting]
Main.lower_vector <-Main.lower_vector[sorting]
Main.upper_vector<-Main.upper_vector[sorting]
TI.lower_vector <-TI.lower_vector[sorting]
TI.upper_vector<-TI.upper_vector[sorting]

Main.L.t <- ifelse(Main.L.t <= 0, 0, Main.L.t)
tempC    <- t(cbind(Main.L.t, TI.L.t,Main.lower_vector,Main.upper_vector,TI.lower_vector,TI.upper_vector))
tempC<- as.data.frame(tempC[,c(47:57)]) #top 10 variables with the highest total effect 



Global_sa_top_ten<-as.data.frame(melt(tempC[1,]))
Global_sa_top_ten[,3:4]<-as.data.frame(melt(tempC[2,]))
Global_sa_top_ten <- Global_sa_top_ten[ ,-c(3) ]
colnames(Global_sa_top_ten)<-c("variable","main","total")
Global_sa_top_ten<-melt(Global_sa_top_ten, id="variable")
colnames(Global_sa_top_ten)<-c("variable","indices","value")


#organising the upper and lower Ci for the First order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_main_lower<-as.data.frame(melt(tempC[3,]))
Global_sa_main_lower <- Global_sa_main_lower[ ,-c(3) ]
colnames(Global_sa_main_lower)<-c("variable","lower")
Global_sa_main_lower<-melt(Global_sa_main_lower, id="variable")
colnames(Global_sa_main_lower)<-c("variable","indices","lower")

Global_sa_main_upper<-as.data.frame(melt(tempC[4,]))
Global_sa_main_upper <- Global_sa_main_upper[ ,-c(3) ]
colnames(Global_sa_main_upper)<-c("variable","upper")
Global_sa_main_upper<-melt(Global_sa_main_upper, id="variable")
colnames(Global_sa_main_upper)<-c("variable","indices","upper")

Global_sa_main<-cbind(Global_sa_main_upper,Global_sa_main_lower)
Global_sa_main <- Global_sa_main[ ,-c(1,2,4,5) ]

#organising the upper and lower Ci for the total order indices into the correct format so that it can be added to the global_sa_top_ten file
Global_sa_total_lower<-as.data.frame(melt(tempC[5,]))
Global_sa_total_lower <- Global_sa_total_lower[ ,-c(3) ]
colnames(Global_sa_total_lower)<-c("variable","lower")
Global_sa_total_lower<-melt(Global_sa_total_lower, id="variable")
colnames(Global_sa_total_lower)<-c("variable","indices","lower")

Global_sa_total_upper<-as.data.frame(melt(tempC[6,]))
Global_sa_total_upper <- Global_sa_total_upper[ ,-c(3) ]
colnames(Global_sa_total_upper)<-c("variable","upper")
Global_sa_total_upper<-melt(Global_sa_total_upper, id="variable")
colnames(Global_sa_total_upper)<-c("variable","indices","upper")

Global_sa_total<-cbind(Global_sa_total_upper,Global_sa_total_lower)
Global_sa_total <- Global_sa_total[ ,-c(1,2,4,5) ]

Global_sa_CI<-rbind(Global_sa_main,Global_sa_total)

Global_sa_top_ten<-cbind(Global_sa_top_ten,Global_sa_CI)

p_8<-ggplot(Global_sa_top_ten,(aes(fill=indices,x=variable, y=value )))+
  geom_bar(position="dodge", stat="identity")+
  ylim(-0.1,1)+
  theme_classic()+
  scale_fill_brewer(palette="Set1")+
  coord_flip()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  ggtitle('8 hours')


grid.arrange( p_0.5, p_1, p_1.5, p_2, p_4, p_8, ncol=3, nrow =2, top="Concentration in Lung tissue Rat 250mg oral exposure")


