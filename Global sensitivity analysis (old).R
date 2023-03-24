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

#Simulations
set.seed(15204)                       #to ensure a reproducible output if random input is used
MW                         <-132.16   #The molecular weight of Cinnamaldehyde
BW                         <-70      #Body weight in Kg
Volume_exposure_chamber    <-10       #volume exposure chamber in L


#--Physio-chemical parameters--#
#Values calculated in QSAR calculations file
#-Cinnamaldehyde-#

P_F      <-  47.75  #Fat/Blood partition coefficient
P_L      <-  1.83   #liver/Blood partition coefficient
P_SI     <-  1.81   #Small intestine/Blood partition coefficients
P_RP     <-  1.81   #Richly perfused tissues/Blood partition coefficients
P_SP     <-  1.50   #Slowly perfused tissues/Blood partition coefficients
P_B      <-  0.29  #Blood/Air Partition Coefficient 
P_Pu     <-  1.81   #lung/Blood partition coefficient

#-Cinnamyl Alcohol-#
P_OH_F    <-  49.26 #Fat/Blood partition coefficient
P_OH_L    <-  1.18  #liver/Blood partition coefficient
P_OH_SI   <-  1.18  #Small intestine/Blood partition coefficients
P_OH_RP   <-  1.18  #Richly perfused tissues/Blood partition coefficients
P_OH_SP   <-  1.53  #Slowly perfused tissues/Blood partition coefficients
P_OH_Pu   <-  1.18  #Lung/Blood partition coefficients

#-Tissues volumes in L #
V_F      <- 14.03 #Fat
V_L      <- 1.68 #Liver
V_SI     <- 0.6 #Small intestine
V_A      <- 1.31 #Arterial Blood
V_V      <- 3.87 #Venous Blood
V_RP     <- 3.14 #Richly perfused (RP)
V_SP     <- 40.39 #Slowly perfused (SP)
V_Pu     <- 4.98 #Lung

#-Cardiac parameters-#

Q_C       <-   390   #Cardiac output in L/h

#-Blood flow to tissues in L #
Q_F      <- 20.28  #Fat
Q_L      <- 54.99  #Liver
Q_SI     <- 33.54  #Small intestine
Q_RP     <- 184.47  #Richly perfused (RP)
Q_SP     <- 96.72  #Slowly perfused (SP)
Q_Pu     <- Q_C    #Lung


#inhalation parameters
P_V     <-   540             #Alveolar ventilation L/h

#----GSH parameters----#
#--GSH synthesis in μmol/kg tissue/h--#

G_SYN_L     <- 1122 * V_L        #Liver 
G_SYN_SI    <- 27   * V_SI        #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
k_L_GLOS    <- 0.142        #Liver
k_SI_GLOS   <- 0.044        #Small intestine

#--Initial GSH concentration--#
init_GSH_L  <- 7111 * V_L   #initial amount of GSH in the liver in μmol
init_GSH_SI <- 1555  * V_SI #initial amount of GSH in the small intestine in μmol

k_GSH <- 6.6 * 10^(-4)      #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
C_PRO_L     <- 5319   * V_L     #Liver
C_PRO_SI    <- 245    * V_SI      #Small intestine

#--Chemical parameters--#
Ka <- 0.62                   #Absorption rate constant for uptake in the Small intestine in per H

#----Liver----#
S9_scaling_L <- 143  * (V_L * 1000) #scaling factor for S9 fraction per g tissue

#----Liver----#

#-first rate order constants-#
k_L_OH  <- 4.2*10^-2 *  60 / 1000 * S9_scaling_L        #Scaled first rate order constant for the enzymatic oxidation of cinnamyl alcohol in the liver in umol/h

#--Michaelis menten constants--#
Km_L_CA     <-  8.5         #Km for enzymatic oxidation of Cinnamaldehyde into Cinnamic acid in the liver in μM
Km_L_AO     <-  330         #Km for enzymatic reduction of Cinnamaldehyde into cinnamyl alcOHol in the liver in μM
Km_L_GST    <-  1.7*10^3    #Km for enzymatic conjugation of Cinnamaldehyde with GST in the liver in μM  
Km_L_GST_G  <-  100         #Km toward GSH for enzymatic conjugation of Cinnamaldehyde in the liver (μM)

#--Vmax values--#
Vsmax_L_CA    <-  9.7  * 60 / 1000 * S9_scaling_L  #Scaled Vmax for enzymatic oxidation of Cinnamaldehyde in the liver in μmol/h 
Vsmax_L_AO    <-  73   * 60 / 1000 * S9_scaling_L      #Scaled Vmax for enzymatic reduction of Cinnamaldehyde in the liver in μmol/h
Vsmax_L_GST   <-  37    * 60 / 1000 * S9_scaling_L      #Scaled Vmax for enzymatic conjugation of Cinnamaldehyde with GSH in the liver in μmol/h

#----Small intestines----#
S9_scaling_SI <-  11.4 * (V_SI *1000) #scaling factor fraction S9 protein per g tissue

#--Michaelis Menten constants--#
Km_SI_CA    <- 70          #Km for enzymatic oxidation of Cinnamaldehyde into Cinnamic acid in the Small Intestine in μM
Km_SI_AO    <- 90          #Km for enzymatic reduction of Cinnamaldehyde into Cinnamyl alcOHol in the Small Intestine in μM
Km_SI_OH    <- 290         #Km for enzymatic oxidation of Cinnamly alcOHol into Cinnamaldehyde in the Small Intestine in μM
Km_SI_GST   <- 0           #Km for enzymatic conjugation of Cinnamaldehye with GST in the Small Intestine in μM (RAT value)
Km_SI_GST_G <- 100         #Km toward GSH for enzymatic conjugation of Cinnamaldehyde in the small intestine (μM)

#-Vmax values-#
Vsmax_SI_CA    <- 21    * 60/1000 * S9_scaling_SI   #Scaled Vmax for enzymatic oxidation of Cinnamaldehyde into Cinnamic acid in the Small Intestine in μmol/h 
Vsmax_SI_AO    <- 30   * 60/1000 * S9_scaling_SI    #Scaled Vmax for enzymatic reduction of Cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
Vsmax_SI_OH    <- 5.0  * 60/1000 * S9_scaling_SI   #Scaled Vmax for enzymatic Oxidation of Cinnamyl alcohol into Cinnamaldehyde in the Small Intestine in μmol/h 
Vsmax_SI_GST   <- 0    * 60/1000 * S9_scaling_SI   #Scaled Vmax for enzymatic Conjugation of Cinnamaldehyde with GSH in the in the small intestine in μmol/h (RAT value)

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
                 BW,
                 V_F,
                 V_L,
                 V_SI,
                 V_A,
                 V_V,
                 V_RP,
                 V_SP,
                 V_Pu,
                 Q_C,
                 Q_Pu,
                 Q_F,
                 Q_L,
                 Q_SI,
                 Q_RP,
                 Q_SP,
                 P_V,
                 G_SYN_L,
                 G_SYN_SI,
                 k_L_GLOS,
                 k_SI_GLOS,
                 init_GSH_L,
                 init_GSH_SI,
                 k_GSH,
                 C_PRO_L,
                 C_PRO_SI,
                 Ka,
                 k_L_OH,
                 Km_L_CA,
                 Km_L_AO,
                 Km_L_GST,
                 Km_L_GST_G,
                 Vsmax_L_CA,
                 Vsmax_L_AO,
                 Vsmax_L_GST,
                 Km_SI_CA,
                 Km_SI_AO,
                 Km_SI_OH,
                 Km_SI_GST,
                 Km_SI_GST_G,
                 Vsmax_SI_CA,
                 Vsmax_SI_AO,
                 Vsmax_SI_OH,
                 Vsmax_SI_GST,
                 Volume_exposure_chamber,
                 S9_scaling_SI,
                 S9_scaling_L)


#Creating empty vectors for use later
colnames <- c(colnames(dist_para))
par_var <- length(colnames)

#Generating a upper and lower bound for use in the analysis 
Mean <- dist_para

Lower <- Mean - 0.1 * Mean
Upper <- Mean + 0.1 * Mean

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


#Removing variation in volume exposure chamber as this is supposed to be fixed
X1[,58]<- 10
X2[,58]<- 10


#the number of bootstrap replicates
n_boot <- 1000


#Sobol design
sa <- soboljansen(model=NULL, X1, X2, nboot = n_boot, conf = 0.95, events = ex)


phys <- sa$X


#Writing the result into a file so that the environment can be cleaned to conserve memory
write.csv(phys,"D:/PBK/Cinnamaldehyde-pbk\\GSA_phys_human_250mg_10%", row.names = TRUE)



#Loading extracted simulation data. 
solve.pbk.sa <- read.csv("D:/PBK/Cinnamaldehyde-pbk\\SA_human_250mg_oral_CV", row.names=1)

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


SimRes = as.data.frame(matrix(NA,124000,8))
SimRes[,1]=tab1[,2]
SimRes[,2]=tab2[,2]
SimRes[,3]=tab3[,2]
SimRes[,4]=tab4[,2]
SimRes[,5]=tab5[,2]
SimRes[,6]=tab6[,2]
SimRes[,7]=tab7[,2]
SimRes[,8]=tab8[,2]


write.csv(SimRes,"D:/PBK/Cinnamaldehyde-pbk\\SimRes_oral_Human2_250g.csv", row.names = TRUE)

SimRes <- read.csv("D:/PBK/Cinnamaldehyde-pbk\\SimRes.inhalation_C_Pu_popgen-100mg_4000.csv", row.names=1)

#Redefining these two variables as these are also used with dist_parm creation but not all of thet variables in dist_parm are used in the SA calculation
#so using them here would create an error.
colnames <- colnames(X1)
par_var <- length(X1)


#Sobol analysis plot blood Nrow is the number of paramters in the model
t_A<-(c(0.2,0.5,1,1.5,2,3,4,8))
FOI          = TI          = TI.borninf           = TI.bornsup          = matrix(NA, nrow = par_var, ncol = length(t_A))  
rownames(FOI)= rownames(TI)= rownames(TI.borninf) = rownames(TI.bornsup)= colnames

t_SA <-0.2


for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    FOI[,i]       <- sa$S[,1]    #First order indices
    TI[,i]        <- sa$T[,1]    #Total indices
    TI.borninf[,i] <- sa$T[,4]   #Lower CL total indices
    TI.bornsup[,i] <- sa$T[,5]   #Upper CL total indices
   
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


FOI.L = as.matrix(FOI[,1:length(t_A)])       # as.matrix
TI.L  = as.matrix(TI[,1:length(t_A)])

FOI.L.t <- apply(FOI.L, 1, mean, na.rm=TRUE)
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)

sorting = order(TI.L.t, decreasing = F)
TI.L.t  = TI.L.t[sorting]
FOI.L.t = FOI.L.t[sorting]

FOI.L.t = ifelse(FOI.L.t <= 0, 0, FOI.L.t)
tempC    = t(cbind(FOI.L.t, TI.L.t))

tempC2 <- as.data.frame(tempC[,c(50:60)])

#t_SA = 0.2
sa.plot_0.2 <-as.data.frame(tempC[,c(50:60)])
rownames(sa.plot_0.2 <- c("0.2h main","0.2h total"))

par(mfrow=c(1,1), las=1, mai=c(0.35,1,0.35,0.1), mgp = c(3.5,0.5,0))
plot_Pu_0.2 <- barplot(as.matrix(tempC2), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="", cex.lab=1.5 , xlim=c(0,1.1) )

t_SA <-0.5


for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    FOI[,i]       <- sa$S[,1]    #First order indices
    TI[,i]        <- sa$T[,1]    #Total indices
    TI.borninf[,i] <- sa$T[,4]   #Lower CL total indices
    TI.bornsup[,i] <- sa$T[,5]   #Upper CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


FOI.L = as.matrix(FOI[,1:length(t_A)])       # as.matrix
TI.L  = as.matrix(TI[,1:length(t_A)])

FOI.L.t <- apply(FOI.L, 1, mean, na.rm=TRUE)
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)

sorting = order(TI.L.t, decreasing = F)
TI.L.t  = TI.L.t[sorting]
FOI.L.t = FOI.L.t[sorting]

FOI.L.t = ifelse(FOI.L.t <= 0, 0, FOI.L.t)
tempC    = t(cbind(FOI.L.t, TI.L.t))

tempC2 <- as.data.frame(tempC[,c(50:60)])

#t_SA = 0.5
sa.plot_0.5 <-as.data.frame(tempC[,c(50:60)])
rownames(sa.plot_0.5) <- c("0.5h main","0.5h total")

par(mfrow=c(1,1), las=1, mai=c(0.35,1,0.35,0.1), mgp = c(3.5,0.5,0))
#colnames(tempC2) <- c("Q_SI", "Ka", "V_SP", "Q_RP", "P_SP", "Q_SP", "k_GSH", "C_PRO_L", "VL", "QC")
plot_Pu_0.5 <- barplot(as.matrix(tempC2), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="", cex.lab=1.5 , xlim=c(0,1.1) )

t_SA <-1


for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    FOI[,i]       <- sa$S[,1]    #First order indices
    TI[,i]        <- sa$T[,1]    #Total indices
    TI.borninf[,i] <- sa$T[,4]   #Lower CL total indices
    TI.bornsup[,i] <- sa$T[,5]   #Upper CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


FOI.L = as.matrix(FOI[,1:length(t_A)])       # as.matrix
TI.L  = as.matrix(TI[,1:length(t_A)])

FOI.L.t <- apply(FOI.L, 1, mean, na.rm=TRUE)
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)

sorting = order(TI.L.t, decreasing = F)
TI.L.t  = TI.L.t[sorting]
FOI.L.t = FOI.L.t[sorting]

FOI.L.t = ifelse(FOI.L.t <= 0, 0, FOI.L.t)
tempC    = t(cbind(FOI.L.t, TI.L.t))

tempC2 <- as.data.frame(tempC[,c(50:60)])
#t_SA = 1
sa.plot_1 <-as.data.frame(tempC[,c(50:60)])
rownames(sa.plot_1) <- c("1h main","1h total")


par(mfrow=c(1,1), las=1, mai=c(0.35,1,0.35,0.1), mgp = c(3.5,0.5,0))
#colnames(tempC2) <- c("Q_SI", "Ka", "V_SP", "Q_RP", "P_SP", "Q_SP", "k_GSH", "C_PRO_L", "VL", "QC")
plot_Pu_1 <- barplot(as.matrix(tempC2), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="", cex.lab=1.5 , xlim=c(0,1.1) )

t_SA <-1.5


for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    FOI[,i]       <- sa$S[,1]    #First order indices
    TI[,i]        <- sa$T[,1]    #Total indices
    TI.borninf[,i] <- sa$T[,4]   #Lower CL total indices
    TI.bornsup[,i] <- sa$T[,5]   #Upper CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


FOI.L = as.matrix(FOI[,1:length(t_A)])       # as.matrix
TI.L  = as.matrix(TI[,1:length(t_A)])

FOI.L.t <- apply(FOI.L, 1, mean, na.rm=TRUE)
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)

sorting = order(TI.L.t, decreasing = F)
TI.L.t  = TI.L.t[sorting]
FOI.L.t = FOI.L.t[sorting]

FOI.L.t = ifelse(FOI.L.t <= 0, 0, FOI.L.t)
tempC    = t(cbind(FOI.L.t, TI.L.t))

tempC2 <- as.data.frame(tempC[,c(50:60)])

#t_SA = 1.5
sa.plot_1.5 <-as.data.frame(tempC[,c(50:60)])
rownames(sa.plot_1.5) <- c("1.5h main","1.5h total")

par(mfrow=c(1,1), las=1, mai=c(0.35,1,0.35,0.1), mgp = c(3.5,0.5,0))
#colnames(tempC2) <- c("Q_SI", "Ka", "V_SP", "Q_RP", "P_SP", "Q_SP", "k_GSH", "C_PRO_L", "VL", "QC")
plot_Pu_1.5 <- barplot(as.matrix(tempC2), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="", cex.lab=1.5 , xlim=c(0,1.1) )

t_SA <-2


for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    FOI[,i]       <- sa$S[,1]    #First order indices
    TI[,i]        <- sa$T[,1]    #Total indices
    TI.borninf[,i] <- sa$T[,4]   #Lower CL total indices
    TI.bornsup[,i] <- sa$T[,5]   #Upper CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


FOI.L = as.matrix(FOI[,1:length(t_A)])       # as.matrix
TI.L  = as.matrix(TI[,1:length(t_A)])

FOI.L.t <- apply(FOI.L, 1, mean, na.rm=TRUE)
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)

sorting = order(TI.L.t, decreasing = F)
TI.L.t  = TI.L.t[sorting]
FOI.L.t = FOI.L.t[sorting]

FOI.L.t = ifelse(FOI.L.t <= 0, 0, FOI.L.t)
tempC    = t(cbind(FOI.L.t, TI.L.t))

tempC2 <- as.data.frame(tempC[,c(50:60)])

par(mfrow=c(1,1), las=1, mai=c(0.35,1,0.35,0.1), mgp = c(3.5,0.5,0))
#colnames(tempC2) <- c("Q_SI", "Ka", "V_SP", "Q_RP", "P_SP", "Q_SP", "k_GSH", "C_PRO_L", "VL", "QC")
plot_Pu_2 <- barplot(as.matrix(tempC2), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="", cex.lab=1.5 , xlim=c(0,1.1) )

#t_SA = 2
sa.plot_2 <-as.data.frame(tempC[,c(50:60)])
rownames(sa.plot_2) <- c("2h main","2h total")

t_SA <-4


for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    FOI[,i]       <- sa$S[,1]    #First order indices
    TI[,i]        <- sa$T[,1]    #Total indices
    TI.borninf[,i] <- sa$T[,4]   #Lower CL total indices
    TI.bornsup[,i] <- sa$T[,5]   #Upper CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


FOI.L = as.matrix(FOI[,1:length(t_A)])       # as.matrix
TI.L  = as.matrix(TI[,1:length(t_A)])

FOI.L.t <- apply(FOI.L, 1, mean, na.rm=TRUE)
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)

sorting = order(TI.L.t, decreasing = F)
TI.L.t  = TI.L.t[sorting]
FOI.L.t = FOI.L.t[sorting]

FOI.L.t = ifelse(FOI.L.t <= 0, 0, FOI.L.t)
tempC    = t(cbind(FOI.L.t, TI.L.t))

tempC2 <- as.data.frame(tempC[,c(50:60)])

par(mfrow=c(1,1), las=1, mai=c(0.35,1,0.35,0.1), mgp = c(3.5,0.5,0))
#colnames(tempC2) <- c("Q_SI", "Ka", "V_SP", "Q_RP", "P_SP", "Q_SP", "k_GSH", "C_PRO_L", "VL", "QC")
plot_Pu_4 <- barplot(as.matrix(tempC2), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="", cex.lab=1.5 , xlim=c(0,1.1) )

#t_SA = 4
sa.plot_4 <-as.data.frame(tempC[,c(50:60)])
rownames(sa.plot_4) <- c("4h main","4h total")

t_SA <-8


for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    sa=tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    FOI[,i]       <- sa$S[,1]    #First order indices
    TI[,i]        <- sa$T[,1]    #Total indices
    TI.borninf[,i] <- sa$T[,4]   #Lower CL total indices
    TI.bornsup[,i] <- sa$T[,5]   #Upper CL total indices
    
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}


FOI.L = as.matrix(FOI[,1:length(t_A)])       # as.matrix
TI.L  = as.matrix(TI[,1:length(t_A)])

FOI.L.t <- apply(FOI.L, 1, mean, na.rm=TRUE)
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)

sorting = order(TI.L.t, decreasing = F)
TI.L.t  = TI.L.t[sorting]
FOI.L.t = FOI.L.t[sorting]

FOI.L.t = ifelse(FOI.L.t <= 0, 0, FOI.L.t)
tempC    = t(cbind(FOI.L.t,TI.L.t))

tempC2 <- as.data.frame(tempC[,c(50:60)])

#t_SA = 8
sa.plot_8 <-as.data.frame(tempC[,c(50:60)])
rownames(sa.plot_8) <- c("8h total","8h main")

par(mfrow=c(2,3),las=1, mar=c(3,5,3,3.5), mgp = c(3.5,0.5,0)) 
mtext("Main and total sensitivity indexes", side=2, line=1)
#barplot(as.matrix(sa.plot_0.2), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="20 min", cex.lab=1.5 , xlim=c(0,1.1) )
barplot(as.matrix(sa.plot_0.5), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="30 min", cex.lab=1.5 , xlim=c(0,1.1) )
barplot(as.matrix(sa.plot_1), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="1 hour", cex.lab=1.5 , xlim=c(0,1.1) )
barplot(as.matrix(sa.plot_1,5), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="1.5 hours", cex.lab=1.5 , xlim=c(0,1.1) )
barplot(as.matrix(sa.plot_2), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="2 hours", cex.lab=1.5 , xlim=c(0,1.1) )
barplot(as.matrix(sa.plot_4), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="4 hours", cex.lab=1.5 , xlim=c(0,1.1) )
barplot(as.matrix(sa.plot_8), col=c("firebrick1","firebrick4"), horiz =T , beside =T , main="8 hours", cex.lab=1.5 , xlim=c(0,1.1))





