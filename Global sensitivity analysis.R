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

P_F      <-  1.62   #Fat/Blood partition coefficient
P_L      <-  0.59   #liver/Blood partition coefficient
P_SI     <-  0.59   #Small intestine/Blood partition coefficients
P_RP     <-  0.59   #Richly perfused tissues/Blood partition coefficients
P_SP     <-  0.78   #Slowly perfused tissues/Blood partition coefficients
P_B      <-  1.25E5 #Blood/Air Partition Coefficient 
P_Pu     <-  0.59   #lung/Blood partition coefficient

#-Cinnamyl Alcohol-#
P_OH_F    <-  1.64 #Fat/Blood partition coefficient
P_OH_L    <-  0.59  #liver/Blood partition coefficient
P_OH_SI   <-  0.59  #Small intestine/Blood partition coefficients
P_OH_RP   <-  0.59  #Richly perfused tissues/Blood partition coefficients
P_OH_SP   <-  0.78 #Slowly perfused tissues/Blood partition coefficients
P_OH_Pu   <-  0.59 #Lung/Blood partition coefficients

#--Physiological Parameters--#

fV_fat   <- 0.2142  #Unit less fraction of total tissue volume that consist of fat
fV_liver <- 0.0257  #Unit less fraction of total tissue volume that consist of  liver
fV_si    <- 0.0091  #Unit less fraction of total tissue volume that consist of Small intestine
fV_a     <- 0.02  #Unit less fraction of total tissue volume that consist of arterial blood
fV_v     <- 0.059  #Unit less fraction of total tissue volume that consist of venous blood
fV_rp    <- 0.048  #Unit less fraction of total tissue volume that consist of Richly perfused tissue
fV_sp    <- 0.6164  #Unit less fraction of total tissue volume that consist of Slowly perfused tissue
fV_pu    <- 0.076  #Unit less fraction of total tissue volume that consist of

V_adjust <- fV_fat + fV_liver + fV_si + fV_a + fV_v + fV_rp + fV_sp +fV_pu

#-Tissues volumes in L #
V_F      <- (fV_fat / V_adjust)* BW    #Fat
V_L      <- (fV_liver / V_adjust)* BW  #Liver
V_SI     <- (fV_si / V_adjust)* BW     #Small intestine
V_A      <- (fV_a / V_adjust)* BW  #Arterial Blood
V_V      <- (fV_v / V_adjust)* BW  #Venous Blood
V_RP     <- (fV_rp / V_adjust)* BW  #Richly perfused (RP)
V_SP     <- (fV_sp / V_adjust)* BW   #Slowly perfused (SP)
V_Pu     <- (fV_pu / V_adjust)* BW    #Lung

#-Cardiac parameters-#

Q_C       <- 390   #Cardiac output in L/h

fQC_fat   <- 0.052  #Unit less fraction of total Cardiac output that goes to fat
fQC_liver <- 0.141  #Unit less fraction of total Cardiac output that goes to the liver
fQC_si    <- 0.086  #Unit less fraction of total Cardiac output that goes to the Small intestine
fQC_rp    <- 0.473  #Unit less fraction of total Cardiac output that goes to the Richly perfused tissue
fQC_sp    <- 0.248  #Unit less fraction of total Cardiac output that goes to the Slowly perfused tissue

Q_adjust <- fQC_fat + fQC_liver + fQC_si + fQC_rp + fQC_sp

#-Blood flow to tissues in % cardiac output-#
Q_F      <- (fQC_fat / Q_adjust)* Q_C    #Fat
Q_L      <- (fQC_liver / Q_adjust)* Q_C  #Liver
Q_SI     <- (fQC_si / Q_adjust)* Q_C     #Small intestine
Q_RP     <- (fQC_rp / Q_adjust)* Q_C  #Richly perfused (RP)
Q_SP     <- (fQC_sp / Q_adjust)* Q_C   #Slowly perfused (SP)
Q_Pu     <- Q_C    #Lung

#inhalation parameters
P_V     <- 540             #Alveolar ventilation L/h

#----GSH parameters----#
#--GSH synthesis in μmol/kg tissue/h--#

G_SYN_L     <- 1122         #Liver 
G_SYN_SI    <- 27           #Small intestine

#-Apparent first order rate constant GSH turn over(RAT?) per h-#
k_L_GLOS    <- 0.142        #Liver
k_SI_GLOS   <- 0.044        #Small intestine

#--Initial GSH concentration--#
init_GSH_L  <- 5639 * V_L   #initial amount of GSH in the liver in μmol/kg
init_GSH_SI <- 1250 * V_SI  #initial amount of GSH in the small intestine in μmol/kg

k_GSH <- 6.6 * 10^(-4)      #The second-order rate constant of the chemical reaction of cinnamaldehyde with GSH in μmol/h
k_DNA <- 1.6 * 10^(-8)      #The second-order rate constant of the reaction between cinnamaldehyde and 2ʹ-dG in μmol/h

#----Protein reactive sites in μmol/kg tissue----#
C_PRO_L     <- 5319         #Liver
C_PRO_SI    <- 245          #Small intestine

#----DNA parameters----#
C_L_dG     <-  1.36         #Concentration of 2ʹ-dG in the liver μmol/kg liver
T_0.5      <-  38.5         #Half-life of DNA adduct in the liver in hours


#--Chemical parameters--#
Ka <- 5.0                   #Absorption rate constant for uptake in the Small intestine in per H

#----Liver----#

#-first rate order constants-#
k_L_OH  <- 4.2e-02          #Scaled first rate order constant for the enzymatic oxidation of cinnamyl alcohol in the liver in umol/h

#--Michaelis menten constants--#
Km_L_CA     <-  8.5         #Km for enzymatic oxidation of Cinnamaldehyde into Cinnamic acid in the liver in μM
Km_L_AO     <-  330         #Km for enzymatic reduction of Cinnamaldehyde into cinnamyl alcOHol in the liver in μM
Km_L_GST    <-  100         #Km for enzymatic conjugation of Cinnamaldehyde with GST in the liver in μM  
Km_L_GST_G  <-  1.7*10^3    #Km toward GSH for enzymatic conjugation of Cinnamaldehyde in the liver (μM)

#--Vmax values--#
Vsmax_L_CA    <-  9.7       #Scaled Vmax for enzymatic oxidation of Cinnamaldehyde in the liver in μmol/h 
Vsmax_L_AO    <-  73        #Scaled Vmax for enzymatic reduction of Cinnamaldehyde in the liver in μmol/h
Vsmax_L_GST   <-  37        #Scaled Vmax for enzymatic conjugation of Cinnamaldehyde with GSH in the liver in μmol/h

#----Small intestines----#
#--Michaelis Menten constants--#
Km_SI_CA    <- 70          #Km for enzymatic oxidation of Cinnamaldehyde into Cinnamic acid in the Small Intestine in μM
Km_SI_AO    <- 90          #Km for enzymatic reduction of Cinnamaldehyde into Cinnamyl alcOHol in the Small Intestine in μM
Km_SI_OH    <- 290         #Km for enzymatic oxidation of Cinnamly alcOHol into Cinnamaldehyde in the Small Intestine in μM
Km_SI_GST   <- 600         #Km for enzymatic conjugation of Cinnamaldehye with GST in the Small Intestine in μM (RAT value)
Km_SI_GST_G <- 100         #Km toward GSH for enzymatic conjugation of Cinnamaldehyde in the small intestine (μM)

#-Vmax values-#
Vsmax_SI_CA    <- 21       #Scaled Vmax for enzymatic oxidation of Cinnamaldehyde into Cinnamic acid in the Small Intestine in μmol/h 
Vsmax_SI_AO    <- 30       #Scaled Vmax for enzymatic reduction of Cinnamaldehyde into Cinnamyl alcOHol in  the Small Intestine in μmol/h 
Vsmax_SI_OH    <- 5.0      #Scaled Vmax for enzymatic Oxidation of Cinnamyl alcohol into Cinnamaldehyde in the Small Intestine in μmol/h 
Vsmax_SI_GST   <- 63       #Scaled Vmax for enzymatic Conjugation of Cinnamaldehyde with GSH in the in the small intestine in μmol/h (RAT value)

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
                 fV_fat,   
                 fV_liver, 
                 fV_si,    
                 fV_a,     
                 fV_v,     
                 fV_rp,    
                 fV_sp,    
                 fV_pu,   
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
                 fQC_fat,   
                 fQC_liver,
                 fQC_si,    
                 fQC_rp,    
                 fQC_sp,    
                 Q_adjust, 
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
                 k_DNA,
                 C_PRO_L,
                 C_PRO_SI,
                 C_L_dG,
                 T_0.5,
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
                 Volume_exposure_chamber)


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

#Defining empty vectors to be used in the following calculations #Human
Flow_volumes <-c()
Tissue_volumes<-c()
Total_body_weight<-c()
Total_QC<-c()
#Based on the Changes in the volume or blood fractions recalculating blood flows and tissue volumes to ensure mass balance
for(i in 1:nrow(X1)){
  #Tissue volumes calculations
    X1[i,23]<- sum(X1[i,15:22])
    X1[i,24]<- (X1[i,15]/X1[i,23])* X1[i,14]
    X1[i,25]<- (X1[i,16]/X1[i,23])* X1[i,14]
    X1[i,26]<- (X1[i,17]/X1[i,23])* X1[i,14]
    X1[i,27]<- (X1[i,18]/X1[i,23])* X1[i,14]
    X1[i,28]<- (X1[i,19]/X1[i,23])* X1[i,14]
    X1[i,29]<- (X1[i,20]/X1[i,23])* X1[i,14]
    X1[i,30]<- (X1[i,21]/X1[i,23])* X1[i,14]
    X1[i,31]<- (X1[i,22]/X1[i,23])* X1[i,14]
    #Blood flow calculations
    X1[i,38]<- sum(X1[i,33:37])
    X1[i,39]<- X1[i,32]
    X1[i,40]<- (X1[i,33]/X1[i,38])* X1[i,32]
    X1[i,41]<- (X1[i,34]/X1[i,38])* X1[i,32]
    X1[i,42]<- (X1[i,35]/X1[i,38])* X1[i,32]
    X1[i,43]<- (X1[i,36]/X1[i,38])* X1[i,32]
    X1[i,44]<- (X1[i,37]/X1[i,38])* X1[i,32]
    Flow_volumes <-append(Flow_volumes,sum(X1[i,40:44]))
    Tissue_volumes<-append(Tissue_volumes ,sum(X1[i,24:31]))
    Total_body_weight <- append(Total_body_weight ,X1[i,14])
    Total_QC<- append(Total_QC ,X1[i,32])
    names<- c("total flow","QC","Total volume","BW")
    Mass_check <- as.data.frame(cbind(Flow_volumes,Total_QC,Tissue_volumes,Total_body_weight), col.names = names)
  }
#Calculations for X2
for(i in 1:nrow(X2)){
  #Tissue volumes calculations
  X2[i,23]<- sum(X2[i,15:22])
  X2[i,24]<- (X2[i,15]/X2[i,23])* X2[i,14]
  X2[i,25]<- (X2[i,16]/X2[i,23])* X2[i,14]
  X2[i,26]<- (X2[i,17]/X2[i,23])* X2[i,14]
  X2[i,27]<- (X2[i,18]/X2[i,23])* X2[i,14]
  X2[i,28]<- (X2[i,19]/X2[i,23])* X2[i,14]
  X2[i,29]<- (X2[i,20]/X2[i,23])* X2[i,14]
  X2[i,30]<- (X2[i,21]/X2[i,23])* X2[i,14]
  X2[i,31]<- (X2[i,22]/X2[i,23])* X2[i,14]
  #Blood flow calculations
  X2[i,38]<- sum(X2[i,33:37])
  X2[i,39]<- X2[i,32]
  X2[i,40]<- (X2[i,33]/X2[i,38])* X2[i,32]
  X2[i,41]<- (X2[i,34]/X2[i,38])* X2[i,32]
  X2[i,42]<- (X2[i,35]/X2[i,38])* X2[i,32]
  X2[i,43]<- (X2[i,36]/X2[i,38])* X2[i,32]
  X2[i,44]<- (X2[i,37]/X2[i,38])* X2[i,32]
}

#Removing variation in volume exposure chamber as this is supposed to be fixed
X1[,76]<- 10
X2[,76]<- 10

#Removing unnecessary variables 
X1<- X1[-c(15:23,33:38)]
X2<- X2[-c(15:23,33:38)]
#the number of bootstrap replicates
n_boot <- 1000


X1 <-GSA_pop[1:2000,]
X2 <-GSA_pop[2001:4000,]

X1<- X1[-c(1:2,33)]
X2<- X2[-c(1:2,33)]


#Sobol design
sa <- soboljansen(model=NULL, X1, X2, nboot = n_boot, conf = 0.95, events = ex)


phys <- sa$X


#Writing the result into a file so that the environment can be cleaned to conserve memory
write.csv(phys,"D:/PBK/Cinnamaldehyde-pbk\\GSA_phys", row.names = TRUE)



#Loading extracted simulation data. 
solve.pbk.sa <- read.csv("D:/PBK/Cinnamaldehyde-pbk\\solve.pbk.sa", row.names=1)

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


SimRes = as.data.frame(matrix(NA,126000,8))
SimRes[,1]=tab1[,2]
SimRes[,2]=tab2[,2]
SimRes[,3]=tab3[,2]
SimRes[,4]=tab4[,2]
SimRes[,5]=tab5[,2]
SimRes[,6]=tab6[,2]
SimRes[,7]=tab7[,2]
SimRes[,8]=tab8[,2]


write.csv(SimRes,"D:/PBK/Cinnamaldehyde-pbk\\SimRes.inhalation_C_Pu_popgen-100mg_4000.csv", row.names = TRUE)

SimRes <- read.csv("D:/PBK/Cinnamaldehyde-pbk\\SimRes.inhalation_C_Pu_popgen-100mg_4000.csv", row.names=1)

#Redefining these two variables as these are also used with dist_parm creation but not all of thet variables in dist_parm are used in the SA calculation
#so using them here woul create an error.
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

tempC2 <- as.data.frame(tempC[,c(51:61)])

#t_SA = 0.2
sa.plot_0.2 <-as.data.frame(tempC[,c(51:61)])
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

tempC2 <- as.data.frame(tempC[,c(51:61)])

#t_SA = 0.5
sa.plot_0.5 <-as.data.frame(tempC[,c(51:61)])
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

tempC2 <- as.data.frame(tempC[,c(51:61)])
#t_SA = 1
sa.plot_1 <-as.data.frame(tempC[,c(51:61)])
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

tempC2 <- as.data.frame(tempC[,c(51:61)])

#t_SA = 1.5
sa.plot_1.5 <-as.data.frame(tempC[,c(51:61)])
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

tempC2 <- as.data.frame(tempC[,c(51:61)])

par(mfrow=c(1,1), las=1, mai=c(0.35,1,0.35,0.1), mgp = c(3.5,0.5,0))
#colnames(tempC2) <- c("Q_SI", "Ka", "V_SP", "Q_RP", "P_SP", "Q_SP", "k_GSH", "C_PRO_L", "VL", "QC")
plot_Pu_2 <- barplot(as.matrix(tempC2), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="", cex.lab=1.5 , xlim=c(0,1.1) )

#t_SA = 2
sa.plot_2 <-as.data.frame(tempC[,c(51:61)])
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

tempC2 <- as.data.frame(tempC[,c(51:61)])

par(mfrow=c(1,1), las=1, mai=c(0.35,1,0.35,0.1), mgp = c(3.5,0.5,0))
#colnames(tempC2) <- c("Q_SI", "Ka", "V_SP", "Q_RP", "P_SP", "Q_SP", "k_GSH", "C_PRO_L", "VL", "QC")
plot_Pu_4 <- barplot(as.matrix(tempC2), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="", cex.lab=1.5 , xlim=c(0,1.1) )

#t_SA = 4
sa.plot_4 <-as.data.frame(tempC[,c(51:61)])
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

tempC2 <- as.data.frame(tempC[,c(51:61)])

#t_SA = 8
sa.plot_8 <-as.data.frame(tempC[,c(51:61)])
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





