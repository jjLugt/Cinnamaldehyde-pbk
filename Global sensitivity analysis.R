#author: Joris Jean van der Lugt
#date: 20-05-2021
#Global sensitivity analysis


library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(shiny)
library(truncnorm)
library(reshape2)
library(sensitivity)

#Generate a population using a population data set either using Popgen or not

#Generating a upper and lower bound for use in the analysis 
colnames <- c(colnames(phys))
par_var <- length(colnames)

Mean <- phys[1,]

Lower <- Mean - 0.1 * Mean
Upper <- Mean + 0.1 * Mean

Lower <- Mean - 0.2 * Mean
Upper <- Mean + 0.2 * Mean


#create data frames for population
n_sim  <- 2000               #number of iterations
X1 <- matrix(NA, nrow = n_sim, ncol = par_var)
colnames(X1) <- colnames
X1 <- as.data.frame(X1)
var <- X1

X2 <- matrix(NA, nrow = n_sim, ncol = par_var)
colnames(X2) <- colnames
X2 <- as.data.frame(X2)
var <- X2

#create distribution population
for(i in 1:par_var){
  X1[,i] <- runif(n_sim, min = Lower[,i], max = Upper[,i])
  X2[,i] <- runif(n_sim, min = Lower[,i], max = Upper[,i])
}

n_boot <- 1000

#Sobol design
sa <- soboljansen(model= NULL , X1, X2, nboot = n_boot, conf = 0.95, events = ex)
phys <- sa$X


P_F<-phys$P_F
P_L<-phys$P_L
P_SI<-phys$P_SI
P_RP<-phys$P_RP
P_SP<-phys$P_SP
P_OH_F<-phys$P_OH_F
P_OH_L<-phys$P_OH_L
P_OH_SI<-phys$P_OH_SI
P_OH_RP<-phys$P_OH_RP
P_OH_SP<-phys$P_OH_SP
Age <- phys$Age
Height_start <- phys$Height_start
Height_cv <- phys$Height_cv
Height <- phys$Height
BW_start <- phys$BW_start
BW_cv <- phys$BW_cv
BW<-phys$BW
BSA <- phys$BSA
V_L<-phys$V_L
V_F<-phys$V_F
V_F_min <- phys$V_F_min
V_B <- phys$V_B
V_A<-phys$V_A
V_V<-phys$V_V
V_SI<-phys$V_SI
V_RP<-phys$V_RP
V_SP<-phys$V_SP
Q_C<-phys$Q_C
Q_SI<-phys$Q_SI
Q_F<-phys$Q_F
Q_L<-phys$Q_L
Q_RP<-phys$Q_RP
Q_SP<-phys$Q_SP
G_SYN_L<-phys$G_SYN_L
G_SYN_SI<-phys$G_SYN_SI
k_L_GLOS<-phys$k_L_GLOS
k_SI_GLOS<-phys$k_SI_GLOS
init_GSH_L<-phys$init_GSH_L
init_GSH_SI<-phys$init_GSH_SI
k_GSH<-phys$k_GSH
k_DNA<-phys$k_DNA
C_PRO_L<-phys$C_PRO_L
C_PRO_SI<-phys$C_PRO_SI
C_L_dG<-phys$C_L_dG
T_0.5<-phys$T_0.5
Ka<-phys$Ka
k_L_OH <- phys$k_L_OH
Km_L_CA<-phys$Km_L_CA
Km_L_AO<-phys$Km_L_AO
Km_L_GST<-phys$Km_L_GST
Km_L_GST_G<-phys$Km_L_GST_G
Vsmax_L_CA<-phys$Vsmax_L_CA
Vsmax_L_AO<-phys$Vsmax_L_AO
Vsmax_L_GST<-phys$Vsmax_L_GST
Km_SI_CA<-phys$Km_SI_CA
Km_SI_AO<-phys$Km_SI_AO
Km_SI_OH<-phys$Km_SI_OH
Km_SI_GST<-phys$Km_SI_GST
Km_SI_GST_G<-phys$Km_SI_GST_G
Vsmax_SI_CA<-phys$Vsmax_SI_CA
Vsmax_SI_AO<-phys$Vsmax_SI_AO
Vsmax_SI_OH<-phys$Vsmax_SI_OH
Vsmax_SI_GST<-phys$Vsmax_SI_GST
DOSE<- phys$DOSE 



#Run the model after assigning the sobol dataset to the variables 


#Analysing the generated data set 
solve.pbk$vec_t=rep(seq(0,8,0.01),times=66000)
solve.pbk.sa=as.data.frame(matrix(NA,1602000,2))
colnames(solve.pbk.sa)=c("time","CV")
solve.pbk.sa[,1]=solve.pbk$time
solve.pbk.sa[,2]=solve.pbk$C_V
solve.pbk.sa=solve.pbk.sa[which(solve.pbk.sa[,"time"]==0.2|solve.pbk.sa[,"time"]==0.5|solve.pbk.sa[,"time"]==1|solve.pbk.sa[,"time"]==1.5| 
                                  solve.pbk.sa[,"time"]==2|solve.pbk.sa[,"time"]==3|solve.pbk.sa[,"time"]==4|
                                  solve.pbk.sa[,"time"]==8),]
SimRes = as.data.frame(matrix(NA,66000,8))

tab1=solve.pbk.sa[which(solve.pbk.sa[,"time"]==0.2),]
tab2=solve.pbk.sa[which(solve.pbk.sa[,"time"]==0.5),]
tab3=solve.pbk.sa[which(solve.pbk.sa[,"time"]==1),]
tab4=solve.pbk.sa[which(solve.pbk.sa[,"time"]==1.5),]
tab5=solve.pbk.sa[which(solve.pbk.sa[,"time"]==2),]
tab6=solve.pbk.sa[which(solve.pbk.sa[,"time"]==3),]
tab7=solve.pbk.sa[which(solve.pbk.sa[,"time"]==4),]
tab8=solve.pbk.sa[which(solve.pbk.sa[,"time"]==8),]

SimRes[,1]=tab1[,2]
SimRes[,2]=tab2[,2]
SimRes[,3]=tab3[,2]
SimRes[,4]=tab4[,2]
SimRes[,5]=tab5[,2]
SimRes[,6]=tab6[,2]
SimRes[,7]=tab7[,2]
SimRes[,8]=tab8[,2]

#Sobol analysis plot blood
t_A<-(c(0.2,0.5,1,1.5,2,3,4,8))
par(mfrow=c(1,1), las=3, cex=0.7)
FOI          = TI          = TI.borninf           = TI.bornsup          = matrix(NA, nrow = par_var, ncol = length(t_A))  
rownames(FOI)= rownames(TI)= rownames(TI.borninf) = rownames(TI.bornsup)= colnames

t_SA <- 4

for(i in 1:length(t_A)){
  print(i)
  if (t_A[i] %in% t_SA) {
    tell(sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
    FOI[,i]       = sa$S[,1]    #First order indices
    TI[,i]        = sa$T[,1]    #Total indices
    TI.borninf[,i] = sa$T[,4]   #Lower CL total indices
    TI.bornsup[,i] = sa$T[,5]   #Upper CL total indices
   
    plot(sa, main=colnames(SimRes)[i],las=3, cex=0.7)
  }
}



