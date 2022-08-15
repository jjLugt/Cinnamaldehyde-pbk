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
colnames <- c(colnames(parameters))
par_var <- length(colnames)

#this is stupid work on this!! mean for the data sets is trivial to calculate
Mean <- parameters

Lower <- Mean - 0.1 * Mean
Upper <- Mean + 0.1 * Mean

#Lower <- Mean - 0.2 * Mean
#Upper <- Mean + 0.2 * Mean


#create data frames for population
n_sim  <- 100               #number of iterations
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
sa <- soboljansen(model=NULL, X1, X2, nboot = n_boot, conf = 0.95, events = ex)
phys <- sa$X


P_F<-phys$P_F
P_L<-phys$P_L
P_SI<-phys$P_SI
P_RP<-phys$P_RP
P_B<-phys$P_B
P_SP<-phys$P_SP
P_Pu<-phys$P_Pu
P_OH_F<-phys$P_OH_F
P_OH_L<-phys$P_OH_L
P_OH_SI<-phys$P_OH_SI
P_OH_RP<-phys$P_OH_RP
P_OH_SP<-phys$P_OH_SP
P_OH_Pu<-phys$P_OH_Pu
Age <- phys$Age
Height <- phys$Height
BW<-phys$BW
V_L<-phys$V_L
V_F<-phys$V_F
V_B <- phys$V_B
V_A<-phys$V_A
V_V<-phys$V_V
V_SI<-phys$V_SI
V_Pu<-phys$V_Pu
V_RP<-phys$V_RP
V_SP<-phys$V_SP
Q_C<-phys$Q_C
Q_SI<-phys$Q_SI
Q_F<-phys$Q_F
Q_L<-phys$Q_L
Q_Pu<-phys$Q_Pu
Q_RP<-phys$Q_RP
Q_SP<-phys$Q_SP
P_V<-phys$P_V
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
Oral_Dose<- phys$Oral_Dose 
Inhalation_Dose<-phys$Inhalation_Dose
Volume_exposure_chamber<-phys$Volume_exposure_chamber

parameters=cbind(P_F,
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
                 Age,
                 Height,
                 BW,
                 V_F,
                 V_L,
                 V_SI,
                 V_B,
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
                 Oral_Dose,
                 Inhalation_Dose,
                 Volume_exposure_chamber)

#defining the begin situation of the model Inhalation variation 
inits <- c("A_GI"         =0,
           "A_P_Art"      =0,
           "A_Inhalation" =0,
           "A_Exhalation" =0,
           "A_Pu"         =0,
           "A_OH_Pu"      =0,
           "A_V"          =0,
           "A_OH_V"       =0,
           "A_F"          =0,
           "A_OH_F"       =0,
           "AM_L_CA"      =0,
           "AM_L_AO"      =0,
           "AM_L_AG_GST"  =0,
           "AM_L_AG_CHEM" =0,
           "AM_L_AP"      =0,
           "AM_L_DA_FORM" =0,
           "AM_L_DA"      =0,
           "A_OH_M_L_C_A" =0,
           "A_OH_L"       =0,
           "A_L"          =0,
           "AM_Lc_GSH"    =init_GSH_L, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =init_GSH_SI,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);


#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop <- solve(PBK_Cinnamaldehyde, parameters, events = ex, inits, cores=6) #Solve the PBPK model


#Analysing the generated data set 
solve.pbk_nonpop$vec_t=rep(seq(0,8,0.1),times=6500)
solve.pbk.sa=as.data.frame(matrix(NA,526500,2))
colnames(solve.pbk.sa)=c("time","C_V")
solve.pbk.sa[,1]=solve.pbk_nonpop$time
solve.pbk.sa[,2]=solve.pbk_nonpop$C_V
solve.pbk.sa=solve.pbk.sa[which(solve.pbk.sa[,"time"]==0.2|solve.pbk.sa[,"time"]==0.5|solve.pbk.sa[,"time"]==1|solve.pbk.sa[,"time"]==1.5| 
                                  solve.pbk.sa[,"time"]==2|solve.pbk.sa[,"time"]==3|solve.pbk.sa[,"time"]==4|
                                  solve.pbk.sa[,"time"]==8),]
SimRes = as.data.frame(matrix(NA,6500,8))

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

dev.off()
FOI.L = as.matrix(FOI[,1:length(t_A)])       # as.matrix
TI.L  = as.matrix(TI[,1:length(t_A)])

FOI.L.t <- apply(FOI.L, 1, mean, na.rm=TRUE)
TI.L.t <- apply(TI.L, 1, mean, na.rm=TRUE)

sorting = order(TI.L.t, decreasing = F)
TI.L.t  = TI.L.t[sorting]
FOI.L.t = FOI.L.t[sorting]

FOI.L.t = ifelse(FOI.L.t <= 0, 0, FOI.L.t)
tempC    = t(cbind(FOI.L.t, TI.L.t))

tempC2 <- as.data.frame(tempC[,c(54:63)])

par(mfrow=c(1,1), las=1, mai=c(0.35,1,0.35,0.1), mgp = c(3.5,0.5,0))
colnames(tempC2) <- c("Q_SI", "Ka", "V_SP", "Q_RP", "P_SP", "Q_SP", "k_GSH", "C_PRO_L", "VL", "QC")
O_CV_0.2 <- barplot(as.matrix(tempC2), col=c("firebrick1","firebrick4"), horiz = T, beside =T , main="", cex.lab=1.5 , xlim=c(0,1) )

write.csv(SimRes, file = "SimRes_08-15-2022_inhalation_parameters_no_inhalation.csv")
