phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[1:2500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=1:2500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=1:2500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=1:2500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=1:2500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=1:2500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\\\solve.pbk_nonpop1", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[2501:5000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=2501:5000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=2501:5000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=2501:5000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=2501:5000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=2501:5000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\\\solve.pbk_nonpop2", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[5001:7500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=5001:7500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=5001:7500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=5001:7500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=5001:7500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=5001:7500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop3", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[7501:10000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=7501:10000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=7501:10000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=7501:10000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=7501:10000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=7501:10000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop4", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[10001:12500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde
Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw

ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=10001:12500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=10001:12500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=10001:12500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=10001:12500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=10001:12500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop5", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[12501:15000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=12501:15000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=12501:15000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=12501:15000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=12501:15000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=12501:15000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop6", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[15001:17500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=15001:17500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=15001:17500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=15001:17500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=15001:17500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=15001:17500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop7", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[17501:20000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=17501:20000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=17501:20000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=17501:20000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=17501:20000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=17501:20000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop8", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[20001:22500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=20001:22500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=20001:22500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=20001:22500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=20001:22500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=20001:22500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop9", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[22501:25000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=22501:25000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=22501:25000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=22501:25000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=22501:25000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=22501:25000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop10", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[25001:27500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=25001:27500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=25001:27500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=25001:27500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=25001:27500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=25001:27500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop11", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[27501:30000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=27501:30000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=27501:30000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=27501:30000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=27501:30000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=27501:30000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop12", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[30001:32500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=30001:32500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=30001:32500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=30001:32500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=30001:32500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=30001:32500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop13", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[32501:35000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=32501:35000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=32501:35000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=32501:35000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=32501:35000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=32501:35000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop14", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[35001:37500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=35001:37500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=35001:37500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=35001:37500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=35001:37500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=35001:37500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop15", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[37501:40000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=37501:40000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=37501:40000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=37501:40000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=37501:40000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=37501:40000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop16", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[40001:42500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=40001:42500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=40001:42500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=40001:42500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=40001:42500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=40001:42500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop17", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[42501:45000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=42501:45000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=42501:45000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=42501:45000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=42501:45000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=42501:45000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop18", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[45001:47500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=45001:47500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=45001:47500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=45001:47500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=45001:47500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=45001:47500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop19", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[47501:50000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=47501:50000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=47501:50000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=47501:50000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=47501:50000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=47501:50000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop20", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[50001:52500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=50001:52500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=50001:52500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=50001:52500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=50001:52500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=50001:52500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop21", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[52501:55000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=52501:55000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=52501:55000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=52501:55000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=52501:55000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=52501:55000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop22", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[55001:57500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=55001:57500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=55001:57500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=55001:57500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=55001:57500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=55001:57500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop23", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[57501:60000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=57501:60000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=57501:60000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=57501:60000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=57501:60000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=57501:60000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop24", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[60001:62500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=60001:62500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=60001:62500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=60001:62500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=60001:62500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=60001:62500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model

write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop25", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[62501:65000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=62501:65000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=62501:65000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=62501:65000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=62501:65000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=62501:65000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop26", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[65001:67500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=65001:67500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=65001:67500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=65001:67500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=65001:67500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=65001:67500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop27", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[67501:70000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=67501:70000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=67501:70000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=67501:70000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=67501:70000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=67501:70000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop28", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[70001:72500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=70001:72500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=70001:72500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=70001:72500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=70001:72500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=70001:72500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop29", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[72501:75000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=72501:75000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=72501:75000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=72501:75000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=72501:75000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=72501:75000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop30", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[75001:77500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=75001:77500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=75001:77500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=75001:77500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=75001:77500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=75001:77500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop31", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[77501:80000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=77501:80000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=77501:80000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=77501:80000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=77501:80000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=77501:80000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop32", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[80001:82500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=80001:82500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=80001:82500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=80001:82500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=80001:82500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=80001:82500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop33", row.names = TRUE)



phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[82501:85000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=82501:85000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=82501:85000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=82501:85000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=82501:85000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=82501:85000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop34", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[85001:87500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=85001:87500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=85001:87500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=85001:87500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=85001:87500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=85001:87500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop35", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[87501:90000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=87501:90000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=87501:90000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=87501:90000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=87501:90000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=87501:90000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop36", row.names = TRUE)



phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[90001:92500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=90001:92500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=90001:92500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=90001:92500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=90001:92500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=90001:92500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop37", row.names = TRUE)




phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[92501:95000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=92501:95000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=92501:95000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=92501:95000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=92501:95000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=92501:95000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop38", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[95001:97500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=95001:97500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=95001:97500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=95001:97500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=95001:97500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=95001:97500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop39", row.names = TRUE)




phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[97501:100000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=97501:100000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=97501:100000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=97501:100000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=97501:100000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=97501:100000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop40", row.names = TRUE)




phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[100001:102500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=100001:102500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=100001:102500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=100001:102500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=100001:102500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=100001:102500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop41", row.names = TRUE)

phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[102501:105000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=102501:105000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=102501:105000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=102501:105000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=102501:105000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=102501:105000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop42", row.names = TRUE)



phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[105001:107500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=105001:107500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=105001:107500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=105001:107500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=105001:107500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=105001:107500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop43", row.names = TRUE)



phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[107501:110000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=107501:110000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=107501:110000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=107501:110000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=107501:110000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=107501:110000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop44", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[110001:112500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=110001:112500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=110001:112500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=110001:112500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=110001:112500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=110001:112500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop45", row.names = TRUE)



phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[112501:115000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=112501:115000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=112501:115000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=112501:115000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=112501:115000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=112501:115000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop46", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[115001:117500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=115001:117500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=115001:117500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=115001:117500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=115001:117500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=115001:117500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop47", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[117501:120000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=117501:120000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=117501:120000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=117501:120000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=117501:120000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=117501:120000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop48", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[120001:122500,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=120001:122500,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=120001:122500,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=120001:122500,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=120001:122500,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=120001:122500,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop49", row.names = TRUE)



phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[122501:125000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=122501:125000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=122501:125000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=122501:125000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=122501:125000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=122501:125000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop50", row.names = TRUE)


phys <- read_csv("phys_inhalation_pop")

#Running the Global SA directly takes to much memory so it is necessary to split up the data set in parts
phys1<-phys[125001:126000,]
P_F<-phys1$P_F
P_L<-phys1$P_L
P_SI<-phys1$P_SI
P_RP<-phys1$P_RP
P_B<-phys1$P_B
P_SP<-phys1$P_SP
P_Pu<-phys1$P_Pu
P_OH_F<-phys1$P_OH_F
P_OH_L<-phys1$P_OH_L
P_OH_SI<-phys1$P_OH_SI
P_OH_RP<-phys1$P_OH_RP
P_OH_SP<-phys1$P_OH_SP
P_OH_Pu<-phys1$P_OH_Pu
BW<-phys1$BW
V_L<-phys1$V_L
V_F<-phys1$V_F
V_A<-phys1$V_A
V_V<-phys1$V_V
V_SI<-phys1$V_SI
V_Pu<-phys1$V_Pu
V_RP<-phys1$V_RP
V_SP<-phys1$V_SP
Q_C<-phys1$Q_C
Q_SI<-phys1$Q_SI
Q_F<-phys1$Q_F
Q_L<-phys1$Q_L
Q_Pu<-phys1$Q_Pu
Q_RP<-phys1$Q_RP
Q_SP<-phys1$Q_SP
P_V<-phys1$P_V
G_SYN_L<-phys1$G_SYN_L
G_SYN_SI<-phys1$G_SYN_SI
k_L_GLOS<-phys1$k_L_GLOS
k_SI_GLOS<-phys1$k_SI_GLOS
init_GSH_L<-phys1$init_GSH_L
init_GSH_SI<-phys1$init_GSH_SI
k_GSH<-phys1$k_GSH
k_DNA<-phys1$k_DNA
C_PRO_L<-phys1$C_PRO_L
C_PRO_SI<-phys1$C_PRO_SI
C_L_dG<-phys1$C_L_dG
T_0.5<-phys1$T_0.5
Ka<-phys1$Ka
k_L_OH <- phys1$k_L_OH
Km_L_CA<-phys1$Km_L_CA
Km_L_AO<-phys1$Km_L_AO
Km_L_GST<-phys1$Km_L_GST
Km_L_GST_G<-phys1$Km_L_GST_G
Vsmax_L_CA<-phys1$Vsmax_L_CA
Vsmax_L_AO<-phys1$Vsmax_L_AO
Vsmax_L_GST<-phys1$Vsmax_L_GST
Km_SI_CA<-phys1$Km_SI_CA
Km_SI_AO<-phys1$Km_SI_AO
Km_SI_OH<-phys1$Km_SI_OH
Km_SI_GST<-phys1$Km_SI_GST
Km_SI_GST_G<-phys1$Km_SI_GST_G
Vsmax_SI_CA<-phys1$Vsmax_SI_CA
Vsmax_SI_AO<-phys1$Vsmax_SI_AO
Vsmax_SI_OH<-phys1$Vsmax_SI_OH
Vsmax_SI_GST<-phys1$Vsmax_SI_GST
Volume_exposure_chamber=phys1$Volume_exposure_chamber

parameters1 <- cbind(P_F,
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

#exposure
amount.units               <-"umol"
time.units                 <-"h"
nbr.doses                  <-1        #number of doses
time.0                     <-0        #time start dosing
time.end                   <-8        #time end of simulation
time.frame                 <-0.1     #time steps of simulation
MW                         <-132.16   #The molecular weight of Cinnamaldehyde

Inhalation_Dose_in_mg_bw   <-100        #The inhaled dose in mg/kg-bw
Oral_Dose_in_mg_bw         <-0      #Dose in mg/kg-bw




ex1 <- eventTable(amount.units = amount.units, time.units = time.units) %>%
  et(id=125001:126000,seq(from = time.0, to = time.end, by = time.frame))%>%
  et(id=125001:126000,amt=(Oral_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3  , dur=0.01, cmt="A_GI", nbr.doses=nbr.doses)%>%
  et(id=125001:126000,amt=(Inhalation_Dose_in_mg_bw) * phys1$BW/ MW  * 1e+3 , dur=0.01, cmt="A_Inhalation", nbr.doses=nbr.doses)%>%
  et(id=125001:126000,amt=phys1$init_GSH_SI, dur=0.01, cmt="AM_SIc_GSH", nbr.doses=1)%>%
  et(id=125001:126000,amt=phys1$init_GSH_L, dur=0.01, cmt="AM_Lc_GSH", nbr.doses=1)

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
           "AM_Lc_GSH"    =0, 
           "AM_SI_CA"     =0,
           "AM_SI_AO"     =0,
           "AM_SI_AG_GST" =0,
           "AM_SI_AG_CHEM"=0,
           "AM_SI_AP"     =0,
           "A_OH_M_SI_C_A"=0,
           "A_OH_SI"      =0,
           "A_SI"         =0,
           "AM_SIc_GSH"   =0,
           "A_RP"         =0,
           "A_OH_RP"      =0,
           "A_SP"         =0,
           "A_OH_SP"      =0
);



#Run the model after assigning the sobol dataset to the variables 
solve.pbk_nonpop1 <- solve(PBK_Cinnamaldehyde, parameters1, events = ex1, inits) #Solve the PBPK model


write.csv(solve.pbk_nonpop1,"D:/PBK/Cinnamaldehyde-pbk\\solve.pbk_nonpop51", row.names = TRUE)

