#Combined data processing of Human Cinnamaldehyde pbk models
#Before running code first run all model associated with this data processing
#Human desolve model
#human rxode model
#Human population model


#Rxode
#cinnamaldehyde model 
#Mass balance calculation rxode inhalation complete
mass_df <- solve.pbk_nonpop/BW * MW /1e+3
mass_df <- mass_df[,c(70:82,84:89,91:95,97:103)]
mass_at_t <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])

#cinnamaldehyde model 
#Mass balance calculation rxode inhalation incomplete
mass_df <- solve.pbk_nonpop/BW * MW /1e+3
mass_df <- mass_df[,c(65:76,78,79,80,82, 83,85:89,92:97)]
mass_at_t <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])

#Cinnamylalcohol model

#Mass balance calculation rxode non metabolism
Mass_df <- solve.pbk_nonpop/BW * MW /1e+3
mass_df <-Mass_df[,c(45:53,55,56,58:63)]
mass_at_t <- data.frame(mass=as.numeric())


for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])

#Rxode 
#Normal model
mass_df <- solve.pbk_nonpop/BW * MW /1e+3
mass_df <- mass_df[,c(59:68,70:74,77:81,83:89)]
mass_at_t <- data.frame(mass=as.numeric())


for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])
#Rxode data visualisation 
pL_GSH = ggplot(solve.pbk_nonpop, aes(time, AM_Lc_GSH)) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of GSH in the liver")
pL_GSH + scale_y_log10()

pA_L = ggplot(solve.pbk_nonpop, aes(time, A_L )) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of Cinnamaldehyde in the liver")
pA_L 

pA_SP = ggplot(solve.pbk_nonpop, aes(time, A_SP )) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of Cinnamaldehyde in Slowely perfused tissue")
pA_SP

pA_GST = ggplot(solve.pbk_nonpop, aes(time, AM_L_AG_GST )) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of Cinnamaldehyde metabolised with GSH")
pA_GST

#Rxode
#Tissue results 
tissue_results <- solve.pbk_nonpop[,c(1,3,4,7,13,19,25,39)]
pC_V = ggplot(tissue_results, aes(time, C_V )) + 
  geom_line() +
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Cinnamaldehyde concentration in venous blood")
pC_V

pC_A = ggplot(tissue_results, aes(time, C_A )) + 
  geom_line() +
  labs(x = "Time in hours", y = "umol/L") +
  ggtitle("Cinnamaldehyde concentration arterial blood")
pC_A

pC_L = ggplot(tissue_results, aes(time, C_L )) + 
  geom_line() +
labs(x = "Time in hours", y = "umol/L") +
  ggtitle("Cinnamaldehyde concentration in the liver")
pC_L

pC_SI = ggplot(tissue_results, aes(time, C_SI )) + 
  geom_line() +
labs(x = "Time in hours", y = "umol/L") +
  ggtitle("Cinnamaldehyde concentration in the small intestine")
pC_SI
pC_F = ggplot(tissue_results, aes(time, C_F )) + 
  geom_line() +
labs(x = "Time in hours", y = "umol/L") +
  ggtitle("Cinnamaldehyde concentration in fat")
pC_F

#Comparison between human desolve and rxode model

combined_liver <- data.frame(solve.pbk_nonpop["time"],solve.pbk_nonpop["A_L"],df_pbk_results['A_L'])
colnames(combined_liver)=c("time","Human_desolve","Human_rxode")       #Add column names
c_liver <- melt(data = combined_liver, id.vars="time", variable.name= "lever", variable.value = "concentratie")



Pcl <- ggplot(c_liver, aes(time, y = value, col=lever )) +
  geom_line()+
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of Cinnamaldehyde in the liver")
Pcl

#combined figure of all three models liver outputs
combined_liver <- data.frame(tab_C_L["time"],tab_C_L[2],tab_C_L[3],tab_C_L[4],tab_C_L_popgen[2],tab_C_L_popgen[3],tab_C_L_popgen[4], df_pbk_results['A_L'],solve.pbk_nonpop["A_L"])
colnames(combined_liver)=c("time","C_L_P2.5","C_L_P50","C_L_P97.5","C_L_P2.5_pop","C_L_P50_pop","C_L_P97.5_pop","Human_desolve","Human_rxode")       #Add column names

ggL <- ggplot(combined_liver)+
  geom_line(aes(x=time, y=C_L_P2.5), linetype = "dashed")+
  geom_line(aes(x=time, y=C_L_P50), color = "red", size = 1)+
  geom_line(aes(x=time, y=C_L_P97.5), linetype = "dashed")+
  geom_line(aes(x=time, y=C_L_P2.5_pop), linetype = "dashed")+
  geom_line(aes(x=time, y=C_L_P50_pop), color = "purple", size = 1)+
  geom_line(aes(x=time, y=C_L_P97.5_pop), linetype = "dashed")+
  geom_line(aes(x=time, y=Human_desolve),color="green", size = 1)+
  geom_line(aes(x=time, y=Human_rxode),color="blue",size=1)+
  labs(y = "liver concentration ",
       x = "Time (h)")  +
  theme_classic()

ggL+ scale_y_log10()

#combined figure of all three models venous blood outputs
combined_venous <- data.frame(tab_C_V["time"],tab_C_V[2],tab_C_V[3],tab_C_V[4],df_pbk_results['A_V'],solve.pbk_nonpop["A_V"])
colnames(combined_venous)=c("time","C_L_P2.5","C_L_P50","C_L_P97.5","Human_desolve","Human_rxode")       #Add column names

ggV <- ggplot(combined_venous)+
  geom_line(aes(x=time, y=C_L_P2.5), linetype = "dashed")+
  geom_line(aes(x=time, y=C_L_P50), color = "red", size = 1)+
  geom_line(aes(x=time, y=C_L_P97.5), linetype = "dashed")+
  geom_line(aes(x=time, y=Human_desolve),color="green", size = 1)+
  geom_line(aes(x=time, y=Human_rxode),color="blue",size=1)+
  labs(y = "Venous blood concentration ",
       x = "Time (h)")  +
  theme_classic()

ggV+ scale_y_log10()


#----paramter histograms-----#

#bw
hist(phys$BW)
hist(var_m$BW)
hist(var_f$BW)

mean(phys$BW)
average(var_f$BW)


#height

hist(phys$Height,breaks = 100)
hist(var_f$Height, breaks = 100)
hist(var_m$Height,breaks = 100)

hist(phys$Height)
hist(var_f$Height)
hist(var_m$Height)


hist(phys$Height,breaks = 100)
hist(var_f$Height, breaks = 100)
hist(var_m$Height,breaks = 100)

hist(phys$Height)
hist(var_f$Height)
hist(var_m$Height)


mean(phys$Height)
mean(var_f$Height)
mean(var_m$Height)




mean(phys$Q_C)
mean(var_m$Q_C)
mean(var_f$Q_C)

hist(phys$Q_C,breaks = 100)
hist(var_m$Q_C,breaks = 100)
hist(var_f$Q_C,breaks = 100)
#loading data file popgen
popgen.data <- read.csv("~/Biologie/Major internship/popgen data.csv")
popgen.data$Age = as.numeric(popgen.data$Age)
popgen.data$Cardiac.Output = as.numeric(popgen.data$Cardiac.Output)


#Body weight comparison
hist(phys$BW,breaks = 100)
hist(popgen.data$Body.Mass[0:2000],breaks = 100)

#Cardiac output comparison
hist(popgen.data$Cardiac.Output,breaks = 100)
hist(phys$Q_C,breaks = 100)

mean(popgen.data$Cardiac.Output[0:2000])
mean(phys$Q_C)

#Liver tissue blood flow
hist(phys$Q_L,breaks = 100)
hist(popgen.data$Liver.flow[0:2000] , breaks = 100)
mean(popgen.data$Liver.flow[0:2000])
mean(phys$Q_L)

#liver mass
hist(phys$V_L,breaks = 100)
hist(popgen.data$Liver.mass[0:2000] , breaks = 100)
mean(popgen.data$Liver.mass[0:2000])
mean(phys$V_L)


#RP tissue blood flow  
hist(phys$Q_RP,breaks = 100)
hist(popgen.data$Richly.Perfused.flow[0:2000] , breaks = 100)
mean(popgen.data$Richly.Perfused.flow[0:2000])
mean(phys$Q_RP)


#RP tissue mass   
hist(phys$V_RP,breaks = 100)
hist(popgen.data$Richly.Perfused.mass[0:2000] , breaks = 100)
mean(popgen.data$Richly.Perfused.mass[0:2000])
mean(phys$V_RP)


#Population based model data visualisation and analysis
tab_solve_C_V=as.data.frame(matrix(NA,time.end/time.frame+1,(N+NF)))    #Create an empty data frame with amount of timepoints=amount of rows and amount of individuals=amount of columns
for (i in 1:(N+NF)) {
  tab.i=solve.pbk[which(solve.pbk[,"sim.id"]==i),]                  #Put all individuals in data frame
  tab.i=as.data.frame(tab.i)
  tab_solve_C_V[,i]=tab.i$C_V
}

tab_C_V=as.data.frame(matrix(NA,time.end/time.frame+1,4))       #Create an empty data frame with amount of timepoints=amount of rows and 4 columns
tab_C_V[,1]=c(seq(time.0,time.end,by=time.frame))               #Timepoints in first column
for (i in 1:(time.end/time.frame+1)) {
  tab_C_V[i,2]=quantile(tab_solve_C_V[i,],0.025, na.rm = TRUE)      #Lower bound of confidence interval in second column
  tab_C_V[i,3]=quantile(tab_solve_C_V[i,],0.5, na.rm = TRUE)        #Median in third column
  tab_C_V[i,4]=quantile(tab_solve_C_V[i,],0.975, na.rm = TRUE)      #Upper bound of confidence interval in fourth column
}
colnames(tab_C_V)=c("time","CV_P2.5","CV_P50","CV_P97.5")       #Add column names

gg <- ggplot(tab_C_V)+
  geom_line(aes(x=time, y=CV_P2.5), linetype = "dashed")+
  geom_line(aes(x=time, y=CV_P50), color = "red", size = 1)+
  geom_line(aes(x=time, y=CV_P97.5), linetype = "dashed")+
  labs(y = "Blood concentration ",
       x = "Time (h)")  +
  theme_classic()

gg

tab_solve_C_L=as.data.frame(matrix(NA,time.end/time.frame+1,(N+NF)))    #Create an empty data frame with amount of timepoints=amount of rows and amount of individuals=amount of columns
for (i in 1:(N+NF)) {
  tab.i=solve.pbk[which(solve.pbk[,"sim.id"]==i),]                  #Put all individuals in data frame
  tab.i=as.data.frame(tab.i)
  tab_solve_C_L[,i]=tab.i$C_L
}


tab_C_L=as.data.frame(matrix(NA,time.end/time.frame+1,4))       #Create an empty data frame with amount of timepoints=amount of rows and 4 columns
tab_C_L[,1]=c(seq(time.0,time.end,by=time.frame))               #Timepoints in first column
for (i in 1:(time.end/time.frame+1)) {
  tab_C_L[i,2]=quantile(tab_solve_C_L[i,],0.025, na.rm = TRUE)      #Lower bound of confidence interval in second column
  tab_C_L[i,3]=quantile(tab_solve_C_L[i,],0.5, na.rm = TRUE)        #Median in third column
  tab_C_L[i,4]=quantile(tab_solve_C_L[i,],0.975, na.rm = TRUE)      #Upper bound of confidence interval in fourth column
}

colnames(tab_C_L)=c("time","C_L_P2.5","C_L_P50","C_L_P97.5")


#Rat model calc
pL_GSH = ggplot(solve.pbk, aes(time, AM_Lc_GSH)) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of GSH in the liver")
pL_GSH + scale_y_log10()

pA_L = ggplot(solve.pbk, aes(time, A_L )) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of Cinnamaldehyde in the liver")
pA_L 

pA_V = ggplot(solve.pbk, aes(time, A_V )) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of Cinnamaldehyde blood")
pA_V 


data_250mg <- read_excel("D:/Joris/Toxicology and Environmental Health/Master stage/Comparison data/data cnma in blood ug 250mg-kg dose .xlsx")

blood_amount_total <- solve.pbk[,1]
blood_amount_total <- cbind(blood_amount_total,(rowSums (solve.pbk[,44:45])))
colnames(blood_amount_total) <- c("time","concentration")
blood_amount_total <- as.data.frame(blood_amount_total)

comparison_data <- merge.data.frame(blood_amount_total, data_250mg, by="time", all= "TRUE")

pBlood = ggplot(blood_amount_total, aes(time, concentration )) + 
  geom_line() + 
  labs(x = "Time in hours", y = "umol") +
  ggtitle("Amount of Cinnamaldehyde blood")
pBlood  + scale_y_log10()

fig_blood <- plot_ly(blood_amount_total, x=~time, y =~concentration , type = 'scatter', mode = 'lines')%>%
  add_trace(blood_amount_total, x=~time, y =~"rat 1" , type = 'scatter', mode = 'markers')


fig_blood 


