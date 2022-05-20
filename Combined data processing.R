#Combined data processing of Human Cinnamaldehyde pbk models
#Before running code first run all model associated with this data processing
#Human desolve model
#human rxode model
#Human population model


#Rxode data manipulation
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

