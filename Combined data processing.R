#Combined data processing of Human Cinnamaldehyde pbk models
#Before running code first run all model associated with this data processing
#Human desolve model
#human rxode model
#Human population model

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
combined_liver <- data.frame(tab_C_L["time"],tab_C_L[2],tab_C_L[3],tab_C_L[4],df_pbk_results['A_L'],solve.pbk_nonpop["A_L"])
colnames(combined_liver)=c("time","C_L_P2.5","C_L_P50","C_L_P97.5","Human_desolve","Human_rxode")       #Add column names

ggL <- ggplot(combined_liver)+
  geom_line(aes(x=time, y=C_L_P2.5), linetype = "dashed")+
  geom_line(aes(x=time, y=C_L_P50), color = "red", size = 1)+
  geom_line(aes(x=time, y=C_L_P97.5), linetype = "dashed")+
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
hist(phys$Height)
hist(var_f$Height)
hist(var_m$Height)

mean(phys$Height)
mean(var_f$Height)
mean(var_m$Height)
