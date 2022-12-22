#Comparison to in vivo data 
#author: Joris Jean van der Lugt
#date: 17-11-2022

library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(truncnorm)
library(reshape2)
library(plotly)
library(PKNCA)
library(bio3d)

#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_V + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]


#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


Combined_data_file_for_graph_500mg <- read_delim("Combined data file for graph 500mg.csv", 
                                                 delim = ";", escape_double = FALSE, trim_ws = TRUE)
colnames(Combined_data_file_for_graph_500mg) <- c("time","umol.l","ID")

p1 <- plot_ly(Combined_data_file_for_graph_500mg, x=~time, y=Combined_data_file_for_graph_500mg$'umol.l', 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 500mg/kg-bw',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in μmol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1

g <- ggplot(Combined_data_file_for_graph_500mg,aes(time,Combined_data_file_for_graph_500mg$'umol.l',color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+
  labs(subtitle="Oral dose 500mg/kg-BW", 
       y="Cinnamaldehyde concentration (μmol/L)", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")
ggsave(plot=g,"CNMA in blood 500mg oral comparison.png",
       width= 11.69, height= 8.3, dpi= 250)

RAT_data_obs_1 <- Combined_data_file_for_graph_500mg[1:15,]
RAT_data_obs_2 <- Combined_data_file_for_graph_500mg[16:30,]
RAT_data_obs_3 <- Combined_data_file_for_graph_500mg[31:46,]
RAT_data_Zao   <- Combined_data_file_for_graph_500mg[47:57,]
SIM_data_pred    <-Combined_data_file_for_graph_500mg[c(79,89,100,114,124,132,150,169,190,210,231,252,273,302,322),]
SIM_data_pred_kiwa <- Combined_data_file_for_graph_500mg[c(380,390,401,415,425,433,451,470,491,511,532,553,574,603,623),]
SIM_data_ka    <-Combined_data_file_for_graph_500mg[c(681,691,702,716,726,734,752,771,792,812,833,854,875,904,924),]

SIM_data_pred[,4]<- as.data.frame(RAT_data_obs_1[,2])
SIM_data_pred[,5]<- as.data.frame(RAT_data_obs_2[,2])
SIM_data_pred[,6]<- as.data.frame(RAT_data_obs_3[1:15,2])
SIM_data_pred[,7]<- as.data.frame(SIM_data_pred_kiwa[,2])
SIM_data_pred[,8]<- as.data.frame(SIM_data_ka[,2])

colnames(SIM_data_pred)<- c("Time","sim","ID","rat_1","rat_2","rat_3","SIM-Kiwa","ka")


#sim rat 1
rsmesim_rat1<-sqrt(mean((SIM_data_pred$sim - SIM_data_pred$rat_1)^2))

#sim rat 2
rsmesim_rat2<- sqrt(mean((SIM_data_pred$sim - SIM_data_pred$rat_2)^2))

#sim rat 3
rsmesim_rat3<- sqrt(mean((SIM_data_pred$sim - SIM_data_pred$rat_3)^2))

mean_RSME<-(rsmesim_rat1+rsmesim_rat2+rsmesim_rat3)/3

RSMD_500mg_1<- rmsd(a=SIM_data_pred$sim, b=SIM_data_pred$rat_1)
RSMD_500mg_2<- rmsd(a=SIM_data_pred$sim, b=SIM_data_pred$rat_2)
RSMD_500mg_3<- rmsd(a=SIM_data_pred$sim, b=SIM_data_pred$rat_3)

p<-ggplot(SIM_data_pred, aes(x=SIM_data_pred$sim, y=SIM_data_pred$rat_1)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  geom_line(data=tibble(x=c(0, 60),y=c(0,60)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 60),y=c(0,60)),linetype = "dashed",
            aes(x=x,y=y*2),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 60),y=c(0,60)),linetype = "dashed",
            aes(x=x,y=y*0.5),
            inherit.aes = FALSE)+
  annotate("text", x = 10, y = 75, label = "RSMD: 21.9", size= 8)+
  labs(x='Predicted Values(μmol/L)', y='Actual Values(μmol/L)', title='Predicted vs. Actual Values 500mg/kg-BW oral dose',subtitle="Yuan et al 1992")+
  geom_point(aes(x=SIM_data_pred$sim,y=SIM_data_pred$rat_2))+
  geom_point(aes(x=SIM_data_pred$sim,y=SIM_data_pred$rat_3))+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
p
ggsave(plot=p,"Pred vs actual 500mg oral data.png",
       width= 11.69, height= 8.3, dpi= 250)

#sim rat 1
rsmesimka_rat1<-sqrt(mean((SIM_data_pred$ka - SIM_data_pred$rat_1)^2))

#sim rat 2
rsmesimka_rat2<- sqrt(mean((SIM_data_pred$ka - SIM_data_pred$rat_2)^2))

#sim rat 3
rsmesimka_rat3<- sqrt(mean((SIM_data_pred$ka - SIM_data_pred$rat_3)^2))

mean_RSME_ka<-(rsmesimka_rat1+rsmesimka_rat2+rsmesimka_rat3)/3

RSMD_500mg_ka_1<- rmsd(a=SIM_data_pred$ka, b=SIM_data_pred$rat_1)
RSMD_500mg_ka_2<- rmsd(a=SIM_data_pred$ka, b=SIM_data_pred$rat_2)
RSMD_500mg_ka_3<- rmsd(a=SIM_data_pred$ka, b=SIM_data_pred$rat_3)

#plotting new plot with adjusted parameters for ka
p<-ggplot(SIM_data_pred, aes(x=SIM_data_pred$ka, y=SIM_data_pred$rat_1)) +
  geom_point() +
  geom_abline(intercept=0, slope=1)+
  geom_line(data=tibble(x=c(0, 15),y=c(0,15)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0,15),y=c(0,15)),linetype = "dashed",
            aes(x=x,y=y*2),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0,15),y=c(0,15)),linetype = "dashed",
            aes(x=x,y=y*0.5),
            inherit.aes = FALSE)+
  labs(x='Predicted Values(μmol/L)', y='Actual Values(μmol/L)', title='Predicted vs. Actual Values 500mg/kg-BW oral dose',subtitle="Yuan et al 1992, Ka = 0.2")+
  annotate("text", x = 4, y = 20, label = "RSMD: 5.5", size= 8)+
  geom_point(aes(x=SIM_data_pred$ka,y=SIM_data_pred$rat_2))+
  geom_point(aes(x=SIM_data_pred$ka,y=SIM_data_pred$rat_3))+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
p
ggsave(plot=p,"Pred vs actual 500mg oral ka .png",
       width= 11.69, height= 8.3, dpi= 250)


#Making a pred vs pred plot

Kiwa_data      <-Combined_data_file_for_graph_500mg[c(58,60:72),]
SIM_pred_pred  <-Combined_data_file_for_graph_500mg[c(374,375,376,379,384,401,427,468,505,538,574,607,640,667),]
SIM_pred_pred[,3]<- as.data.frame(Kiwa_data[,2])
SIM_pred_pred[,4]<-SIM_pred_pred[,2]/SIM_pred_pred[,3]
SIM_pred_pred[,5]<-SIM_pred_pred[,3]/SIM_pred_pred[,2]
colnames(SIM_pred_pred)<- c("Time","sim","kiwa","sim/kiwa","kiwa/sim")


#sim ka vs kiwa
rsmesim_kiwa<-sqrt(mean((SIM_pred_pred$sim - SIM_pred_pred$kiwa)^2))


RSME_sim_kiwa<-(rsmesim_kiwa)


p<-ggplot(SIM_pred_pred, aes(x=sim, y=kiwa)) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  ylim(0,1000)+
  xlim(0,1000)+
  annotate("text", x = 200, y = 750, label = "Mean RSME: 167.03", size= 8)+
  labs(x='Predicted Values R', y='Predicted values Kiwa', title='Predicted vs. Predicted Values 500mg oral dose')+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
ggsave(plot=p,"Pred vs pred 500mg oral.png",
       width= 11.69, height= 8.3, dpi= 250)


#Calculating AUC values vor the comparison graphs made above  

AUC_data <-PKNCAconc(Combined_data_file_for_graph_500mg, umol.l~time|ID)

#Dosing data per Paper/simulation run 
#Loading dosing file
d_dose <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/500m_comparison_doses.csv", sep=";")                       

dose_obj <- PKNCAdose(d_dose, dose~time|ID)

#letting pknc chose the end time of the auc calc
data_obj_automatic <- PKNCAdata(AUC_data, dose_obj)

#Computing the data both manual and automatic
results_obj_automatic <- pk.nca(data_obj_automatic)

#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_automatic)))
summary(results_obj_automatic)

#saving data 
write.csv(results_obj_automatic$result,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//500mg_oral_results.csv")


#RAT_1 calculation
x <- c(Combined_data_file_for_graph_500mg[1:15,1])
y <- c(Combined_data_file_for_graph_500mg[1:15,2])
RAT_1_AUC <- auc(x[["time"]],y[["umol.l"]], type=c("spline"))


#RAT_2
x <- c(Combined_data_file_for_graph_500mg[16:30,1])
y <- c(Combined_data_file_for_graph_500mg[16:30,2])
RAT_2_AUC <- auc(x[["time"]],y[["umol.l"]], type=c("spline"))

#RAT_3
x <- c(Combined_data_file_for_graph_500mg[31:46,1])
y <- c(Combined_data_file_for_graph_500mg[31:46,2])
RAT_3_AUC <- auc(x[["time"]],y[["umol.l"]], type=c("spline"))

#Zhao
x <- c(Combined_data_file_for_graph_500mg[47:57,1])
y <- c(Combined_data_file_for_graph_500mg[47:57,2])
ZHAO_AUC <- auc(x[["time"]],y[["umol.l"]], type=c("spline"))

#Generating a file for the comparison with in vivo data for 250 mg exposure  
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]

#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


Combined_data_file_for_graph_250mg <- read_delim("Combined data file for graph 250mg.csv", 
                                                 delim = ";", escape_double = FALSE, trim_ws = TRUE)

p1 <- plot_ly(Combined_data_file_for_graph_250mg, x=~time, y=Combined_data_file_for_graph_250mg$`ug.ml`, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 250mg/kg-bw',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in umol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1

g <- ggplot(Combined_data_file_for_graph_250mg,aes(time,Combined_data_file_for_graph_250mg$`ug.ml`,color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+
  labs(subtitle="Oral dose 250mg/kg-BW", 
       y="Cinnamaldehyde concentration (μmol/L)", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="Cinnamaldehyde PBK model")
ggsave(plot=g,"CNMA in blood 250mg oral comparison.png",
       width= 11.69, height= 8.3, dpi= 250)

RAT_data_obs_1 <- Combined_data_file_for_graph_250mg[21:36,]
RAT_data_obs_2 <- Combined_data_file_for_graph_250mg[37:52,]
RAT_data_obs_3 <- Combined_data_file_for_graph_250mg[53:68,]
RAT_data_Zao   <- Combined_data_file_for_graph_250mg[1:8,]
SIM_data_pred  <-Combined_data_file_for_graph_250mg[c(70,71,74,80,84,89,94,104,114,124,134,144,154,164,179,189),]
SIM_data_Ka    <-Combined_data_file_for_graph_250mg[c(312,313,316,322,326,331,336,346,356,366,376,387,396,411,421,431),]


SIM_data_pred[,4]<- as.data.frame(RAT_data_obs_1[,2])
SIM_data_pred[,5]<- as.data.frame(RAT_data_obs_2[,2])
SIM_data_pred[,6]<- as.data.frame(RAT_data_obs_3[,2])
SIM_data_pred[,7]<- as.data.frame(SIM_data_Ka[,2])


#Calculating residual values
SIM_data_pred[,8] <- SIM_data_pred[,2]-SIM_data_pred[,4]
SIM_data_pred[,9] <- SIM_data_pred[,4]-SIM_data_pred[,7]
SIM_data_pred[,10] <- SIM_data_pred[,5]-SIM_data_pred[,7]


colnames(SIM_data_pred)<- c("Time","sim","ID","rat_1", "rat_2","rat_3","SIM_ka","Residual_Rat_1","Residual_Rat1_ka","Residual_Rat2_ka")



#sim rat 1
rsmesim_rat1_250<-sqrt(mean((SIM_data_pred$sim - SIM_data_pred$rat_1)^2))

#sim rat 2
rsmesim_rat2_250<- sqrt(mean((SIM_data_pred$sim - SIM_data_pred$rat_2)^2))

#sim rat 3
rsmesim_rat3_250<- sqrt(mean((SIM_data_pred$sim - SIM_data_pred$rat_3)^2))

mean_RSME_250<-mean(c(rsmesim_rat1_250,rsmesim_rat2_250,rsmesim_rat3_250))

RSMD_250mg_1<- rmsd(a=SIM_data_pred$sim, b=SIM_data_pred$rat_1)
RSMD_250mg_2<- rmsd(a=SIM_data_pred$sim, b=SIM_data_pred$rat_2)
RSMD_250mg_3<- rmsd(a=SIM_data_pred$sim, b=SIM_data_pred$rat_3)

predvsout_250mg <- ggplot(SIM_data_pred, aes(x=SIM_data_pred$sim, y=SIM_data_pred$rat_1)) +
  geom_point() +
  #geom_abline(intercept=0, slope=1)+
  geom_line(data=tibble(x=c(0, 30),y=c(0,30)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 30),y=c(0,30)),linetype = "dashed",
            aes(x=x,y=y*2),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 30),y=c(0,30)),linetype = "dashed",
            aes(x=x,y=y*0.5),
            inherit.aes = FALSE)+
  labs(x='Predicted Values(μmol/L)', y='Actual Values(μmol/L)', title='Predicted vs. Actual Values 250mg/kg-BW oral',subtitle="Yuan et al 1992")+
  annotate("text", x = 10, y = 40, label = "RSMD: 18.9", size= 8)+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))+
  geom_point(aes(x=SIM_data_pred$sim,y=SIM_data_pred$rat_2))+
geom_point(aes(x=SIM_data_pred$sim,y=SIM_data_pred$rat_3))
predvsout_250mg

ggsave(plot=predvsout_250mg,"Pred vs data 250mg oral.png",
       width= 11.69, height= 8.3, dpi= 250)

predvsout_250mg


#sim rat 1
rsmesim_rat1_250ka<-sqrt(mean((SIM_data_pred$SIM_ka - SIM_data_pred$rat_1)^2))

#sim rat 2
rsmesim_rat2_250ka<- sqrt(mean((SIM_data_pred$SIM_ka - SIM_data_pred$rat_2)^2))

#sim rat 3
rsmesim_rat3_250ka<- sqrt(mean((SIM_data_pred$SIM_ka - SIM_data_pred$rat_3)^2))

mean_RSME_250ka<-mean(c(rsmesim_rat1_250ka,rsmesim_rat2_250ka,rsmesim_rat3_250ka))
RSMD_250mg_1_ka<- rmsd(a=SIM_data_pred$SIM_ka, b=SIM_data_pred$rat_1)
RSMD_250mg_2_ka<- rmsd(a=SIM_data_pred$SIM_ka, b=SIM_data_pred$rat_2)
RSMD_250mg_3_ka<- rmsd(a=SIM_data_pred$SIM_ka, b=SIM_data_pred$rat_3)

predvsout_250mg_ka <- ggplot(SIM_data_pred, aes(x=SIM_data_pred$SIM_ka, y=SIM_data_pred$rat_1)) +
  geom_point() +
  geom_abline(intercept=0, slope=1)+
  geom_line(data=tibble(x=c(0, 7.5),y=c(0,7.5)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 7.5),y=c(0,7.5)),linetype = "dashed",
            aes(x=x,y=y*2),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 7.5),y=c(0,7.5)),linetype = "dashed",
            aes(x=x,y=y*0.5),
            inherit.aes = FALSE)+
  labs(x='Predicted Values(μmol/L)', y='Actual Values(μmol/L)', title='Predicted vs. Actual Values 250mg/kg-BW oral',subtitle="Yuan et al 1992, Ka= 0.2")+
  annotate("text", x = 2.5, y = 17, label = "RSMD: 4.45", size= 8)+
  geom_point(aes(x=SIM_data_pred$SIM_ka,y=SIM_data_pred$rat_2))+
  geom_point(aes(x=SIM_data_pred$SIM_ka,y=SIM_data_pred$rat_3))+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
predvsout_250mg_ka

ggsave(plot=predvsout_250mg_ka,"Pred vs actual 250mg oral ka data.png",
       width= 11.69, height= 8.3, dpi= 250)

predvsout_250mg_ka


#plotting residual rat 1
ggplot(SIM_data_pred, aes(x=SIM_data_pred$Time, y=SIM_data_pred$Residual_Rat_1)) +
  geom_point() +
  geom_abline(slope = 0,intercept = 0)+ 
  labs(x= "Time", y='Residual values', title='Residual Values 250mg oral Yuan data')

#plotting residual rat 1  vs ka
ggplot(SIM_data_pred, aes(x=SIM_data_pred$Time, y=SIM_data_pred$Residual_Rat1_ka)) +
  geom_point() +
  geom_abline(slope = 0,intercept = 0)+ 
  labs(x= "Time", y='Residual values', title='Residual Values 250mg oral Yuan data ka')

#plotting residual rat 2  vs ka
ggplot(SIM_data_pred, aes(x=SIM_data_pred$Time, y=SIM_data_pred$Residual_Rat2_ka)) +
  geom_point() +
  geom_abline(slope = 0,intercept = 0)+ 
  labs(x= "Time", y='Residual values', title='Residual Values 250mg oral Yuan data')


#Calculating AUC values for the comparison graphs made above  

AUC_data <-PKNCAconc(Combined_data_file_for_graph_250mg, ug.ml~time|ID)

#Dosing data per Paper/simulation run 
#Loading dosing file
d_dose <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/250mg_comparison_doses.csv", sep=";")                       

dose_obj <- PKNCAdose(d_dose, dose~time|ID)


#letting pknc chose the end time of the auc calc
data_obj_automatic <- PKNCAdata(AUC_data, dose_obj)


#Computing the data both manual and automatic
results_obj_automatic <- pk.nca(data_obj_automatic)



#look at the data to get an impression
knitr::kable(head(as.data.frame(results_obj_automatic)))
summary(results_obj_automatic)

#saving data 
write.csv(results_obj_automatic$result,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//250mg_oral_results.csv")


#RAT_1 calculation
x <- c(Combined_data_file_for_graph_250mg[21:36,1])
y <- c(Combined_data_file_for_graph_250mg[21:36,2])
RAT_1_AUC <- auc(x[["time"]],y[["ug.ml"]], type=c("spline"))


#RAT_2
x <- c(Combined_data_file_for_graph_250mg[37:52,1])
y <- c(Combined_data_file_for_graph_250mg[37:52,2])
RAT_2_AUC <- auc(x[["time"]],y[["ug.ml"]], type=c("spline"))

#RAT_3
x <- c(Combined_data_file_for_graph_250mg[53:68,1])
y <- c(Combined_data_file_for_graph_250mg[53:68,2])
RAT_3_AUC <- auc(x[["time"]],y[["ug.ml"]], type=c("spline"))

#Zhao
x <- c(Combined_data_file_for_graph_250mg[1:8,1])
y <- c(Combined_data_file_for_graph_250mg[1:8,2])
ZHAO_AUC <- auc(x[["time"]],y[["ug.ml"]], type=c("spline"))

#Kiwa
x <- c(Combined_data_file_for_graph_250mg[9:20,1])
y <- c(Combined_data_file_for_graph_250mg[9:20,2])
Kiwa_AUC <- auc(x[["time"]],y[["ug.ml"]], type=c("spline"))


#iv concentration
#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]

#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


Combined_data_file_for_graph <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/IV_blood_Concentration.csv", sep=";")

p1 <- plot_ly(Combined_data_file_for_graph, x=~Time, y=Combined_data_file_for_graph$umol.L, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 20mg/kg-bw IV',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in umol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1


g <- ggplot(Combined_data_file_for_graph,aes(Time,umol.L,color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+
  labs(subtitle="IV dose 20mg/kg/bw", 
       y="Cinnamaldehyde concentration in umol/l", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")


KIWA_data <- Combined_data_file_for_graph[9:27,]
SIM_data_pred <-Combined_data_file_for_graph[c(29,30,31,32,33,34,35,36,37,38,40,43,46,49,52,55,58,61,63),]
SIM_data_kiwa <-Combined_data_file_for_graph[c(70,71,72,73,74,75,76,77,78,80,83,86,89,92,95,98,101,103,104),]

SIM_data_pred[,4]<- as.data.frame(KIWA_data[,2])
SIM_data_pred[,5]<- as.data.frame(SIM_data_kiwa[,2])
colnames(SIM_data_pred)<- c("Time","sim","ID","KIWA","sim kiwa")



#sim ka vs kiwa
RMSE_sim_kiwa_iv<-sqrt(mean((SIM_data_pred$sim - SIM_data_pred$KIWA)^2))


RSME_sim_kiwa<-(rsme_sim_kiwa)


p <-ggplot(SIM_data_pred, aes(x=SIM_data_pred$"sim kiwa", y=SIM_data_pred$KIWA)) +
  geom_point(size=3) +
  geom_abline(intercept=0, slope=1) +
  labs(x='Predicted Values R', y='Predicted values Kiwa', title='Predicted vs. predicted Values 20mg IV dose Kiwa')+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
ggsave(plot=p,"Pred vs pred 20mg iv.png",
       width= 11.69, height= 8.3, dpi= 250)



ZAO_data_2014 <- Combined_data_file_for_graph[1:8,]
SIM_data    <-Combined_data_file_for_graph[c(29,30,31,33,38,43,48,58),]

SIM_data[,3]<- as.data.frame(ZAO_data_2014[,2])


colnames(SIM_data)<- c("Time","sim","ZAO")

RSME_sim_IV<-sqrt(mean((SIM_data$sim - SIM_data$ZAO)^2))

RSMD<- rmsd(a=SIM_data$sim, b=SIM_data$ZAO)

p<-ggplot(SIM_data, aes(x=SIM_data$sim, y=SIM_data$ZAO)) +
  geom_point(size=3) +
  geom_abline(intercept=0, slope=1) +
  geom_line(data=tibble(x=c(0, 11),y=c(0,11)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 11),y=c(0,11)),linetype = "dashed",
            aes(x=x,y=y*2),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 11),y=c(0,11)),linetype = "dashed",
            aes(x=x,y=y*0.5),
            inherit.aes = FALSE)+
  labs(x='Predicted Values(μmol/L)', y='Actual Values(μmol/L)', title='Predicted vs. Actual Values 20mg IV dose',subtitle="Zao et al 2014")+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
p
ggsave(plot=p,"Pred vs actual 20mg iv.png",
       width= 11.69, height= 8.3, dpi= 250)
p


#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]

#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)



write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")

Combined_data_file_for_graph <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/IV_10mg_comparison.csv", sep=";")

p1 <- plot_ly(Combined_data_file_for_graph, x=~Time, y=Combined_data_file_for_graph$ug.ml, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 10mg/kg-BW IV',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in umol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1


g <- ggplot(Combined_data_file_for_graph,aes(Time,ug.ml,color=ID))

g + geom_point()+ 
  labs(subtitle="IV dose 10mg/kg/bw", 
       y="Cinnamaldehyde concentration in umol/L", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")

Shetty_data <- Combined_data_file_for_graph[1:8,]
SIM_data<- Combined_data_file_for_graph[c(11,14,19,24,29,40,89,129),]

SIM_data[,3]<- as.data.frame(Shetty_data[,2])


colnames(SIM_data)<- c("Time","sim","Shetty")

RSME_sim_IV_shetty<-sqrt(mean((SIM_data$sim - SIM_data$Shetty)^2))

RSMD<- rmsd(a=SIM_data$sim, b=SIM_data$Shetty)

p <-ggplot(SIM_data, aes(x=SIM_data$sim, y=SIM_data$Shetty)) +
  geom_point(size=3) +
  geom_abline(intercept=0, slope=1) +
  geom_line(data=tibble(x=c(0, 75),y=c(0,75)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 75),y=c(0,75)),linetype = "dashed",
            aes(x=x,y=y*2),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 75),y=c(0,75)),linetype = "dashed",
            aes(x=x,y=y*0.5),
            inherit.aes = FALSE)+
  labs(x='Predicted Values(μmol/L) ', y='Actual Values(μmol/L)', title='Predicted vs. Actual Values 10mg/kg-BW IV dose',subtitle="Shetty et al 2020")+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
p
ggsave(plot=p,"Pred vs actual 10mg iv.png",
       width= 11.69, height= 8.3, dpi= 250)
p




####-------------DONG et al 2022 data comparison--------------####
#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]

#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")



Combined_data_file_for_graph <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/375 oral dose comparison DONG et al 2022.csv", sep=";")

p1 <- plot_ly(Combined_data_file_for_graph, x=~Time, y=Combined_data_file_for_graph$umol.l, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 375mg/kg/bw',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in umol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1


g <- ggplot(Combined_data_file_for_graph,aes(Time,umol.l,color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+ 
  labs(subtitle="oral375mg/kg/bw dose", 
       y="Cinnamaldehyde concentration in umol/l", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")



Dong_data <- Combined_data_file_for_graph[1:10,]
SIM_data <-Combined_data_file_for_graph[c(12,13,14,16,19,21,31,51,71,91),]
SIM_data_kiwa<-Combined_data_file_for_graph[c(93,94,95,97,100,102,112,132,152,172),]

SIM_data[,4]<- as.data.frame(Dong_data[,2])
SIM_data[,5]<- as.data.frame(SIM_data_kiwa[,2])
colnames(SIM_data)<- c("Time","sim","ID","Dong","simka")


#sim ka vs kiwa
RMSE_sim_dong<-sqrt(mean((SIM_data$sim - SIM_data$Dong)^2))


RSME_sim_kiwa<-(rsme_sim_kiwa)

RSMD<- rmsd(a=SIM_data$sim, b=SIM_data$Dong)

p <-ggplot(SIM_data, aes(x=SIM_data$"sim", y=SIM_data$Dong)) +
  geom_point(size=3) +
  geom_abline(intercept=0, slope=1) +
  geom_line(data=tibble(x=c(0, 5000),y=c(0,5000)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 5000),y=c(0,5000)),linetype = "dashed",
            aes(x=x,y=y*2),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 5000),y=c(0,5000)),linetype = "dashed",
            aes(x=x,y=y*0.5),
            inherit.aes = FALSE)+
  labs(x='Predicted Values(μmol/L)', y='Predicted values(μmol/L)', title='Predicted vs. predicted Values 375mg/kg-BW dose',subtitle="Dong et al 2022")+
  theme_classic()+
  annotate("text", x = 2000, y = 7500, label = "RSMD: 10201", size= 8)+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
p
ggsave(plot=p,"Pred vs actual 375mg oral.png",
       width= 11.69, height= 8.3, dpi= 250)

RMSE_sim_dong_ka<-sqrt(mean((SIM_data$simka - SIM_data$Dong)^2))

RSMD<- rmsd(a=SIM_data$simka, b=SIM_data$Dong)
p <-ggplot(SIM_data, aes(x=SIM_data$"simka", y=SIM_data$Dong)) +
  geom_point(size=3) +
  geom_abline(intercept=0, slope=1) +
  geom_line(data=tibble(x=c(0, 5000),y=c(0,5000)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 5000),y=c(0,5000)),linetype = "dashed",
            aes(x=x,y=y*2),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 5000),y=c(0,5000)),linetype = "dashed",
            aes(x=x,y=y*0.5),
            inherit.aes = FALSE)+
  labs(x='Predicted Values(μmol/L)', y='Predicted values(μmol/L)', title='Predicted vs. predicted Values 375mg/kg-BW dose',subtitle="Dong et al 2022 Ka= 0.2")+
  theme_classic()+
  annotate("text", x = 2000, y = 7500, label = "RSMD: 10259", size= 8)+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
p
ggsave(plot=p,"Pred vs actual ka = 0.2 375mg oral.png",
      width= 11.69, height= 8.3, dpi= 250)



####-------------JI et al 2015 data comparison--------------####
#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]

#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")



Combined_data_file_for_graph <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/15 mg oral dose .csv", sep=";")

p1 <- plot_ly(Combined_data_file_for_graph, x=~Time, y=Combined_data_file_for_graph$umol.l, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 15mg/kg-bw IV',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in umol/l'),
         legend  =list(title= list(text='Type of Data')))
p1


Ji_data <- Combined_data_file_for_graph[1:11,]
SIM_data<- Combined_data_file_for_graph[c(13,15,17,22,27,32,42,52,92,132,252),]
SIM_data_ka<- Combined_data_file_for_graph[c(264,266,268,273,278,283,293,303,343,383,503),]
SIM_data[,3]<- as.data.frame(Ji_data[,2])
SIM_data[,4]<- as.data.frame(SIM_data_ka[,2])

colnames(SIM_data)<- c("Time","sim","Ji","simka")

RMSE_sim_ji<-sqrt(mean((SIM_data$sim - SIM_data$Ji)^2))

RMSD_sim_ji<- rmsd(a=SIM_data$sim, b=SIM_data$Ji )

p <-ggplot(SIM_data, aes(x=SIM_data$sim, y=SIM_data$Ji)) +
  geom_point(size=3) +
  geom_abline(intercept=0, slope=1) +
  geom_line(data=tibble(x=c(0, 1.5),y=c(0,1.5)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 1.5),y=c(0,1.5)),linetype = "dashed",
            aes(x=x,y=y*10),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 1.5),y=c(0,1.5)),linetype = "dashed",
            aes(x=x,y=y*0.1),
            inherit.aes = FALSE)+
  labs(x='Predicted Values(μmol/L)', y='Actual Values(μmol/L)', title='Predicted vs. Actual Values 15mg/kg-BW oral dose',subtitle="Ji et al 2015" )+
  theme_classic()+
  annotate("text", x = 0.25, y = 2, label = "RSMD: 0.76", size= 8)+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
p
ggsave(plot=p,"Pred vs actual 15mg 10 fold difference oral.png",
       width= 11.69, height= 8.3, dpi= 250)

g <- ggplot(Combined_data_file_for_graph,aes(Time,umol.l,color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+ 
  labs(subtitle="IV dose 15mg/kg/bw", 
       y="Cinnamaldehyde concentration in umol/l", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")

RMSE_sim_ka<-sqrt(mean((SIM_data$simka - SIM_data$Ji)^2))

RMSD_sim_ka<- rmsd(a=SIM_data$simka, b=SIM_data$Ji )

p <-ggplot(SIM_data, aes(x=SIM_data$simka, y=SIM_data$Ji)) +
  geom_point(size=3) +
  geom_abline(intercept=0, slope=1) +
  geom_line(data=tibble(x=c(0, 0.3),y=c(0,0.3)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 0.3),y=c(0,0.3)),linetype = "dashed",
            aes(x=x,y=y*5),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 0.3),y=c(0,0.3)),linetype = "dashed",
            aes(x=x,y=y*0.2),
            inherit.aes = FALSE)+
  labs(x='Predicted Values(μmol/L)', y='Actual Values(μmol/L)', title='Predicted vs. Actual Values 15mg/kg-BW oral dose',subtitle="Ji et al 2015 ka= 0.2" )+
  theme_classic()+
  annotate("text", x = 0.05, y = 0.40, label = "RSMD: 0.102", size= 8)+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
p
ggsave(plot=p,"Pred vs actual ka 15mg 5 fold difference oral.png",
       width= 11.69, height= 8.3, dpi= 250)


####-------------YOUNG et al 2020 data comparison--------------####
#Generating a file for the comparison with in vivo data 
blood_data <- as.data.frame(solve.pbk_rat[,c(1,3,4)])
#write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")


#going from umol/l to umol
blood_data[2]<- blood_data[2] * V_V
blood_data[3]<-blood_data[3]  * V_A

#combining C_v + C_A to create a combined amount of cinnamaldehyde
blood_data[3]<-blood_data[2]+ blood_data[3]

#Dropping the second colum 
blood_data<-blood_data[-c(2)]

#going back an amount to an amount per L
blood_data[2]<- blood_data[2]/ (V_V + V_A)


write.csv(blood_data,"D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK//Blood_Data.csv")



Combined_data_file_for_graph <- read.csv("D:/Joris/Toxicology and Environmental Health/Master stage/R/Cinnamaldehyde PBK/Yong et al 2020 50mg-kg-bw.csv", sep=";")

p1 <- plot_ly(Combined_data_file_for_graph, x=~Time, y=Combined_data_file_for_graph$umol.l, 
              color = ~ ID , 
              colors = "Set2",
              type= "scatter",
              mode= "markers",
              hovertext= ~sample)%>%
  layout(title= 'Blood concentration comparison dose = 50mg/kg/bw',
         xaxis= list(title= 'Time (hours)'),
         yaxis= list(title= 'Cinnamaldehyde concentration in umol/l', type="log"),
         legend  =list(title= list(text='Type of Data')))
p1


g <- ggplot(Combined_data_file_for_graph,aes(Time,umol.l,color=ID))

g + geom_point()+ scale_y_continuous(trans='log10')+ 
  labs(subtitle="oral 50 mg/kg/bw dose", 
       y="Cinnamaldehyde concentration in umol/l", 
       x="Time in Hours", 
       title="Cinnamaldehyde concentration in blood", 
       caption="PBK model")

Yong_data <- Combined_data_file_for_graph[1:13,]
SIM_data<- Combined_data_file_for_graph[c(15,16,17,19,22,24,29,34,44,54,74,94,134),]
SIM_data_ka<- Combined_data_file_for_graph[c(136,137,138,140,143,145,150,155,165,175,195,215,255),]
SIM_data[,3]<- as.data.frame(Yong_data[,2])
SIM_data[,4]<- as.data.frame(SIM_data_ka[,2])



colnames(SIM_data)<- c("Time","sim","Yong","simka")


RMSE_sim_Yong<-sqrt(mean((SIM_data$sim - SIM_data$Yong)^2))

RMSD_sim_Yong<-rmsd(a=SIM_data$sim, b=SIM_data$Yong)
p <-ggplot(SIM_data, aes(x=SIM_data$sim, y=SIM_data$Yong)) +
  geom_point(size=3) +
  geom_abline(intercept=0, slope=1) +
  geom_line(data=tibble(x=c(0, 5),y=c(0,5)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 5),y=c(0,5)),linetype = "dashed",
            aes(x=x,y=y*2),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 5),y=c(0,5)),linetype = "dashed",
            aes(x=x,y=y*0.5),
            inherit.aes = FALSE)+
  labs(x='Predicted Values(μmol/L)', y='Actual Values(μmol/L)', title='Predicted vs. Actual Values 50mg/kg-BW oral dose',subtitle="Yong et al 2020" )+
  theme_classic()+
  annotate("text", x = 1, y = 6, label = "RSMD: 1.74", size= 8)+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
p
ggsave(plot=p,"Pred vs actual 50mg oral.png",
       width= 11.69, height= 8.3, dpi= 250)


RMSE_sim_ka<-sqrt(mean((SIM_data$simka - SIM_data$Yong)^2))

RMSD_sim_ka<-rmsd(a=SIM_data$simka, b=SIM_data$Yong)

p <-ggplot(SIM_data, aes(x=SIM_data$simka, y=SIM_data$Yong)) +
  geom_point(size=3) +
  geom_abline(intercept=0, slope=1) +
  geom_line(data=tibble(x=c(0, 2),y=c(0,2)),
            aes(x=x,y=y),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 2),y=c(0,2)),linetype = "dashed",
            aes(x=x,y=y*2),
            inherit.aes = FALSE)+
  geom_line(data=tibble(x=c(0, 2),y=c(0,2)),linetype = "dashed",
            aes(x=x,y=y*0.5),
            inherit.aes = FALSE)+
  labs(x='Predicted Values(μmol/L)', y='Actual Values(μmol/L)', title='Predicted vs. Actual Values 50mg/kg-BW oral dose',subtitle="Yong et al 2020 ka= 0.2" )+
  theme_classic()+
  annotate("text", x = 0.5, y = 3, label = "RSMD: 1.82", size= 8)+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
p
ggsave(plot=p,"Pred vs actual ka 5 fold difference 50mg oral.png",
       width= 11.69, height= 8.3, dpi= 250)


