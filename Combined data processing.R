#author: Joris Jean van der Lugt
#date: 05-08-2022
#Combined data processing of Human Cinnamaldehyde pbk models
#Before running code first run all model associated with this data processing
#human rxode model


library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(shiny)
library(truncnorm)
library(reshape2)
library(plotly)
library(PKNCA)


#Rxode
#cinnamaldehyde model Human
#Mass balance calculation rxode inhalation complete
mass_df_single <- solve.pbk/BW * MW /1e+3
mass_df_single <- mass_df_single[,c(66:79,81,82,85:89,93:97)]
mass_at_t_single <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df_single)){
  mass_at_t_single[nrow(mass_at_t_single) + 1,] <- rowSums(mass_df_single[i,])
}
plot(mass_at_t_single[,1])


#cinnamaldehyde model Rat
#Mass balance calculation rxode inhalation complete
mass_df <- solve.pbk_rat/BW * MW /1e+3
mass_df <- mass_df[,c(66:84,86,88:92,94,96)]
mass_at_t <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])


#Rxode 
#population mass balance
mass_df <- solve.pbk/phys[2,3] * MW /1e+3
mass_df <- mass_df[242:482,c(67:80,82,83,86:90,94:98)]
mass_at_t <- data.frame(mass=as.numeric())


for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])

results_human<-solve.pbk[c(1:241),c(1:98)]
plot(results_human$AM_L_CA)



#Non inhalation mass balance 
mass_df <- solve.pbk/BW * MW /1e+3
mass_df <- mass_df[,c(59:68,70,71,74:81,85:89)]
mass_at_t <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])

#non inhalation Rat
mass_df <- solve.pbk_rat/BW * MW /1e+3
mass_df <- mass_df[,c(57:71,73,75,77:80,82,84)]
mass_at_t <- data.frame(mass=as.numeric())

for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])


#Metabolic output checking
#human model
#Checking proportions of urinary metabolites
#Urinary metabolites are all ultimately derived from the oxidation of Cinnamaldehyde.

data_urinary_metabolites <- solve.pbk[,c(1,78,87)]

total_dose <- Oral_Dose
total_dose<- Inhalation_Dose

total_urinary_metabolites_24h <- data_urinary_metabolites[241,2] + data_urinary_metabolites[241,3]

percentage_metabolised <- total_urinary_metabolites_24h/total_dose *100

#-------------population model
# data extraction  

meta_extraction <- unique(solve.pbk[solve.pbk$time == 24,c("time", "id","AM_L_CA","AM_SI_CA")])

dose_extraction <- unique(ex[ex$cmt == "A_GI",c("id","amt")])

percentage_metabolised_popgen <-c()
for (i in 1:2000){
  total <- meta_extraction[i,3] + meta_extraction[i,4]
  dose  <- dose_extraction[i,2]
  metabolised <- total/dose *100
  percentage_metabolised_popgen <-append(percentage_metabolised_popgen, metabolised)
}


meta_extraction <- meta_extraction[-c(1,3,4)]
meta_extraction[,2] <- as.data.frame(percentage_metabolised_popgen)


colnames(meta_extraction)<-c("id","percentage metabolised")
melt_meta_extraction <- melt(meta_extraction,id=c("id")) 



rownames(single_meta) <-c("2001")

melt_meta_extraction$id[meta_extraction$id == 1:1000] <- "Male"  
melt_meta_extraction$id[meta_extraction$id == 1001:2000] <- "Female" 

male_melt_meta_extraction<-melt_meta_extraction[1:1000,1:3]
female_melt_meta_extraction<-melt_meta_extraction[1001:2000,1:3]

p_metabolism<-ggplot(melt_meta_extraction ,aes(x=variable,y=value))+
  geom_boxplot(notch=TRUE)+
  geom_jitter(aes(col=id),alpha=0.6)+
  scale_color_manual(values = c( "Male" = "blue",
                                 "Female" = "red"),
                     labels= c( "Male", "Female"),
                     name= "Sex")+
  labs(x='', y='Percentage metabolised', title='Carboxylic acid metabolism Human')+
  theme_classic()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size = 15),
        title = element_text(size=20))+
  theme(legend.text = element_text(size=15, color="black"),
        legend.position ="top")
ggplotly(p_metabolism)


p_metabolism<-ggplot(male_melt_meta_extraction ,aes(x=variable,y=value))+
  geom_boxplot(notch=TRUE)+
  geom_jitter(aes(col=id),alpha=0.6)+
  labs(x='', y='Percentage metabolised', title='Carboxylic acid metabolism Human')+
  theme_classic()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size = 15),
        title = element_text(size=20))+
  theme(legend.text = element_text(size=15, color="black"),
        legend.position ="top")
ggplotly(p_metabolism)

p_metabolism<-ggplot(female_melt_meta_extraction ,aes(x=variable,y=value))+
  geom_boxplot(notch=TRUE)+
  geom_jitter(aes(col=id),alpha=0.6)+
  labs(x='', y='Percentage metabolised', title='Carboxylic acid metabolism Human')+
  theme_classic()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size = 15),
        title = element_text(size=20))+
  theme(legend.text = element_text(size=15, color="black"),
        legend.position ="top")
ggplotly(p_metabolism)


boxplot(male_metabolised_data,female_metabolised_data,
        ylab= "percentage metabolised into carboxylic acid metabolites",
        names= c("Male", "Female"),
        col="orange",
        las=2,
        title="Amount metabolized to carboxylic acid")+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
hist(male_metabolised_data)
hist(female_metabolised_data)


#significance analysis 
#data analysis for Joris
#info: Levene test is less sensitive than the Bartlett test to departures from normality.

library(ggpubr)
library(plyr)
library(multcompView)
library(car)


setwd("C:/Users/OrfeasPetropoulos/OneDrive - Wageningen University & Research/joris")

data<- read.csv("melt_boxplot_250mg_inhalation_cmax")
unique(data$id)
data$id <- factor(data$id)
unique(data$variable)

#subset the data for analysis
{
  data_Lung <- subset(data, variable %in% c("Lung"))
  data_Blood <- subset(data, variable %in% c("Blood"))
  data_Fat <- subset(data, variable %in% c("Fat"))
  data_SlowlyPerfused <- subset(data, variable %in% c("Slowly Perfused"))
  data_Liver <- subset(data, variable %in% c("Liver"))
  data_RichlyPerfused <- subset(data, variable %in% c("Richly Perfused"))
  data_SmallIntestine <- subset(data, variable %in% c("Small Intestine"))
  data_MALE <- subset(data, id %in% c("Male"))
  data_FEMALE <- subset(data, id %in% c("Female"))
} 


###Data analysis

#1. Is there difference between the organs?
{
  #1a--> within Female
  results <- aov(value ~ variable , data_FEMALE)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1b--> within Male
  results <- aov(value ~ variable , data_MALE)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
} 


#2. Is there difference between males and females for the same organ?
#2.a. LUNG DATA
{
  shapiro.test(data_Lung$value) # we can reject the null hypothesis that the data tested follows a normal distribution
  plot(density(na.omit(data_Lung$value)))
  ggqqplot(data_Lung$value) #it is indeed right skewed
  leveneTest(value ~ id, data = data_Lung) #sign.we can reject the null hypothesis that the data tested has equal variances
  #because of big sample size, we are allowed to perform parametric test
  results <- aov(value ~ id , data_Lung)
  # t.test(value ~ id , data = data_Lung, var.equal= F) #is the same as the anova here
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld) #Significant differences between blocks
  
} 

#2.b. BLOOD DATA
{
  shapiro.test(data_Blood$value) # we can reject the null hypothesis that the data tested follows a normal distribution
  plot(density(na.omit(data_Blood$value)))
  ggqqplot(data_Blood$value) #it is indeed right skewed
  leveneTest(value ~ id, data = data_Blood) #sign.we can reject the null hypothesis that the data tested has equal variances
  #because of big sample size, we are allowed to perform parametric test
  results <- aov(value ~ id , data_Blood)
  # t.test(value ~ id , data = data_Lung, var.equal= F) #is the same as the anova here
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld) #Significant differences between blocks
  
} 

#2.c. FAT DATA
{
  shapiro.test(data_Fat$value) # we can reject the null hypothesis that the data tested follows a normal distribution
  plot(density(na.omit(data_Fat$value)))
  ggqqplot(data_Fat$value) #it is indeed right skewed
  leveneTest(value ~ id, data = data_Fat) #sign.we can reject the null hypothesis that the data tested has equal variances
  #because of big sample size, we are allowed to perform parametric test
  results <- aov(value ~ id , data_Fat)
  # t.test(value ~ id , data = data_Lung, var.equal= F) #is the same as the anova here
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld) #Significant differences between blocks
  
} 

#2.d. SlowlyPerfused DATA
{
  shapiro.test(data_SlowlyPerfused$value) # we can reject the null hypothesis that the data tested follows a normal distribution
  plot(density(na.omit(data_SlowlyPerfused$value)))
  ggqqplot(data_SlowlyPerfused$value) #it is indeed right skewed
  leveneTest(value ~ id, data = data_SlowlyPerfused) #sign.we can reject the null hypothesis that the data tested has equal variances
  #because of big sample size, we are allowed to perform parametric test
  results <- aov(value ~ id , data_SlowlyPerfused)
  # t.test(value ~ id , data = data_Lung, var.equal= F) #is the same as the anova here
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld) #Significant differences between blocks
  
} 

#2.e. LIVER DATA
{
  shapiro.test(data_Liver$value) # we can reject the null hypothesis that the data tested follows a normal distribution
  plot(density(na.omit(data_Liver$value)))
  ggqqplot(data_Liver$value) #it is indeed right skewed
  leveneTest(value ~ id, data = data_Liver) #sign.we can reject the null hypothesis that the data tested has equal variances
  #because of big sample size, we are allowed to perform parametric test
  results <- aov(value ~ id , data_Liver)
  # t.test(value ~ id , data = data_Lung, var.equal= F) #is the same as the anova here
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld) #Significant differences between blocks
  
} 

#2.f. Richly Perfused DATA
{
  shapiro.test(data_RichlyPerfused$value) # we can reject the null hypothesis that the data tested follows a normal distribution
  plot(density(na.omit(data_RichlyPerfused$value)))
  ggqqplot(data_RichlyPerfused$value) #it is indeed right skewed
  leveneTest(value ~ id, data = data_RichlyPerfused) #sign.we can reject the null hypothesis that the data tested has equal variances
  #because of big sample size, we are allowed to perform parametric test
  results <- aov(value ~ id , data_RichlyPerfused)
  # t.test(value ~ id , data = data_Lung, var.equal= F) #is the same as the anova here
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld) #Significant differences between blocks
  
} 

#2.g. Small Intestine DATA
{
  shapiro.test(data_SmallIntestine$value) # we can reject the null hypothesis that the data tested follows a normal distribution
  plot(density(na.omit(data_SmallIntestine$value)))
  ggqqplot(data_SmallIntestine$value) #it is indeed right skewed
  leveneTest(value ~ id, data = data_SmallIntestine) #NOT sign.we cannot reject the null hypothesis that the data tested has equal variances
  #because of big sample size, we are allowed to perform parametric test
  results <- aov(value ~ id , data_SmallIntestine)
  # t.test(value ~ id , data = data_Lung, var.equal= F) #is the same as the anova here
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld) #Significant differences between blocks
  
} 


#because almost all variables have UNEQUAL VARIANCE, I peformed 
# a welch test instead of the one-way anova, just to be sure 
# but the results remain the same!

library(rstatix)
welch_anova_test(data_Lung, value ~ id)
welch_anova_test(data_Blood, value ~ id)
welch_anova_test(data_Fat, value ~ id)
welch_anova_test(data_RichlyPerfused, value ~ id)
welch_anova_test(data_SlowlyPerfused, value ~ id)
welch_anova_test(data_Liver, value ~ id)
welch_anova_test(data_SmallIntestine, value ~ id) #for small intestine not needed, as it had equal variances


fig <- ggplot(subset(data, variable %in% c("Lung")),
              mapping = aes(x = id, y= value))+
  geom_boxplot(aes(fill= id),position="identity", size=0.5)+
  geom_dotplot(fill="black", binaxis='y', stackdir='center',position=position_dodge(0.5), binwidth= 0.03)+
  
  #adding the significance letters manualy
  geom_text(aes(x, y, label=lab),
            data=data.frame(x=c("Female","Male"),
                            y=c(6,8.5),
                            lab=c("a", "b"),
                            size=4))

fig



#visualisation of popgen results
#Venous blood concentration
tab_solve_C_V=as.data.frame(matrix(NA,time.end/time.frame+1,(N+NF)))    #Create an empty data frame with amount of timepoints=amount of rows and amount of individuals=amount of columns
for (i in 1:(N+NF)) {
  tab.i=solve.pbk[which(solve.pbk[,"id"]==i),]                  #Put all individuals in data frame
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
  labs(y = "Venous Blood concentration ",
       x = "Time (h)")  +
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))

gg

#arterial blood concentration
tab_solve_C_A=as.data.frame(matrix(NA,time.end/time.frame+1,(N+NF)))    #Create an empty data frame with amount of timepoints=amount of rows and amount of individuals=amount of columns
for (i in 1:(N+NF)) {
  tab.i=solve.pbk[which(solve.pbk[,"id"]==i),]                  #Put all individuals in data frame
  tab.i=as.data.frame(tab.i)
  tab_solve_C_A[,i]=tab.i$C_A
}

tab_C_A=as.data.frame(matrix(NA,time.end/time.frame+1,4))       #Create an empty data frame with amount of timepoints=amount of rows and 4 columns
tab_C_A[,1]=c(seq(time.0,time.end,by=time.frame))               #Timepoints in first column
for (i in 1:(time.end/time.frame+1)) {
  tab_C_A[i,2]=quantile(tab_solve_C_A[i,],0.025, na.rm = TRUE)      #Lower bound of confidence interval in second column
  tab_C_A[i,3]=quantile(tab_solve_C_A[i,],0.5, na.rm = TRUE)        #Median in third column
  tab_C_A[i,4]=quantile(tab_solve_C_A[i,],0.975, na.rm = TRUE)      #Upper bound of confidence interval in fourth column
}
colnames(tab_C_A)=c("time","CA_P2.5","CA_P50","CA_P97.5")       #Add column names

gg <- ggplot(tab_C_A)+
  geom_line(aes(x=time, y=CA_P2.5), linetype = "dashed")+
  geom_line(aes(x=time, y=CA_P50), color = "red", size = 1)+
  geom_line(aes(x=time, y=CA_P97.5), linetype = "dashed")+
  labs(y = "Arterial Blood concentration ",
       x = "Time (h)")  +
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))

gg





