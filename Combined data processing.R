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



mass_df_single <- solve.pbk/BW * MW /1e+3
mass_df_single <- mass_df_single[,c(65:69,71:77,79,80,83:87,91:95)]
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
mass_df <- mass_df[1:481,c(67:80,82,83,86:90,94:98)]
mass_at_t <- data.frame(mass=as.numeric())


for (i in 1:nrow(mass_df)){
  mass_at_t[nrow(mass_at_t) + 1,] <- rowSums(mass_df[i,])
}
plot(mass_at_t[,1])

results_human<-solve.pbk[c(1:241),c(1:98)]



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


#setwd("C:/Users/OrfeasPetropoulos/OneDrive - Wageningen University & Research/joris")

data_2.8<- read.csv("melt_boxplot_2.8mg_inhalation_auc")
unique(data_2.8$id)
data_2.8$id <- factor(data_2.8$id)
unique(data_2.8$variable)

#subset the data_2.8 for analysis
{
  data_2.8_Lung <- subset(data_2.8, variable %in% c("Lung"))
  data_2.8_Blood <- subset(data_2.8, variable %in% c("Blood"))
  data_2.8_Fat <- subset(data_2.8, variable %in% c("Fat"))
  data_2.8_SlowlyPerfused <- subset(data_2.8, variable %in% c("Slowly Perfused"))
  data_2.8_Liver <- subset(data_2.8, variable %in% c("Liver"))
  data_2.8_RichlyPerfused <- subset(data_2.8, variable %in% c("Richly Perfused"))
  data_2.8_SmallIntestine <- subset(data_2.8, variable %in% c("Small Intestine"))
  data_2.8_MALE <- subset(data_2.8, id %in% c("Male"))
  data_2.8_FEMALE <- subset(data_2.8, id %in% c("Female"))
} 


data_250<- read.csv("melt_boxplot_250mg_inhalation_auc")
unique(data_250$id)
data_250$id <- factor(data_250$id)
unique(data_250$variable)

#subset the data_250 for analysis
{
  data_Lung <- subset(data_250, variable %in% c("Lung"))
  data_Blood <- subset(data_250, variable %in% c("Blood"))
  data_Fat <- subset(data_250, variable %in% c("Fat"))
  data_SlowlyPerfused <- subset(data_250, variable %in% c("Slowly Perfused"))
  data_Liver <- subset(data_250, variable %in% c("Liver"))
  data_RichlyPerfused <- subset(data_250, variable %in% c("Richly Perfused"))
  data_SmallIntestine <- subset(data_250, variable %in% c("Small Intestine"))
  data_MALE <- subset(data_250, id %in% c("Male"))
  data_FEMALE <- subset(data_250, id %in% c("Female"))
} 

data_250_oral<- read.csv("melt_boxplot_250mg_oral_auc")
unique(data_250_oral$id)
data_250_oral$id <- factor(data_250_oral$id)
unique(data_250_oral$variable)

#subset the data_250_oral for analysis
{
  data_250_oral_Lung <- subset(data_250_oral, variable %in% c("Lung"))
  data_250_oral_Blood <- subset(data_250_oral, variable %in% c("Blood"))
  data_250_oral_Fat <- subset(data_250_oral, variable %in% c("Fat"))
  data_250_oral_SlowlyPerfused <- subset(data_250_oral, variable %in% c("Slowly Perfused"))
  data_250_oral_Liver <- subset(data_250_oral, variable %in% c("Liver"))
  data_250_oral_RichlyPerfused <- subset(data_250_oral, variable %in% c("Richly Perfused"))
  data_250_oral_SmallIntestine <- subset(data_250_oral, variable %in% c("Small Intestine"))
  data_250_oral_MALE <- subset(data_250_oral, id %in% c("Male"))
  data_250_oral_FEMALE <- subset(data_250_oral, id %in% c("Female"))
} 


#Creation of combined 2.8 and 250 inhalation testing-----------
data_250_Lung$id<-as.factor(250)
data_2.8_Lung$id<-as.factor(2.8)

combined_lung<-rbind(data_250_Lung,data_2.8_Lung)

data_250_Blood$id<-as.factor(250)
data_2.8_Blood$id<-as.factor(2.8)

combined_blood<-rbind(data_250_Blood,data_2.8_Blood)

data_250_Fat$id<-as.factor(250)
data_2.8_Fat$id<-as.factor(2.8)

combined_Fat<-rbind(data_250_Fat,data_2.8_Fat)



data_250_SlowlyPerfused$id<-as.factor(250)
data_2.8_SlowlyPerfused$id<-as.factor(2.8)

combined_SlowlyPerfused<-rbind(data_250_SlowlyPerfused,data_2.8_SlowlyPerfused)



data_250_RichlyPerfused$id<-as.factor(250)
data_2.8_RichlyPerfused$id<-as.factor(2.8)

combined_RichlyPerfused<-rbind(data_250_RichlyPerfused,data_2.8_RichlyPerfused)

data_250_SmallIntestine$id<-as.factor(250)
data_2.8_SmallIntestine$id<-as.factor(2.8)

combined_SmallIntestine<-rbind(data_250_SmallIntestine,data_2.8_SmallIntestine)

data_250_Liver$id<-as.factor(250)
data_2.8_Liver$id<-as.factor(2.8)

combined_liver<-rbind(data_250_Liver,data_2.8_Liver)


###Data analysis--------

#1. Is there difference between the organs? for combined------------
{
  #1a--> between exposures lung 
  results <- aov(value ~ id , combined_lung)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1b--> between exposures blood 
  results <- aov(value ~ id , combined_blood)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1c--> between exposures Fat 
  results <- aov(value ~ id , combined_Fat)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1d--> between exposures Slowlyperfused 
  results <- aov(value ~ id , combined_SlowlyPerfused)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1e--> between exposures Richlyperfused 
  results <- aov(value ~ id , combined_RichlyPerfused)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1e--> between exposures SmallIntestine 
  results <- aov(value ~ id , combined_SmallIntestine)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1e--> between exposures SmallIntestine 
  results <- aov(value ~ id , combined_SmallIntestine)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1e--> between exposures Liver
  results <- aov(value ~ id , combined_liver)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  
} 
#Combined 250 oral and inhalation---------------
data_250_Lung$id<-as.factor(250)
data_250_oral_Lung$id<-as.factor(2.8)

combined_lung_oral<-rbind(data_250_Lung,data_250_oral_Lung)

data_250_Blood$id<-as.factor(250)
data_250_oral_Blood$id<-as.factor(2.8)

combined_blood_oral<-rbind(data_250_Blood,data_250_oral_Blood)

data_250_Fat$id<-as.factor(250)
data_250_oral_Fat$id<-as.factor(2.8)

combined_Fat_oral<-rbind(data_250_Fat,data_250_oral_Fat)



data_250_SlowlyPerfused$id<-as.factor(250)
data_250_oral_SlowlyPerfused$id<-as.factor(2.8)

combined_SlowlyPerfused_oral<-rbind(data_250_SlowlyPerfused,data_250_oral_SlowlyPerfused)



data_250_RichlyPerfused$id<-as.factor(250)
data_250_oral_RichlyPerfused$id<-as.factor(2.8)

combined_RichlyPerfused_oral<-rbind(data_250_RichlyPerfused,data_250_oral_RichlyPerfused)

data_250_SmallIntestine$id<-as.factor(250)
data_250_oral_SmallIntestine$id<-as.factor(2.8)

combined_SmallIntestine_oral<-rbind(data_250_SmallIntestine,data_250_oral_SmallIntestine)

data_250_Liver$id<-as.factor(250)
data_250_oral_Liver$id<-as.factor(2.8)

combined_liver_oral<-rbind(data_250_Liver,data_250_oral_Liver)

#1. Is there difference between the organs? combined-------------
{
  #1a--> between exposures lung 
  results <- aov(value ~ id , combined_lung_oral)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1b--> between exposures blood 
  results <- aov(value ~ id , combined_blood_oral)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1c--> between exposures Fat 
  results <- aov(value ~ id , combined_Fat_oral)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1d--> between exposures Slowlyperfused 
  results <- aov(value ~ id , combined_SlowlyPerfused_oral)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1e--> between exposures Richlyperfused 
  results <- aov(value ~ id , combined_RichlyPerfused_oral)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1e--> between exposures SmallIntestine 
  results <- aov(value ~ id , combined_SmallIntestine_oral)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1e--> between exposures SmallIntestine 
  results <- aov(value ~ id , combined_SmallIntestine_oral)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
  #1e--> between exposures Liver
  results <- aov(value ~ id , combined_liver_oral)
  summary(results)
  tukey <- TukeyHSD(results, conf.level = 0.95)
  tukey.cld <- multcompLetters4(results, tukey)
  print(tukey.cld)
  
} 


#1. Is there difference between the organs?---------
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

#male female differences in a organ
#because almost all variables have UNEQUAL VARIANCE, I performed 
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


welch_anova_test(combined_lung, value ~ id)
welch_anova_test(combined_blood, value ~ id)
welch_anova_test(combined_Fat, value ~ id)
welch_anova_test(combined_RichlyPerfused, value ~ id)
welch_anova_test(combined_SlowlyPerfused, value ~ id)
welch_anova_test(combined_liver, value ~ id)
welch_anova_test(combined_SmallIntestine, value ~ id) #for small intestine not needed, as it had equal variances


welch_anova_test(combined_lung_oral, value ~ id)
welch_anova_test(combined_blood_oral, value ~ id)
welch_anova_test(combined_Fat_oral, value ~ id)
welch_anova_test(combined_RichlyPerfused_oral, value ~ id)
welch_anova_test(combined_SlowlyPerfused_oral, value ~ id)
welch_anova_test(combined_liver_oral, value ~ id)
welch_anova_test(combined_SmallIntestine_oral, value ~ id) #for small intestine not needed, as it had equal variances



#visualization of popgen results
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
  labs(y = "Blood concentration (μmol/l) ",
       x = "Time (h)",
       title='CNMA blood concentration')  +
  xlim(0,12)+
theme_classic()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size = 15),
        title = element_text(size=20))+
  theme(legend.text = element_text(size=15, color="black"),
        legend.position ="top")

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
  labs(y = "Arterial Blood concentration (μmol/l) ",
       x = "Time (h)")  +
  theme_classic()+
  xlim(0,4)+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))

gg


#fat concentration
tab_solve_C_F=as.data.frame(matrix(NA,time.end/time.frame+1,(N+NF)))    #Create an empty data frame with amount of timepoints=amount of rows and amount of individuals=amount of columns
for (i in 1:(N+NF)) {
  tab.i=solve.pbk[which(solve.pbk[,"id"]==i),]                  #Put all individuals in data frame
  tab.i=as.data.frame(tab.i)
  tab_solve_C_F[,i]=tab.i$C_F
}

tab_C_F=as.data.frame(matrix(NA,time.end/time.frame+1,4))       #Create an empty data frame with amount of timepoints=amount of rows and 4 columns
tab_C_F[,1]=c(seq(time.0,time.end,by=time.frame))               #Timepoints in first column
for (i in 1:(time.end/time.frame+1)) {
  tab_C_F[i,2]=quantile(tab_solve_C_F[i,],0.025, na.rm = TRUE)      #Lower bound of confidence interval in second column
  tab_C_F[i,3]=quantile(tab_solve_C_F[i,],0.5, na.rm = TRUE)        #Median in third column
  tab_C_F[i,4]=quantile(tab_solve_C_F[i,],0.975, na.rm = TRUE)      #Upper bound of confidence interval in fourth column
}
colnames(tab_C_F)=c("time","CF_P2.5","CF_P50","CF_P97.5")       #Add column names

gg <- ggplot(tab_C_F)+
  geom_line(aes(x=time, y=CF_P2.5), linetype = "dashed")+
  geom_line(aes(x=time, y=CF_P50), color = "red", size = 1)+
  geom_line(aes(x=time, y=CF_P97.5), linetype = "dashed")+
  labs(y = "Cinnamaldehyde concentratin in Fat μmol/l ",
       x = "Time (h)")  +
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
ggsave(plot=gg,"CNMA concentration in fat tissue 2.8mg/kg-BW",
       width= 11.69, height= 8.3, dpi= 250)
gg

#Lung concentration
tab_solve_C_Pu=as.data.frame(matrix(NA,time.end/time.frame+1,(N+NF)))    #Create an empty data frame with amount of timepoints=amount of rows and amount of individuals=amount of columns
for (i in 1:(N+NF)) {
  tab.i=solve.pbk[which(solve.pbk[,"id"]==i),]                  #Put all individuals in data frame
  tab.i=as.data.frame(tab.i)
  tab_solve_C_Pu[,i]=tab.i$C_Pu
}

tab_C_Pu=as.data.frame(matrix(NA,time.end/time.frame+1,4))       #Create an empty data frame with amount of timepoints=amount of rows and 4 columns
tab_C_Pu[,1]=c(seq(time.0,time.end,by=time.frame))               #Timepoints in first column
for (i in 1:(time.end/time.frame+1)) {
  tab_C_Pu[i,2]=quantile(tab_solve_C_Pu[i,],0.025, na.rm = TRUE)      #Lower bound of confidence interval in second column
  tab_C_Pu[i,3]=quantile(tab_solve_C_Pu[i,],0.5, na.rm = TRUE)        #Median in third column
  tab_C_Pu[i,4]=quantile(tab_solve_C_Pu[i,],0.975, na.rm = TRUE)      #Upper bound of confidence interval in fourth column
}
colnames(tab_C_Pu)=c("time","C_P2.5","C_P50","C_P97.5")       #Add column names

gg <- ggplot(tab_C_Pu)+
  geom_line(aes(x=time, y=C_P2.5), linetype = "dashed")+
  geom_line(aes(x=time, y=C_P50), color = "red", size = 1)+
  geom_line(aes(x=time, y=C_P97.5), linetype = "dashed")+
  labs(y = "Cinnamaldehyde concentratin in Lung μmol/l ",
       x = "Time(h)",  
  title='CNMA Concentratin in Lung Tissue')  +
  theme_classic()+
  xlim(0,12)+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
ggsave(plot=gg,"CNMA concentration in Lung tissue 2.8mg/kg-BW",
       width= 11.69, height= 8.3, dpi= 250)

gg



#Lung concentration
tab_solve_C_L=as.data.frame(matrix(NA,time.end/time.frame+1,(N+NF)))    #Create an empty data frame with amount of timepoints=amount of rows and amount of individuals=amount of columns
for (i in 1:(N+NF)) {
  tab.i=solve.pbk[which(solve.pbk[,"id"]==i),]                  #Put all individuals in data frame
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
colnames(tab_C_L)=c("time","C_P2.5","C_P50","C_P97.5")       #Add column names

gg <- ggplot(tab_C_L)+
  geom_line(aes(x=time, y=C_P2.5), linetype = "dashed")+
  geom_line(aes(x=time, y=C_P50), color = "red", size = 1)+
  geom_line(aes(x=time, y=C_P97.5), linetype = "dashed")+
  labs(y = "Cinnamaldehyde concentratin in Lung μmol/l ",
       x = "Time (h)")  +
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
ggsave(plot=gg,"CNMA concentration in Liver tissue 2.8mg/kg-BW",
       width= 11.69, height= 8.3, dpi= 250)
gg
