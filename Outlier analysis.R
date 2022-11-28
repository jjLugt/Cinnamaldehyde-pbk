

Liver_check <-as.data.frame(parameters[,18])
Liver_check[,2]<- 1:2000
Liver_check[,3]<-metabolised_data_frame
#Liver_check[,4]<- parameters[,49]
#Liver_check[,5]<- parameters[,17]
colnames(Liver_check) <-c("V_L","id","percentage")#,"vsmax","fat")


melt_liver_check<-melt(Liver_check,id=c("percentage"))


melt_liver_check<-melt_liver_check[c(1:2000),c(1:3)]
melt_liver_check[,4]<-1:2000
colnames(melt_liver_check) <-c("percentage","variable","value","id")


melt_liver_check$id[melt_liver_check$id == 1:1000] <- "Male"  
melt_liver_check$id[melt_liver_check$id == 1001:2000] <- "Female" 

p<-ggplot(melt_liver_check, aes(x=value, y=percentage,colour=id)) +
  geom_point(alpha=0.6)+
  scale_color_manual(values = c( "Male" = "blue",
                                 "Female" = "red"),
                     labels= c( "Male", "Female"),
                     name= "Sex")+
  labs(x='Volume liver in L', y='Percentage metabolised', title='intersex differences human')+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "top",
        title = element_text(size=20))
  
ggplotly(p)
  
Q_C_check <-metabolised_data_frame
Q_C_check[,2] <- parameters[,26]
colnames(Q_C_check) <-c("percentage","Q_C")

melt_Q_C_check<-melt(Q_C_check,id=c("percentage"))
melt_Q_C_check[,4]<-1:2000
colnames(melt_Q_C_check) <-c("percentage","variable","value","id")

melt_Q_C_check$id[melt_Q_C_check$id == 1:1000] <- "Male"  
melt_Q_C_check$id[melt_Q_C_check$id == 1001:2000] <- "Female" 

p<-ggplot(melt_Q_C_check, aes(x=value, y=percentage,colour=id)) +
  geom_point()+
  scale_color_manual(values = c( "Male" = "blue",
                                 "Female" = "red"),
                     labels= c( "Male", "Female"),
                     name= "Sex")+
  labs(x='Cardiac output', y='Percentage metabolised', title='Intrasex differences human')+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "top",
        title = element_text(size=20))

ggplotly(p)

