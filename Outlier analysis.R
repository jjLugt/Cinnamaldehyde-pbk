

Liver_check <-as.data.frame(parameters[,18])
Liver_check[,2]<- 1:2000
Liver_check[,3]<-metabolised_data_frame
Liver_check[,4]<- parameters[,49]
Liver_check[,5]<- parameters[,17]
colnames(Liver_check) <-c("V_L","id","percentage","vsmax","fat")




orderd_liver <-arrange(Liver_check,desc(percentage))

#Liver_check[,3] <- orderd_auc_cb$PPORRES


p<-ggplot(Liver_check, aes(x=V_L, y=percentage)) +
  geom_point() +
  labs(x='Liver', y='Percentage metabolised', title='Liver weight/percentage metabolized to carboxylic acid metabolites')+
  theme_classic()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        title = element_text(size=20))
p
unique(orderd_liver[orderd_liver$id == "1018",c("id","V_L")])
unique(orderd_auc_cb[orderd_auc_cb$id == "1018",c("id","PPORRES")])


hist(var_f_pop$V_L)
hist(var_m_pop$V_L)

