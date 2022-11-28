#author: Joris Jean van der Lugt
#date: 26-08-2021
#SA data extraction
library(RxODE)
library(tidyverse)
library(readxl)
library(readr)
library(truncnorm)
library(reshape2)
library(sensitivity)
library(PKNCA)



#importing created PBK results files and stripping both the time and C_V colum from them
solve <- read_csv("solve.pbk_nonpop1")

solve.pbk.sa <-as.data.frame(solve$time)
solve.pbk.sa <-cbind(solve.pbk.sa,solve$C_Pu)
colnames(solve.pbk.sa) <- c("time","C_Pu")

solve2 <- read_csv("solve.pbk_nonpop2")
solve.pbk.sa2 <-as.data.frame(solve2$time)
solve.pbk.sa2 <-cbind(solve.pbk.sa2,solve2$C_Pu)
colnames(solve.pbk.sa2) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa2)

solve3 <- read_csv("solve.pbk_nonpop3")
solve.pbk.sa3 <-as.data.frame(solve3$time)
solve.pbk.sa3 <-cbind(solve.pbk.sa3,solve3$C_Pu)
colnames(solve.pbk.sa3) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa3)


solve4 <- read_csv("solve.pbk_nonpop4")
solve.pbk.sa4 <-as.data.frame(solve4$time)
solve.pbk.sa4 <-cbind(solve.pbk.sa4,solve4$C_Pu)
colnames(solve.pbk.sa4) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa4)

solve5 <- read_csv("solve.pbk_nonpop5")
solve.pbk.sa5 <-as.data.frame(solve5$time)
solve.pbk.sa5 <-cbind(solve.pbk.sa5,solve5$C_Pu)
colnames(solve.pbk.sa5) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa5)

solve6 <- read_csv("solve.pbk_nonpop6")
solve.pbk.sa6 <-as.data.frame(solve6$time)
solve.pbk.sa6 <-cbind(solve.pbk.sa6,solve6$C_Pu)
colnames(solve.pbk.sa6) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa6)


solve7 <- read_csv("solve.pbk_nonpop7")
solve.pbk.sa7 <-as.data.frame(solve7$time)
solve.pbk.sa7 <-cbind(solve.pbk.sa7,solve7$C_Pu)
colnames(solve.pbk.sa7) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa7)


solve8 <- read_csv("solve.pbk_nonpop8")
solve.pbk.sa8 <-as.data.frame(solve8$time)
solve.pbk.sa8 <-cbind(solve.pbk.sa8,solve8$C_Pu)
colnames(solve.pbk.sa8) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa8)

solve9 <- read_csv("solve.pbk_nonpop9")
solve.pbk.sa9 <-as.data.frame(solve9$time)
solve.pbk.sa9 <-cbind(solve.pbk.sa9,solve9$C_Pu)
colnames(solve.pbk.sa9) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa9)


solve10 <- read_csv("solve.pbk_nonpop10")
solve.pbk.sa10 <-as.data.frame(solve10$time)
solve.pbk.sa10 <-cbind(solve.pbk.sa10,solve10$C_Pu)
colnames(solve.pbk.sa10) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa10)


solve11 <- read_csv("solve.pbk_nonpop11")
solve.pbk.sa11 <-as.data.frame(solve11$time)
solve.pbk.sa11 <-cbind(solve.pbk.sa11,solve11$C_Pu)
colnames(solve.pbk.sa11) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa11)

solve12 <- read_csv("solve.pbk_nonpop12")
solve.pbk.sa12 <-as.data.frame(solve12$time)
solve.pbk.sa12 <-cbind(solve.pbk.sa12,solve12$C_Pu)
colnames(solve.pbk.sa12) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa12)

solve13 <- read_csv("solve.pbk_nonpop13")
solve.pbk.sa13 <-as.data.frame(solve13$time)
solve.pbk.sa13 <-cbind(solve.pbk.sa13,solve13$C_Pu)
colnames(solve.pbk.sa13) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa13)

solve14 <- read_csv("solve.pbk_nonpop14")
solve.pbk.sa14 <-as.data.frame(solve14$time)
solve.pbk.sa14 <-cbind(solve.pbk.sa14,solve14$C_Pu)
colnames(solve.pbk.sa14) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa14)

solve15 <- read_csv("solve.pbk_nonpop15")
solve.pbk.sa15 <-as.data.frame(solve15$time)
solve.pbk.sa15 <-cbind(solve.pbk.sa15,solve15$C_Pu)
colnames(solve.pbk.sa15) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa15)

solve16 <- read_csv("solve.pbk_nonpop16")
solve.pbk.sa16 <-as.data.frame(solve16$time)
solve.pbk.sa16 <-cbind(solve.pbk.sa16,solve16$C_Pu)
colnames(solve.pbk.sa16) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa16)

solve17 <- read_csv("solve.pbk_nonpop17")
solve.pbk.sa17 <-as.data.frame(solve17$time)
solve.pbk.sa17 <-cbind(solve.pbk.sa17,solve17$C_Pu)
colnames(solve.pbk.sa17) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa17)


solve18 <- read_csv("solve.pbk_nonpop18")
solve.pbk.sa18 <-as.data.frame(solve18$time)
solve.pbk.sa18 <-cbind(solve.pbk.sa18,solve18$C_Pu)
colnames(solve.pbk.sa18) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa18)

solve19 <- read_csv("solve.pbk_nonpop19")
solve.pbk.sa19 <-as.data.frame(solve19$time)
solve.pbk.sa19 <-cbind(solve.pbk.sa19,solve19$C_Pu)
colnames(solve.pbk.sa19) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa19)

solve20 <- read_csv("solve.pbk_nonpop20")
solve.pbk.sa20 <-as.data.frame(solve20$time)
solve.pbk.sa20 <-cbind(solve.pbk.sa20,solve20$C_Pu)
colnames(solve.pbk.sa20) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa20)

solve21 <- read_csv("solve.pbk_nonpop21")
solve.pbk.sa21 <-as.data.frame(solve21$time)
solve.pbk.sa21 <-cbind(solve.pbk.sa21,solve21$C_Pu)
colnames(solve.pbk.sa21) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa21)

solve22 <- read_csv("solve.pbk_nonpop22")
solve.pbk.sa22 <-as.data.frame(solve22$time)
solve.pbk.sa22 <-cbind(solve.pbk.sa22,solve22$C_Pu)
colnames(solve.pbk.sa22) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa22)

solve23 <- read_csv("solve.pbk_nonpop23")
solve.pbk.sa23 <-as.data.frame(solve23$time)
solve.pbk.sa23 <-cbind(solve.pbk.sa23,solve23$C_Pu)
colnames(solve.pbk.sa23) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa23)

solve24 <- read_csv("solve.pbk_nonpop24")
solve.pbk.sa24 <-as.data.frame(solve24$time)
solve.pbk.sa24 <-cbind(solve.pbk.sa24,solve24$C_Pu)
colnames(solve.pbk.sa24) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa24)

solve25 <- read_csv("solve.pbk_nonpop25")
solve.pbk.sa25 <-as.data.frame(solve25$time)
solve.pbk.sa25 <-cbind(solve.pbk.sa25,solve25$C_Pu)
colnames(solve.pbk.sa25) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa25)

solve26 <- read_csv("solve.pbk_nonpop26")
solve.pbk.sa26 <-as.data.frame(solve26$time)
solve.pbk.sa26 <-cbind(solve.pbk.sa26,solve26$C_Pu)
colnames(solve.pbk.sa26) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa26)

solve27 <- read_csv("solve.pbk_nonpop27")
solve.pbk.sa27 <-as.data.frame(solve27$time)
solve.pbk.sa27 <-cbind(solve.pbk.sa27,solve27$C_Pu)
colnames(solve.pbk.sa27) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa27)

solve28 <- read_csv("solve.pbk_nonpop28")
solve.pbk.sa28 <-as.data.frame(solve28$time)
solve.pbk.sa28 <-cbind(solve.pbk.sa28,solve28$C_Pu)
colnames(solve.pbk.sa28) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa28)


solve29 <- read_csv("solve.pbk_nonpop29")
solve.pbk.sa29 <-as.data.frame(solve29$time)
solve.pbk.sa29 <-cbind(solve.pbk.sa29,solve29$C_Pu)
colnames(solve.pbk.sa29) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa29)

solve30 <- read_csv("solve.pbk_nonpop30")
solve.pbk.sa30 <-as.data.frame(solve30$time)
solve.pbk.sa30 <-cbind(solve.pbk.sa30,solve30$C_Pu)
colnames(solve.pbk.sa30) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa30)

solve31 <- read_csv("solve.pbk_nonpop31")
solve.pbk.sa31 <-as.data.frame(solve31$time)
solve.pbk.sa31 <-cbind(solve.pbk.sa31,solve31$C_Pu)
colnames(solve.pbk.sa31) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa31)


solve32 <- read_csv("solve.pbk_nonpop32")
solve.pbk.sa32 <-as.data.frame(solve32$time)
solve.pbk.sa32 <-cbind(solve.pbk.sa32,solve32$C_Pu)
colnames(solve.pbk.sa32) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa32)


solve33 <- read_csv("solve.pbk_nonpop33")
solve.pbk.sa33 <-as.data.frame(solve33$time)
solve.pbk.sa33 <-cbind(solve.pbk.sa33,solve33$C_Pu)
colnames(solve.pbk.sa33) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa33)

solve34 <- read_csv("solve.pbk_nonpop34")
solve.pbk.sa34 <-as.data.frame(solve34$time)
solve.pbk.sa34 <-cbind(solve.pbk.sa34,solve34$C_Pu)
colnames(solve.pbk.sa34) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa34)


solve35 <- read_csv("solve.pbk_nonpop35")
solve.pbk.sa35 <-as.data.frame(solve35$time)
solve.pbk.sa35 <-cbind(solve.pbk.sa35,solve35$C_Pu)
colnames(solve.pbk.sa35) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa35)


solve36 <- read_csv("solve.pbk_nonpop36")
solve.pbk.sa36 <-as.data.frame(solve36$time)
solve.pbk.sa36 <-cbind(solve.pbk.sa36,solve36$C_Pu)
colnames(solve.pbk.sa36) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa36)

solve37 <- read_csv("solve.pbk_nonpop37")
solve.pbk.sa37 <-as.data.frame(solve37$time)
solve.pbk.sa37 <-cbind(solve.pbk.sa37,solve37$C_Pu)
colnames(solve.pbk.sa37) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa37)

solve38 <- read_csv("solve.pbk_nonpop38")
solve.pbk.sa38 <-as.data.frame(solve38$time)
solve.pbk.sa38 <-cbind(solve.pbk.sa38,solve38$C_Pu)
colnames(solve.pbk.sa38) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa38)

solve39 <- read_csv("solve.pbk_nonpop39")
solve.pbk.sa39 <-as.data.frame(solve39$time)
solve.pbk.sa39 <-cbind(solve.pbk.sa39,solve39$C_Pu)
colnames(solve.pbk.sa39) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa39)

solve40 <- read_csv("solve.pbk_nonpop40")
solve.pbk.sa40 <-as.data.frame(solve40$time)
solve.pbk.sa40 <-cbind(solve.pbk.sa40,solve40$C_Pu)
colnames(solve.pbk.sa40) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa40)

solve41 <- read_csv("solve.pbk_nonpop41")
solve.pbk.sa41 <-as.data.frame(solve41$time)
solve.pbk.sa41 <-cbind(solve.pbk.sa41,solve41$C_Pu)
colnames(solve.pbk.sa41) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa41)


solve42 <- read_csv("solve.pbk_nonpop42")
solve.pbk.sa42 <-as.data.frame(solve42$time)
solve.pbk.sa42 <-cbind(solve.pbk.sa42,solve42$C_Pu)
colnames(solve.pbk.sa42) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa42)

solve43 <- read_csv("solve.pbk_nonpop43")
solve.pbk.sa43 <-as.data.frame(solve43$time)
solve.pbk.sa43 <-cbind(solve.pbk.sa43,solve43$C_Pu)
colnames(solve.pbk.sa43) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa43)

solve44 <- read_csv("solve.pbk_nonpop44")
solve.pbk.sa44 <-as.data.frame(solve44$time)
solve.pbk.sa44 <-cbind(solve.pbk.sa44,solve44$C_Pu)
colnames(solve.pbk.sa44) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa44)

solve45 <- read_csv("solve.pbk_nonpop45")
solve.pbk.sa45 <-as.data.frame(solve45$time)
solve.pbk.sa45 <-cbind(solve.pbk.sa45,solve45$C_Pu)
colnames(solve.pbk.sa45) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa45)

solve46 <- read_csv("solve.pbk_nonpop46")
solve.pbk.sa46 <-as.data.frame(solve46$time)
solve.pbk.sa46 <-cbind(solve.pbk.sa46,solve46$C_Pu)
colnames(solve.pbk.sa46) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa46)

solve47 <- read_csv("solve.pbk_nonpop47")
solve.pbk.sa47 <-as.data.frame(solve47$time)
solve.pbk.sa47 <-cbind(solve.pbk.sa47,solve47$C_Pu)
colnames(solve.pbk.sa47) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa47)

solve48 <- read_csv("solve.pbk_nonpop48")
solve.pbk.sa48 <-as.data.frame(solve48$time)
solve.pbk.sa48 <-cbind(solve.pbk.sa48,solve48$C_Pu)
colnames(solve.pbk.sa48) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa48)

solve49 <- read_csv("solve.pbk_nonpop49")
solve.pbk.sa49 <-as.data.frame(solve49$time)
solve.pbk.sa49 <-cbind(solve.pbk.sa49,solve49$C_Pu)
colnames(solve.pbk.sa49) <- c("time","C_Pu")
solve.pbk.sa <-rbind(solve.pbk.sa,solve.pbk.sa49)


write.csv(solve.pbk.sa,"D:/PBK/Cinnamaldehyde-pbk\\SA_rat_250mg_oral_C_Pu_corrected", row.names = TRUE)

