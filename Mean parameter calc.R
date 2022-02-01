#Simulations
set.seed(15204)         #to ensure a reproducible output
amount.units <-"umol"
time.units   <-"h"
nbr.doses    <-1        #number of doses
time.0       <-0        #time start dosing
time.end     <-8        #time end of simulation
time.frame   <-0.01     #time steps of simulation
N           <-1000     #Number of males
NF          <-1000     #Number of females
Dose_in_mg   <-250      #Dose in mg/kg-bw
MW           <-132.16   #The molecular weight of Cinnamaldehyde

colnames <-c("Age","Height_start","Height_cv","Height","BW_start","BW_cv","BW","BSA","V_L","V_F","V_F_min",
             "V_B","V_A","V_V","V_SI","V_RP","V_SP","Q_C","Q_F","Q_L","Q_SI","Q_RP","Q_SP")

par_var_m <- length(colnames)

#create data frames for population males
var_m <- matrix(NA, nrow = N, ncol = par_var_m)
colnames(var_m) <- colnames
var_m <- as.data.frame(var_m)



Age                    <- 30                                       #Age (years)
var_m$Age              <- Age
var_m$Height_start     <- 180    #Body height baseline (cm)
#var_m$Height_cv        <- rnorm(N,0,0.039)                                     #Variation in body height
var_m$Height           <- 180                                                   #Body height (cm)
var_m$BW_start         <- 70                    #Body weight baseline (kg)
#var_m$BW_cv            <- rnorm(N,0,0.15)                                      #Variation in body weight
var_m$BW               <- 70                   #Body weight (kg)
var_m$BSA              <- 0.007184 * var_m$Height^0.725 * var_m$BW^0.425       #Body surface area (m2)

#-Tissues volumes-#

var_m$V_L    <- (1072.8 * (var_m$BSA) - 345.7) / 1000                             #Volume liver tissue (l)
var_m$V_F  <- (1.36 * var_m$BW)/(var_m$Height/100)-42                         #Volume adipose tissue (L)
var_m$V_F_min  <- 0.05 * var_m$BW                                           #Minimum of adipose tissue should be at least 5% of body weight
var_m$V_F      <- ifelse(var_m$V_F < var_m$V_F_min, var_m$V_F_min, var_m$V_F)      #To ensure that adipose tissue is at least 5% of body weight
var_m$V_B       <-(((13.1 * var_m$Height + 18.05 * var_m$BW - 480) / 0.5723) / 1000)         #Volume blood (L)
var_m$V_A       <-var_m$V_B / 3                                                   #Volume arterial blood (L)
var_m$V_V       <-var_m$V_B * (2/3)                                               #Volume venous blood (L) 
var_m$V_SI     <-0.021 * (var_m$BW - var_m$V_F * 0.92) / 1.05                           #Volume gut tissue (L)
var_m$V_RP     <-(2.331 * 10^-3 * var_m$Age + 0.1253 * var_m$BW^0.8477 + var_m$Height^0.3821 - 4.725) - var_m$V_SI - var_m$V_L   #Volume richly perfused tissue (L)
var_m$V_SP     <-var_m$BW - var_m$V_B - var_m$V_RP -var_m$V_SI - var_m$V_L - var_m$V_F  #Volume slowly perfused tissue (L)

#-Cardiac parameters-#

var_m$Q_C           <- var_m$BSA * 60 * (3 - 0.01 * (var_m$Age - 20))           #Cardiac output (L/h)
var_m$Q_SI          <- var_m$Q_C * 0.15                                         #Blood flow to the gut (L/h)
var_m$Q_F           <- var_m$Q_C * 0.05                                         #Blood flow to adipose tissue (L/h)
var_m$Q_L           <- var_m$Q_C * 0.065                                        #Blood flow to liver via hepatic artery (L/h)
var_m$Q_RP          <- 0.626 * var_m$Q_C - var_m$Q_SI - var_m$Q_L               #Blood flow to richly perfused tissue (L/h)
var_m$Q_SP          <- 0.374 * var_m$Q_C - var_m$Q_F                            #Blood flow to slowly perfused tissue (L/h)



mean_BW <- mean(var_m$BW)
mean_height <-mean(var_m$Height)
mean_Q_C <- mean(var_m$Q_C)
mean_Q_F <- mean(var_m$Q_F)
mean_Q_L <- mean(var_m$Q_L)
mean_Q_SI <- mean(var_m$Q_SI)
mean_Q_RP <- mean(var_m$Q_RP)
mean_Q_SP <- mean(var_m$Q_SP)
mean_V_L <- mean(var_m$V_L)
mean_V_V <-mean(var_m$V_V)
mean_V_A <-mean(var_m$V_A)
mean_V_F <-mean(var_m$V_F)
mean_V_SI <-mean(var_m$V_SI)
mean_V_RP <-mean(var_m$V_RP)
mean_V_SP <-mean(var_m$V_SP)





