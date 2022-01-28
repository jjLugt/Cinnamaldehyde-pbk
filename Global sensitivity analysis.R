#install libraries
library("sensitivity")

#names for the lists
colnames <- c(colnames(phys))
par_var <- length(colnames)

Mean <- phys[1,]
Lower <- Mean - 0.1*Mean
Upper <- Mean + 0.1*Mean

#create data frames for population
n_sim  <- 1000                #number of iterations
X1 <- matrix(NA, nrow = n_sim, ncol = par_var)
colnames(X1) <- colnames
X1 <- as.data.frame(X1)
var <- X1

X2 <- matrix(NA, nrow = n_sim, ncol = par_var)
colnames(X2) <- colnames
X2 <- as.data.frame(X2)
var <- X2

#create distribution population
for(i in 1:par_var){
  X1[,i] <- runif(n_sim, min = Lower[,i], max = Upper[,i])
  X2[,i] <- runif(n_sim, min = Lower[,i], max = Upper[,i])
}

n_boot <- 1000

#Sobol design
sa <- soboljansen(model = NULL, X1, X2, nboot = n_boot, conf = 0.95)
phys <- sa$X

