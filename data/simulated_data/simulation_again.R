


###############  Simulation Study ################################
###############  Continuous Outcome ##############################
###############  Hengshi Yu ######################################
library(mvtnorm)
beta_0 <- 20
beta_trt <- -0.4007
beta_age <- -0.2400
rho <- 0.4278
rho_min <- 0.3
rho_max <- 0.55
rho_std <- 0.1346
n <- 100000
ptrt <- 0.525
mean_age <- 28
var_age <- 55
min_age <- 18
max_age <- 57
beta_1 <- 0.1
beta_2 <- -0.5
beta_3 <- 0.23
mean_1 <- 10
var_1 <- 1
low_2 <- 0
up_2 <- 4
p3 <- 0.34
sigmai <- 344.9

simu <- function(n, m){
  Z <- matrix(NA,n,m)     ## age groups of n*m individuals in n clusters
  mu <- matrix(NA,n,m)    ## means of n*m individuals in n clusters
  y <- matrix(NA,n,m)     ## results of n*m individuals in n clusters
  X <- matrix(NA,n,m)
  X1 <- matrix(NA,n,m)
  X2 <- matrix(NA,n,m)
  X3 <- matrix(NA,n,m)
  R <- matrix(NA, m,m)
  A <- diag(sigmai, m)

for(i in 1:n){
  X[i, ] <- rbinom(1, 1, ptrt)
    for(j in 1:m){
      Z[i,j] <- rnorm(1, mean_age, var_age)
      Z[i,j] <- (Z[i,j] < min_age) * min_age + (Z[i,j] > max_age) * max_age + Z[i,j] * (Z[i,j]>= min_age) * (Z[i,j] <= max_age)  
      X1[i,j] <- rnorm(1, mean_1, var_1)
      X2[i,j] <- runif(1, low_2, up_2)
      X3[i,j] <- rbinom(1,1, p3)
      mu[i,j] <- beta_0  + beta_1 * X1[i,j] + beta_2 * X2[i,j] + beta_3 * X3[i,j] +  beta_trt * X[i,j] + beta_age * Z[i,j]
    }
    rho_value <- rho
    #rho_value <- rho_value * (rho_value <= rho_max) *(rho_value >= rho_min) + rho_min * (rho_value < rho_min) + rho_max * (rho_value > rho_max)
    for(j in 1:m){
      for(k in 1:m){
        R[j,k]<-ifelse(j==k, 1, rho_value)  
      }
    }
    Variance <- A^(1/2)%*% R %*% A^(1/2)
    y[i, ] <- rmvnorm(1,mean=mu[i,],sigma=Variance)
   }
  idmatrix <- matrix(NA,n,m)
  for(i in 1:n){
    idmatrix[i,] <- i
  }


  
  yob<-as.vector(t(y)) ## transpose the outcome matrix to be a vector
  treatment<-as.vector(t(X))
  age <- as.vector(t(Z))
  variable_1 <- as.vector(t(X1))
  variable_2 <- as.vector(t(X2))
  variable_3 <- as.vector(t(X3))
  id <- as.vector(t(idmatrix))
  obs<-data.frame(id, variable_1, variable_2, variable_3, treatment, age, yob)
  return(obs = obs)
}
n <- 1000000
m <- 5
npart = 100
usedata <- simu(n,m)

parts <- n/npart
for(i in 1:parts){
  write.csv(usedata, file = paste0("./output/new_simulated_", i,  "_data.csv"), row.names = FALSE)  
}
