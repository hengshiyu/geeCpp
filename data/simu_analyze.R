library(gee)
library(geeM)
library(geepack)
load("new_simulated_1_I.RData")
n <- 10000
m <- 5
nm = n*m
data_use <- usedata[1:nm, ]

gee_result <- gee(formula = yob ~variable_1 + variable_2 + variable_3 + treatment + age, id = id, data = data_use, family =gaussian, corstr="exchangeable")
geeM_result <- geem(formula = yob ~variable_1 + variable_2 + variable_3 + treatment + age, id = id, data = data_use, family =gaussian, corstr="exchangeable")
geeglm_result <-  geeglm(formula = yob ~variable_1 + variable_2 + variable_3 + treatment + age, id = id, data = data_use, family =gaussian, corstr="exchangeable")

# gee
summary(gee_result)
# geeM
summary(geeM_result)
#geeglm
summary(geeglm_result)

time_gee <- c()
time_geeM <- c()
time_geeglm <- c()
steps <- c(500, 1500, 2500, 10000, 300000)
for(i in 1:length(steps)){
  n <- steps[i]
  m <- 5
  nm = n*m
  data_use <- usedata[1:nm, ]
  
  time_gee <-  c(time_gee, as.numeric(system.time(gee_result <- gee(formula = yob ~variable_1 + variable_2 + variable_3 + treatment + age, id = id, data = data_use, family =gaussian, corstr="exchangeable"))[3]) )
  time_geeM <- c(time_geeM, as.numeric(system.time(geeM_result <- geem(formula = yob ~variable_1 + variable_2 + variable_3 + treatment + age, id = id, data = data_use, family =gaussian, corstr="exchangeable"))[3]))
  time_geeglm <- c(time_geeglm, as.numeric(system.time(geeglm_result <-  geeglm(formula = yob ~variable_1 + variable_2 + variable_3 + treatment + age, id = id, data = data_use, family =gaussian, corstr="exchangeable"))[3]))
}

datacollect <- cbind(steps, time_gee, time_geeM, time_geeglm)
colnames(datacollect) <- c("n", "gee", "geeM (geem)", "geepack (geeglm)")