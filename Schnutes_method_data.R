library(dplyr)
library(deSolve)

setwd("~/github/fisheries_models")

create_schnute_data <- function(path, r,k,B0,q,E,sigma_obs_C,sigma_obs_E){
  
  
  # define model 
  dBdt <- function(t,x,p){
    dB = (r*x[1]*(1-x[1]/K) - q*E(t)*x[1])
    dC =  q*E(t)*x[1]
    dE = E(t)
    list(B = dB, C = dC, E = dE)
  } 
  
  # run ode
  interval <- 0.01
  times <- seq(0, 30, by = interval)
  out <- ode(y = B0, times = times, func = dBdt, parms = c())

  ## integrate over time steps of one year and add observaiton errors
  C_ls <- c()
  E_ls <- c()
  for(i in 1:30){
    C <- out[,"C"][times <= i & times > i-1 ]
    E <- out[,"E"][times <= i & times > i-1 ]
    C_ls <- append(C_ls, sum(C)*interval + rnorm(1,0,sigma_obs_C))
    E_ls <- append(E_ls, sum(E)*interval + rnorm(1,0,sigma_obs_E))
  }
  
  data = data.frame(E = E_ls[2:30], C = C_ls[2:30])
  plot(data$E, data$C/data$E)
  write.csv(data, path)
}



r <- 3
K <- 1
B0 <- K
q <- 0.25
E <- function(t) t*4.0/30

## define parametes 
sigma_obs_C <- 0.01
sigma_obs_E <- 0.01
create_schnute_data("data/sim_schnute_1.csv",r,K,B0,q,E,sigma_obs_C,sigma_obs_E)




r <- 1
K <- 1
B0 <- K
q <- 0.25
E <- function(t) t*1.0/30
sigma_obs_C <- 0.01
sigma_obs_E <- 0.001
create_schnute_data("data/sim_schnute_2.csv",r,K,B0,q,E,sigma_obs_C,sigma_obs_E)






