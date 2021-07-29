library(dplyr)
library(deSolve)

setwd("~/github/fisheries_models")

create_schnute_data <- function(path, N, r,k,B0,q,E,sigma_obs_C,sigma_obs_E){
  
  
  # define model 
  dBdt <- function(t,x,p){
    dB = (r*x[1]*(1-x[1]/K) - q*E(t)*x[1])
    dC =  q*E(t)*x[1]
    dE = E(t)
    list(B = dB, C = dC, E = dE)
  } 
  
  # run ode
  interval <- 0.01
  times <- seq(0, N, by = interval)
  out <- ode(y = B0, times = times, func = dBdt, parms = c())

  ## integrate over time steps of one year and add observaiton errors
  C_ls <- c()
  E_ls <- c()
  for(i in 1:N){
    C <- out[,"C"][times <= i & times > i-1 ]
    E <- out[,"E"][times <= i & times > i-1 ]
    C_ls <- append(C_ls, sum(C)*interval + rnorm(1,0,sigma_obs_C))
    E_ls <- append(E_ls, sum(E)*interval + rnorm(1,0,sigma_obs_E))
  }
  
  data = data.frame(E = E_ls[2:N], C = C_ls[2:N])
  plot(data$E, data$C/data$E)
  write.csv(data, path)
}



r <- 3
K <- 0.95
B0 <- K
q <- 0.5
E <- function(t) t*4.0/30

## define parametes 
sigma_obs_C <- 0.005
sigma_obs_E <- 0.005
N <- 30
create_schnute_data("data/sim_schnute_1.csv",N,r,K,B0,q,E,sigma_obs_C,sigma_obs_E)





# up and down the isocline 
r <- 1.0
K <- 2.0
B0 <- K
q <- 0.25
E <- function(t){ 
  if(t < 15){
    return(t*4.0/30)
  }else{
    return(15*4.0/30-(t-15)*4.0/30 )
  }
}

## define parametes 
sigma_obs_C <- 0.001
sigma_obs_E <- 0.001
N <- 30
create_schnute_data("data/sim_schnute_2.csv",N,r,K,B0,q,E,sigma_obs_C,sigma_obs_E)





# up and down the isocline 
r <- 1.0
K <- 2.0
B0 <- K
q <- 0.25
E <- function(t){ 
  if(t < 15){
    return(t*4.0/30)
  }else{
    return(15*4.0/30-(t-15)*4.0/30 )
  }
}

## define parametes 
sigma_obs_C <- 0.005
sigma_obs_E <- 0.005
N <- 30
create_schnute_data("data/sim_schnute_3.csv",N,r,K,B0,q,E,sigma_obs_C,sigma_obs_E)






# up and down the isocline 
r <- 1.0
K <- 3
B0 <- K
q <- 0.2

E <- function(t) 7.0*sin(t)^2
## define parametes 
sigma_obs_C <- 0.0005
sigma_obs_E <- 0.00001
N <- 100
create_schnute_data("data/sim_schnute_4.csv",N,r,K,B0,q,E,sigma_obs_C,sigma_obs_E)


