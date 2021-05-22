## JAGS Code for Implementing One-Step Ahead Forecast With RS-ZIMLP Models
## Author: Yanling Li
## Date: May 21th, 2021
library(rjags)
source("postcalc.R")

# set-up 
r1 = 1  # seed
N = 200 # number of persons 
O = 60 # number of time points
condition = "moderate" # simulation conditions
#condition = "high"
K = 5 # length of the forecast window

# load simulated data set and prepare a data list for JAGS
load(paste0("SimulatedData_RSZIMLP_",condition,"_N",N,"T",O,"_",r1,".Rdata"))


jags_data <- list(Y = Y_obs, 
                  P = dim(Y)[1],
                  maxT = dim(Y)[2],
                  X = X_obs,
                  K = K)

modelString = "
model{
  for (pp in 1:P) {
    
    # random intercepts 
    intercept[pp] ~ dnorm(int0, tau_int)
    intercept_pred[pp] <- intercept[pp]

    # 1st observation
    odds[pp,1,1] <- 1
    odds[pp,1,2] <- exp(pi0)
    
    midx[pp,1]~dcat(odds[pp,1,])
    
    eta[pp,1] ~ dnorm(0,0.01)
    
    Y[pp,1] ~ dpois((midx[pp,1]-1)*exp(eta[pp,1])+0.00001)
    
    X[pp,1,1] ~ dnorm(0,0.01)
    X[pp,1,2] ~ dnorm(0,0.01)
    
    # model fitting based on data from up to time maxT-K
    for(tt in 2:(maxT-K)){        	    
      
      # from regime midx[pp,tt-1] to regime 1 (i.e., ZI regime)
      odds[pp,tt,1] <- exp(alpha[1,1,midx[pp,tt-1]]+alpha[2,1,midx[pp,tt-1]]*X[pp,tt-1,1])
      # from regime midx[pp,tt-1] to regime 2 (i.e., drinking regime)
      odds[pp,tt,2] <- exp(alpha[1,2,midx[pp,tt-1]]+alpha[2,2,midx[pp,tt-1]]*X[pp,tt-1,1])
      
      midx[pp,tt]~dcat(odds[pp,tt,])

      mu_eta[pp,tt] <- intercept[pp] + AR * (eta[pp,tt-1] - intercept[pp]) + coeffX * X[pp,tt-1,2]
      eta[pp,tt] ~ dnorm(mu_eta[pp,tt], tau_noise)
      Y[pp,tt] ~ dpois((midx[pp,tt]-1)*exp(eta[pp,tt])+0.00001)
      
      X[pp,tt,1] ~ dnorm(ARX[1]*X[pp,tt-1,1], tauX[1])
      X[pp,tt,2] ~ dnorm(ARX[2]*X[pp,tt-1,2], tauX[2])
    } # close loop over time points
    
    
    # implement one-step ahead forecast
    for(tt in (maxT-K+1):maxT){

      # from regime midx[pp,tt-1] to regime 1 (i.e., ZI regime)
      odds[pp,tt,1] <- exp(alpha_pred[1,1,midx[pp,tt-1]]+alpha_pred[2,1,midx[pp,tt-1]]*X[pp,tt-1,1])
      # from regime midx[pp,tt-1] to regime 2 (i.e., drinking regime)
      odds[pp,tt,2] <- exp(alpha_pred[1,2,midx[pp,tt-1]]+alpha_pred[2,2,midx[pp,tt-1]]*X[pp,tt-1,1])

      midx[pp,tt]~dcat(odds[pp,tt,])  # estimated regime indicator
      midx_pred[pp,tt]~dcat(odds[pp,tt,]) # predictive regime indicator

      mu_eta[pp,tt] <- intercept_pred[pp] + AR_pred * (eta[pp,tt-1] - intercept_pred[pp]) + coeffX_pred * X[pp,tt-1,2]

      eta[pp,tt] ~ dnorm(mu_eta[pp,tt], tau_noise_pred) # estimated eta
      eta_pred[pp,tt] ~ dnorm(mu_eta[pp,tt], tau_noise_pred) # predictive eta

      Y[pp,tt] ~ dpois((midx[pp,tt]-1)*exp(eta[pp,tt])+0.00001) # observed Y
      Y_pred[pp,tt] ~ dpois((midx_pred[pp,tt]-1)*exp(eta_pred[pp,tt])+0.00001) # predictive Y

      X[pp,tt,1] ~ dnorm(ARX_pred[1]*X[pp,tt-1,1], tauX_pred[1])
      X[pp,tt,2] ~ dnorm(ARX_pred[2]*X[pp,tt-1,2], tauX_pred[2])
    } # close loop over time points

  } # close loop over persons
  
  # priors
  # parameters in the multilevel AR-X model 
  int0 ~ dnorm(0, 0.01) 
  tau_int ~ dgamma(0.001,0.001)
  int_sd <- pow(tau_int,-1/2)
  AR ~ dnorm(0, 1)
  coeffX ~ dnorm(0, 0.01)
  tau_noise ~ dgamma(0.001,0.001)
  sd_noise <- pow(tau_noise,-1/2)

  # parameters in the initial regime model
  pi0 ~ dnorm(0, 0.01)
  # parameters in the RS model
  for(i in 1:2){
    alpha[i,1,1] <- 0
    alpha[i,1,2] ~ dnorm(0,0.01)
    alpha[i,2,1] ~ dnorm(0,0.01)
    alpha[i,2,2] <- 0
  }

  # paramters in the AR model for covariates
  for(i in 1:2){
  ARX[i] ~ dnorm(0,1)
  tauX[i] ~ dgamma(0.001,0.001)
  sd_noiseX[i] <- pow(tauX[i], -1/2)
  }
  
  # stop updating hyper-parameters before the forecast window by 
  # saving hyper-parameters estimated based on data up to time maxT-K into new variables
  AR_pred <- AR
  coeffX_pred <- coeffX
  tau_noise_pred <- tau_noise

  for(i in 1:2){
    alpha_pred[i,1,1] <- 0
    alpha_pred[i,1,2] <- alpha[i,1,2]
    alpha_pred[i,2,1] <- alpha[i,2,1]
    alpha_pred[i,2,2] <- 0
  }

  for(i in 1:2){
    ARX_pred[i] <- ARX[i]
    tauX_pred[i] <- tauX[i]
  }

} # end of model
"

writeLines(modelString, con = "rszimlp.txt")


inits1 <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = r1)
inits2 <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = r1+500)

jagsModel <- jags.model(file = "rszimlp.txt", data = jags_data,
                        inits = list(inits1,inits2),
                        n.chains = 2, n.adapt = 4000) 
update(jagsModel, n.iter = 1000)

parameterlist <-c("int0","int_sd","AR","coeffX","sd_noise",
                  "pi0","alpha","ARX","sd_noiseX", 
                  "midx","Y","Y_pred")

codaSamples <- coda.samples(jagsModel, 
                            variable.names = parameterlist,
                            n.iter = 20000)
resulttable <- zcalc(codaSamples)
warns = warnings()

save(warns, resulttable, file=paste0("Result_RSZIMLP_",condition,"_N",N,"T",O,"_",r1,".Rdata"))
