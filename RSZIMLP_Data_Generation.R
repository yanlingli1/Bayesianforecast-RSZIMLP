## Data Generation Code
## Author: Yanling Li
## Date: May 21th, 2021

# set-up 
r1 = 1  # seed
N = 200 # number of persons 
O = 60 # number of time points
condition = "moderate" # simulation conditions
#condition = "high"

set.seed(r1)

# paramters in the multilevel AR-X model
int0 = 2 # population mean of random intercepts
int_sd = 0.5 # standard deviation of random intercepts
AR = 0.3 # auto-regression parameter
coeffX = 0.5 # coefficient of the covariate (i.e., x2)
sd_noise = 0.5 # standard deviation of the process noise

# parameters in the initial regime model
pi0 = -2
# parameters in the RS model
if(condition == "moderate"){
  alpha21 = -2.5 
}
if(condition == "high"){
  alpha21 = -3.5
}

alpha = array(rep(NA, 2*2*2), dim=c(2,2,2))
alpha[1,1,1] = 0 
alpha[1,1,2] = -2.5 
alpha[1,2,1] = alpha21
alpha[1,2,2] = 0

alpha[2,1,1] = 0 
alpha[2,1,2] = 0.2
alpha[2,2,1] = 0.2 
alpha[2,2,2] = 0

# parameters in the AR model for covariates
# the first covariate is the covariate in the RS model 
# the second covariate is the covariate in the multilevel AR-X model
ARX = c(0.9, 0.6) # AR parameters
sd_noiseX = c(0.5, 0.5) # standard deviations of process noises


# simulate data
odds = array(rep(NA,N*O*2),dim = c(N,O,2)) # odds of regime transitions
midx = matrix(NA,N,O) # regime indicator 
intercept=rnorm(N,int0,int_sd) # random intercepts
eta=matrix(NA,N,O) # latent log mean of Y
Y=matrix(NA,N,O) # observed variable, Y~Poisson(exp(eta))
X=array(rep(NA,N*O*2),dim = c(N,O,2)) # covariates

for(pp in 1:N){
  # the first time point
  odds[pp,1,1] = 1
  odds[pp,1,2] = exp(pi0)
  midx[pp,1] = sample.int(n=2,size=1,prob=odds[pp,1,])
  eta[pp,1] = rnorm(1, 0, 0.5)
  Y[pp,1] = (midx[pp,1]-1)*rpois(1, exp(eta[pp,1]))
  X[pp,1,] = rnorm(2, 0, 0.5)
  for(oo in 2:O){
    # from regime midx[pp,oo-1] to regime 1 (i.e., ZI regime)
    odds[pp,oo,1] = exp(alpha[1,1,midx[pp,oo-1]]+alpha[2,1,midx[pp,oo-1]]*X[pp,oo-1,1]) 
    # from regime midx[pp,oo-1] to regime 2 (i.e., drinking regime)
    odds[pp,oo,2] = exp(alpha[1,2,midx[pp,oo-1]]+alpha[2,2,midx[pp,oo-1]]*X[pp,oo-1,1])
    # regime indicator (1: ZI regime; 2: drinking regime)
    midx[pp,oo]=sample.int(n=2,size=1,prob=odds[pp,oo,])
    # multilevel AR-X model
    eta[pp,oo] = intercept[pp] + AR * (eta[pp,oo-1]- intercept[pp])  + coeffX * X[pp,oo-1,2] + rnorm(1,0,sd_noise)
    # ZIP model
    Y[pp,oo] = (midx[pp,oo]-1)*rpois(1, exp(eta[pp,oo])) 
    # model for covariates
    X[pp,oo,1] = rnorm(1, ARX[1]*X[pp,oo-1,1], sd_noiseX[1])
    X[pp,oo,2] = rnorm(1, ARX[2]*X[pp,oo-1,2], sd_noiseX[2])
  } # close loop over time points
} # close loop over persons

# generate missing data
missing_gen = function(lambda10,lambda11,lambda12,lambda20,lambda21,lambda22){
  n=N
  nt=O
  
  logit1=matrix(rep(NA,(nt*n)),nrow=n)
  pr1=matrix(rep(NA,(nt*n)),nrow=n)
  ry=matrix(rep(0,(nt*n)),nrow=n)
  
  logit2=matrix(rep(NA,(nt*n)),nrow=n)
  pr2=matrix(rep(NA,(nt*n)),nrow=n)
  rx1=matrix(rep(0,(nt*n)),nrow=n)
  rx2=matrix(rep(0,(nt*n)),nrow=n)
  
  for (i in 1:n){
    for (t in 2:(nt-1)){ # not generate missingness in the first and last observations
      logit1[i,t]=lambda10+lambda11*x1[i,t]+lambda12*x2[i,t]
      pr1[i,t]=exp(logit1[i,t])/(1+exp(logit1[i,t]))
      ry[i,t]=rbinom(1,1,pr1[i,t])
      
      logit2[i,t]=lambda20+lambda21*x3[i,t]+lambda22*x4[i,t]
      pr2[i,t]=exp(logit2[i,t])/(1+exp(logit2[i,t]))
      rx1[i,t]=rbinom(1,1,pr2[i,t])
      rx2[i,t]=rbinom(1,1,pr2[i,t])
    }#end of t loop
  }#end of i loop
  return(list(ry=ry,rx1=rx1,rx2=rx2))
}

x1=matrix(runif(N*O,-3,3),nrow=N,ncol=O)
x2=matrix(runif(N*O,-3,3),nrow=N,ncol=O)
x3=matrix(runif(N*O,-3,3),nrow=N,ncol=O)
x4=matrix(runif(N*O,-3,3),nrow=N,ncol=O)

miss=missing_gen(lambda10=-1.1,lambda11=.6,lambda12=.6,
                 lambda20=-1.1,lambda21=.6,lambda22=.6)

print(sum(miss[[1]])/(N*O))
print(sum(miss[[2]])/(N*O))
print(sum(miss[[3]])/(N*O))

Y_obs=Y
Y_obs[miss[[1]]==1]=NA

X1_obs=X[,,1]
X1_obs[miss[[2]]==1]=NA

X2_obs=X[,,2]
X2_obs[miss[[3]]==1]=NA

X_obs = array(rep(NA,N*O*2),dim = c(N,O,2))
X_obs[,,1] = X1_obs
X_obs[,,2] = X2_obs


# save data
save(Y_obs, X_obs, N, O, 
     int0, int_sd, AR, coeffX, sd_noise, 
     alpha, pi0,
     ARX, sd_noiseX,
     midx, file=paste0("SimulatedData_RSZIMLP_",condition,"_N",N,"T",O,"_",r1,".Rdata"))



