# Simulations from  "Instrument Validity Tests with Causal Trees: With an Application to the Same-sex Instrument", Raphael Guber (2018), September 25th, 2018
library(causalTree)  # see Susan Athey's github page on how to install this package: https://github.com/susanathey/causalTree
library(treeClust)  # for getting the leaf of estimation sample
library(foreign)
library(lattice)
library(fBasics)
library(mvtnorm)
library(AER)

set.seed(2)

# run the simulations.R file before running the simulations
#### SIMULATIONS

N_sim <- 1000
# N_sim <- 3000

VCOV <- matrix(c(5,0.5,0.5,5),2,2) # covariance matrix

repetitions<-10000

reject <- reject_hm <- pval_hm <- c()


set.seed(2)
while(length(reject)<repetitions){
  
  #set.seed(13032018)
  
  errors <-(rmvnorm(N_sim,rep(0,2),VCOV))
  
  # IV:
  Z<-as.numeric(rnorm(N_sim)>0)
  
  # covariates:
  u <- 3 # number of covariates
  X <- matrix(rnorm(N_sim*u), ncol=u)
  colnames(X) <- paste("Xvar", 1:u, sep="")
  betaXY = c(runif(u, 0, 1)) # rep(0,3) 
  
  # violation of exclusion restriction:
  gamma <- 0 # DGP 1: no violation
  # gamma <- 1 # DGP 2: global violation
  # gamma <- ifelse(X[,1]<0 & X[,2]>0 & X[,3]>0,5,0) # DGP 3: local violation
  
  # first stage coefficient/violation of monotonicity:
  alpha <- 0.2
  # alpha <- 0.6
  # alpha <- ifelse(X[,1]<0 & X[,2]>0,-10,3) # DGP4
  
  # treatment
  D<-as.numeric((alpha*Z+errors[,1])>0)
  Ylatent<-as.vector(2*D+gamma*Z+X%*%betaXY+errors[,2]) # covariates predict Y
  Y <- as.numeric(cut(Ylatent, breaks=c(-Inf,quantile(Ylatent,seq(0.25,0.75,0.25)),Inf))) # cut Ylatent into 4 equally spaced values
  
  data <- as.data.frame(cbind(Y,D,Z,X))
  
  
  # Huber and Mellace test. Use partially recentered p-value
  gridval_hm=sort(unique(Y))
  test_hm <-ivtestprob(y=Y,d=D,z=Z,gridval=gridval_hm, n_boot1=999, n_boot2=999)
  reject_hm <- c(reject_hm, test_hm$partialr.pval<=.05)
  pval_hm <- c(pval_hm, test_hm$partialr.pval)
  
  # My test
  test <- ctlatetest(data=data, minsize=50,share=1/2,form="Xvar1+Xvar2+Xvar3")
  reject <- c(reject,test$`test statistic`>test$`critical value`[[1]])
  
}


mean(reject) # rejection rate of my test
mean(reject_hm) # rejection rate of Huber and Mellace test


