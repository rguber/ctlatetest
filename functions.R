# Programs for "Instrument Validity Tests with Causal Trees: With an Application to the Same-sex Instrument", Raphael Guber (2018), September 25th, 2018
library(causalTree)  # see Susan Athey's github page on how to install this package: https://github.com/susanathey/causalTree
library(treeClust)  # for getting the leaf of estimation sample
library(foreign)
library(lattice)
library(fBasics)
library(mvtnorm)
library(AER)

set.seed(2)

#### PROGRAMS

# lines 24 to 401 contain code for the implementation of the IV validity tests by Toru Kitagawa's "A test for instrument validity."
# Econometrica 83 (5), 2043 - 2063) and Huber and Mellace's "Testing instrument validity for late identication based on inequality
# moment constraints." Review of Economics and Statistics 97 (2), 398 - 411.

# lines 405 to 415 contain helper functions for the ctlatetest function

# lines 426 to 721 contain the code for the ctlatetest function



#### Kitagawa's (2015) test (I thank Toru Kitagawa for providing his code):
{
  generate.grids.coarse<-function(y,d,z,subsets=5){
    
    
    #dataz1<-cbind("y"=y[z==1],"d"=d[z==1],"z"=z[z==1])
    #dataz0<-cbind("y"=y[z==0],"d"=d[z==0],"z"=z[z==0])
    
    #data<-as.matrix(rbind(dataz1,dataz0))
    
    #Yd1z0<-sort(subset(data[,1], (data[,2]==1)&(data[,3]==0)))
    #Yd0z1<-sort(subset(data[,1], (data[,2]==0)&(data[,3]==1)))
    #Y1grid<-unique(Yd1z0)
    #Y0grid<-unique(Yd0z1)
    #L1<-length(Y1grid)
    #L0<-length(Y0grid)
    #Y1grid<-Y1grid[floor(seq.int(from=1, to=L1,length.out=subsets+1))]
    #Y0grid<-Y0grid[floor(seq.int(from=1, to=L0,length.out=subsets+1))]
    
    Y1grid <- Y0grid <- seq(from=min(y),to=max(y),length.out = subsets+1) # from Huber Mellace, much easier
    
    
    return(list("Y1grid"=Y1grid, "Y0grid"=Y0grid))
  }
  
  # This function computes the empirical probabilities over the interval classes constructed by the function "generate.grids.coarse".
  
  compute.PQ.grids<-function(y,d,z,subsets=5){
    
    dataz1<-cbind("y"=y[z==1],"d"=d[z==1],"z"=z[z==1])
    dataz0<-cbind("y"=y[z==0],"d"=d[z==0],"z"=z[z==0])
    
    data<-as.matrix(rbind(dataz1,dataz0))
    
    m <- length(dataz1[,1])
    n <- length(dataz0[,1])
    N <- m+n
    
    grids<-generate.grids.coarse(y,d,z,subsets)
    
    p1<-sum(data[1:m,2])/m
    p0<-sum(data[(m+1):(m+n),2])/n
    
    Yd1z1<-sort(subset(data[1:m,1], data[1:m,2]==1))
    Yd1z0<-sort(subset(data[(m+1):(m+n),1], data[(m+1):(m+n),2]==1))
    Yd0z1<-sort(subset(data[1:m,1], data[1:m,2]==0))
    Yd0z0<-sort(subset(data[(m+1):(m+n),1], data[(m+1):(m+n),2]==0))
    
    Qcdf.d1<-findInterval(grids$Y1grid,Yd1z0)/n
    Pcdf.d1<-findInterval(grids$Y1grid,Yd1z1)/m
    Qcdf.d0<-findInterval(grids$Y0grid,Yd0z0)/n
    Pcdf.d0<-findInterval(grids$Y0grid,Yd0z1)/m
    
    g1length<-length(grids$Y1grid)
    g0length<-length(grids$Y0grid)
    
    QV.d1<-triang(Qcdf.d1 - matrix(kronecker(rep(1,g1length), c(0,Qcdf.d1[1:(g1length-1)])), nrow=g1length, byrow=T))
    PV.d1<-triang(Pcdf.d1 - matrix(kronecker(rep(1,g1length), c(0,Pcdf.d1[1:(g1length-1)])), nrow=g1length, byrow=T))
    QV.d0<-triang(Qcdf.d0 - matrix(kronecker(rep(1,g0length), c(0,Qcdf.d0[1:(g0length-1)])), nrow=g0length, byrow=T))
    PV.d0<-triang(Pcdf.d0 - matrix(kronecker(rep(1,g0length), c(0,Pcdf.d0[1:(g0length-1)])), nrow=g0length, byrow=T))
    return(list("Q1"=QV.d1,"P1"=PV.d1, "Q0"=QV.d0, "P0"=PV.d0))
  }
  
  compute.stat.grids.H<-function(y,d,z,xi=0.1,subsets=5){
    
    dataz1<-cbind("y"=y[z==1],"d"=d[z==1],"z"=z[z==1])
    dataz0<-cbind("y"=y[z==0],"d"=d[z==0],"z"=z[z==0])
    
    data<-as.matrix(rbind(dataz1,dataz0))
    
    m <- length(dataz1[,1])
    n <- length(dataz0[,1])
    N <- m+n
    lambda<-m/N
    
    # grids<-generate.grids.coarse(y,d,z)
    PQ<-compute.PQ.grids(y,d,z,subsets=subsets)
    weight.d1<-pmax(xi,sqrt(lambda*PQ$Q1*(1-PQ$Q1)+(1-lambda)*PQ$P1*(1-PQ$P1)))
    weight.d0<-pmax(xi,sqrt(lambda*PQ$Q0*(1-PQ$Q0)+(1-lambda)*PQ$P0*(1-PQ$P0)))
    
    Q_P.d1<-c(PQ$Q1-PQ$P1)
    P_Q.d0<-c(PQ$P0-PQ$Q0)
    
    weightedQ_P.d1<-Q_P.d1/weight.d1
    weightedP_Q.d0<-P_Q.d0/weight.d0
    
    Tw<-(as.numeric(n)*as.numeric(m)/(n+m))^(0.5)*max(weightedQ_P.d1,weightedP_Q.d0)
    return(Tw)
  }
  
  boot.pval<-function(y,d,z,B=500,xi=0.1,subsets=5){
    
    
    
    dataz1<-cbind("y"=y[z==1],"d"=d[z==1],"z"=z[z==1])
    dataz0<-cbind("y"=y[z==0],"d"=d[z==0],"z"=z[z==0])
    
    data<-as.matrix(rbind(dataz1,dataz0))
    
    m <- length(dataz1[,1])
    n <- length(dataz0[,1])
    N <- m+n
    lambda<-m/N
    
    teststat <- compute.stat.grids.H(y,d,z,subsets=subsets,xi=xi)
    
    #grids<-generate.grids.coarse(y,d,z,subsets=subsets)
    #PQ<-compute.PQ.grids(y,d,z,subsets=subsets)
    
    bootstats <- rep(NA,B)
    
    for(b in 1:B) {
      index <- sample(c(1:(m+n)), size=(m+n), replace=TRUE)
      yboot <- y[index]
      dboot <- d[index]
      zboot <- z[index]
      bootstats[b] <- compute.stat.grids.H(yboot,dboot,zboot,subsets=subsets,xi=xi)
    }                                                                               
    pval<-sum(teststat <= bootstats)/B
    return(pval)
  }
}
##### Huber and Mellace's (2015) probability test based on Bennett (2009) and Chen and Schroeter (2012):
# (I thank Martin Huber for providing his code)

# IV probability test based on Bennett (2009) and Chen and Schroeter (2012)
ivtestprob<-function(y,d,z,n_boot1=999, n_boot2=999, gridval=NULL, subsets=5){
  n=length(y)
  p1_1<-mean(d[z==1])
  p1_0<-mean(d[z==0])
  p0_1<-1-mean(d[z==1])
  p0_0<-1-mean(d[z==0])
  z1<-z[d==1]
  z0<-z[d==0]
  y1<-y[d==1]
  y0<-y[d==0]
  y11=y1[z1==1]
  y01=y1[z1==0]
  y10=y0[z0==1]
  y00=y0[z0==0]
  
  if ((is.null(gridval)*1)==1){
    gridval=seq(from=min(y),to=max(y),length.out = subsets+1)
    grid=subsets
  }
  if ((is.null(gridval)*1)==0){
    grid=length(gridval)-1
  }
  p11<-c()
  p01<-c()
  p10<-c()
  p00<-c()
  
  gridval[length(gridval)]=gridval[length(gridval)]+0.00001
  for (i in 1:(grid)){
    p11<-c(p11,mean(y11>=gridval[i] & y11<gridval[i+1]))
    p01<-c(p01,mean(y01>=gridval[i] & y01<gridval[i+1]))
    p10<-c(p10,mean(y10>=gridval[i] & y10<gridval[i+1]))
    p00<-c(p00,mean(y00>=gridval[i] & y00<gridval[i+1]))
  }
  lengthp=length(p11)
  
  q<-p1_0/p1_1
  r<-p0_1/p0_0
  
  p11min<-(p11-(1-q))/q
  p11max<-p11/q
  p00min<-(p00-(1-r))/r
  p00max<-p00/r
  
  
  #first bootstrap 
  Tn1b = c()
  Tn2b = c()
  Tn3b = c()
  Tn4b = c()
  
  while(length(Tn1b)<n_boot1*lengthp){
    sboot<-sample(1:length(y),length(y),TRUE)
    zb<-z[sboot]
    db<-d[sboot]
    yb<-y[sboot]
    
    obsb=length(yb)
    obs1b<-sum(db)
    obs0b<-sum(1-db)
    
    p1_1b<-mean(db[zb==1])
    p1_0b<-mean(db[zb==0])
    p0_1b<-1-mean(db[zb==1])
    p0_0b<-1-mean(db[zb==0])
    
    if ((p1_0b/p1_1b<=1) & (p0_1b/p0_0b<=1)){
      z1b<-zb[db==1]
      z0b<-zb[db==0]
      y1b<-yb[db==1]
      y0b<-yb[db==0]
      y11b=y1b[z1b==1]
      y01b=y1b[z1b==0]
      y10b=y0b[z0b==1]
      y00b=y0b[z0b==0]
      
      p11b<-c()
      p01b<-c()
      p10b<-c()
      p00b<-c()
      for (i in 1:(grid)){
        p11b<-c(p11b,mean(y11b>=gridval[i] & y11b<gridval[i+1]))
        p01b<-c(p01b,mean(y01b>=gridval[i] & y01b<gridval[i+1]))
        p10b<-c(p10b,mean(y10b>=gridval[i] & y10b<gridval[i+1]))
        p00b<-c(p00b,mean(y00b>=gridval[i] & y00b<gridval[i+1]))
      }
      
      qb<-p1_0b/p1_1b
      rb<-p0_1b/p0_0b
      p11minb<-(p11b-(1-qb))/qb
      p11maxb<-p11b/qb
      p00minb<-(p00b-(1-rb))/rb
      p00maxb<-p00b/rb
      
      Tn1b = cbind(Tn1b,(p11minb-p01b)) 
      Tn2b = cbind(Tn2b,(p01b-p11maxb)) 
      Tn3b = cbind(Tn3b,(p00minb-p10b)) 
      Tn4b = cbind(Tn4b,(p10b-p00maxb)) 
    }
  }
  
  # variances
  varTn1b=diag(var(t(Tn1b)))
  varTn2b=diag(var(t(Tn2b)))
  varTn3b=diag(var(t(Tn3b)))
  varTn4b=diag(var(t(Tn4b)))
  
  # create objects for partial and full recentering
  FCbsample1 <-matrix(0,nrow=lengthp,ncol=n_boot2)
  FCbsample2 <-matrix(0,nrow=lengthp,ncol=n_boot2)
  FCbsample3 <-matrix(0,nrow=lengthp,ncol=n_boot2)
  FCbsample4 <-matrix(0,nrow=lengthp,ncol=n_boot2)
  
  PCbsample1 <-matrix(0,nrow=lengthp,ncol=n_boot2)
  PCbsample2 <-matrix(0,nrow=lengthp,ncol=n_boot2)
  PCbsample3 <-matrix(0,nrow=lengthp,ncol=n_boot2)
  PCbsample4 <-matrix(0,nrow=lengthp,ncol=n_boot2)
  
  FCp_vals1 <- matrix(0, nrow = n_boot2, ncol = lengthp)
  FCp_vals2 <- matrix(0, nrow = n_boot2, ncol = lengthp)
  FCp_vals3 <- matrix(0, nrow = n_boot2, ncol = lengthp)
  FCp_vals4 <- matrix(0, nrow = n_boot2, ncol = lengthp)
  
  PCp_vals1 <- matrix(0, nrow = n_boot2, ncol = lengthp)
  PCp_vals2 <- matrix(0, nrow = n_boot2, ncol = lengthp)
  PCp_vals3 <- matrix(0, nrow = n_boot2, ncol = lengthp)
  PCp_vals4 <- matrix(0, nrow = n_boot2, ncol = lengthp)
  
  #Fix value of sequence delta_n
  delta_n = sqrt(2*log(log(n))/n)
  
  # compute test statistics
  Tn1 = (p11min-p01) 
  Tn2 = (p01-p11max)
  Tn3 = (p00min-p10) 
  Tn4 = (p10-p00max)  
  Tn=c(Tn1, Tn2, Tn3, Tn4) 
  
  # compute recentering terms
  maxVec1=( Tn1 >= (-delta_n*sqrt(varTn1b))  )*Tn1    + ( Tn1 <  (-delta_n*sqrt(varTn1b))  )* (-delta_n*sqrt(varTn1b) )
  maxVec2=( Tn2 >= (-delta_n*sqrt(varTn2b))  )*Tn2   +( Tn2 < (-delta_n*sqrt(varTn2b))  )* (-delta_n*sqrt(varTn2b) )
  maxVec3=( Tn3 >= (-delta_n*sqrt(varTn3b))  )*Tn3    + ( Tn3 <  (-delta_n*sqrt(varTn3b))  )* (-delta_n*sqrt(varTn3b) )
  maxVec4=( Tn4 >= (-delta_n*sqrt(varTn4b))  )*Tn4   +( Tn4 < (-delta_n*sqrt(varTn4b))  )* (-delta_n*sqrt(varTn4b) )
  
  Tn1 = sqrt(n)*Tn1 
  Tn2 = sqrt(n)*Tn2
  Tn3 = sqrt(n)*Tn3 
  Tn4 = sqrt(n)*Tn4
  
  Tn1b=sqrt(n)*Tn1b
  Tn2b=sqrt(n)*Tn2b
  Tn3b=sqrt(n)*Tn3b
  Tn4b=sqrt(n)*Tn4b
  
  # Generate bootstrap samples (using first bootstrap)- full recentering
  FCBootstrap1<-matrix(Tn1b-rep(Tn1,n_boot1),lengthp,n_boot1)
  FCBootstrap2<-matrix(Tn2b-rep(Tn2,n_boot1),lengthp,n_boot1)
  FCBootstrap3<-matrix(Tn3b-rep(Tn3,n_boot1),lengthp,n_boot1)
  FCBootstrap4<-matrix(Tn4b-rep(Tn4,n_boot1),lengthp,n_boot1)
  
  # Generate bootstrap samples (using first bootstrap)- partial recentering
  PCBootstrap1<-matrix(Tn1b-rep((sqrt(n)* maxVec1),n_boot1),lengthp,n_boot1)
  PCBootstrap2<-matrix(Tn2b-rep((sqrt(n)* maxVec2),n_boot1),lengthp,n_boot1)
  PCBootstrap3<-matrix(Tn3b-rep((sqrt(n)* maxVec3),n_boot1),lengthp,n_boot1)
  PCBootstrap4<-matrix(Tn4b-rep((sqrt(n)* maxVec4),n_boot1),lengthp,n_boot1)
  
  # Generate empirical distribution functions - full recentering
  FC_ECDF1 <- apply(FCBootstrap1, 1, ecdf)
  FC_ECDF2 <- apply(FCBootstrap2, 1, ecdf)
  FC_ECDF3 <- apply(FCBootstrap3, 1, ecdf)
  FC_ECDF4 <- apply(FCBootstrap4, 1, ecdf)
  
  # Generate empirical distribution functions - partial recentering
  PC_ECDF1 <- apply(PCBootstrap1, 1, ecdf)
  PC_ECDF2 <- apply(PCBootstrap2, 1, ecdf)
  PC_ECDF3 <- apply(PCBootstrap3, 1, ecdf)
  PC_ECDF4 <- apply(PCBootstrap4, 1, ecdf)
  
  # second bootstrap - draw observations from first bootstrap
  bsample2 <- sample(n_boot1, n_boot2, replace=TRUE)
  
  # draw samples of recentered test statistics out of first bootstrap - full recentering
  FCbsample1 = FCBootstrap1[1:lengthp,bsample2]
  FCbsample2 = FCBootstrap2[1:lengthp,bsample2]
  FCbsample3 = FCBootstrap3[1:lengthp,bsample2]
  FCbsample4 = FCBootstrap4[1:lengthp,bsample2]
  
  # draw samples of recentered test statistics out of first bootstrap - partial recentering
  PCbsample1 = PCBootstrap1[1:lengthp,bsample2]
  PCbsample2 = PCBootstrap2[1:lengthp,bsample2]
  PCbsample3 = PCBootstrap3[1:lengthp,bsample2]
  PCbsample4 = PCBootstrap4[1:lengthp,bsample2]
  
  # empirical distributions of fully recentered p-values
  
  p_val_1 =1
  p_val_2 =1
  p_val_3 =1
  p_val_4 =1
  
  for (k in 1:lengthp) {
    FCp_vals1[1:n_boot2,k] = 1-FC_ECDF1[[k]](FCbsample1[k,1:n_boot2])	
    FCp_vals2[1:n_boot2,k] = 1-FC_ECDF2[[k]](FCbsample2[k,1:n_boot2])
    FCp_vals3[1:n_boot2,k] = 1-FC_ECDF3[[k]](FCbsample3[k,1:n_boot2])	
    FCp_vals4[1:n_boot2,k] = 1-FC_ECDF4[[k]](FCbsample4[k,1:n_boot2])
    FCp_vals=cbind(FCp_vals1,FCp_vals2,FCp_vals3,FCp_vals4)
    
    # empirical distributions of partially recentered p-values
    PCp_vals1[1:n_boot2,k] = 1-FC_ECDF1[[k]](PCbsample1[k,1:n_boot2])
    PCp_vals2[1:n_boot2,k] = 1-FC_ECDF2[[k]](PCbsample2[k,1:n_boot2])
    PCp_vals3[1:n_boot2,k] = 1-FC_ECDF3[[k]](PCbsample3[k,1:n_boot2])
    PCp_vals4[1:n_boot2,k] = 1-FC_ECDF4[[k]](PCbsample4[k,1:n_boot2])
    PCp_vals=cbind(PCp_vals1,PCp_vals2,PCp_vals3,PCp_vals4)
    
    # compute the minimum p-value based on the original sample
    p_val_1 =  min( 1- FC_ECDF1[[k]](Tn1[k]), p_val_1 )
    p_val_2 =  min( 1- FC_ECDF2[[k]](Tn2[k]), p_val_2 )
    p_val_3 =  min( 1- FC_ECDF3[[k]](Tn3[k]), p_val_3 )
    p_val_4 =  min( 1- FC_ECDF4[[k]](Tn4[k]), p_val_4 )
    p_val =  min(p_val_1,p_val_2,p_val_3,p_val_4)
  }  # end of loop of second bootstrap
  
  # compute adjusted p-values - full recentering 
  min_p_vec <- apply(FCp_vals,1,min) # minimum p-value vector
  edfF_n <- ecdf(min_p_vec)    # empirical distribution function
  FCpval = edfF_n(p_val)  
  
  # compute adjusted p-values - partial recentering 
  min_p_vec <- apply(PCp_vals,1,min)
  edfF_n <- ecdf(min_p_vec)
  PCpval = edfF_n(p_val) 
  
  # CHEN AND SZROETER TEST
  J=c(n*varTn1b,n*varTn2b,n*varTn3b, n*varTn4b)
  J[J==0]=0.000000000001
  theta=1/sqrt(J)
  #Tuning parameter for the sample size
  K=1/delta_n
  # Smooth the indicator function
  psi=pnorm(K*theta*Tn);
  Lambda=(dnorm(K*theta*Tn)*K)/sqrt(n)
  # Create a diagonal matrix with the variances of the elements of Tn
  delta=diag(theta);
  # Compute the P-value
  Q1=sqrt(n)*t(psi)%*%(delta%*%Tn)-matrix(1,1,length(Tn))%*%(Lambda);
  Q2=sqrt(t(psi)%*%delta%*%t(t(J)%*%delta*psi))
  Pc=min(1-pnorm(Q1/max(Q2,0.00000001)),1)
  
  list(ChSz.pval=Pc, fullr.pval=FCpval, partialr.pval=PCpval)
  
}


#### Causal Tree IV validity test

# helper functions
row2leaf <- function(row,rp) # converts row number from tree$frame into leaf number
{
  as.numeric(rownames(rp$frame)[row])
}

parent <- function(x) { # gives the number of the parent leaf
  if (x[1] != 1)
    c(Recall(if (x %% 2 == 0L) x / 2 else (x - 1) / 2), x) else x
}


# The Causal Tree IV validity test procedure. It is essentially a wrapper for Athey and Imbens (2016) causalTree function/package
# minsize is the minimum number of obvservations with either D=1 or D=0 in a leaf
# share denotes the share of observations used for sample A (remainder is used for sample B). Sample A (B) is used for building the
# tree while sample B is used for estimating the inequalities. The roles are then swapped
# form is a string of type "x1+x2+x3" where x1, x2, x3 are covariates to search for violations
# beta_N is a parameter of the Chernozhukov, Chetverikov and Kato (2016) (CCK) many moment inequalities test
# alpha is the significance level
# B is the number of bootstrap replications for the CCK test

ctlatetest <- function(data, minsize=1000, share=1/2,form, beta_N=0.0001, alpha=0.05, B=2000) {
  
  Y <- data$Y
  Z <- data$Z
  D <- data$D
  
  N  <- nrow(data)
  
  P <- Q <-  matrix(NA,nrow=N, ncol=length(unique(Y)))
  
  for (x in 1:length(unique(Y))) {
    P[,x] <- Q[,x] <- ifelse(Y==sort(unique(Y))[x],1,0)
  }
  
  # pseudo outcomes:
  
  P <- P*D
  Q <- -Q*(1-D)
  
  # store zetas
  zeta1A <-zeta1B <- zeta0A <- zeta0B <- c()
  
  # split sample:
  
  trID <- which(data$Z==1)
  conID <- which(data$Z==0)
  
  
  index <- c(sample(trID,length(trID)*share),sample(conID,length(conID)*share))
  
  
  sampleA <- data[index,]
  sampleB <- data[-index,]
  
  
  Nb <- nrow(sampleB)
  Na <- nrow(sampleA)
  
  ## built trees and calculate zetas
  # first for P (d=1):
  for (x in 1:ncol(P)) {
    sampleA$P <- P[index,x]
    sampleB$P <- P[-index,x]
    
    
    # first build tree:
    formula <- paste("P~",  form) 
    
    # no honest splitting here, we do that manually
    invisible(capture.output(tree <- causalTree(formula, data = sampleA, treatment = sampleA$Z,
                                                split.Rule = "CT", cv.option = "CT", split.Honest = F, cv.Honest = F, split.Bucket = F, 
                                                xval = 5, minsize = minsize)))
    
    # Pruning (ignored now)
    opcp <- tree$cptable[,1][which.min(tree$cptable[,4])] 
    finaltree <- prune(tree, opcp)
    
    tree <- finaltree
    
    # now save information about in which leaf observations from sample B fall:
    leaves <- row2leaf(rpart.predict.leaves(tree, sampleB, type = "where"),tree)
    allleaves <- lapply(leaves,parent)
    
    size <- nrow(tree$frame)
    
    
    leafPdummy <- zeta1Btemp<-  matrix(NA, nrow=Nb, ncol=size)
    
    
    for (z in 1:size) {
      l <- sort(as.numeric(rownames(tree$frame)))[z]
      leafPdummy[,z] <- ifelse(lapply(lapply(allleaves,is.element,l),max)==1,1,0)
      
      # compute P(Z=1 & leaf=x)
      pZ <- mean(sampleB$Z*leafPdummy[,z])
      pNZ <- mean((1-sampleB$Z)*leafPdummy[,z])
      zeta1Btemp[,z] <-leafPdummy[,z]*sampleB$P*((1-sampleB$Z)/pNZ-sampleB$Z/pZ)
    }
    
    zeta1B <- cbind(zeta1B, zeta1Btemp)
    
    rm(tree)
    # now switch roles of samples
    invisible(capture.output(tree <- causalTree(formula, data = sampleB, treatment = sampleB$Z,
                                                split.Rule = "CT", cv.option = "CT", split.Honest = F, cv.Honest = F, split.Bucket = F, 
                                                xval = 5, minsize = minsize)))
    
    # Pruning
    opcp <- tree$cptable[,1][which.min(tree$cptable[,4])] 
    finaltree <- prune(tree, opcp)
    
    tree <- finaltree
    
    # now save information about in which leaf observations from sample A fall:
    leaves <- row2leaf(rpart.predict.leaves(tree, sampleA, type = "where"),tree)
    allleaves <- lapply(leaves,parent)
    
    size <- nrow(tree$frame)
    
    
    leafPdummy <- zeta1Atemp<-  matrix(NA, nrow=Na, ncol=size)
    
    
    for (z in 1:size) {
      l <- sort(as.numeric(rownames(tree$frame)))[z]
      leafPdummy[,z] <- ifelse(lapply(lapply(allleaves,is.element,l),max)==1,1,0)
      
      # compute P(Z|leaf=x)
      pZ <- mean(sampleA$Z*leafPdummy[,z])
      pNZ <- mean((1-sampleA$Z)*leafPdummy[,z])
      zeta1Atemp[,z] <-leafPdummy[,z]*sampleA$P*((1-sampleA$Z)/pNZ-sampleA$Z/pZ)
    }
    zeta1A <- cbind(zeta1A, zeta1Atemp)
  }
  
  
  # now for Q (d=0):
  for (x in 1:ncol(Q)) {
    sampleA$Q <- Q[index,x]
    sampleB$Q <- Q[-index,x]
    
    
    # first build tree:
    formula <- paste("Q~",  form) 
    
    # no honest splitting here, we do that manually
    invisible(capture.output(tree <- causalTree(formula, data = sampleA, treatment = sampleA$Z,
                                                split.Rule = "CT", cv.option = "CT", split.Honest = F, cv.Honest = F, split.Bucket = F, 
                                                xval = 5, minsize = minsize)))
    
    # Pruning 
    opcp <- tree$cptable[,1][which.min(tree$cptable[,4])] 
    finaltree <- prune(tree, opcp)
    
    tree <- finaltree
    
    # now save information about in which leaf observations from sample B fall:
    leaves <- row2leaf(rpart.predict.leaves(tree, sampleB, type = "where"),tree)
    allleaves <- lapply(leaves,parent)
    
    size <- nrow(tree$frame)
    
    
    leafQdummy <- zeta0Btemp<-  matrix(NA, nrow=Nb, ncol=size)
    
    
    for (z in 1:size) {
      l <- sort(as.numeric(rownames(tree$frame)))[z]
      leafQdummy[,z] <- ifelse(lapply(lapply(allleaves,is.element,l),max)==1,1,0)
      
      # compute P(Z|leaf=x)
      pZ <- mean(sampleB$Z*leafQdummy[,z])
      pNZ <- mean((1-sampleB$Z)*leafQdummy[,z])
      
      zeta0Btemp[,z] <-leafQdummy[,z]*sampleB$Q*((1-sampleB$Z)/pNZ-sampleB$Z/pZ)
    }
    
    zeta0B <- cbind(zeta0B, zeta0Btemp)
    
    rm(tree)
    # now switch roles of samples
    invisible(capture.output(tree <- causalTree(formula, data = sampleB, treatment = sampleB$Z,
                                                split.Rule = "CT", cv.option = "CT", split.Honest = F, cv.Honest = F, split.Bucket = F, 
                                                xval = 5, minsize = minsize)))
    
    # Pruning 
    opcp <- tree$cptable[,1][which.min(tree$cptable[,4])] 
    finaltree <- prune(tree, opcp)
    
    tree <- finaltree
    
    # now save information about in which leaf observations from sample A fall:
    leaves <- row2leaf(rpart.predict.leaves(tree, sampleA, type = "where"),tree)
    allleaves <- lapply(leaves,parent)
    
    size <- nrow(tree$frame)
    
    
    
    leafQdummy <- zeta0Atemp<-  matrix(NA, nrow=Na, ncol=size)
    
    
    for (z in 1:size) {
      l <- sort(as.numeric(rownames(tree$frame)))[z]
      leafQdummy[,z] <- ifelse(lapply(lapply(allleaves,is.element,l),max)==1,1,0)
      
      # compute P(Z|leaf=x)
      pZ <- mean(sampleA$Z*leafQdummy[,z])
      pNZ <- mean((1-sampleA$Z)*leafQdummy[,z])
      
      zeta0Atemp[,z] <-leafQdummy[,z]*sampleA$Q*((1-sampleA$Z)/pNZ-sampleA$Z/pZ)
    }
    zeta0A <- cbind(zeta0A, zeta0Atemp)
  }
  
  
  zetaA<- cbind(zeta1A,zeta0A)
  zetaB<- cbind(zeta1B,zeta0B)
  
  p <- ncol(zetaA)+ncol(zetaB)
  
  zeta <- matrix(NA, nrow=N,ncol=p)
  zeta[1:Na,1:ncol(zetaA)] <- zetaA
  zeta[(Na+1):N,(ncol(zetaA)+1):p] <- zetaB
  
  #### CCK test procedure
  # critical value function
  
  N <- (Nb+Na)/2
  
  
  # compute c(beta_N)
  
  cbeta_N <- qnorm(1-beta_N/p)/sqrt(1-qnorm(1-beta_N/p)/N) 
  
  # generate pre-selected zetas
  zeta_sel_no <-  c()
  
  
  for (x in 1:ncol(zeta)) {
    
    # calculate standard dev
    zeta_sd <- sqrt(var(zeta[,x],na.rm=TRUE))
    
    # preliminary T stat
    T_prelim <- sqrt(N)*mean(zeta[,x],na.rm=T)/zeta_sd
    
    # Do not select a moment if its T stat is smaller than some value, or its std is zero
    if (T_prelim>-2*cbeta_N & zeta_sd!=0 & T_prelim!="NaN") {
      zeta_sel_no<-c(zeta_sel_no, x)
    }
  }
  
  # selected zetas
  zeta_sel <- zeta[,zeta_sel_no]
  dim(zeta_sel)
  
  # compute means and stds
  zeta_sel_m <- zeta_sel_sd <-c()
  for (x in 1:ncol(zeta_sel)) {
    
    # calculate mean
    zeta_sel_m <- c(zeta_sel_m,mean(zeta_sel[,x],na.rm=TRUE))
    
    # calculate standard dev
    zeta_sel_sd <- c(zeta_sel_sd,sd(zeta_sel[,x],na.rm=TRUE))
    
  }
  
  
  # test stat:
  T <- c()
  for (x in 1:ncol(zeta_sel)) {
    
    T <-c(T,sqrt(N)*zeta_sel_m[x]/zeta_sel_sd[x] )
    
  }
  
  T <- max(T)
  
  # critical value procedure:
  W <-c()
  
  # now perform bootstrap
  
  # number of reps:
  
  
  while(length(W)<B) {
    sboot<-c(sample(1:Na,Na,TRUE),sample((Na+1):(Na+Nb),Nb),TRUE)
    zeta_sel_boot <-zeta_sel[sboot,] 
    
    W_pre <-c()
    
    for (x in 1:ncol(zeta_sel_boot)) {
      
      W_pre <-c(W_pre,sqrt(N)*mean(zeta_sel_boot[,x]-zeta_sel_m[x], na.rm=TRUE)/zeta_sel_sd[x] )
      
    }
    
    W <-c(W,max(W_pre))
    
  }
  
  # critical value
  
  q <- 1-alpha+2*beta_N
  
  c_alpha <- quantile(W,q, ra.rm=TRUE)
  
  # p-value:
  pseudo_pval <- mean(as.numeric(W>T))+2*beta_N
  
  return(list("no inequalities"=p, "test statistic" =T, "critical value"=c_alpha, "Pseudo p-value"=pseudo_pval, "inequalities"=dim(zeta_sel)))
  
}

