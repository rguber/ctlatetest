###########################################
# Application to US census data from "Instrument Validity Tests with Causal Trees: With an Application to the Same-sex Instrument", 
# Raphael Guber (2018), September 25th, 2018
##########################################

library(causalTree)  
library(treeClust)  
library(foreign)
library(readstata13)

set.seed(2)

samesexiv<-read.dta("~/ae98.dta") # this imports the Stata file ae98.dta

samesexiv$hispm <- as.factor(samesexiv$hispm)
samesexiv$blackm <- as.factor(samesexiv$blackm)
samesexiv$othracem <- as.factor(samesexiv$othracem)
samesexiv$sexk <- as.factor(samesexiv$sexk)

# take random 50% sample for speed:
set.seed(04052018) 
samesexiv <- samesexiv[sample(1:nrow(samesexiv),size=0.50* nrow(samesexiv)),]

samesexiv$Y <- samesexiv$weeksm
samesexiv$D <- samesexiv$morekids
samesexiv$Z <- samesexiv$samesex

####### run the functions.R file before executing the code below #####

#### Kitagawa's (2015) test:
boot.pval(y=samesexiv$weeksm,d=samesexiv$morekids,z=samesexiv$samesex,xi=1,subsets=53,B=2000)

##### Huber and Mellace's (2015, ReStat) probability test:
gridval=sort(unique(samesexiv$weeksm))
ivtestprob(y=samesexiv$weeksm,d=samesexiv$morekids,z=samesexiv$samesex,gridval=gridval, n_boot1=2000, n_boot2=2000)

# My test

test <- ctlatetest(data=samesexiv,minsize=2000,form="sexk + sex2nd + agem +schoolk +blackm+ hispm+ othracem+ agefstm ",B=2000)



