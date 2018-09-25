###########################################
# Application to Chinese census data from "Instrument Validity Tests with Causal Trees: With an Application to the Same-sex Instrument", 
# Raphael Guber (2018), September 25th, 2018
##########################################

library(causalTree)  
library(treeClust)  
library(foreign)
library(readstata13)

set.seed(2)

samesexiv <-read.dta13("~/china_fulldata.dta") # this imports the Stata file china_fulldata.dta

# take random 50% sample for speed:
set.seed(04052018) 

samesexiv$dayswrk <- as.numeric(samesexiv$dayswrk)-1
samesexiv$age <- as.numeric(samesexiv$age) 

samesexiv$Y <- samesexiv$dayswrk
samesexiv$D <- samesexiv$morekids
samesexiv$Z <- samesexiv$samesex


samesexiv <- samesexiv[sample(1:nrow(samesexiv),size=0.50* nrow(samesexiv)),]

samesexiv$Y <- samesexiv$dayswrk
samesexiv$D <- samesexiv$morekids
samesexiv$Z <- samesexiv$samesex

####### run the functions.R file before executing the code below #####

# Kitagawas (2015) test:
boot.pval(y=samesexiv$dayswrk,d=samesexiv$morekids,z=samesexiv$samesex,xi=1,subsets=7,B=2000)

##### Huber and Mellace's (2015, ReStat) probability test:

gridval=sort(unique(samesexiv$dayswrk))
ivtestprob(y=samesexiv$dayswrk,d=samesexiv$morekids,z=samesexiv$samesex,gridval=gridval, n_boot1=2000, n_boot2=2000)

# My test
test <- ctlatetest(data=samesexiv,minsize=2000,form="age+ ethniccn +lit+  sex1+ sex2+ edattain + agefirstbirth",B=2000)



