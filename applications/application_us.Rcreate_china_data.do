clear all

* take raw IPUMS data. The IPUMS data have to be downloaded by yourself at https://international.ipums.org/international/
use "~\ipums_china_2000.dta"
egen famid=group(sample serial)
keep if famunit==1 // keep only households with one family

compress
save rawdata_china.dta, replace


* create sample of children

keep if related>=3000 & related<4000 // keep all kinds of children
drop if age>=18

gen nage=-age
bysort famid (nage): gen childid=_n
bysort famid: gen nchildren=_N
drop nage

keep famid childid age sex   nchildren

xtset famid childid

reshape wide age sex  , i(famid) j(childid)

compress
save child_sample.dta, replace

* create sample of mothers
clear all

use rawdata_china.dta

keep if relate==1 | relate==2 // keep only spouses or household heads

keep if age>=21 & age<=35  // keep only persons between age 21 and 35
keep if sex==2 // keep only women
* keep if marst==2 // keep only married women
bysort famid: gen nmom=_n
keep if nmom==1
drop nmom

compress
save mother_sample.dta, replace


* merge children to mothers
merge 1:1 famid using child_sample
keep if _m==3
drop _m

* keep only if at least 2 children
keep if nchildren>1

* age at first birth
gen agefirstbirth=age-age1
drop if agefirstbirth<15

* keep only mothers whose oldest child is 18
keep if age1<18

* outcome
recode dayswrk (9 0=0) 
							
gen emp=dayswrk>0

recode edattain (4=3)

* treatment
gen morekids=nchildren>2

* samesex instrument
gen samesex=sex1==sex2


* lump ethnicities below .5% of sample together
bysort ethniccn: gen obs=_N
qui count
replace ethniccn=99 if obs<r(N)*0.005

keep age famid ethniccn lit edattain sex1 sex2 agefirstbirth morekids ///
 samesex emp dayswrk


compress
saveold china_fulldata.dta, replace version(13)
erase mother_sample.dta
erase child_sample.dta
