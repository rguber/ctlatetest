clear
set more off
cap log close


use AngristEvans1980.dta // I use the original Angrist and Evans (1998) data of 
						 // the 1980 US census. These can be downloaded at 
						 // Joshua Angrist's data archive at 
						 // https://economics.mit.edu/faculty/angrist/data1/data/angev98

gen id=_n

* construct the sexes of 1st two kids;
destring sex2nd, replace
destring sex3rd, replace
destring ageq2nd, replace
destring ageq3rd, replace
destring agem, replace
destring aged, replace
destring qtrbthd, replace
destring asex2nd, replace
destring aage2nd, replace

//X and Z
gen boy1st=(sexk==0)
gen boy2nd=(sex2nd==0)
gen boys2=((sexk==0) & (sex2nd==0))
gen girls2=((sexk==1) & (sex2nd==1))
gen samesex=((boys2==1) | (girls2==1))
gen morekids=(kidcount>2)
tab samesex if kidcount>1
gen twin2nd=(ageq2nd==ageq3rd)
gen samesextwin2nd=(sex2nd==sex3rd & ageq2nd==ageq3rd)

//Y
gen workedm=(weeksm>0)
gen incomem=income1m+max(0,income2m)
** deflate wages **
replace incomem=incomem*2.099173554

//further things
gen blackm=(racem==2)
gen hispm=(racem==12)
gen whitem=(racem==1)
gen othracem=1-blackm-hispm-whitem

gen yobd=80-aged if qtrbthd==0
replace yobd=79-aged if qtrbthd!=0

gen ageqm=4*(80-yobm)-qtrbthm-1
gen agefstm=int((ageqm-ageqk)/4)
gen ageqd=4*(80-yobd)-qtrbthd
gen agefstd=int((ageqd-ageqk)/4)

replace qtrmar=qtrmar-1
gen yom=yobm+agemar if (qtrbthm<=qtrmar)
replace yom=yobm+agemar+1 if (qtrbthm>qtrmar)
gen dom_q = yom+(qtrmar/4)
gen do1b_q = yobk+((qtrbkid)/4)
gen illegit=0
replace illegit=1 if (dom_q-do1b_q)>0

gen lincomem=log(incomem)

//sample restrictions:
//keep if marital==0		//sample restriction for married sample
drop if agem==.
//drop if aged==.			//sample restriction for married sample
drop if agem<21
drop if agem>35
keep if kidcount>=2
//keep if timesmar==1		//sample restriction for married sample
//drop if illegit==1		//sample restriction for married sample
drop if ageq2nd<=4
drop if agefstm<15
//drop if agefstd<15		//sample restriction for married sample
drop if aage==1
drop if aqtrbrth==1
drop if aage2nd==1
drop if asex==1
drop if asex2nd==1

assert twin2nd!=.
assert morekids!=.
assert worked!=.



gen agefstm2=agefstm^2

label def sexk 1 "male" 0 "female"

label val  sexk sex2nd sexk


compress



keep samesex morekids agem agefstm schoolk sexk sex2nd blackm hispm othracem hoursm incomem ///
weeksm 

saveold "H:\Ideas\Sharpen bounds\application\samesex\ae98\ae98.dta", replace ///
version(12)


