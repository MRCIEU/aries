varnames <- c(
    ## misc
    "c755", ## maternal occupation
    "b032", ## parity
    "c645a", ## maternal education
    "bestgest", ## gestational age
    "c804", ## ethnicity
    ## child birthweight
    "kz030", ## birthweight
    ## mother bmi
    "dw042", ## pregnancy weight and BMI at 12 weeks gest
    "fm1ms111",## 18 years postnatal
    "fm2ms112",## 20 years postnatal
    ## father bmi
    "ff1ms111",##20 years postnatal
    ## children bmi
    "f7ms026a",##7.5 years
    "fdms026a",##10.5 years
    "fjmr022a",##17.5 years
    "fkms1040",##24 years
    ## mother prenatal alcohol
    "a260",
    "b721",
    "c373",
    ## mother prenatal smoking
    "b669",#Smoking before pregnancy,antenatal,mother
    "b670",#Smoking during first timester,antenatal,mother
    "b671",#Smoking during second trimester,antenatal,mother
    "c482",#Smoking during third trimester,antenatal,mother                  
    "cotinine_Preg_Trim1_avg", ## blood cotinine, first trimester
    "cdadj_preg", ## cadmium adjusted
    ## child smoking 14y
    "ccr700",#Smoked,child
    "ccr740",#Frequency (per day groups),child
    ## child smoking, tf3, 15.5y
    "fh8410",#Smoked any (by mother),child
    "Cotinine_TF3", ## cotinine levels
    ##  child smoking 16.5y
    "ccs4000",#Smoked,F17,child
    "ccs4040", ## G6: Number of cigarettes YP smokes on a daily basis
    "cotinine_TF4_avg", ## cotinine
    ## child smoking 24+
    "YPD7000", ## ever smoked a cig
    "YPD7001", ## number of cigs in their lifetime
    "YPD7020", ## Respondent smokes every day
    "YPD7030", ## Respondent smokes every week
    ## mother smoking
    "a200",
    "b650",
    "e179",#Smoked since birth,0.15,mother
    "e171",
    "f620",
    "g820",
    "h720",
    "j735",
    "k6180",#Frequency (per day groups),5,mother
    "l5050",
    "l5051",
    "m5160",
    "n5000",
    "r6010",               
    "s1300",
    "s1301",
    "t5520", ## currently a smoker
    "t5560", ## ever smoked in the past
    ## mother stopped smoking
    "b659",#(Y=1, N=2, missing=-1) Mother has stopped smoking (pregnancy)
    "n5006",#(Y=1, N=2, missing=-1) Mother has now stopped smoking  (8y)
    "r6016",#(Y=1, N=2, missing=-1) Respondent has now stopped smoking (11y)
    "t5582",
    ## partner smoking
    "b685",#Frequency (per day; by mother),antenatal,partner
    "pb077",#Smokes,antenatal,partner              
    "pb074", #Stopped smoking now
    "pb077", #Smoked REG in last 9MTHS
    "pc260", #NO smoked daily in last 2MTHS of PREG
    "pd620", #Cigarettes Smoked Per Day
    "pe450", #Number of Cigarettes Smoked Per Day
    "pf7090", #Number of cigarettes smokes by partner per day 
    "ph6180", #Number of cigarettes respondent smokes per day
    "pj5050", #Number of cigarettes smoked by respondent per day on weekdays
    "pj5051", #Number of cigarettes smoked by respondent per day on weekends 
    "pk5160", #E9: Number of cigarettes partner smokes per day              
    "pl5000", #J1a: Respondent has ever been a smoker
    "pl5006", #J1d: Respondent has stopped smoking
    "pp6010", #G2a: Respondent has ever been a smoker
    "pp6016", #G2d: Respondent has now stopped smoking
    "pq1300", #A8a1: Number of cigarettes partner smokes per weekday
    "pq1301", #A8a2: Number of cigarettes partner smokes per weekend day
    "pq1302", #A8b1: Partner smokes a pipe
    "pq1303", #A8b2: Partner smokes cigars/cigarillos              
    "fa5520", #E21: Currently a smoker (cigarettes or tobacco): FoF1
    "fa5560", #E31: Ever smoked in past: FoF1
    "fa5582", # age when stopped 
    ## partner drinking
    "pc281", ## Amount of alcohol consumed per day
    "pe452", ## Amount of alcohol consumed per week
    "pf6100", ## Amount of alcohol consumed per week
    "pg5050", ## Amount of alcohol consumed per week
    "ph6190", ## Amount of alcohol consumed per week
    ## mother drinking
    "t5500", ##E60a: Frequency of alcohol consumption
    "k6190", ##F12a: Amount of alcohol mother drinks
    "h723", ##Quantity of alcohol mother drinks
    "g822", ##Quantity of alcohol mother drinks
    "f625", ##Alcohol consumption
    "e221", ##FREQ of alcohol use since birth
    ## mother: Alcohol Use Disorders Identification Test (AUDIT)
    "t5510")

varnames <- tolower(varnames)
variable.index <- alspac::findVars(varnames)
variable.index <- variable.index[which(tolower(variable.index$name) %in% varnames),]
stopifnot(all(varnames %in% tolower(variable.index$name)))
variable.index <- alspac::filterVars(variable.index,
    b032=c(cat1="Current"),
    c645a=c(cat1="Current"),
    c804=c(cat1="Current"),
    kz030=c(cat1="Current"),
    t5510=c(cat1="Current"),
    bestgest=c(obj="^mz"))

stopifnot(
    length(setdiff(tolower(varnames),
                   tolower(variable.index$name)))
    ==0)

alspac <- alspac::extractVars(variable.index)


set.cond <- function(x,condition,value) {
    original <- x
    x[which(condition)] <- value
    #if (is.numeric(x) & length(unique(x)) > 10)
    #    print(quantile(x, na.rm=T))
    #else
    #    print(table(original=original, new=x, useNA="always"))
    x
}
map <- function(x, ...) {
    key <- list(...)
    x[which(!x %in% names(key))] <- NA
    xp <- unlist(key)[match(x, names(key))]
    names(xp) <- names(x)
    #print(table(original=x, mapped=xp, useNA="ifany"))
    xp
}
or.na <- function(...) {
    x <- sapply(list(...), function(x) x)
    is.na <- rowSums(!is.na(x)) == 0
    ret <- rowSums(x, na.rm=T) > 0
    ret[is.na] <- NA
    ret
}
not.false <- function(x) {
    idx <- which(is.na(x))
    x[idx] <- TRUE
    x
}
is.true <- function(x) {
    idx <- which(is.na(x))
    x[idx] <- FALSE
    x
}
and.na <- function(...) {
    x <- do.call(cbind, list(...))
    rowSums(x, na.rm=T) == ncol(x)
}


alspac$occupation <- set.cond(alspac$c755,!alspac$c755 %in% 1:6, NA)
alspac$parity <- set.cond(alspac$b032,alspac$b032<0,NA)
alspac$parity <- set.cond(alspac$parity,alspac$b032>2,3)
alspac$education <- set.cond(alspac$c645a,alspac$c645a<0,NA)
alspac$gestational.age <- set.cond(alspac$bestgest, alspac$bestgest<0,NA)
alspac$ethnicity <- set.cond(alspac$c804,alspac$c804<0,NA)==1
alspac$birthweight <- set.cond(alspac$kz030,alspac$kz030<0,NA)

alspac$bmi.pre <- with(alspac, set.cond(dw042,dw042<0,NA))
alspac$bmi.middle.18 <- with(alspac,set.cond(fm1ms111,fm1ms111<0,NA))
alspac$bmi.middle.20 <- with(alspac,set.cond(fm2ms112,fm2ms112<0,NA))
alspac$bmi.middle <- with(alspac, rowMeans(cbind(bmi.middle.18, bmi.middle.20), na.rm=T))

alspac$bmi.father <- with(alspac, set.cond(ff1ms111,ff1ms111<0,NA))
alspac$bmi.child.7 <- with(alspac, set.cond(f7ms026a,f7ms026a<0,NA))
alspac$bmi.child.10 <- with(alspac, set.cond(fdms026a,fdms026a<0,NA))
alspac$bmi.child.17 <- with(alspac, set.cond(FJMR022a,FJMR022a<0,NA))
alspac$bmi.child.24 <- with(alspac, set.cond(FKMS1040,FKMS1040<0,NA))

alspac$smoking.pre <- with(alspac, set.cond(b669, b669<0, NA)) > 0

alspac$smoking.prenatal <- with(alspac, {
    sustained <- rowSums(cbind(set.cond(b670,b670 < 0,NA) > 0,
                               set.cond(b671,b671 < 0,NA) > 0),
                         na.rm=T)
    sustained[which(sustained == 1)] <- NA
    sustained > 0
})

alspac$cotinine.prenatal <- with(alspac,
    set.cond(cotinine_Preg_Trim1_avg, cotinine_Preg_Trim1_avg<0,NA))

alspac$cadmium.prenatal <- with(alspac,
    set.cond(Cdadj_preg, Cdadj_preg <0, NA))


never.smoking.recent <- with(alspac,
    cbind(and.na(r6010==2, is.na(r6016)), ## has never been a smoker
          or.na(s1300==0, s1301==0), ## smokes no cigarettes
          and.na(t5520 == 2, t5560==2))) ## is not and never been a smoker

never.smoking <- with(alspac,
    cbind(a200==0, ## no cigarettes
          and.na(b650==2, is.na(b659)), ## never smoked nor quit
          ## no cigarettes
          c482==0, e171==2, e179==0,
          f620==0, g820==0, h720==0, j735==0, k6180==0,
          l5050==0, l5051==0, m5160==0,
          and.na(n5000==2, is.na(n5006)), ## has never been a smoker
          never.smoking.recent))

ever.smoking <- with(alspac,
    cbind(a200>0, ## number of cigarettes
          or.na(b650==1, b659 == 1), ## ever smoked or has quit
          ## number of cigarettes
          c482>0, e179>0, f620>0, g820>0, h720>0, j735>0, k6180>0,
          l5050>0, l5051>0, m5160>0,
          ## has ever been a smoker or has quit
          or.na(n5000==1, !is.na(n5006)),
          or.na(r6010==1, !is.na(r6016)), ## has been a smoker and not quit
          or.na(s1300>0, s1301>0), ## number of cigarettes
          or.na(t5520 == 1, t5560==1))) ## is or has been a smoker

ever.quit <- with(alspac,
    (is.true(b659==1) | is.true(n5006==1) | is.true(r6016==1) ## has quit
     | !is.na(t5582))) ## age when stopped

ever.smoker <- rowSums(ever.smoking, na.rm=T) > 0
current.smoker <- with(alspac, ever.smoker & t5520==1)
never.smoker <- rowSums(ever.smoking, na.rm=T) == 0 & rowSums(never.smoking.recent, na.rm=T) > 0
former.smoker <- ever.smoker & ever.quit & alspac$t5520 == 2

alspac$smoking.middle <- NA
alspac$smoking.middle[which(never.smoker)] <- "never"
alspac$smoking.middle[which(current.smoker)] <- "current"
alspac$smoking.middle[which(former.smoker)] <- "former"


alspac$smoking.mother.test <- NA
alspac$smoking.mother.test[with(alspac, which(and.na(t5520==2, t5560==2)))] <- 0 ## not currently and never been a smoker
alspac$smoking.mother.test[with(alspac, which(or.na(t5520==1, t5560==1)))] <- 1 ## is or has been a smoker
with(alspac, table(smoking.middle, smoking.mother.test, useNA="always"))
alspac$smoking.middle[which(alspac$smoking.mother.test == 0)] <- "never"


never.smoking <- with(alspac,
    cbind(b685 == 0, ## number of times smokes per day
          pb077 == 1, ## smoked regularly in the last 9 months
          pc260 == 0, ## number of times smokes per day last 2 months preg
          pd620 == 0, ## number of cigarettes smoked per day
          pe450 == 0, ## number of cigarettes smoked per day
          pf7090 == 0, ## number of cigarettes smoked per day
          ph6180 == 0, ## number of cigarettes smoked per day
          pj5050 + pj5051 == 0, ## number of cigarettes per day
          pk5160 == 0, ## number of cigarettes per day
          pl5000 == 2, ## ever smoker
          pp6010 == 2, ## ever smoker
          !or.na(pq1300 + pq1301 > 0, pq1302 != 3, pq1303 != 3),
          ## number of cigarettes per day, pipe, cigar
          and.na(fa5560==2, fa5520==2))) ## never and not currently a smoker

ever.smoking <- with(alspac,
    cbind(b685 > 0, ## number of times smokes per day
          pb077 > 1, ## smoked regularly in the last 9 months
          pc260 > 0, ## number of times smokes per day last 2 months preg
          pd620 > 0, ## number of cigarettes smoked per day
          pe450 > 0, ## number of cigarettes smoked per day
          pf7090 > 0, ## number of cigarettes smoked per day
          ph6180 > 0, ## number of cigarettes smoked per day
          pj5050 + pj5051 > 0, ## number of cigarettes per day
          pk5160 > 0, ## number of cigarettes per day
          pl5000 == 1, ## ever smoker
          pp6010 == 1, ## ever smoker
          or.na(pq1300 + pq1301 > 0, pq1302 != 3, pq1303 != 3),
          ## number of cigarettes per day, pipe, cigar
          or.na(fa5560==1, fa5520==1))) ## ever or currently a smoker

ever.quit <- with(alspac,
    (is.true(pb074==1) | is.true(pl5006==1) | is.true(pp6016==1) ## has quit
     | !is.na(fa5582))) ## age when stopped

ever.smoker <- rowSums(ever.smoking, na.rm=T) > 0
current.smoker <- with(alspac, ever.smoker & fa5520==1)
never.smoker <- (rowSums(ever.smoking, na.rm=T) == 0
                 & rowSums(never.smoking, na.rm=T) > 5)
former.smoker <- with(alspac, ever.smoker & ever.quit & fa5520==2)

alspac$smoking.father <- NA
alspac$smoking.father[which(never.smoker)] <- "never"
alspac$smoking.father[which(current.smoker)] <- "current"
alspac$smoking.father[which(former.smoker)] <- "former"

alspac$smoking.father.test <- NA
alspac$smoking.father.test[with(alspac, which(and.na(fa5560==2, fa5520==2)))] <- 0 ## never and not currently a smoker
alspac$smoking.father.test[with(alspac, which(or.na(fa5560==1, fa5520==1)))] <- 1 ## ever or currently a smoker
with(alspac, table(smoking.father, smoking.father.test, useNA="always"))

alspac$smoking.child.17 <- NA
alspac$smoking.child.17[which(or.na(alspac$ccr700 == 2,
                                    alspac$ccs4000 == 2,
                                    alspac$fh8410 == 2))] <- 0
alspac$smoking.child.17[which(or.na(alspac$ccr740 >= 1,
                                    alspac$ccs4040 >= 1))] <- 1


alspac$smoking.child.24 <- NA
alspac$smoking.child.24[which(or.na(alspac$YPD7000 == 0,
                                    alspac$YPD7001 %in% 1:4))] <- 0
alspac$smoking.child.24[which(or.na(alspac$YPD7001 == 5,
                                    alspac$YPD7020 == 1,
                                    alspac$YPD7030 == 1))] <- 1

alspac$cotinine.child.15 <- with(alspac,
    set.cond(Cotinine_TF3, Cotinine_TF3<0,NA))
alspac$cotinine.child.17 <- with(alspac,
    set.cond(cotinine_TF4_avg, cotinine_TF4_avg<0,NA))

alspac$alcohol.prenatal <- with(alspac,
    rowSums(cbind(set.cond(a260, a260 < 0, NA) > 0,
                  set.cond(b721, b721 < 0, NA) > 1,
                  set.cond(c373, c373 < 0, NA) > 0),
            na.rm=T) > 0)

alspac$alcohol.middle <- rowMeans(sapply(alspac[,c("t5500","k6190","h723","g822","f625","e221")], function(variable) {
    variable <- set.cond(variable, variable<0,NA)
    as.vector(scale(variable))
}), na.rm=T)

alspac$alcohol.audit <- with(alspac, set.cond(t5510, t5510<0,NA))

alspac$alcohol.father <- rowMeans(sapply(alspac[,c("pc281","pe452","pf6100","pg5050","ph6190")], function(variable) {
    variable <- set.cond(variable, variable <0 | variable > 6, NA)
    as.vector(scale(variable))
}), na.rm=T)
