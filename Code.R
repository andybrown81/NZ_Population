########################################
###     Load Data and Packages       ###
########################################

### Set working directory   
setwd("\\\\bournemouth.ac.uk/data/staff/home/aabrown/Profile/Desktop"

### Install and load packages 
# R version: 3.4.0
install.packages("doParallel")
install.packages("devtools")
install_github("ahb108/rcarbon")

library(doParallel)
library(devtools)
library(rcarbon)


### Read in 14C dates        
nzdates <- read.csv("NZall.csv")

### Subset dates based on CRA. Pre 800 CRA is likely to be an 'old wood' date
nzdates <- subset(nzdates,CRA>0)
nzdates <- subset(nzdates, CRA<800)

### Subset by region 
northern_dates <- subset(nzdates, Region=="Northern")
central_dates <- subset(nzdates, Region=="Central")
southern_dates <- subset(nzdates, Region=="Southern")

### Check region lengths
length(northern_dates)
length(central_dates)
length(southern_dates)



#################################
###   Calibrate and Binning   ###
#################################


### Calibrate dates using southern hemisphere curve
caldates_northern <- calibrate(ages=northern_dates$CRA,errors=northern_dates$Error,calCurves="shcal13")
caldates_central <- calibrate(ages=central_dates$CRA,errors=central_dates$Error,calCurves="shcal13")
caldates_southern <- calibrate(ages=southern_dates$CRA,errors=southern_dates$Error,calCurves="shcal13")

### Make bins (100year cut-off)
bins_N <- binPrep(sites=northern_dates$Site.Number,ages=northern_dates$CRA,h=100)
bins_C <- binPrep(sites=central_dates$Site.Number,ages=central_dates$CRA,h=100)
bins_S <- binPrep(sites=southern_dates$Site.Number,ages=southern_dates$CRA,h=100)

### Check number of unique bins (Northern - 102 | Central - 55 | Southern - 58)
NorthernBins <- unique(bins_N)
length(NorthernBins)
CentralBins <- unique(bins_C)
length(CentralBins)
SouthernBins <- unique(bins_S)
length(SouthernBins)



#######################
###   Create SPDs   ###
#######################
 


### Make regional SPDs using 50 year rolling mean (Northern, Central, Southern in decending order)
spd.n <- spd(x=caldates_northern,bins=bins_N,timeRange=c(800,0),runm=50)
spd.c <- spd(x=caldates_central,bins=bins_C,timeRange=c(800,0),runm=50)
spd.s <- spd(x=caldates_southern,bins=bins_S,timeRange=c(800,0),runm=50)

### Plot SPDs
par(mfrow = c(3, 1))
plot(spd.n)
plot(spd.c)
plot(spd.s)



#########################################
###   Model Testing Northern Region   ###
#########################################


### Exponential Model Test (5000 simulations)
expMod_northern <- modelTest(x=caldates_northern,bins=bins_N,errors=northern_dates$Error,runm=50,timeRange=c(800,0),model="explog",calCurves="shcal13",nsim=5000,ncores=2)

### Logistic Model (5000 simulations)
grd_N <- spd.n$grid
x <- grd_N$calBP
y <- grd_N$PrDens
log.ss <- nls(y~SSlogis(x, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.2,xmid=400,scale=-100))
logisticFit_N <- data.frame(calBP=x,PrDens=SSlogis(x,coefficients(log.ss)[1],coefficients(log.ss)[2],coefficients(log.ss)[3]))

### Check fit of logistic model to observed data
plot(spd.n)
lines(logisticFit_N,col="red",lty=2)

### Run model
logMod_northern <- modelTest(x=caldates_northern,bins=bins_N,errors=northern_dates$Error,predgrid=logisticFit_N,runm=50,timeRange=c(800,0),model="custom",calCurves="shcal13",nsim=5000,ncores=2)

### Plot regional models
par(mfrow = c(2, 1))
plot(expMod_northern, ylim = c(0.0, 0.5)) 
plot(logMod_northern, ylim = c(0.0, 0.5))

### Get global p values
expMod_northern$pval
logMod_northern$pval

########################################
###   Model Testing Central Region   ###
########################################


### Exponential Model (5000 simulations)
expMod_central <- modelTest(x=caldates_central,bins=bins_C,errors=central_dates$Error,runm=50,timeRange=c(800,0),model="explog",calCurves="shcal13",nsim=5000,ncores=2)

### Logistic Model
grd_C <- spd.c$grid
x <- grd_C$calBP
y <- grd_C$PrDens
log.ss <- nls(y~SSlogis(x, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.8,xmid=400,scale=-100))
logisticFit_C <- data.frame(calBP=x,PrDens=SSlogis(x,coefficients(log.ss)[1],coefficients(log.ss)[2],coefficients(log.ss)[3]))

### Check fit of model to observed data and run model
plot(spd.c)
lines(logisticFit_C, col = "red", lty = 2)
logMod_central <- modelTest(x=caldates_central, bins=bins_C, errors=central_dates$Error, predgrid=logisticFit_C, runm=50, timeRange=c(800,0), model="custom", calCurves="shcal13", nsim=5000,ncores=2)

### Plot regional models
plot(expMod_central, ylim = c(0.0,0.2))
plot(logMod_central, ylim = c(0.0,0.2))

### Get global p values
expMod_central$pval
logMod_central$pval


#########################################
###   Model Testing Southern Region   ###
#########################################


### Exponential Model (5000 simulations)
expMod_southern <- modelTest(x=caldates_southern,bins=bins_S,errors=southern_dates$Error,runm=50,timeRange=c(800,0),model="explog",calCurves="shcal13",nsim=5000,ncores=2)   

### Logistic Model
grd_S <- spd.s$grid
x <- grd_S$calBP
y <- grd_S$PrDens
log.ss <- nls(y~SSlogis(x, Asym, xmid, scale),control=nls.control(maxiter=500),start=list(Asym=0.2,xmid=400,scale=-100))
logisticFit_S <- data.frame(calBP=x,PrDens=SSlogis(x,coefficients(log.ss)[1],coefficients(log.ss)[2], coefficients(log.ss)[3]))

### Check fit of model to observed data and run model
plot(spd.s)
lines(logisticFit_S,col="red",lty=2)
logMod_southern <- modelTest(x=caldates_southern,bins=bins_S,errors=southern_dates$Error,predgrid=logisticFit_S,runm=50,timeRange=c(800,0),model="custom",calCurves="shcal13",nsim=5000,ncores=2)

### Plot regional models
plot(expMod_southern, ylim = c(0,0.25))
plot(logMod_southern, ylim = c(0,0.25))

### Get global p values
expMod_southern$pval 
logMod_southern$pval



############################
###   Plot All Regions   ###
############################


par(mfrow = c(2, 3))

plot(expMod_northern, ylim = c(0.0,0.5))
plot(logMod_northern, ylim = c(0.0,0.5))
lines(logisticFit_N, col = "red", lty = 2)

plot(expMod_central, ylim = c(0.0,0.3))
plot(logMod_central, ylim = c(0.0,0.3))
lines(logisticFit_C, col = "red", lty = 2)

plot(expMod_southern, ylim = c(0,0.4))
plot(logMod, ylim = c(0,0.4))
lines(logisticFit_S,col="red",lty=2)



###########################
###  Permutation Tests	###
###########################

### Run permutations
Perm <- permTest(x=caldates, marks = nzdates$Region, timeRange = c(800,0), nsim=5000, runm=50)

### Plot results for each region
par(mfrow = c(1, 3))
plot.SpdPermTest(Perm, focalm='1')
plot.SpdPermTest(Perm, focalm='2')
plot.SpdPermTest(Perm, focalm='3')




