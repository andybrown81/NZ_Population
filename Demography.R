########################################
###     Load Data and Packages       ###
########################################

### Set working directory   
setwd("")

### Install and load packages 
# R version: 3.4.0
install.packages("doParallel")
install.packages("devtools")
install.packages("zoo")
install.packages("sp")
install.packages("R.utils")
install.packages("rworldmap")
install.packages("rcarbon")

library(doParallel)
library(devtools)
library(zoo)
library(sp)
library(R.utils)
library(rworldmap)
library(rcarbon)


### Read in 14C dates        
nzdates <- read.csv("SPDdata.csv")

### Subset dates based on CRA. Pre 800 CRA is likely to be an 'old wood' date
nzdates <- subset(nzdates,C14Age>0)
nzdates <- subset(nzdates, C14Age<800)

### Subset by region 
northern_dates <- subset(nzdates, Region=="Northern")
central_dates <- subset(nzdates, Region=="Central")
southern_dates <- subset(nzdates, Region=="Southern")

### Check region lengths (ND - 154 | CD - 72 | SD - 89)
length(unique(northern_dates$LabID))
length(unique(central_dates$LabID))
length(unique(southern_dates$LabID))

setInternet2(TRUE)

###################################
###   Calibration and Binning   ###
###################################


### Calibrate dates using southern hemisphere curve
caldates_northern <- calibrate(northern_dates$C14Age,northern_dates$C14SD,calCurves="shcal13", normalised=FALSE, verbose=FALSE)
caldates_central <- calibrate(central_dates$C14Age,central_dates$C14SD,calCurves="shcal13", normalised=FALSE, verbose=FALSE)  
caldates_southern <- calibrate(southern_dates$C14Age,southern_dates$C14SD,calCurves="shcal13", normalised=FALSE, verbose=FALSE)  


### Make bins (100 year cut-off)
bins_N <- binPrep(sites=northern_dates$SiteID,ages=northern_dates$C14Age,h=100)
bins_C <- binPrep(sites=central_dates$SiteID,ages=central_dates$C14Age,h=100)
bins_S <- binPrep(sites=southern_dates$SiteID,ages=southern_dates$C14Age,h=100)

### Check number of unique bins (Northern - 110 | Central - 46 | Southern - 54)
length(unique(bins_N))
length(unique(bins_C))
length(unique(bins_S))



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
set.seed(12345)
expMod_northern <- modelTest(x=caldates_northern,bins=bins_N,errors=northern_dates$C14SD,runm=50,timeRange=c(800,0),model="exponential",calCurves="shcal13",nsim=5000,ncores=2)
par(mfrow = c(1, 1))
plot(expMod_northern)
lines(expMod_northern$fit,lwd=2,lty=2)


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
set.seed(12345)
logMod_northern <- modelTest(x=caldates_northern,bins=bins_N,errors=northern_dates$C14SD,predgrid=logisticFit_N,runm=50,timeRange=c(800,0),model="custom",calCurves="shcal13",nsim=5000,ncores=2)

### Plot regional models
par(mfrow = c(1, 2))
plot(expMod_northern, ylim = c(0.0, 0.5)) 
plot(logMod_northern, ylim = c(0.0, 0.5))

### Check global p values
expMod_northern$pval
logMod_northern$pval

########################################
###   Model Testing Central Region   ###
########################################


### Exponential Model (5000 simulations)
set.seed(12345)
expMod_central <- modelTest(x=caldates_central,bins=bins_C,errors=central_dates$C14SD,runm=50,timeRange=c(800,0),model="exponential",calCurves="shcal13",nsim=5000,ncores=2)
plot(expMod_central)
lines(expMod_central$fit,lwd=2,lty=2)

### Logistic Model
grd_C <- spd.c$grid
x <- grd_C$calBP
y <- grd_C$PrDens
log.ss <- nls(y~SSlogis(x, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.2,xmid=400,scale=-100))
logisticFit_C <- data.frame(calBP=x,PrDens=SSlogis(x,coefficients(log.ss)[1],coefficients(log.ss)[2],coefficients(log.ss)[3]))

### Check fit of model to observed data 
plot(spd.c)
lines(logisticFit_C, col = "red", lty = 2)

### Run model
set.seed(12345)
logMod_central <- modelTest(x=caldates_central, bins=bins_C, errors=central_dates$C14SD, predgrid=logisticFit_C, runm=50, timeRange=c(800,0), model="custom", calCurves="shcal13", nsim=5000,ncores=2)

### Plot regional models
plot(expMod_central, ylim = c(0.0,0.25))
plot(logMod_central, ylim = c(0.0,0.25))

### Check global p values
expMod_central$pval
logMod_central$pval


#########################################
###   Model Testing Southern Region   ###
#########################################


### Exponential Model (5000 simulations)
set.seed(12345)
expMod_southern <- modelTest(x=caldates_southern,bins=bins_S,errors=southern_dates$C14SD,runm=50,timeRange=c(800,0),model="exponential",calCurves="shcal13",nsim=5000,ncores=2)   
plot(expMod_southern)
lines(expMod_southern$fit,lwd=2,lty=2)

### Logistic Model
grd_S <- spd.s$grid
x <- grd_S$calBP
y <- grd_S$PrDens
log.ss <- nls(y~SSlogis(x, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.2,xmid=400,scale=-100))
logisticFit_S <- data.frame(calBP=x,PrDens=SSlogis(x,coefficients(log.ss)[1],coefficients(log.ss)[2], coefficients(log.ss)[3]))

### Check fit of model to observed data 
plot(spd.s)
lines(logisticFit_S,col="red",lty=2)

### Run model
set.seed(12345)
logMod_southern <- modelTest(x=caldates_southern,bins=bins_S,errors=southern_dates$C14SD,predgrid=logisticFit_S,runm=50,timeRange=c(800,0),model="custom",calCurves="shcal13",nsim=10,ncores=2)

### Plot regional models
plot(expMod_southern, ylim = c(0,0.25))
plot(logMod_southern, ylim = c(0,0.25))

### Check global p values
expMod_southern$pval 
logMod_southern$pval



############################
###   Plot All Regions   ###
############################


par(mfrow = c(2, 3))

plot(expMod_northern, ylim = c(0.0,0.5))
plot(logMod_northern, ylim = c(0.0,0.5))
lines(logisticFit_N, col = "red", lty = 2)

plot(expMod_central, ylim = c(0.0,0.25))
plot(logMod_central, ylim = c(0.0,0.25))
lines(logisticFit_C, col = "red", lty = 2)

plot(expMod_southern, ylim = c(0,0.2))
plot(logMod_southern, ylim = c(0,0.2))
lines(logisticFit_S,col="red",lty=2)

Pvalues <- matrix( c(expMod_northern$pval, logMod_northern$pval, expMod_central$pval, logMod_central$pval, expMod_southern$pval, logMod_southern$pval), nrow=6, ncol=1, byrow=FALSE)

###########################
###  Permutation Tests	###
###########################

### Calibrate NZ dates
caldates <- calibrate(ages=nzdates$CRA,errors=nzdates$Error,calCurves="shcal13")

### Run permutations
set.seed(12345)
Perm <- permTest(x=caldates, marks = nzdates$Region, timeRange = c(800,0), nsim=5000, runm=50)

### Plot results for each region
par(mfrow = c(1, 3))
plot.SpdPermTest(Perm, focalm='2')
plot.SpdPermTest(Perm, focalm='1')
plot.SpdPermTest(Perm, focalm='3')

### Get global p values (Note 1 = Central, 2 = Northern, 3 = Southern)
Perm$pValueList


#########################################
###     Spatial Permutation Test      ###
#########################################


###   Calibrate and Binning   ###


caldates <- calibrate(nzdates$C14Age,nzdates$C14SD,calCurves="shcal13", normalised=FALSE, verbose=FALSE)
bins <- binPrep(sites=nzdates$SiteID,ages=nzdates$C14Age,h=100)


###   Prepare Spatial Data    ###

sites <- unique(data.frame(id=nzdates$SiteID,lat=nzdates$Latitude,lon=nzdates$Longitude))
row.names(sites) <- sites$id
sites <- sites[,-1]


coordinates(sites) <- c("lon","lat")
proj4string(sites) <- CRS("+proj=longlat +datum=WGS84")

###   Compute distance matrix
d <- spDists(sites,sites,longlat=TRUE)

###   Compute spatial weights
w <- spweights(d,h=50)



################################
###    Execute SPpermTest    ###
################################

breaks <- seq(700,100,-300) #500 year blocks
timeRange <- c(700,100) #set the timerange of analysis in calBP, older date first

res <- SPpermTest(calDates=caldates,bins=bins,timeRange=timeRange,
                  locations=sites,permute="locations",nsim=1000,
                  breaks=breaks,spatialweights=w,ncores=3,verbose=FALSE)



###   Plot Results             

library(rworldmap)
base <- getMap(resolution="low") #extract basemap
#extract bounding coordinates of the site distribution
xrange <- bbox(sites)[1,]
yrange <- bbox(sites)[2,]





