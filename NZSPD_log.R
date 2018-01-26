########################################
###     Load Data and Packages       ###
########################################


### Install and load packages 
# R version: 3.4.0
install.packages("doParallel")
install.packages("devtools")
install.packages("zoo")
install.packages("sp")
install.packages("R.utils")
install.packages("rworldmap")
install_github("ahb108/rcarbon")
install.packages("RCurl")

library(doParallel)
library(devtools)
library(zoo)
library(sp)
library(R.utils)
library(rworldmap)
library(rcarbon)
library(RCurl)


require(RCurl)

### Read in 14C dates 
nzdates <-read.csv(text=getURL("https://raw.githubusercontent.com/andybrown81/NZ_Population/master/SPDdata.csv"))

### Subset dates based on CRA. Pre 800 CRA is likely to be an 'old wood' date
nzdates <- subset(nzdates,C14Age>0)
nzdates <- subset(nzdates, C14Age<800)

### Subset by region 
northern_dates <- subset(nzdates, Region=="Northern")
central_dates <- subset(nzdates, Region=="Central")
southern_dates <- subset(nzdates, Region=="Southern")

### Check region lengths (ND - 153 | CD - 73 | SD - 89)
length(unique(northern_dates$LabID))
length(unique(central_dates$LabID))
length(unique(southern_dates$LabID))


###################################
###   Calibration and Binning   ###
###################################

### Calibrate dates using southern hemisphere curve
caldates_northern <- calibrate(northern_dates$C14Age,northern_dates$C14SD,calCurves="shcal13", normalised=FALSE, verbose=FALSE)
caldates_central <- calibrate(central_dates$C14Age,central_dates$C14SD,calCurves="shcal13", normalised=FALSE, verbose=FALSE)  
caldates_southern <- calibrate(southern_dates$C14Age,southern_dates$C14SD,calCurves="shcal13", normalised=FALSE, verbose=FALSE)  
caldates <- calibrate(nzdates$C14Age,nzdates$C14SD,calCurves="shcal13", normalised=FALSE, verbose=FALSE)

### Make bins (100 year cut-off)
bins_N <- binPrep(sites=northern_dates$SiteID,ages=northern_dates$C14Age,h=100)
bins_C <- binPrep(sites=central_dates$SiteID,ages=central_dates$C14Age,h=100)
bins_S <- binPrep(sites=southern_dates$SiteID,ages=southern_dates$C14Age,h=100)
bins <- binPrep(sites=nzdates$SiteID,ages=nzdates$C14Age,h=100)

### Check number of unique bins (Northern - 112 | Central - 47 | Southern - 55)
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
spd.nz <- spd(x=caldates,bins=bins,timeRange=c(800,0),runm=50)

### Plot SPDs
par(mfrow = c(2, 2))
plot(spd.n)
plot(spd.c)
plot(spd.s)
plot(spd.nz)


#########################################
###         Model Testing NZ          ###
#########################################

### Logistic Model
grd_nz <- spd.nz$grid
x <- grd_nz$calBP
y <- grd_nz$PrDens
log.ss <- nls(y~SSlogis(x, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.8,xmid=400,scale=-100))
logisticFit_nz <- data.frame(calBP=x,PrDens=SSlogis(x,coefficients(log.ss)[1],coefficients(log.ss)[2], coefficients(log.ss)[3]))

### Check fit of model to observed data 
plot(spd.nz)
lines(logisticFit_nz,col="red",lty=2)

### Run model
set.seed(12345)
logMod_nz <- modelTest(x=caldates, bins=bins, errors=nzdates$C14SD,predgrid=logisticFit_nz,runm=50,timeRange=c(800,0),model="custom",calCurves="shcal13",nsim=5000,ncores=2, raw=TRUE)

### Plot model
plot(logMod_nz, ylim = c(0,0.5))

### Check global p values
logMod_nz$pval


#########################################
###   Model Testing Northern Region   ###
#########################################


### Logistic Model (5000 simulations)
grd_N <- spd.n$grid
x <- grd_N$calBP
y <- grd_N$PrDens
log.ss <- nls(y~SSlogis(x, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.6,xmid=400,scale=-100))
logisticFit_N <- data.frame(calBP=x,PrDens=SSlogis(x,coefficients(log.ss)[1],coefficients(log.ss)[2],coefficients(log.ss)[3]))

### Check fit of logistic model to observed data
plot(spd.n)
lines(logisticFit_N,col="red",lty=2)

### Run model
set.seed(12345)
logMod_northern <- modelTest(x=caldates_northern,bins=bins_N,errors=northern_dates$C14SD,predgrid=logisticFit_N,runm=50,timeRange=c(800,0),model="custom",calCurves="shcal13",nsim=5000,ncores=2, raw=TRUE)

### Plot regional models
par(mfrow = c(1, 1))
plot(logMod_northern, ylim = c(0.0, 0.5))

### Check global p values
logMod_northern$pval


########################################
###   Model Testing Central Region   ###
########################################

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
logMod_central <- modelTest(x=caldates_central, bins=bins_C, errors=central_dates$C14SD, predgrid=logisticFit_C, runm=50, timeRange=c(800,0), model="custom", calCurves="shcal13", nsim=5000,ncores=2, raw=TRUE)

### Plot regional models
plot(logMod_central, ylim = c(0.0,0.25))

### Check global p values
logMod_central$pval


#########################################
###   Model Testing Southern Region   ###
#########################################

### Logistic Model
grd_S <- spd.s$grid
x <- grd_S$calBP
y <- grd_S$PrDens
log.ss <- nls(y~SSlogis(x, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.8,xmid=400,scale=-100))
logisticFit_S <- data.frame(calBP=x,PrDens=SSlogis(x,coefficients(log.ss)[1],coefficients(log.ss)[2], coefficients(log.ss)[3]))

### Check fit of model to observed data 
plot(spd.s)
lines(logisticFit_S,col="red",lty=2)

### Run model
set.seed(12345)
logMod_southern <- modelTest(x=caldates_southern,bins=bins_S,errors=southern_dates$C14SD,predgrid=logisticFit_S,runm=50,timeRange=c(800,0),model="custom",calCurves="shcal13",nsim=5000,ncores=2,raw=TRUE)

### Plot regional models
plot(logMod_southern, ylim = c(0,0.25))

### Check global p values
logMod_southern$pval



############################
###   Plot All Regions   ###
############################

par(mfrow = c(2, 2))

plot(logMod_nz, ylim = c(0.0,0.6))
lines(logisticFit_nz, col = "red", lty = 2)

plot(logMod_northern, ylim = c(0.0, 0.5))
lines(logisticFit_N, col = "red", lty = 2)

plot(logMod_central, ylim = c(0.0,0.25))
lines(logisticFit_C, col = "red", lty = 2)

plot(logMod_southern, ylim = c(0,0.2))
lines(logisticFit_S,col="red",lty=2)

Pvalues <- matrix(c(logMod_nz$pval, logMod_northern$pval, logMod_central$pval, logMod_southern$pval),  nrow=4, ncol=1, byrow=FALSE)
rownames(Pvalues) <- c("New Zealand", "Northern", "Central", "Southern")
colnames(Pvalues) <- c("Global Pvalue")
Pvalues



###########################
###  Permutation Tests	###
###########################


### Run permutations
set.seed(12345)
Perm <- permTest(x=caldates, marks = nzdates$Region, timeRange = c(800,0), nsim=5000, runm=50)

### Plot results for each region
par(mfrow = c(3, 1))
plot.SpdPermTest(Perm, focalm='2')
plot.SpdPermTest(Perm, focalm='1')
plot.SpdPermTest(Perm, focalm='3')

?plot.SpdPermTest

### Get global p values (Note 1 = Central, 2 = Northern, 3 = Southern)
Perm$pValueList


#########################################
###     Spatial Permutation Test      ###
#########################################

###   Prepare Spatial Data    ###
sites <- unique(data.frame(id=nzdates$SiteID,lat=nzdates$Latitude,lon=nzdates$Longitude))
row.names(sites) <- sites$id
sites <- sites[,-1]
coordinates(sites) <- c("lon","lat")
proj4string(sites) <- CRS("+proj=longlat +datum=WGS84")

###   Compute distance matrix
d <- spDists(sites,sites,longlat=TRUE)

###   Compute spatial weights
w <- spweights(d,h=100)


### Compute SPD and geometric growth rate ###

breaksize=100
nBreaks=length(breaks)-1
breaks <- seq(750,150,-100) 
timeRange <- c(750,150) 

spd=spd(x=caldates,timeRange=timeRange,bins=bins,datenormalised=FALSE,spdnormalised=TRUE)
nBreaks=length(breaks)-1
spd.blocksum=numeric()
spd.roc=numeric()

for (i in 1:nBreaks)
{
  spd.blocksum[i]=sum(spd$grid$PrDens[spd$grid$calBP<=breaks[i]&spd$grid$calBP>=breaks[i+1]])
}

for (i in 1:c(nBreaks-1))
{
  spd.roc[i]=(spd.blocksum[i+1]/spd.blocksum[i])^(1/breaksize)-1
}


###    Execute SPpermTest    ###


res <- SPpermTest(calDates=caldates,bins=bins,timeRange=timeRange,locations=sites,permute="locations",nsim=1000,
                  breaks=breaks,spatialweights=w,ncores=3,verbose=FALSE)


###   Plot Results             

library(rworldmap)
base <- getMap(resolution="low") #extract basemap
#extract bounding coordinates of the site distribution
xrange <- bbox(sites)[1,]
yrange <- bbox(sites)[2,]


par(mfrow=c(2,3))

for (i in 1:5)
{
  par(mar=c(0.1,0.1,0,0.5))
  plot(base,col="antiquewhite3",border="antiquewhite3",xlim=xrange,ylim=yrange)
  plot(res,index=i,add=TRUE,option="raw",breakRange=c(-0.01,0.01),breakLength=9,baseSize=1) 
  legend("topleft",legend=c(NA),border=NA,title=as.roman(i),cex=2,bty="n")
}


## Figure 5
plot(res,option="rawlegend",breakRange=c(-0.01,0.01),breakLength=9,rd=3,legSize=1.6)


#### Significance Testing ####


for (i in 1:5)
{
  par(mar=c(0.1,0.1,0,0.5))	
  plot(base,col="antiquewhite3",border="antiquewhite3",xlim=xrange,ylim=yrange)
  plot(res,index=i,add=TRUE,option="test",baseSize=1)
  legend("topleft",legend=c(NA),border=NA,title=as.roman(i),cex=2,bty="n")
}

## Figure 6
plot(res,option="testlegend",legSize=2)


#####################################
###         Figure 4              ###
#####################################

par(mfrow=c(1,2))

plot(spd.nz$grid$calBP, spd.nz$grid$PrDens,col=rgb(0,0,0,0.5),lwd=0.5,type="n",xlim=c(800,0),main="a",ylim=c(0,0.5),xlab="cal BP",ylab="SPD")

for (x in c(1,3,5,7))
{
  rect(xleft=breaks[x],xright=breaks[x+1],ybottom=-100,ytop=100,border=NA,col=rgb(0,0,0,0.2))
}

for (x in c(2,4,6))
{
  rect(xleft=breaks[x],xright=breaks[x+1],ybottom=-100,ytop=100,border=NA,col=rgb(0,0,0,0.05))
}


lines(spd.nz$grid$calBP, spd.nz$grid$PrDens,col=rgb(0,0,0,0.5),lwd=0.5)
lines(spd.nz$grid$calBP,rollmean(spd.nz$grid$PrDens,k=200,fill=NA),lwd=2)


plot(spd.roc,ylim=range(-0.05,0.05),xlim=c(0.5,5.5),xlab="Transitions",type="l",ylab="Geometric Growth Rate",axes=F,main="b")
points(spd.roc,pch=20)
abline(h=0,lty=2,col="red")
axis(side=1,at=1:5,labels=as.roman(1:5))
axis(side=2)





