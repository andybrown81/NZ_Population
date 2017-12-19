library(RCurl)
require(RCurl)

nzdates <-read.csv(text=getURL("https://raw.githubusercontent.com/andybrown81/NZ_Population/master/SPDdata.csv"))

nzdates <- subset(nzdates,C14Age>0)
nzdates <- subset(nzdates, C14Age<800)

northern_dates <- subset(nzdates, Region=="Northern")
caldates_northern <- calibrate(northern_dates$C14Age,northern_dates$C14SD,calCurves="shcal13", normalised=FALSE, verbose=FALSE)
bins_N <- binPrep(sites=northern_dates$SiteID,ages=northern_dates$C14Age,h=100)
spd.n <- spd(x=caldates_northern,bins=bins_N,timeRange=c(800,0),runm=50)

### Logistic Model (Toy 50 sims)
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
logMod_northern <- modelTest(caldates_northern,bins=bins_N,errors=northern_dates$C14SD,predgrid=logisticFit_N,runm=50,timeRange=c(800,0),model="custom",calCurves="shcal13",nsim=50,ncores=2, raw=TRUE, datenormalised=FALSE)

### not recognising logMod_northern as a true SpdModelTest class.
P2P_Test <- p2pTest(x=logMod_northern,p1=400,p2=100)


