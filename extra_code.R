spd2gg <- function(spd,breaks)
{
require(rcarbon)	
nBreaks = length(breaks)-1
timeRange = eval(parse(text=spd$metadata[2]))
timeSequence = timeRange[1]:timeRange[2]
obs=numeric(length=nBreaks)    
for (x in 1:nBreaks)
{
	index=which(timeSequence<=breaks[x]&timeSequence>breaks[x+1])
	obs[x]=sum(spd$grid[index,2])
}

res=numeric(length=nBreaks-1)
for (i in 1:(nBreaks-1))
{
d=abs(breaks[i+1]-breaks[i]) 	
res[i]=(obs[i+1]/obs[i])^(1/d)-1
}

return(list(sumblock=obs,geomg=res,breaks=breaks))
}




plot.spd2gg<- function(x,...)
{
require(rcarbon)
breaks=x$breaks
obs=x$sumblock
res-x$geomg
par(mar=c(4,4,4,4))
nn = paste(breaks[-length(breaks)],breaks[-1],sep="-")
barplot(x$sumblock,names.arg=nn,ylab="Summed Probability",,space=0,col="bisque3",border=NA,...)
par(new=T)
xx = 1:c(length(nn)-1)
plot(0,0,xlim=c(0,length(nn)),ylim=range(res),axes=FALSE,xlab="",ylab="")
lines(xx,res,lwd=2,col="darkgreen")
points(xx,res,pch=20,col="darkgreen")
axis(4,col="darkgreen", col.axis="darkgreen")
mtext(side=4,"geometric growth rate",col="darkgreen",line=2)
abline(h=0,lty=2,col="blue")
}

searchSSPD(x)
{
require(rcarbon)	




}

