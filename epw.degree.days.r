##Calculate degree days from the existing and morphed weather files
library(scales)
source('/storage/home/ssobie/code/repos/crcm5/read.write.epw.r',chdir=T)


##--------------------------------
##Degree Day Values
dd <- function(tmean,tbase) {
  g <- tmean - tbase
  days <- sum(g[g>0], na.rm=T)
  return(round(days))
}

gdd<-function(data,fac){tapply(data,fac, dd, tbase=5)}   ##Growing degree days
cdd<-function(data,fac){tapply(data,fac, dd, tbase=18)}  ##Cooling degree days
hdd<-function(data,fac){tapply(-data,fac,dd, tbase=-18)} ##Heating degree days
fdd<-function(data,fac){tapply(-data,fac,dd, tbase=0)} ##Freezing degree days
ffd<-function(data,fac){tapply(data,fac,fd)} ##Frost Free days
s30<-function(data,fac){tapply(data,fac,s3)} ##Days over 30 degC

find.epw.dd <- function(epw.files,epw.dir) {

  dd.matrix <- matrix(NA,nrow=length(epw.files),ncol=3)
  colnames(dd.matrix) <- c('Name','CDD','HDD')

  for (i in seq_along(epw.files)) { 
    pef.split <- strsplit(epw.files[i],'_')[[1]]
    if (grepl('MORPHED',epw.files[i])) {
      print(pef.split[4])
      dd.matrix[i,1] <- pef.split[4]
    } else {
      print(pef.split[3])
      dd.matrix[i,1] <- pef.split[3]
    }
    epw.present <- read.epw.file(epw.dir,epw.files[i])
    dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))
    day.fac <- as.factor(format(dates,'%m-%d'))
    epw.tas <- epw.present$data[,7]
    tas.day <- tapply(epw.tas,day.fac,mean)
    dd.matrix[i,2] <- dd(tas.day,18)
    dd.matrix[i,3] <- dd(-tas.day,-18)
  }
  return(dd.matrix)
}

make.chull <- function(x,y) {
  test <- chull(x,y)
  hpts <- c(test,test[1])
  xp <- x[hpts]
  yp <- y[hpts]
  rv <- list(xp=xp,yp=yp)
  return(rv)
}

epw.dir <- '/storage/data/projects/rci/weather_files/wx_files/'
morphed.dir <- '/storage/data/projects/rci/weather_files/wx_files/morphed_files/tas_only/'
if (1==1) {
print('Past')
past.epw.files <- list.files(path=epw.dir,pattern='_CWEC.epw')
dd.past <- find.epw.dd(past.epw.files,epw.dir)

print('2020s')
epw.files.2020s <- list.files(path=morphed.dir,pattern='_2011-2040_CWEC.epw')
dd.2020s <- find.epw.dd(epw.files.2020s,morphed.dir)

print('2050s')
epw.files.2050s <- list.files(path=morphed.dir,pattern='_2041-2070_CWEC.epw')
dd.2050s <- find.epw.dd(epw.files.2050s,morphed.dir)

print('2080s')
epw.files.2080s <- list.files(path=morphed.dir,pattern='_2071-2100_CWEC.epw')
dd.2080s <- find.epw.dd(epw.files.2080s,morphed.dir)

dd.ix <- dd.past[,1] %in% dd.2020s[,1]
dd.past.match <- dd.past[dd.ix,]
}


png('/storage/data/projects/rci/weather_files/epw.cdd.hdd.png',width=1000,height=600)
par(mar=c(5.1,5.1,4.1,2.1))
plot(as.numeric(dd.past.match[,2]),as.numeric(dd.past.match[,3]),
        xlim=c(-10,1100),ylim=c(1000,7000),pch=18,yaxs='i',xaxs='i',
        cex=1.5,cex.axis=1.75,cex.lab=1.75,cex.main=2,
        xlab='Cooling Degree Days',ylab='Heating Degree Days',
        main='BC Weather Files Heating and Cooling Degree Days')
abline(h=seq(0,10000,1000),v=seq(0,2000,100),col='lightgray',lty=2)

ch.2080s <- make.chull(as.numeric(dd.2080s[,2]),as.numeric(dd.2080s[,3]))
polygon(ch.2080s$xp,ch.2080s$yp,col=alpha('red',0.3),border=alpha('red',0.3))

ch.2050s <- make.chull(as.numeric(dd.2050s[,2]),as.numeric(dd.2050s[,3]))
polygon(ch.2050s$xp,ch.2050s$yp,col=alpha('orange',0.3),border=alpha('orange',0.3))

ch.2020s <- make.chull(as.numeric(dd.2020s[,2]),as.numeric(dd.2020s[,3]))
polygon(ch.2020s$xp,ch.2020s$yp,col=alpha('gold',0.3),border=alpha('gold',0.3))

ch.past <- make.chull(as.numeric(dd.past.match[,2]),as.numeric(dd.past.match[,3]))
polygon(ch.past$xp,ch.past$yp,col=alpha('gray',0.3),border=alpha('gray',0.3))

points(as.numeric(dd.past.match[,2]),as.numeric(dd.past.match[,3]),col='black',pch=18,cex=1.5)
points(as.numeric(dd.2020s[,2]),as.numeric(dd.2020s[,3]),col='gold',pch=18,cex=1.5)
points(as.numeric(dd.2050s[,2]),as.numeric(dd.2050s[,3]),col='orange',pch=18,cex=1.5)
points(as.numeric(dd.2080s[,2]),as.numeric(dd.2080s[,3]),col='red',pch=18,cex=1.5)
box(which='plot')
legend('topright',leg=c('Past','2020s','2050s','2080s'),col=c('Black','Gold','Orange','Red'),pch=18,cex=1.75)
dev.off()