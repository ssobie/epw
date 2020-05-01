##Script to evaluate YVR temperature stats to deal with
##the EPW file year selection issues

##-------------------------------------------------------------------
##Read EC downloaded data

read_ec_obs_file <- function(var.name,interval,dir) {

   prefix <- 'eng-daily-'
   ec.file <- paste0(dir,prefix,interval,'.csv')
   ec.raw <- read.csv(ec.file,skip=24,header=T)
  
   ix <- grep(var.name,names(ec.raw))
  
   data <- ec.raw[,ix]
   dates <- as.Date(paste(ec.raw[,2],sprintf('%02d',ec.raw[,3]),sprintf('%02d',ec.raw[,4]),sep='-'))

   rv <- list(data=data,dates=dates)
   return(rv)
}

all_years_from_ec_obs <- function(var.title,dir) {
    
    ec.2011 <- read_ec_obs_file(var.title,'01012011-12312011',yvr.dir)
    ec.2012 <- read_ec_obs_file(var.title,'01012012-12312012',yvr.dir)
    ec.2013_1 <- read_ec_obs_file(var.title,'01012013-06122013',yvr.dir)
    ec.2013_2 <- read_ec_obs_file(var.title,'06132013-12312013',yvr.dir)
    ec.2014 <- read_ec_obs_file(var.title,'01012014-12312014',yvr.dir)
    ec.2015 <- read_ec_obs_file(var.title,'01012015-12312015',yvr.dir)
    ec.2016 <- read_ec_obs_file(var.title,'01012016-12312016',yvr.dir)
    ec.2017 <- read_ec_obs_file(var.title,'01012017-12312017',yvr.dir)
    ec.2018 <- read_ec_obs_file(var.title,'01012018-12312018',yvr.dir)

    ec.data <- c(ec.2011$data,ec.2012$data,ec.2013_1$data,ec.2013_2$data,
                 ec.2014$data,ec.2015$data,ec.2016$data,ec.2017$data,
                 ec.2018$data)
    ec.dates <- c(ec.2011$dates,ec.2012$dates,ec.2013_1$dates,ec.2013_2$dates,
                 ec.2014$dates,ec.2015$dates,ec.2016$dates,ec.2017$dates,
                 ec.2018$dates)
    rv <- list(data=ec.data,dates=ec.dates)
    return(rv)
}

##-------------------------------------------------------------------
##Read CDCD file

read_cdcd_file <- function(var.name,dir) {
   file <- paste0(dir,'1108447_',var.name,'_EC.csv')
   cdcd <- read.csv(file,header=T)
   years <- cdcd[,1]
   jdays <- cdcd[,2]
   dates <- as.Date(paste(years,jdays,sep='-'),'%Y-%j')
   
   rv <- list(data=cdcd[,3],dates=dates)
   return(rv)
}

##-------------------------------------------------------------------
##Create time series of variable


##-------------------------------------------------------------------
##

yvr.dir <- '/storage/data/projects/rci/weather_files/YVR/YVR_EC_obs/'

tasmax.cdcd <- read_cdcd_file('MAX_TEMP',yvr.dir)
tasmin.cdcd <- read_cdcd_file('MIN_TEMP',yvr.dir)
tas.cdcd <- read_cdcd_file('MEAN_TEMP',yvr.dir)

tasmax.ec <- all_years_from_ec_obs('Max.Temp...C.',yvr.dir)
tasmin.ec <- all_years_from_ec_obs('Min.Temp...C.',yvr.dir)
tas.ec <- all_years_from_ec_obs('Mean.Temp...C.',yvr.dir)

tasmax.series <- c(tasmax.cdcd$data,tasmax.ec$data)
tasmax.dates <- c(tasmax.cdcd$dates,tasmax.ec$dates)

tasmin.series <- c(tasmin.cdcd$data,tasmin.ec$data)
tasmin.dates <- c(tasmin.cdcd$dates,tasmin.ec$dates)

tas.series <- c(tas.cdcd$data,tas.ec$data)
tas.dates <- c(tas.cdcd$dates,tas.ec$dates)

tas.yr.fac <- as.factor(format(tas.dates,'%Y'))
tas.mn.fac <- as.factor(format(tas.dates,'%Y-%m'))

txtn.series <- (tasmax.series+tasmin.series)/2

x <- 1937:2018

tmy.years <- c(2011,2000,2000,2013,2003,2005,2002,2011,2003,2009,2007,2011)

##--------------------------------------------------------------------------------

for (mn in 1:12) {
   print(month.abb[mn])
   plot.file <- paste0('/storage/data/projects/rci/weather_files/YVR/',
                       'yvr.obs.',tolower(month.abb[mn]),'.tmy.temperatures.png')

   png(file=plot.file,width=5,height=5,units='in',res=600,pointsize=6,bg='white')

   mon.tasmax <- tapply(tasmax.series,tas.mn.fac,mean,na.rm=T)
   mon.tasmin <- tapply(tasmin.series,tas.mn.fac,mean,na.rm=T)
   mon.tas <- tapply(tas.series,tas.mn.fac,mean,na.rm=T)

   mon.ix <- grep(paste0('*-',sprintf('%02d',mn)),names(mon.tasmax))
   mon.tasmax <- mon.tasmax[mon.ix]

   par(mfrow=c(3,1))
   par(mar=c(4,6,3,16),mgp=c(4,1.5,0))
   plot(x,mon.tasmax,col='red',pch=18,xlab='',ylab='Max. Temperature (\u00B0C)',
       cex.axis=2,cex.lab=2,cex=2,cex.main=2,
       main=paste0(month.abb[mn],' Average Maximum Temperature at YVR'),
       xlim=c(1935,2020),ylim=range(mon.tasmax),axes=F)
   abline(h=seq(-40,40,1),col='gray',lty=3,lwd=0.5)
   axis(1,at=seq(1940,2020,20),label=seq(1940,2020,20),cex.axis=2)
   axis(2,at=seq(-30,40,1),label=seq(-30,40,1),cex.axis=2)
   month <- which(names(mon.tasmax) == paste0(tmy.years[mn],'-',sprintf('%02d',mn)))
   abline(h=mon.tasmax[month],col='gray',lwd=1.5)
   abline(h=mean(mon.tasmax[35:64]),lwd=1.5)
   abline(h=mean(mon.tasmax[50:79]),lty=2,lwd=1.5)
   abline(h=mean(mon.tasmax[62:78]),lty=3,lwd=1.5)
   points(x,mon.tasmax,col='red',pch=18,cex=1)
   #text(x=2020,y=23,'2002',cex=0.95)
   #text(x=2020,y=21,'1971-\n2000',cex=0.95)
   par(xpd=NA)
   legend('topright',inset=c(-0.28,0),leg=c(tmy.years[mn],'1971-2000','1986-2015','1998-2014'),
           col=c('gray','black','black','black'),lty=c(1,1,2,3),cex=2,lwd=2)
   par(xpd=FALSE)
   box(which='plot')

##--------------------------------------------------------------------
   mon.ix <- grep(paste0('*-',sprintf('%02d',mn)),names(mon.tas))
   mon.tas <- mon.tas[mon.ix]

   plot(x,mon.tas,col='red',pch=18,xlab='',ylab='Avg. Temperature (\u00B0C)',
       cex.axis=2,cex.lab=2,cex=2,cex.main=2,
       main=paste0(month.abb[mn],' Daily Average Temperature at YVR'),
       xlim=c(1935,2020),ylim=range(mon.tas),axes=F)
   abline(h=seq(-40,40,1),col='gray',lty=3,lwd=0.5)
   axis(1,at=seq(1940,2020,20),label=seq(1940,2020,20),cex.axis=2)
   axis(2,at=seq(-30,40,1),label=seq(-30,40,1),cex=0.95,cex.axis=2)
   month <- which(names(mon.tas) == paste0(tmy.years[mn],'-',sprintf('%02d',mn)))
   abline(h=mon.tas[month],col='gray',lwd=1.5)
   abline(h=mean(mon.tas[35:64]),lwd=1.5)
   abline(h=mean(mon.tas[50:79]),lty=2,lwd=1.5)
   abline(h=mean(mon.tas[62:78]),lty=3,lwd=1.5)
   points(x,mon.tas,col='red',pch=18,cex=1) 
   #text(x=2020,y=23,'2002',cex=0.95)
   #text(x=2020,y=21,'1971-\n2000',cex=0.95)
   par(xpd=NA)
   legend('topright',inset=c(-0.28,0),leg=c(tmy.years[mn],'1971-2000','1986-2015','1998-2014'),
           col=c('gray','black','black','black'),lty=c(1,1,2,3),cex=2,lwd=2)
   par(xpd=FALSE)

   box(which='plot')

##--------------------------------------------------------------------

   mon.ix <- grep(paste0('*-',sprintf('%02d',mn)),names(mon.tasmin))
   mon.tasmin <- mon.tasmin[mon.ix]

   plot(x,mon.tasmin,col='blue',pch=18,xlab='Year',ylab='Min. Temperature (\u00B0C)',
       cex.axis=2,cex.lab=2,cex=2,cex.main=2,
       xlim=c(1935,2020),ylim=range(mon.tasmin),
       main=paste0(month.abb[mn],' Average Minimum Temperature at YVR'),axes=F)
   abline(h=seq(-40,40,1),col='gray',lty=3,lwd=0.5)
   ##abline(h=0,col='gray',lwd=0.5)
   axis(1,at=seq(1940,2020,20),label=seq(1940,2020,20),cex.axis=2)
   axis(2,at=seq(-30,40,1),label=seq(-30,40,1),cex.axis=2)

   month <- which(names(mon.tasmin) == paste0(tmy.years[mn],'-',sprintf('%02d',mn)))
   abline(h=mon.tasmin[month],col='gray',lwd=1.5)
   abline(h=mean(mon.tasmin[35:64]),lwd=1.5)
   abline(h=mean(mon.tasmin[50:79]),lty=2,lwd=1.5)
   abline(h=mean(mon.tasmin[62:78]),lty=3,lwd=1.5)
   ##text(x=2015,y=-9,'2011',cex=0.95)
   par(xpd=NA)
   legend('topright',inset=c(-0.29,0),leg=c(tmy.years[mn],'1971-2000','1986-2015','1998-2014'),
          col=c('gray','black','black','black'),lty=c(1,1,2,3),cex=2,lwd=2)
   par(xpd=FALSE)
   points(x,mon.tasmin,col='blue',pch=18,cex=1)
   box(which='plot')

   dev.off()

}