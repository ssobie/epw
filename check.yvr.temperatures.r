##Script to evaluate YVR temperature stats to deal with
##the EPW file year selection issues


##Degree Day Values
dd <- function(tmean,tbase) {
  g <- tmean - tbase
  days <- sum(g[g>0], na.rm=T)
  return(round(days))
}
cdd<-function(data,fac){tapply(data,fac, dd, tbase=18)}  ##Cooling degree days
hdd<-function(data,fac){tapply(-data,fac,dd, tbase=-18)} ##Heating degree days



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
cooling.deg.days <- cdd(tas.series,tas.yr.fac)

txtn.series <- (tasmax.series+tasmin.series)/2
cooling.txtn <- cdd(txtn.series,tas.yr.fac)

x <- 1937:2018

##--------------------------------------------------------------------------------
##Cooling DD plot

cdd.file <- '/storage/data/projects/rci/weather_files/yvr.obs.cdd.png'
png(file=cdd.file,width=3,height=1.5,units='in',res=600,pointsize=6,bg='white')
par(mar=c(4,4,2,1))
plot(x,cooling.txtn,col='red',pch=18,xlab='Year',ylab='Cooling Degree Days',
       cex.axis=0.85,cex.lab=0.85,cex=0.85,main='Cooling Degree Days at YVR')

abline(v=2004,lwd=0.5,col='gray')
text(x=2002,y=1,'2004',cex=0.5)

abline(v=2009,lwd=0.5,col='gray')
text(x=2007,y=1,'2009',cex=0.5)
abline(v=2015,lwd=0.5,col='gray')
text(x=2013,y=1,'2015',cex=0.5)
abline(v=2017,lwd=0.5,col='gray')
text(x=2019,y=1,'2017',cex=0.5)

points(x,cooling.txtn,col='red',pch=18,cex=0.85)
dev.off()

##--------------------------------------------------------------------------------
ann.tasmin <- tapply(tasmin.series,tas.yr.fac,min,na.rm=T)
ann.tasmax <- tapply(tasmax.series,tas.yr.fac,max,na.rm=T)

ann.tx.file <- '/storage/data/projects/rci/weather_files/yvr.obs.ann.tasmax.png'
png(file=ann.tx.file,width=3,height=1.5,units='in',res=600,pointsize=6,bg='white')
par(mar=c(4,4,2,1))
plot(x,ann.tasmax,col='red',pch=18,xlab='Year',ylab='Annual Hottest Temperature',
       cex.axis=0.85,cex.lab=0.85,cex=0.85,main='Annual Hottest Temperature at YVR')
abline(v=2013,lwd=0.5,col='gray')
text(x=2011,y=24.8,'2013',cex=0.5)
abline(v=2009,lwd=0.5,col='gray')
text(x=2007,y=24.3,'2009',cex=0.5)
abline(v=2015,lwd=0.5,col='gray')
text(x=2013,y=24.3,'2015',cex=0.5)
abline(v=2017,lwd=0.5,col='gray')
text(x=2019,y=24.3,'2017',cex=0.5)

points(x,ann.tasmax,col='red',pch=18,cex=0.85)
dev.off()

##--------------------------------------------------------------------------------
su30 <- tapply(tasmax.series,tas.yr.fac,function(x){sum(x>30,na.rm=T)})
su25 <- tapply(tasmax.series,tas.yr.fac,function(x){sum(x>25,na.rm=T)})

su25.file <- '/storage/data/projects/rci/weather_files/yvr.obs.ann.su25.png'
png(file=su25.file,width=3,height=1.5,units='in',res=600,pointsize=6,bg='white')
par(mar=c(4,4,2,1))
plot(x,su25,col='red',pch=18,xlab='Year',ylab='No. Days above 25 degC',
       cex.axis=0.85,cex.lab=0.85,cex=0.85,main='No. of Days above 25 degC at YVR')
abline(v=2018,lwd=0.5,col='gray')
text(x=2019,y=1.5,'2018',cex=0.5)
abline(v=2009,lwd=0.5,col='gray')
text(x=2007,y=0,'2009',cex=0.5)
abline(v=2015,lwd=0.5,col='gray')
text(x=2013,y=0,'2015',cex=0.5)
abline(v=2017,lwd=0.5,col='gray')
text(x=2017,y=0,'2017',cex=0.5)

points(x,su25,col='red',pch=18,cex=0.85)
dev.off()


##---------------------------------------------------------------
mon.tasmax <- tapply(tasmax.series,tas.mn.fac,mean,na.rm=T)
mon.tasmin <- tapply(tasmin.series,tas.mn.fac,mean,na.rm=T)

july.ix <- grep('*-07',names(mon.tasmax))
july.tasmax <- mon.tasmax[july.ix]

july.tx.file <- '/storage/data/projects/rci/weather_files/yvr.obs.july.tasmax.png'
png(file=july.tx.file,width=3,height=1.5,units='in',res=600,pointsize=6,bg='white')
par(mar=c(4,4,2,1))
plot(x,july.tasmax,col='red',pch=18,xlab='Year',ylab='July Average Max Temperature (C)',
       cex.axis=0.85,cex.lab=0.85,cex=0.85,main='July Average Maximum Temperature at YVR')
abline(v=2018,lwd=0.5,col='gray')
text(x=2019,y=20.5,'2018',cex=0.5)
abline(v=2009,lwd=0.5,col='gray')
text(x=2007,y=20,'2009',cex=0.5)
abline(v=2015,lwd=0.5,col='gray')
text(x=2013,y=20,'2015',cex=0.5)
abline(v=2017,lwd=0.5,col='gray')
text(x=2017,y=20,'2017',cex=0.5)

points(x,july.tasmax,col='red',pch=18,cex=0.85)
dev.off()


##--------------------------------------------------------------------
july.tasmin <- mon.tasmin[july.ix]

july.tn.file <- '/storage/data/projects/rci/weather_files/yvr.obs.july.tasmin.png'
png(file=july.tn.file,width=3,height=1.5,units='in',res=600,pointsize=6,bg='white')
par(mar=c(4,4,2,1))
plot(x,july.tasmin,col='red',pch=18,xlab='Year',ylab='July Average Min Temperature (C)',
       cex.axis=0.85,cex.lab=0.85,cex=0.85,main='July Average Minimum Temperature at YVR')


abline(v=2009,lwd=0.5,col='gray')
text(x=2007,y=11.5,'2009',cex=0.5)
abline(v=2015,lwd=0.5,col='gray')
text(x=2013,y=11.5,'2015',cex=0.5)
abline(v=2017,lwd=0.5,col='gray')
text(x=2017,y=11.5,'2017',cex=0.5)

points(x,july.tasmin,col='red',pch=18,cex=0.85)
dev.off()
