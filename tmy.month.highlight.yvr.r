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

##--------------------------------------------------------------------------------
plot.file <- '/storage/data/projects/rci/weather_files/yvr.obs.jan.july.tmy.temperatures.png'

png(file=plot.file,width=3,height=3,units='in',res=600,pointsize=6,bg='white')

mon.tasmax <- tapply(tasmax.series,tas.mn.fac,mean,na.rm=T)
mon.tasmin <- tapply(tasmin.series,tas.mn.fac,mean,na.rm=T)

july.ix <- grep('*-07',names(mon.tasmax))
july.tasmax <- mon.tasmax[july.ix]

par(mfrow=c(2,1))
par(mar=c(2,3,2,1),mgp=c(2,1,0))
plot(x,july.tasmax,col='red',pch=18,xlab='',ylab='Max. Temperature (\u00B0C)',
       cex.axis=0.85,cex.lab=0.85,cex=0.85,cex.main=0.95,
       main='July Average Maximum Temperature at YVR',
       xlim=c(1935,2020),ylim=c(18,26),axes=F)
abline(h=seq(18,26,2),col='gray',lty=3,lwd=0.5)
axis(1,at=seq(1940,2020,20),label=seq(1940,2020,20),cex=0.95)
axis(2,at=seq(18,26,2),label=seq(18,26,2),cex=0.95,cex=0.95)
month <- which(names(july.tasmax) == '2002-07')
abline(h=july.tasmax[month],col='gray')
abline(h=mean(july.tasmax[35:64]))
abline(h=mean(july.tasmax[50:79]),lty=2)
points(x,july.tasmax,col='red',pch=18,cex=0.95)
#text(x=2020,y=23,'2002',cex=0.95)
#text(x=2020,y=21,'1971-\n2000',cex=0.95)
legend('bottomright',leg=c('2002','1971-2000','1986-2015'),
       col=c('gray','black','black'),lty=c(1,1,2))

box(which='plot')

##--------------------------------------------------------------------

jan.ix <- grep('*-01',names(mon.tasmin))
jan.tasmin <- mon.tasmin[jan.ix]

plot(x,jan.tasmin,col='blue',pch=18,xlab='Year',ylab='Min. Temperature (\u00B0C)',
       cex.axis=0.85,cex.lab=0.85,cex=0.85,cex.main=0.95,
       xlim=c(1935,2020),ylim=c(-10,6),
       main='January Average Minimum Temperature at YVR',axes=F)
abline(h=seq(-12,6,4),col='gray',lty=3,lwd=0.5)
##abline(h=0,col='gray',lwd=0.5)
axis(1,at=seq(1940,2020,20),label=seq(1940,2020,20),cex=0.95)
axis(2,at=seq(-12,6,4),label=seq(-12,6,4),cex=0.95,cex=0.95)

month <- which(names(jan.tasmin) == '2011-01')
abline(h=jan.tasmin[month],col='gray')
abline(h=mean(jan.tasmin[35:64]))
abline(h=mean(jan.tasmin[50:79]),lty=2)
##text(x=2015,y=-9,'2011',cex=0.95)

legend('bottomright',leg=c('2002','1971-2000','1986-2015'),
       col=c('gray','black','black'),lty=c(1,1,2))

points(x,jan.tasmin,col='blue',pch=18,cex=0.85)
box(which='plot')

dev.off()
