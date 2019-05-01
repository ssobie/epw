##Script to plot the EPW file series

library(ncdf4)
library(PCICt)
library(zoo)
library(scales)
source('/storage/home/ssobie/code/repos/crcm5/read.write.epw.r',chdir=T)


##------------------------------------------------------------------------------
##Match for EPW fields

get_field_index <- function(var.name) {

   field.names <- c('year', 'month', 'day', 'hour', 'minute',
      'data_source_and_uncertainty_flags', 'dry_bulb_temperature',
      'dew_point_temperature', 'relative_humidity',
      'atmospheric_station_pressure', 'extraterrestrial_horizontal_radiation',
      'extraterrestrial_direct_normal_radition',
      'horizontal_infrared_radiation_intensity', 'global_horizontal_radiation',
      'direct_normal_radiation', 'diffuse_horizontal_radiation',
      'global_horizontal_illuminance', 'direct_normal_illuminance',
      'diffuse_horizontal_illuminance', 'zenith_luminance', 'wind_direction',
      'wind_speed', 'total_sky_cover', 'opaque_sky_cover', 'visibility',
      'ceiling_height', 'present_weather_observation', 'present_weather_codes',
      'precipitable_water', 'aerosol_optical_depth', 'snow_depth',
      'days_since_last_snowfall', 'albedo', 'liquid_precipitation_depth',
      'liquid_precipitation_quantity')
   ix <- grep(var.name,field.names)
}

##------------------------------------------------------------------------------

sub_by_time <- function(var.name,lonc,latc,interval,input.file,gcm,read.dir) {

  print(input.file)              
  nc <- nc_open(paste(read.dir,gcm,'/',input.file,sep=''))
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*86400
  years <- format(time.series,'%Y')
  yrs <- strsplit(interval,'-')[[1]]

  feb.flag <- grep('-02-29',time.series)

  new.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal='365_day')
                           
  new.time0 <- (as.PCICt(format(time.series[1],'%Y-%m-%d'),cal='365_day') - new.origin)/86400
  if (length(feb.flag)==0) {
     new.values <- seq(as.numeric(new.time0),by=1,length.out=length(time.values))
  } else {
     new.values <- seq(as.numeric(new.time0),by=1,length.out=length(time.values[-feb.flag]))
  }

  new.series <- new.origin + new.values*86400
  print(range(time.series))
  print(range(new.series))

  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  
  lon.ix <- which.min(abs(lonc-lon))
  lat.ix <- which.min(abs(latc-lat))

  data.raw <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  data <- data.raw
  if (length(feb.flag!=0)) {
    data <- data.raw[-feb.flag]
  }

  years <- format(new.series,'%Y')
  st <- head(grep(yrs[1],years),1)
  en <- tail(grep(yrs[2],years),1)

  nc_close(nc)
  rv <- list(data=data[st:en],time=new.series[st:en])
  return(rv)
}

##------------------------------------------------------------------------------

make_average_series <- function(series,time,method,rlen=5,agg.fxn) {

  factor <- switch(method,
                   daily='%m-%d',
                   monthly='%m',
                   roll='%m-%d')

  fac <- as.factor(format(time,factor))
  agg.data <- tapply(series,fac,agg.fxn)
  if (method=='roll') {
    agg.data <- rollmean(agg.data,rlen,fill='extend')
  }    
  return(agg.data)
}


##------------------------------------------------------------------------------
###Consider splitting this (and the other Belcher methods) into a separate script
###so that modifications can be run elsewhere and not be confused.

morph_dry_bulb_temp <- function(epw.present,lon,lat,gcm.list,gcm.dir,scenario,
                                method,rlen=NULL,agg.fxn=mean) {

   tas.ix <- get_field_index('dry_bulb_temperature')
   epw.tas <- epw.present$data[,tas.ix]
   dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))
   dy.dates <- as.Date(unique(format(dates,'%Y-%m-%d')))

   epw.agg.mean <- epw.daily.mean <- make_average_series(epw.tas,dates,method='daily',agg.fxn=mean)
   epw.agg.max <- epw.daily.max <-  make_average_series(epw.tas,dates,method='daily',agg.fxn=max)
   epw.agg.min <- epw.daily.min <-  make_average_series(epw.tas,dates,method='daily',agg.fxn=min)

   epw.day.anoms <- epw.tas*0
   ##Daily
   for (d in 1:365) {
      ix <- format(dates,'%j') == sprintf('%03d',d)
      epw.day.anoms[ix] <- epw.tas[ix] - epw.daily.mean[d]
   }
   epw.agg.anoms <- epw.day.anoms

   if (method == 'monthly') {
     dy.dates <- as.Date(unique(format(dates,'%Y-%m-%d')))
     epw.agg.mean <- make_average_series(epw.tas,dates,method='monthly',agg.fxn=mean)
     epw.agg.max <-  make_average_series(epw.tas,dates,method='monthly',agg.fxn=max)
     epw.agg.min <-  make_average_series(epw.tas,dates,method='monthly',agg.fxn=min)
     for (d in seq_along(epw.agg.mean)) {
        ix <- format(dates,'%m') == sprintf('%02d',d)
        epw.agg.anoms[ix] <- epw.tas[ix] - epw.agg.mean[d]
     }
   }
   if (method == 'roll') {
     epw.agg.mean <- make_average_series(epw.tas,dates,method='roll',rlen=rlen,agg.fxn=mean)
     epw.agg.max <-  make_average_series(epw.tas,dates,method='roll',rlen=rlen,agg.fxn=max)
     epw.agg.min <-  make_average_series(epw.tas,dates,method='roll',rlen=rlen,agg.fxn=min)
     for (d in 1:365) {
        ix <- format(dates,'%j') == sprintf('%03d',d)
        epw.agg.anoms[ix] <- epw.tas[ix] - epw.agg.mean[d]
     }     
   }

   tlen <- length(epw.agg.mean)
   
   deltas <- matrix(0,nrow=length(gcm.list),ncol=tlen)
   alphas <- matrix(0,nrow=length(gcm.list),ncol=tlen)

   morphed.tas <- matrix(0,nrow=length(gcm.list),ncol=length(epw.tas*0))

   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      scen.files <- list.files(path=paste0(gcm.dir,gcm),pattern=scenario)
      var.files <- scen.files[grep("tasmax",scen.files)]
      past.tx.file <- var.files[grep("1951-2000",var.files)]
      proj.tx.file <- var.files[grep("2001-2100",var.files)]

      past.tx <- sub_by_time(var.name='tasmax',lonc=lon,latc=lat,
                                  interval='1971-2000',
                                  input.file=past.tx.file,gcm=gcm,read.dir=gcm.dir)
      past.tx.agg <- make_average_series(past.tx$data,past.tx$time,
                                         method,rlen,agg.fxn)

      proj.tx <- sub_by_time(var.name='tasmax',lonc=lon,latc=lat,
                                  interval='2041-2070',
                                  input.file=proj.tx.file,gcm=gcm,read.dir=gcm.dir)
      proj.tx.agg <- make_average_series(proj.tx$data,proj.tx$time,
                                         method,rlen,agg.fxn)

      var.files <- scen.files[grep("tasmin",scen.files)]
      past.tn.file <- var.files[grep("1951-2000",var.files)]
      proj.tn.file <- var.files[grep("2001-2100",var.files)]

      past.tn <- sub_by_time(var.name='tasmin',lonc=lon,latc=lat,
                             interval='1971-2000',
                             input.file=past.tn.file,gcm=gcm,read.dir=gcm.dir)
      past.tn.agg <- make_average_series(past.tn$data,past.tn$time,
                                         method,rlen,agg.fxn)

      proj.tn <- sub_by_time(var.name='tasmin',lonc=lon,latc=lat,
                             interval='2041-2070',
                             input.file=proj.tn.file,gcm=gcm,read.dir=gcm.dir)
      proj.tn.agg <- make_average_series(proj.tn$data,proj.tn$time,
                                         method,rlen,agg.fxn)

      ##
      delta_tasmax <- proj.tx.agg - past.tx.agg
      delta_tasmin <- proj.tn.agg - past.tn.agg
      delta_tas <- (proj.tx.agg+proj.tn.agg)/2 - (past.tx.agg+past.tn.agg)/2
      deltas[g,] <- delta_tas
      past_tas <- (past.tx.agg+past.tn.agg)/2
      print('Delta')
      print(delta_tas)

      alpha <- (delta_tasmax - delta_tasmin) / (epw.agg.max - epw.agg.min)
      print('Alpha')
      print(alpha)        
      alphas[g,] <- alpha 


      if (method=='daily' | method =='roll') {
        for (d in 1:tlen) {
           ix <- format(dates,'%j') == sprintf('%03d',d)  
           morphed.tas[g,ix] <- epw.tas[ix] + delta_tas[d] + alpha[d]*epw.agg.anoms[ix]
           ##morphed.tas[g,ix] <- deltas[g,m] + alphas[g,m]*epw.agg.anoms[ix]
        }
      }

      if (method=='monthly') {
        for (m in 1:tlen) {
           ix <- format(dates,'%m') == sprintf('%02d',m)  
           morphed.tas[g,ix] <- epw.tas[ix] + delta_tas[m] + alpha[m]*epw.agg.anoms[ix]
           ##morphed.tas[g,ix] <- deltas[g,m] + alphas[g,m]*epw.agg.anoms[ix]
        }
      }

   }##GCM loop

  hour.dates <- strftime(paste('1999-',sprintf('%02d',epw.present$data[,2]),'-',sprintf('%02d',epw.present$data[,3]),' ', 
                         sprintf('%02d',epw.present$data[,4]),':00:00', sep=''),format='%Y-%m-%d %H:%M:%S')
browser()

   ens.deltas <- apply(deltas,2,mean)
   ens.alphas <- apply(alphas,2,mean)
   ens.morphed.tas <- apply(morphed.tas,2,mean)
   return(morphed.tas)
   ##epw.present$data[,tas.ix] <- round(ens.morphed.tas,1)
   ##return(epw.present)   
}

##------------------------------------------------------------------------------

##**************************************************************************************

scenario <- 'rcp85'

lon <- -122.36
lat <- 49.03

##'/storage/home/ssobie/code/repos/bc-projected-weather/bcweather/tests/data/' 

epw.dir <- '/storage/data/projects/rci/weather_files/wx_files/offsets/'
write.dir <- '/storage/data/projects/rci/weather_files/wx_files/morphed_files/'
present.epw.file <- 'CAN_BC_1st_and_Clark_offset_from_VANCOUVER-INTL-A_1108395_CWEC.epw'
future.epw.file <- 'MORPHED_ROLL21_TAS_CAN_BC_1st_and_Clark_offset_from_VANCOUVER-INTL-A_1108395_CWEC.epw'
plot.dir <- '/storage/data/projects/rci/weather_files/'


epw.present <- read.epw.file(epw.dir,present.epw.file)

##Create one year of daily dates

gcm.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

gcm.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')

##epw.morphed.tas <- morph_dry_bulb_temp(epw.present,lon,lat,gcm.list,gcm.dir,scenario,method='roll',rlen=21)
dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))
hour.dates <- strftime(paste('1999-',sprintf('%02d',epw.present$data[,2]),'-',sprintf('%02d',epw.present$data[,3]),' ', 
                       sprintf('%02d',epw.present$data[,4]),':00:00', sep=''),format='%Y-%m-%d %H:%M:%S')
dy.dates <- as.Date(unique(format(dates,'%Y-%m-%d')))

load('/storage/data/projects/rci/weather_files/morphed.tas.RData')
load('/storage/data/projects/rci/weather_files/alphas.RData')
load('/storage/data/projects/rci/weather_files/deltas.RData')
   epw.tas <- epw.present$data[,7]
   epw.agg.mean <- epw.daily.mean <- make_average_series(epw.tas,dates,method='daily',agg.fxn=mean)
   epw.day.anoms <- epw.tas*0
   ##Daily
   for (d in 1:365) {
      ix <- format(dates,'%j') == sprintf('%03d',d)
      epw.day.anoms[ix] <- epw.tas[ix] - epw.daily.mean[d]
   }
  ##--------------
  mn <- '-06-|-07-|-08-'
  hr.dates <- hour.dates[grep(mn,hour.dates)]
  hr.dy <- dy.dates[grep(mn,dy.dates)]
  hr.epw <- epw.present$data[grep(mn,dates),7]
  hr.morph <- morphed.tas[,grep(mn,dates)]
  hr.anoms <- epw.day.anoms[grep(mn,dates)]
  
  mn.max <- tapply(hr.epw,as.factor(format(as.Date(hr.dates),'%m-%d')),max)
  mn.min <- tapply(hr.epw,as.factor(format(as.Date(hr.dates),'%m-%d')),min)
  mn.mean <- tapply(hr.epw,as.factor(format(as.Date(hr.dates),'%m-%d')),mean)

  plot.file <- paste0(plot.dir,'agu.examples.abbotsford.jul1.diurnal.tas.png')
  png(plot.file,width=1000,height=400)
  par(mar=c(4.5,4.5,4,2))
  days <- format(as.Date(hr.dy),'%b-%d')
  plot(1:length(days),mn.mean,type='l',lwd=4,xlab='Hour',ylab='Temperature (\u00B0C)',xaxs='i',
       cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Summer Temperature',ylim=c(8,35),col='white',axes=F)
  dx <- c(1,10,20,31,40,50,60,71,81,92)
  axis(1,at=dx,labels=days[dx],cex.axis=1.5,cex=1.5)
  axis(2,at=seq(-40,40,2),labels=seq(-40,40,2),cex.axis=1.5,cex=1.5)
  abline(h=seq(-40,40,2),lty=2,col='gray',lwd=1)

  polygon(c(1:length(days),rev(1:length(days))),
          c(mn.min,rev(mn.max)),col=alpha('gray',0.4),border=alpha('gray',0.4))
  lines(1:length(days),mn.mean,lwd=5)
  
##  morph.dy.mn <- t(apply(hr.morph,1,function(x,y){tapply(x,y,mean)},as.factor(format(as.Date(hr.dates),'%m-%d'))))
##  morph.dy.max <- t(apply(hr.morph,1,function(x,y){tapply(x,y,max)},as.factor(format(as.Date(hr.dates),'%m-%d'))))
##  morph.dy.min <- t(apply(hr.morph,1,function(x,y){tapply(x,y,min)},as.factor(format(as.Date(hr.dates),'%m-%d'))))
##  m50 <- apply(morph.dy.mn,2,mean)
##  m10 <- apply(morph.dy.min,2,quantile,0.1)
##  m90 <- apply(morph.dy.max,2,quantile,0.9)
  m50 <- tapply(hr.morph[3,],as.factor(format(as.Date(hr.dates),'%m-%d')),mean)
  m90 <- tapply(hr.morph[3,],as.factor(format(as.Date(hr.dates),'%m-%d')),max)
  m10 <- tapply(hr.morph[3,],as.factor(format(as.Date(hr.dates),'%m-%d')),min)
  polygon(c(1:length(days),rev(1:length(days))),
          c(m10,rev(m90)),col=alpha('red',0.2),border=alpha('red',0.2))
  lines(1:length(days),m50,lwd=5,col='red')
    
  for (i in 1:12) {
  ##  lines(1:length(days),morph.days[i,],lwd=2,col='green')
  }
  box(which='plot')

##for (i in 1:12) {
##  lines(alphas[i,as.numeric(mn)]*hr.anoms[1:24],lwd=2)
##}
  dev.off()

  ##--------------
browser()  

morphed.tas <- epw.morphed.tas$data[,7]
##write.epw.file(epw.morphed.tas$data,epw.morphed.tas$header,write.dir,future.epw.file)

##browser()


daily.morphed.tas <- t(apply(morphed.tas,1,function(x,dates){tapply(x,as.factor(format(dates,'%m-%d')),mean)},dates))

hour.dates <- strftime(paste('1999-',sprintf('%02d',epw.present.data[,2]),'-',sprintf('%02d',epw.present.data[,3]),' ', 
                          sprintf('%02d',epw.present.data[,4]),':00:00', sep=''),format='%Y-%m-%d %H:%M:%S')

mn <- '07'
hr.dates <- hour.dates[grep(paste0('-',mn,'-'),hour.dates)]
hr.epw <- epw.present.data[grep(paste0('-',mn,'-'),dates),7]
hr.morph <- morphed.tas[,grep(paste0('-',mn,'-'),dates)]
hr.anoms <- epw.day.anoms[grep(paste0('-',mn,'-'),dates)]

if (1==0) {
plot.file <- paste0(plot.dir,'abbotsford.jul.hourly.tas.png')

png(plot.file,width=1200,height=400)
par(mar=c(4.5,4.5,4,2))
plot(1:length(hr.morph),hr.morph,type='l',lwd=4,xlab='Day',ylab='Temperature (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford July Hourly Temperature',ylim=c(10,35),col='white',axes=F)
axis(1,at=seq(0,744,24)[-32],labels=format(dy.dates[grep(paste0('-',mn,'-'),dy.dates)],'%d'),cex.axis=1.5,cex=1.5)
axis(2,at=seq(-15,35,5),labels=seq(-15,35,5),cex.axis=1.5,cex=1.5)
abline(h=seq(-15,35,5),lty=2,col='gray',lwd=2)
abline(h=0,col='gray')
lines(1:length(hr.morph),hr.epw,lwd=4,col='orange')
lines(1:length(hr.morph),hr.morph,lwd=4,col='red')

legend('topleft',leg=c('Future','Present'),col=c('red','orange'),pch=15,cex=1.5)

box(which='plot')
dev.off()
}

browser()

if (1==0) {
plot.file <- paste0(plot.dir,'agu.examples.abbotsford.jul1.diurnal.tas.png')
png(plot.file,width=1200,height=800)
par(mar=c(4.5,4.5,4,2))
plot(1:24,hr.epw[1:24],type='l',lwd=4,xlab='Day',ylab='Temperature (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Jul 1st Hourly Temperature',ylim=c(10,26),col='white',axes=F)
axis(1,at=1:24,labels=1:24,cex.axis=1.5,cex=1.5)
axis(2,at=seq(-40,40,2),labels=seq(-40,40,2),cex.axis=1.5,cex=1.5)
abline(h=seq(-40,40,2),lty=2,col='gray',lwd=1)
abline(h=0,col='gray')
lines(1:24,hr.epw[1:24],lwd=4)

for (i in 1:12) {
  lines(hr.morph[i,1:24],lwd=1,col='gray')
}
lines(apply(hr.morph[,1:24],2,mean),lwd=3,col='orange')
box(which='plot')
dev.off()

}


plot.file <- paste0(plot.dir,'examples.abbotsford.jul1.diurnal.tas.png')
png(plot.file,width=1200,height=800)
par(mar=c(4.5,4.5,4,2))
par(mfrow=c(3,1))

plot(1:24,hr.epw[1:24],type='l',lwd=4,xlab='Hour',ylab='Temperature (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Jul 1 Hourly Temperature',ylim=c(10,26),col='white',axes=F)
axis(1,at=1:24,labels=1:24,cex.axis=1.5,cex=1.5)
axis(2,at=seq(-40,40,2),labels=seq(-40,40,2),cex.axis=1.5,cex=1.5)
abline(h=seq(-40,40,2),lty=2,col='gray',lwd=2)
lines(1:24,hr.epw[1:24],lwd=4)
box(which='plot')

plot(1:24,alphas[2,as.numeric(mn)]*hr.anoms[1:24],type='l',lwd=4,xlab='Hour',ylab='Temperature (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Jul Alphas',ylim=c(-0.8,0.6),col='white',axes=F)
axis(1,at=1:24,labels=1:24,cex.axis=1.5,cex=1.5)
axis(2,at=seq(-1,1,0.25),labels=seq(-1,1,0.25),cex.axis=1.5,cex=1.5)
abline(h=seq(-1,1,0.25),lty=2,col='gray',lwd=2)
for (i in 1:12) {
  lines(alphas[i,as.numeric(mn)]*hr.anoms[1:24],lwd=2)
}
box(which='plot')

plot(1:24,hr.morph[2,1:24],type='l',lwd=4,xlab='Hour',ylab='Temperature (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Jul 1 Hourly Morphed',ylim=c(10,26),col='white',axes=F)
axis(1,at=1:24,labels=1:24,cex.axis=1.5,cex=1.5)
axis(2,at=seq(-40,40,2),labels=seq(-40,40,2),cex.axis=1.5,cex=1.5)
abline(h=seq(-40,40,2),lty=2,col='gray',lwd=2)
for (i in 1:12) {
  lines(hr.morph[i,1:24],lwd=2)
}

box(which='plot')
dev.off()





browser()

plot.file <- paste0(plot.dir,'abbotsford.tas.anomalies.smoothed.11.png')


png(plot.file,width=1200,height=400)
par(mar=c(4.5,4.5,4,2))

plot(1:365,apply(daily.morphed.tas,2,mean),type='l',lwd=4,xlab='Julian Day',ylab='Daily Mean Temperature Anomalies (\u00B0C)',
     cex.main=2,cex.lab=1.5,cex.axis=1.5,main='Abbotsford Temperature Anomalies',ylim=c(0,6.2),col='white')
abline(h=0,col='gray')
lines(6:360,rollmean(apply(daily.morphed.tas,2,mean),11),lwd=4,col='orange')
lines(6:360,rollmean(apply(daily.morphed.tas,2,quantile,0.1),11),lwd=4,col='gold')
lines(6:360,rollmean(apply(daily.morphed.tas,2,quantile,0.9),11),lwd=4,col='red')

box(which='plot')
dev.off()

