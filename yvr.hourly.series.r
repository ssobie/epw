##Script to generate an hourly temperature series by combining BCCAQ2-PRISM 
##with the hourly temperatures from the YVR EPW File

library(pracma)

library(ncdf4)
source('/storage/home/ssobie/code/repos/assessments/single.site.bccaq.to.prism.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/read.write.epw.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/beta.disaggregated.hourly.tas.r',chdir=T)

##-----------------------------------------------------------------
##YVR Coordinates
lonc <- -123.1815
latc <- 49.1967

##-----------------------------------------------------------------
##Create hourly temperature series from BCCAQ2-PRISM data
##and the hourly EPW data

create_hourly_tas <- function(tasmax.series,tasmin.series,
                              morphed.tas,
                              epw.dates,
                              epw.tas) {
   tasmax.data <- tasmax.series$data
   tasmin.data <- tasmin.series$data
   bccaq2.dates <- tasmax.series$dates

   cal <- attr(bccaq2.dates,'cal')
   yrs <- as.numeric(range(format(bccaq2.dates,'%Y')))
   ylen <- yrs[2]-yrs[1]+1

   if (grepl('(gregorian|standard)',cal)) {
      feb.flag <- grep('-02-29',bccaq2.dates)
      tasmax.data <- tasmax.series$data[-feb.flag]
      tasmin.data <- tasmin.series$data[-feb.flag]
      bccaq2.dates <- bccaq2.dates[-feb.flag]
   }
   if (grepl('360',cal)) {
      inx <- seq(30,365*ylen,73)
      insrt <- 1:(ylen*365) %in% inx
      tasmax.fix <- rep(NA,ylen*365)     
      tasmin.fix <- rep(NA,ylen*365)     
      tasmax.fix[,!insrt] <- tasmax.data 
      tasmin.fix[,!insrt] <- tasmin.data 
      for (i in seq_along(inx)) {
         ix <- inx[i]
         tasmax.fix[,ix] <- (tasmax.fix[,ix-1] + tasmax.fix[,ix+1])/2
         tasmin.fix[,ix] <- (tasmin.fix[,ix-1] + tasmin.fix[,ix+1])/2
      }
      tasmax.data <- tasmax.fix
      tasmin.data <- tasmin.fix
      bccaq2.dates <- seq(from=as.Date(paste0(yrs[1],'-01-01')),
                          by='day',
                          to=as.Date(paste0(yrs[2],'-12-31')))
      feb.flag <- grep('-02-29',bccaq2.dates)
      bccaq2.dates <- bccaq2.dates[-feb.flag]
   }

   epw.day.fac <- as.factor(format(epw.dates,'%Y-%m-%d'))

   bccaq2.dirurnal.tas <- tasmax.data-tasmin.data
   epw.tas.max <- tapply(epw.tas,epw.day.fac,max)
   epw.tas.min <- tapply(epw.tas,epw.day.fac,min)
   epw.tas.range <- epw.tas.max-epw.tas.min

   daily.ratio <- bccaq2.dirurnal.tas / rep(epw.tas.range,ylen)
  
   epw.tas.anoms <- tapply(epw.tas,epw.day.fac,function(x){x-min(x)})
   epw.tas.anoms.rep <- rep(epw.tas.anoms,ylen)

   scaled.anoms <- mapply('*',epw.tas.anoms.rep, as.list(daily.ratio),SIMPLIFY=FALSE)
   new.tas <- mapply('+',scaled.anoms,as.list(tasmin.data),SIMPLIFY=FALSE)
   names(new.tas) <- bccaq2.dates
   ##pdf(file='/storage/data/projects/rci/data/cas/wx_files/yvr.hourly.tas.pdf',width=7,height=5)
   png(file='/storage/data/projects/rci/data/cas/wx_files/yvr.hourly.tas.jul.png',width=7,height=3,units='in',res=300)
   par(mar=c(4,4,3,2))
   plot(unlist(new.tas[1:31]),ylim=c(8,32),axes=F,xaxs='i',type='l',col='white',
        ylab='Temperature (\u00B0C)',xlab='Day',main='YVR Hourly Temperature from EPW File (July)') ##BCCAQ2+PRISM and EPW')
   abline(v=seq(0,744,24),lty=2,col='gray')
   ##lines(unlist(new.tas[1:20]))
   ###lines(epw.tas[1:744],col='black')
   lines(epw.tas[4345:5088],col='black')
   ###lines(morphed.tas[4345:5088],col='red',lwd=2)
   axis(1,at=seq(0,744,48),seq(0,31,by=2))
   axis(2,at=seq(5,35,5),seq(5,35,5))
   box(which='plot')
   ##legend('topleft',leg=c('EPW','GCM-PRISM'),col=c('red','black'),pch=16,bg='white')
   legend('topleft',leg=c('EPW'),col=c('black'),pch=15,bg='white',cex=0.7)
   ##legend('topleft',leg=c('EPW','Morphed'),col=c('black','red'),pch=15,bg='white',cex=0.7,pt.cex=1)
   dev.off()
   browser()

   return(new.tas)
}

##-----------------------------------------------------------------
##Read YVR cell and save as RData (for 12 GCMs)

create_bccaq2_prism_series <- function(var.name) {
   bccaq2.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
   save.dir <- '/storage/data/projects/rci/weather_files/YVR/'

   gcm.list <- c('ACCESS1-0','CCSM4','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')

   bccaq2.file <- paste0(save.dir,var.name,'_YVR_BCCAQ2-PRISM_daily_series_1951-2100.RData')
   if (!file.exists(bccaq2.file)) {

      bccaq2.series <- retrieve_bccaq2_prism_series(
                                var.name=var.name,scenario='rcp85',
                                gcm.list=gcm.list,
                                lonc=lonc,latc=latc,
                                base.dir=bccaq2.dir)     

      save(bccaq2.series,file=bccaq2.file)
   } else {
      load(bccaq2.file)
   }
   return(bccaq2.series)
}

##-----------------------------------------------------------------
##Read EPW hourly series for YVR
epw.file <- 'CAN_BC_Vancouver.Intl.AP.718920_CWEC2016.epw'
epw.dir <- '/storage/data/projects/rci/weather_files/wx_2016/'
epw.data <- read.epw.file(epw.dir,epw.file)
epw.tas <- epw.data$data[,7]
epw.years <- rep(1995,length(epw.tas)) ##epw.data$data[,1]
epw.mons <- epw.data$data[,2]
epw.days <- epw.data$data[,3]

##browser()

##-----------------------------------------------------------------
##Read EPW hourly series for YVR morphed with 1971-2000
epw.file <- '2050s_CAN_BC_Vancouver.Intl.AP.718920_CWEC2016.epw'
epw.dir <- '/storage/data/projects/rci/weather_files/wx_2016/wx_2016_morphed_files/Vancouver_Intl_AP/'
epw.data <- read.epw.file(epw.dir,epw.file)
morphed.epw.tas.1971 <- epw.data$data[,7]


##-----------------------------------------------------------------
##Read EPW hourly series for YVR morphed with 1998-2014
epw.file <- '2050s_CAN_BC_Vancouver.Intl.AP.718920_CWEC2016.epw'
epw.dir <- '/storage/data/projects/rci/weather_files/wx_2016/wx_2016_morphed_1998_2014/Vancouver_Intl_AP/'
epw.data <- read.epw.file(epw.dir,epw.file)
morphed.epw.tas.1998 <- epw.data$data[,7]



##-----------------------------------------------------------------
##

tasmax.series <- create_bccaq2_prism_series(var.name='tasmax') 
tasmin.series <- create_bccaq2_prism_series(var.name='tasmin') 



##-----------------------------------------------------------------
##

epw.dates <- as.Date(paste(epw.years,epw.mons,epw.days,sep='-'))
day.fac <- as.factor(format(epw.dates,'%Y-%m-%d'))


##-----------------------------------------------------------------
##
##Originally was making the plots using this function
new.series <- create_hourly_tas(tasmax.series[[1]],tasmin.series[[1]],
                                morphed.tas=morphed.epw.tas.1971,
                                epw.dates,
                                epw.tas)
browser()
##-----------------------------------------------------------------




   tasmax.series <- tasmax.series[[1]]
   tasmin.series <- tasmin.series[[1]]
   bccaq2.dates <- tasmax.series$dates
   cal <- attr(bccaq2.dates,'cal')
   yrs <- as.numeric(range(format(bccaq2.dates,'%Y')))
   ylen <- yrs[2]-yrs[1]+1

   if (grepl('(gregorian|standard)',cal)) {
      feb.flag <- grep('-02-29',bccaq2.dates)
      tasmax.data <- tasmax.series$data[-feb.flag]
      tasmin.data <- tasmin.series$data[-feb.flag]
      bccaq2.dates <- bccaq2.dates[-feb.flag]
   }
   if (grepl('360',cal)) {
      inx <- seq(30,365*ylen,73)
      insrt <- 1:(ylen*365) %in% inx
      tasmax.fix <- rep(NA,ylen*365)     
      tasmin.fix <- rep(NA,ylen*365)     
      tasmax.fix[,!insrt] <- tasmax.data 
      tasmin.fix[,!insrt] <- tasmin.data 
      for (i in seq_along(inx)) {
         ix <- inx[i]
         tasmax.fix[,ix] <- (tasmax.fix[,ix-1] + tasmax.fix[,ix+1])/2
         tasmin.fix[,ix] <- (tasmin.fix[,ix-1] + tasmin.fix[,ix+1])/2
      }
      tasmax.data <- tasmax.fix
      tasmin.data <- tasmin.fix
      bccaq2.dates <- seq(from=as.Date(paste0(yrs[1],'-01-01')),
                          by='day',
                          to=as.Date(paste0(yrs[2],'-12-31')))
      feb.flag <- grep('-02-29',bccaq2.dates)
      bccaq2.dates <- bccaq2.dates[-feb.flag]
   }
   
   ##MetSim Polynomial
   tas.max.time <- seq(17,length.out=31,by=24)
   tas.min.time <- seq(5,length.out=31,by=24)
   tas.ix <- order(c(tas.min.time,tas.max.time))
   tas.min.max <- c(tasmin.data[1:31],tasmax.data[1:31])[tas.ix]
   tas.time <- c(tas.min.time,tas.max.time)[tas.ix]

   metsim <- pchip(xi=tas.time,yi=tas.min.max,
                   x=1:744)

   ##Beta series
   beta.series <- get_beta_series(epw.tas,epw.dates,tasmax.data[1:31],tasmin.data[1:31])

   epw.day.fac <- as.factor(format(epw.dates,'%Y-%m-%d'))

   bccaq2.dirurnal.tas <- tasmax.data-tasmin.data
   epw.tas.max <- tapply(epw.tas,epw.day.fac,max)
   epw.tas.min <- tapply(epw.tas,epw.day.fac,min)
   epw.tas.range <- epw.tas.max-epw.tas.min

   daily.ratio <- bccaq2.dirurnal.tas / rep(epw.tas.range,ylen)
  
   epw.tas.anoms <- tapply(epw.tas,epw.day.fac,function(x){x-min(x)})
   epw.tas.anoms.rep <- rep(epw.tas.anoms,ylen)

   scaled.anoms <- mapply('*',epw.tas.anoms.rep, as.list(daily.ratio),SIMPLIFY=FALSE)
   new.tas <- mapply('+',scaled.anoms,as.list(tasmin.data),SIMPLIFY=FALSE)
   names(new.tas) <- bccaq2.dates
   ##pdf(file='/storage/data/projects/rci/data/cas/wx_files/yvr.hourly.tas.pdf',width=7,height=5)
   png(file='/storage/data/projects/rci/data/cas/wx_files/yvr.hourly.tas.beta.png',width=7,height=5,units='in',res=300)
   plot(unlist(new.tas[1:31]),ylim=c(-15,15),axes=F,xaxs='i',type='l',col='white',
        ylab='Temperature (\u00B0C)',xlab='Day',main='YVR Hourly TAS from EPW File (January)') ##BCCAQ2+PRISM and EPW')
   abline(v=seq(0,744,24),lty=2,col='gray')
   ##lines(metsim,col='purple')
   lines(beta.series,col='goldenrod')
   ##lines(unlist(new.tas[1:31]),col='red')
   ##lines(epw.tas[1:744],col='black')

   ##Add daily max and min
   days <- 31
   hours <- rep(1:31,each=24)
   for (d in 1:days) {    
     lines(which(hours==d),rep(tasmax.data[d],24),col='red')
     lines(which(hours==d),rep(tasmin.data[d],24),col='blue')
   }


##   lines(morphed.epw.tas.1971,col='red')
##   lines(morphed.epw.tas.1998,col='orange')



   axis(1,at=seq(0,744,48),seq(0,31,by=2))
   axis(2,at=seq(-15,15,5),seq(-15,15,5))
   box(which='plot')
   ##legend('bottomright',leg=c('EPW','Morphed (71-00)','Morphed (98-14)'),col=c('black','red','orange'),pch=16,bg='white')
   ##legend('topleft',leg=c('GCM-PRISM'),col=c('red'),pch=16,bg='white')
   ##legend('topleft',leg=c('MetSim'),col=c('purple'),pch=16,bg='white')
   legend('topleft',leg=c('Beta'),col=c('goldenrod'),pch=16,bg='white')
   ##legend('topleft',leg='EPW',col='black',pch=16,bg='white')
   dev.off()
