##Script to generate an hourly temperature series by combining BCCAQ2-PRISM 
##with the hourly temperatures from the YVR EPW File

library(ncdf4)
source('/storage/home/ssobie/code/repos/assessments/single.site.bccaq.to.prism.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/read.write.epw.r',chdir=T)

##-----------------------------------------------------------------
##YVR Coordinates
lonc <- -123.1815
latc <- 49.1967

##-----------------------------------------------------------------
##Create hourly temperature series from BCCAQ2-PRISM data
##and the hourly EPW data

create_hourly_tas <- function(tasmax.series,tasmin.series,
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
   png(file='/storage/data/projects/rci/data/cas/wx_files/yvr.hourly.tas.png',width=7,height=5,units='in',res=300)
   plot(unlist(new.tas[1:20]),ylim=c(-15,15),axes=F,xaxs='i',type='l',
        ylab='Temperature (degC)',xlab='Day',main='YVR Hourly from BCCAQ2+PRISM and EPW')
   abline(v=seq(0,480,24),lty=2,col='gray')
   lines(unlist(new.tas[1:20]))
   lines(epw.tas[1:480],col='red')

   axis(1,at=seq(0,480,24),0:20)
   axis(2,at=seq(-15,15,5),seq(-15,15,5))
   box(which='plot')
   legend('topleft',leg=c('EPW','GCM-PRISM'),col=c('red','black'),pch=16,bg='white')
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

new.series <- create_hourly_tas(tasmax.series[[1]],tasmin.series[[1]],
                                epw.dates,
                                epw.tas)
browser()



tasmax.series.sub <- tasmax.series[[1]]$data[1:365]
tasmax.dates.sub <- tasmax.series[[1]]$dates[1:365]
tasmax.date.fac <- as.factor(format(tasmax.dates.sub,'%m-%d'))

tasmin.series.sub <- tasmin.series[[1]]$data[1:365]
tasmin.dates.sub <- tasmin.series[[1]]$dates[1:365]

bccaq2.dirurnal.tas <- tasmax.series.sub-tasmin.series.sub
epw.tas.max <- tapply(epw.tas,day.fac,max)
epw.tas.min <- tapply(epw.tas,day.fac,min)

epw.tas.anoms <- tapply(epw.tas,day.fac,function(x){x-min(x)})
epw.tas.range <- epw.tas.max-epw.tas.min

daily.ratio <- bccaq2.dirurnal.tas / epw.tas.range

scaled.anoms <- mapply('*',epw.tas.anoms, as.list(daily.ratio),SIMPLIFY=FALSE)
new.tas <- mapply('+',scaled.anoms,as.list(tasmin.series.sub),SIMPLIFY=FALSE)

