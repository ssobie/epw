##Script to plot the EPW file series

library(zoo)

##------------------------------------------------------------------------------
make_average_series <- function(series,time,method,rlen=5,agg.fxn) {

  factor <- switch(method,
                   daily='%m-%d',
                   monthly='%m',
                   roll='%m-%d',
                   seasonal='%m')
  fac <- as.factor(format(time,factor))

  if (method=='daily') {
    agg.daily <- tapply(series,fac,agg.fxn)
  }

  if (method=='monthly') {
     dtime <- as.Date(unique(format(time,'%Y-%m-%d')))
     months <- format(dtime,'%m')
     agg.data <- tapply(series,fac,agg.fxn)   
     agg.daily <- rep(0,365)
     for (mn in 1:12) {
       ix <- months == sprintf('%02d',mn)
       agg.daily[ix] <- agg.data[mn]
     }
  }

  if (method=='seasonal') {
     seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
     fac <- factor(seasons[fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))     
     agg.data <- tapply(series,fac,agg.fxn)   
     daily.seas <- c(rep('DJF',59),rep('MAM',92),rep('JJA',92),rep('SON',91),rep('DJF',31))
     agg.daily <- rep(0,365)
     agg.daily[daily.seas=='DJF'] <- agg.data['DJF']
     agg.daily[daily.seas=='MAM'] <- agg.data['MAM']
     agg.daily[daily.seas=='JJA'] <- agg.data['JJA']
     agg.daily[daily.seas=='SON'] <- agg.data['SON']
  }

  if (method=='roll') {
    agg.data <- tapply(series,fac,agg.fxn)   
    agg.daily <- rollmean(agg.data,as.numeric(rlen),fill='extend')
  }    
  
  return(agg.daily)
}

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##Separate stretching morphing functions

##Relative Humidity
morph_relative_humidity <- function(epw.series,alpha) {
   morphed.series <- epw.series * alpha
   morphed.series[morphed.series > 100] <- 100
   return(morphed.series)
}

##Direct Normal Radiation
morph_direct_normal <- function(epw.series,alpha) {
   morphed.series <- epw.series / alpha
   return(morphed.series)
}

##Atmospheric Station Pressure, Wind Speed, Liquid Precip
morph_by_stretch <- function(epw.series,alpha) {
   morphed.series <- epw.series * alpha
   return(morphed.series)
}

## Sky Cover
morph_sky_cover <- function(epw.series,alpha) {
   morphed.series <- round(epw.series * alpha)
   morphed.series[morphed.series > 10] <- 10
   return(morphed.series)
}

get_morph_function <- function(epw.name) {
   fxn <- list(direct_normal_radiation=morph_direct_normal,
               relative_humidity=morph_relative_humidity,
               atmospheric_station_pressure=morph_by_stretch,
               wind_speed=morph_by_stretch,
               total_sky_cover=morph_sky_cover,
               opaque_sky_cover=morph_sky_cover,
               liquid_precip_quantity=morph_by_stretch)
    rv <- fxn[[epw.name]]
   return(rv) 
}

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##Morph Dry Bulb Temperature using Belcher method

morph_dry_bulb_temp <- function(epw.present,lon,lat,gcm.list,gcm.dir,scenario,interval,
                                method,rlen=NULL,agg.fxn=mean) {

   tas.ix <- get_field_index('dry_bulb_temperature')
   epw.tas <- epw.present$data[,tas.ix]
   dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))

   epw.agg.mean <- epw.daily.mean <- make_average_series(epw.tas,dates,method='daily',agg.fxn=mean)
   epw.agg.max <- epw.daily.max <-  make_average_series(epw.tas,dates,method='daily',agg.fxn=max)
   epw.agg.min <- epw.daily.min <-  make_average_series(epw.tas,dates,method='daily',agg.fxn=min)

   epw.day.anoms <- epw.tas*0
   ##Daily
   for (d in 1:365) {
      ix <- format(dates,'%j') == sprintf('%03d',d)
      epw.day.anoms[ix] <- epw.tas[ix] - epw.daily.mean[d]
   }

   epw.agg.mean <- make_average_series(epw.tas,dates,method=method,agg.fxn=mean)
   epw.agg.max <-  make_average_series(epw.tas,dates,method=method,agg.fxn=max)
   epw.agg.min <-  make_average_series(epw.tas,dates,method=method,agg.fxn=min)

   tlen <- 365
   
   deltas <- matrix(0,nrow=length(gcm.list),ncol=tlen)
   alphas <- matrix(0,nrow=length(gcm.list),ncol=tlen)

   morphed.tas <- matrix(0,nrow=length(gcm.list),ncol=length(epw.tas*0))

   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      alpha.files <- list.files(path=gcm.dir,pattern=paste0('alpha_tasmax_tasmin_',gcm))
      alpha.int.file <- alpha.files[grep(interval,alpha.files)]

      alpha.tx.tn <- read_cell(var.name='alpha_tas',lonc=lon,latc=lat,
                               input.file=alpha.int.file,read.dir=gcm.dir)
      alpha.tx.tn.agg <- make_average_series(alpha.tx.tn$data,alpha.tx.tn$time,
                                          method,rlen,agg.fxn)
 
      delta.ts.files <- list.files(path=gcm.dir,pattern=paste0('delta_tas_',gcm))
      delta.ts.int.file <- delta.ts.files[grep(interval,delta.ts.files)]
      delta.ts <- read_cell(var.name='tas',lonc=lon,latc=lat,
                            input.file=delta.ts.int.file,read.dir=gcm.dir)
      delta.ts.agg <- make_average_series(delta.ts$data,delta.ts$time,
                                          method,rlen,agg.fxn)

      ##
      alpha <- alpha.tx.tn.agg / (epw.agg.max - epw.agg.min)
      for (d in 1:tlen) {
         ix <- format(dates,'%j') == sprintf('%03d',d)  
         morphed.tas[g,ix] <- epw.tas[ix] + delta.ts.agg[d] + alpha[d]*epw.day.anoms[ix]
      }

   }##GCM loop
   ##ens.morphed.tas <- apply(morphed.tas,2,mean)
   ##epw.present$data[,tas.ix] <- round(ens.morphed.tas,1)
   ##browser()
   return(morphed.tas) ##epw.present)   
}

##------------------------------------------------------------------------------
##Morph Dewpoint Temperature using Belcher method

morph_dew_point_temp <- function(epw.present,lon,lat,gcm.list,gcm.dir,scenario,interval,
                                 method,rlen=NULL,agg.fxn=mean) {

   dwpt.ix <- get_field_index('dew_point_temperature')
   epw.dwpt <- epw.present$data[,dwpt.ix]
   dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))

   epw.agg.mean <- epw.daily.mean <- make_average_series(epw.dwpt,dates,method='daily',agg.fxn=mean)
   epw.day.anoms <- epw.dwpt*0
   ##Daily
   for (d in 1:365) {
      ix <- format(dates,'%j') == sprintf('%03d',d)
      epw.day.anoms[ix] <- epw.dwpt[ix] - epw.daily.mean[d]
   }
   tlen <- 365
   
   morphed.dwpt <- matrix(0,nrow=length(gcm.list),ncol=length(epw.dwpt*0))

   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]

      alpha.files <- list.files(path=gcm.dir,pattern=paste0('alpha_dewpoint_',gcm))
      alpha.int.file <- alpha.files[grep(interval,alpha.files)]

      alpha.dwpt <- read_cell(var.name='alpha_dewpoint',lonc=lon,latc=lat,
                              input.file=alpha.int.file,read.dir=gcm.dir)
      alpha.dwpt.agg <- make_average_series(alpha.dwpt$data,alpha.dwpt$time,
                                          method,rlen,agg.fxn)

      delta.dwpt.files <- list.files(path=gcm.dir,pattern=paste0('delta_dewpoint_',gcm))
      delta.dwpt.int.file <- delta.dwpt.files[grep(interval,delta.dwpt.files)]
      delta.dwpt <- read_cell(var.name='dewpoint',lonc=lon,latc=lat,
                            input.file=delta.dwpt.int.file,read.dir=gcm.dir)
      delta.dwpt.agg <- make_average_series(delta.dwpt$data,delta.dwpt$time,
                                          method,rlen,agg.fxn)

      ##
      for (d in 1:tlen) {
         ix <- format(dates,'%j') == sprintf('%03d',d)  
         morphed.dwpt[g,ix] <- epw.agg.mean[d] + delta.dwpt.agg[d] + alpha.dwpt.agg[d]*epw.day.anoms[ix]
      }
   }##GCM loop

   ens.morphed.dwpt <- apply(morphed.dwpt,2,mean)
   epw.present$data[,dwpt.ix] <- round(ens.morphed.dwpt,1)
   return(morphed.dwpt) ##epw.present)   
}

##------------------------------------------------------------------------------
generate_horizontal_radiation <- function(epw.present,lon,lat,gcm.list,gcm.dir,scenario,interval,
                                          method,rlen=NULL,agg.fxn=mean) {

   ghr.ix <- get_field_index('global_horizontal_radiation')
   dhr.ix <- get_field_index('diffuse_horizontal_radiation')
   epw.global <- epw.present$data[,ghr.ix]
   epw.diffuse <- epw.present$data[,dhr.ix]
   flag <- epw.global == 0
   diffuse.to.global.ratio <- epw.diffuse / epw.global
   diffuse.to.global.ratio[flag] <- 0

   dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))
   
   morphed.global <- matrix(0,nrow=length(gcm.list),ncol=length(epw.global*0))
   morphed.diffuse <- matrix(0,nrow=length(gcm.list),ncol=length(epw.diffuse*0))

   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      print(gcm)      
      alpha.files <- list.files(path=gcm.dir,pattern=paste0('alpha_rsds_',gcm))
      alpha.int.file <- alpha.files[grep(interval,alpha.files)]
      print(alpha.int.file)

      alpha.rsds <- read_cell(var.name='rsds',lonc=lon,latc=lat,
                              input.file=alpha.int.file,read.dir=gcm.dir)
      alpha.rsds.agg <- make_average_series(alpha.rsds$data,alpha.rsds$time,
                                          method,rlen,agg.fxn)  
      for (d in 1:365) {
         ix <- format(dates,'%j') == sprintf('%03d',d)  
         morphed.global[g,ix] <- epw.global[ix] * alpha.rsds.agg[d]
         morphed.diffuse[g,ix] <- (epw.global[ix] * alpha.rsds.agg[d]) * diffuse.to.global.ratio[ix]           
      }

   }##GCM loop

   ens.morphed.global <- apply(morphed.global,2,mean)
   ens.morphed.diffuse <- apply(morphed.diffuse,2,mean)
   epw.present$data[,ghr.ix] <- round(ens.morphed.global,0)
   epw.present$data[,dhr.ix] <- round(ens.morphed.diffuse,0)
   rv <- list(global=morphed.global,diffuse=morphed.diffuse)
   return(rv) ##epw.present)   
}

##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
##Morph by stretching using Belcher method

generate_stretched_series <- function(epw.present,epw.var,gcm.var,lon,lat,gcm.list,gcm.dir,scenario,interval,
                                      method,rlen=NULL,agg.fxn=mean) {

   epw.ix <- get_field_index(epw.var)
   morph.fxn <- get_morph_function(epw.var)
   epw.series <- epw.present$data[,epw.ix]
   dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))
   
   morphed.series <- matrix(0,nrow=length(gcm.list),ncol=length(epw.series*0))

   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      alpha.files <- list.files(path=gcm.dir,pattern=paste0('alpha_',gcm.var,'_',gcm))
      alpha.int.file <- alpha.files[grep(interval,alpha.files)]

      alpha <- read_cell(var.name=gcm.var,lonc=lon,latc=lat,
                              input.file=alpha.int.file,read.dir=gcm.dir)
      alpha.agg <- make_average_series(alpha$data,alpha$time,
                                          method,rlen,agg.fxn)  

      for (d in 1:365) {
         ix <- format(dates,'%j') == sprintf('%03d',d)  
         morphed.series[g,ix] <- morph.fxn(epw.series[ix],alpha.agg[d]) ##epw.series[ix] * alpha[d]         
      }
   }##GCM loop

   ens.morphed.series <- apply(morphed.series,2,mean)  
   epw.present$data[,epw.ix] <- round(ens.morphed.series,0)
   return(morphed.series) ##epw.present)   
}

##------------------------------------------------------------------------------

