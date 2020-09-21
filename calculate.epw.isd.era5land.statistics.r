##Calculate statistics from a set of EPW files and store the values
##by location in a csv file

##-----------------------------------------------------------------

source('/storage/home/ssobie/code/repos/epw/read.write.epw.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/epw.support.functions.r',chdir=T)


##----------------------------------------------

##Degree Day Values
dd <- function(tmean,tbase) {
  g <- tmean - tbase
  days <- sum(g[g>0], na.rm=T)
  return(round(days))
}

gdd<-function(data,fac){tapply(data,fac, dd, tbase=5)}   ##Growing degree days
cdd<-function(data,fac){tapply(data,fac, dd, tbase=18)}  ##Cooling degree days
cdd10<-function(data,fac){tapply(data,fac, dd, tbase=10)}  ##Cooling degree days
hdd<-function(data,fac){tapply(-data,fac,dd, tbase=-18)} ##Heating degree days
fdd<-function(data,fac){tapply(-data,fac,dd, tbase=0)} ##Freezing degree days
ffd<-function(data,fac){tapply(data,fac,fd)} ##Frost Free days
s30<-function(data,fac){tapply(data,fac,s3)} ##Days over 30 degC

##----------------------------------------------

series_statistics <- function(input.data,epw.years=NULL,bounds=NULL) {

  dates <- input.data$date
  tas.series <- input.data$series
  hours <- input.data$hour
  tas.len <- round(sum(!is.na(tas.series))/length(tas.series)*100,2)
  print('Fraction of data available')
  print(tas.len)

   if (tas.len < 75) {
      print('Too few data points')
      series.stats <- rep(NA,26)
      return(series.stats)
   }
  
  if (!is.null(bounds)) {

     st <- head(which(format(dates,'%Y') >= bounds[1]),1)          
     en <- tail(which(format(dates,'%Y') <= bounds[2]),1)          
 
     dates <- dates[st:en]
     hours <- hours[st:en]
     tas.series <- tas.series[st:en]

  }

  if (!is.null(epw.years)) {
     isd.mon.series <- 30
     isd.mon.dates <- as.Date('1900-01-01')
     isd.mon.hours <- 24
     for (mn in 1:12) {
        ix <- format(dates,'%Y-%m') == epw.years[mn]
        isd.mon.series <- c(isd.mon.series,tas.series[ix])
        isd.mon.dates <- c(isd.mon.dates,dates[ix])
        isd.mon.hours <- c(isd.mon.hours,hours[ix])
     }

     dates <- isd.mon.dates[-1]
     hours <- isd.mon.hours[-1]
     tas.series <- isd.mon.series[-1]     

     if (length(tas.series) < 8000) {
        return(rep(NA,26))
     }
  }   

  year.mons <- unique(format(dates,'%Y-%m'))
  fac <- as.factor(format(dates,'%Y-%m-%d'))
  day.dates <- as.Date(levels(fac))
  mon.fac <- as.factor(format(day.dates,'%m'))
  seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
  seas.fac <- factor(seasons[mon.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))

  ##Mid-day temperature index
  tas.series.thresh <- tas.series >= 13 & tas.series <= 21
  mid.day.hours <- hours %in% 8:16 ##Hours between 8:00am and 4:00pm
  tas.series.thresh[!mid.day.hours] <- FALSE
  tas.series.midday <- tapply(tas.series.thresh,fac,sum)

  tas.series.daily <- tapply(tas.series,fac,mean,na.rm=T)

  if (length(levels(fac)) < 400) {
     epw.hdd <- dd(-1*tas.series.daily,-18)
     epw.cdd <- dd(tas.series.daily,18)
     epw.cdd.10 <- dd(tas.series.daily,10)
  } else {
     year.fac <- as.factor(format(day.dates,'%Y'))
     epw.hdd <- round(mean(hdd(tas.series.daily,year.fac),na.rm=T),0)
     epw.cdd <- round(mean(cdd(tas.series.daily,year.fac),na.rm=T),0)
     epw.cdd.10 <- round(mean(cdd10(tas.series.daily,year.fac),na.rm=T),0)
  }     

  epw.txx <- round(max(tas.series,na.rm=T),1)
  epw.tnn <- round(min(tas.series,na.rm=T),1)

  epw.975 <- round(quantile(tas.series,0.975,names=F,na.rm=T),1)
  epw.010 <- round(quantile(tas.series,0.01,names=F,na.rm=T),1)
  epw.025 <- round(quantile(tas.series,0.025,names=F,na.rm=T),1)

  tas.series.daily <- tapply(tas.series,fac,mean,na.rm=T)

  tas.monthly <- round(tapply(tas.series.daily,mon.fac,mean,na.rm=T),1)
  tas.seasonal <- round(tapply(tas.series.daily,seas.fac,mean,na.rm=T),1)
  tas.annual <- round(mean(tas.series.daily,na.rm=T),1)


  series.stats <- c(epw.hdd,
                    epw.cdd,
                    epw.cdd.10,
                    epw.txx,
                    epw.tnn,
                    epw.975,
                    epw.010,
                    epw.025,
                    tas.monthly,
                    tas.seasonal,
                    tas.annual,
                    tas.len)

print(series.stats)
  names(series.stats) <- c("hdd",
                           "cdd",
                           "cdd_10",
                           "txxETCCDI",
                           "tnnETCCDI",
                           "tasmax.annual_quantile_975",
                           "tasmin.annual_quantile_010",
                           "tasmin.annual_quantile_025",
                           'tas_jan','tas_feb','tas_mar',
                           'tas_apr','tas_may','tas_jun',
                           'tas_jul','tas_aug','tas_sep',
                           'tas_oct','tas_nov','tas_dec',
                           'tas_win','tas_spr','tas_sum',
                           'tas_fal','tas_ann','tas_len')

   return(series.stats)
}


##----------------------------------------------------------------------
##Fetch NOAA ISD data

get_isd_series <- function(wmo.code,province,isd.var) {

   isd.dir <- paste0("/storage/data/climate/observations/station/noaa_isd/",province,"/")
   isd.file <- list.files(path=isd.dir,pattern=wmo.code)[1]
   print(isd.file)
   isd.csv <- read.csv(paste0(isd.dir,isd.file))
   var.names <- colnames(isd.csv)
   var.ix <- grep(isd.var,var.names)
   isd.dates <- as.Date(paste(isd.csv[,1],
                              sprintf('%02d',isd.csv[,2]),
                              sprintf('%02d',isd.csv[,3]),sep='-'))
   isd.hours <- isd.csv[,4]
   isd.series <- isd.csv[,var.ix]

   rv <- list(date=isd.dates,hour=isd.hours,series=isd.series)
   return(rv)

}

##----------------------------------------------------------------------
##Fetch ERA5-Land Data data

get_era5_land_series <- function(epw.name) {

   save.dir <- "/storage/data/projects/rci/weather_files/era5_land_at_epw_sites/"

   load.file <- paste0(save.dir,epw.name,"_era5-land_tas_hourly_19810101-20191231.RData")
   load(load.file)

   era.hours <- format(as.POSIXct(era.series$date),'%H')

   era.tas <- era.series$series
   rv <- list(date=as.Date(era.series$date),hour=era.hours,series=era.tas)
   return(rv)
}

##----------------------------------------------

get_epw_series <- function(present.epw.file,epw.dir) {

  epw.present <- read.epw.file(epw.dir,present.epw.file)

  tas.ix <- get_field_index('dry_bulb_temperature')
  epw.tas <- epw.present$data[,tas.ix]
  dates <- as.Date(paste(epw.present$data[,1],
                   sprintf('%02d',epw.present$data[,2]),
                   sprintf('%02d',epw.present$data[,3]),sep='-'))

  hours <- sprintf('%02d',epw.present$data[,4])
  rv <- list(date=dates,hour=hours,series=epw.tas)

}

##----------------------------------------------

calc_cwec_values <- function(present.epw.file,epw.dir) {

  epw.present <- read.epw.file(epw.dir,present.epw.file)

  tas.ix <- get_field_index('dry_bulb_temperature')
  epw.tas <- epw.present$data[,tas.ix]
  dates <- as.Date(paste(epw.present$data[,1],
                   sprintf('%02d',epw.present$data[,2]),
                   sprintf('%02d',epw.present$data[,3]),sep='-'))

  hours <- sprintf('%02d',epw.present$data[,4])

  year.mons <- unique(format(dates,'%Y-%m'))
  fac <- as.factor(format(dates,'%Y-%m-%d'))
  day.dates <- as.Date(levels(fac))
  mon.fac <- as.factor(format(day.dates,'%m'))
  seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
  seas.fac <- factor(seasons[mon.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))

  ##Mid-day temperature index
  epw.tas.thresh <- epw.tas >= 13 & epw.tas <= 21
  mid.day.hours <- hours %in% 8:16 ##Hours between 8:00am and 4:00pm
  epw.tas.thresh[!mid.day.hours] <- FALSE
  epw.tas.midday <- tapply(epw.tas.thresh,fac,sum)

  epw.tas.daily <- tapply(epw.tas,fac,mean)

  epw.hdd <- dd(-1*epw.tas.daily,-18)
  epw.cdd <- dd(epw.tas.daily,18)
  epw.cdd.10 <- dd(epw.tas.daily,10)

  epw.txx <- max(epw.tas)
  epw.tnn <- min(epw.tas)

  epw.975 <- quantile(epw.tas,0.975,names=F)
  epw.010 <- quantile(epw.tas,0.01,names=F)
  epw.025 <- quantile(epw.tas,0.025,names=F)

  epw.tas.daily <- tapply(epw.tas,fac,mean)

  tas.monthly <- tapply(epw.tas.daily,mon.fac,mean)
  tas.seasonal <- tapply(epw.tas.daily,seas.fac,mean)
  tas.annual <- mean(epw.tas.daily)

  cwec.col <- c(epw.hdd,
                epw.cdd,
                epw.cdd.10,
                epw.txx,
                epw.tnn,
                epw.975,
                epw.010,
                epw.025,
                tas.monthly,
                tas.seasonal,
                tas.annual)

  names(cwec.col) <- c("hdd",
                       "cdd",
                       "cdd_10",
                       "txxETCCDI",
                       "tnnETCCDI",
                       "tasmax.annual_quantile_975",
                       "tasmin.annual_quantile_010",
                       "tasmin.annual_quantile_025",
                       'tas_jan','tas_feb','tas_mar',
                       'tas_apr','tas_may','tas_jun',
                       'tas_jul','tas_aug','tas_sep',
                       'tas_oct','tas_nov','tas_dec',
                       'tas_win','tas_spr','tas_sum',
                       'tas_fal','tas_ann','tas_len')

   rv <- list(cwec=cwec.col,years=year.mons)

   return(rv)
}


##***********************************************************
##

epw.dir <- '/storage/data/projects/rci/weather_files/canada_cwec_files/canada_epw_files/'
write.dir <- '/storage/data/projects/rci/weather_files/canada_cwec_files/file_statistics/'

stats.names <-  c("Data","lon","lat","hdd","cdd","cdd_10",
                 "txxETCCDI","tnnETCCDI",
                 "tasmax.annual_quantile_975",
                 "tasmin.annual_quantile_010",
                 "tasmin.annual_quantile_025",
                 'tas_jan','tas_feb','tas_mar',
                 'tas_apr','tas_may','tas_jun',
                 'tas_jul','tas_aug','tas_sep',
                 'tas_oct','tas_nov','tas_dec',
                 'tas_win','tas_spr','tas_sum',
                 'tas_fal','tas_ann','tas_len')

province <- 'british_columbia'
   
cwec.2016.files <- list.files(path=paste0(epw.dir,province,'/'),pattern='_CWEC2016.epw')
###cwec.2016.files <- cwec.2016.files[18]

##Cols: epw, isd(epw), era(epw), isd(epw range), era(epw range), isd(full), era(full)
stats.matrix <- matrix(0,nrow=length(stats.names),ncol=8)
stats.matrix[,1] <- stats.names

for (i in seq_along(cwec.2016.files)) {
   file <- cwec.2016.files[i] 
   print(file)
   epw.coords <- get_epw_coordinates(paste0(epw.dir,province,'/'),file)
   names(epw.coords) <- c('lon','lat')
   epw.name <- strsplit(file,'_')[[1]][3]
   wmo.code <-  gsub("[^0-9]","",epw.name)

   ##-----------------------
   ##Stats from EPW file     
   print('EPW Stats')
   epw.series <- get_epw_series(file,paste0(epw.dir,province,'/'))
   year.mons <- unique(format(epw.series$date,'%Y-%m'))
   cwec.stats <- series_statistics(epw.series)
   ###cwec.stats <- calc_cwec_values(file,paste0(epw.dir,province,'/'))

   epw.stats <- c('EPW File',epw.coords,cwec.stats)
   stats.matrix[,2] <- epw.stats     

   ##-----------------------
   ##Stats from ISD file matching EPW months
   print('ISD Stats')
   isd.series <- get_isd_series(wmo.code,province,isd.var='TAS')                
   print(paste0(format(range(isd.series$date),'%Y')[1],' ',format(range(isd.series$date),'%Y')[2]))
   isd.epw.stats <- c('ISD (EPW)','','',series_statistics(isd.series,epw.years=year.mons))
   stats.matrix[,3] <- isd.epw.stats
   st <- head(which(format(isd.series$date,'%Y') >= 1998),1)
   en <- tail(which(format(isd.series$date,'%Y') <= 2014),1)

   isd.range.stats <- c('ISD (Range)',format(isd.series$date,'%Y')[st],
                                      format(isd.series$date,'%Y')[en],
                                      series_statistics(isd.series,bounds=c(1998,2014)))
   stats.matrix[,4] <- isd.range.stats
   isd.all.stats <- c('ISD (All)',format(range(isd.series$date),'%Y')[1],
                                  format(range(isd.series$date),'%Y')[2],
                                  series_statistics(isd.series))
   stats.matrix[,5] <- isd.all.stats

   era.series <- get_era5_land_series(epw.name)
   st <- head(which(format(era.series$date,'%Y') >= 1998),1)
   en <- tail(which(format(era.series$date,'%Y') <= 2014),1)

   era.epw.stats <- c('ERA (EPW)', '','',
                                  series_statistics(era.series,epw.years=year.mons))
   stats.matrix[,6] <- era.epw.stats 
   era.range.stats <- c('ERA (Range)',format(era.series$date,'%Y')[st],
                                      format(era.series$date,'%Y')[en],
                                      series_statistics(era.series,bounds=c(1998,2014)))
   stats.matrix[,7] <- era.range.stats
   era.all.stats <- c('ERA (All)',format(range(era.series$date),'%Y')[1],
                                  format(range(era.series$date),'%Y')[2],
                                  series_statistics(era.series))
   stats.matrix[,8] <- era.all.stats


   write.file <- paste0(epw.name,'_EPW_ISD_ERA5-Land_temperature_statistics.csv')
   write.table(stats.matrix,file=paste0(write.dir,write.file),sep=',',quote=F,row.name=T,col.name=F)

}