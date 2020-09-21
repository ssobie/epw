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

##----------------------------------------------

specific_humidity <- function(dwpt,pas) {
  dwpt <- dwpt + 273
  vape.press <- sat.vape * exp( (lh.vape/R.vapor) * (1/T.zero - 1/dwpt))
  sp.hum <- vape.press * epsilon / pas
}


##----------------------------------------------

calc_cwec_values <- function(present.epw.file,epw.dir) {

  pef.split <- strsplit(present.epw.file,'_')[[1]]
  new.epw.file <- paste0('STATS_CAN_BC_',pef.split[3],'_CWEC.epw')
  epw.present <- read.epw.file(epw.dir,present.epw.file)

  tas.ix <- get_field_index('dry_bulb_temperature')
  epw.tas <- epw.present$data[,tas.ix]
  dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))
  hours <- sprintf('%02d',epw.present$data[,4])
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

  dwpt.ix <- get_field_index('dew_point_temperature')
  epw.dwpt <- epw.present$data[,dwpt.ix]
  pas.ix <- get_field_index('atmospheric_station_pressure')
  epw.pas <- epw.present$data[,pas.ix]

  ##epw.sph <- specific_humidity(epw.dwpt,epw.pas)
  ##epw.twb <- temp.wet.bulb(epw.tas+273,epw.dwpt,epw.pas,epw.sph) - 273

  ##twb.975 <- round(quantile(epw.twb,0.975,names=F),1)
  ##twb.010 <- round(quantile(epw.twb,0.01,names=F),1)
  ##twb.025 <- round(quantile(epw.twb,0.025,names=F),1)

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

##                twb.975,
##                twb.010,
##                twb.025,


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
                       'tas_fal','tas_ann')

##                       "wetbulb.annual_quantile_975",
##                       "wetbulb.annual_quantile_010",
##                       "wetbulb.annual_quantile_025",


   return(cwec.col)
}


##***********************************************************
##

epw.dir <- '/storage/data/projects/rci/weather_files/canada_cwec_files/canada_epw_files/'


stats.names <-  c("lon","lat","hdd","cdd","cdd_10",
                 "txxETCCDI","tnnETCCDI",
                 "tasmax.annual_quantile_975",
                 "tasmin.annual_quantile_010",
                 "tasmin.annual_quantile_025",
                 'tas_jan','tas_feb','tas_mar',
                 'tas_apr','tas_may','tas_jun',
                 'tas_jul','tas_aug','tas_sep',
                 'tas_oct','tas_nov','tas_dec',
                 'tas_win','tas_spr','tas_sum',
                 'tas_fal','tas_ann')

provinces <- 'british_columbia'
provinces <- sort(c('alberta','new_brunswick','nova_scotia','prince_edward_island',
             'yukon','british_columbia','newfoundland_and_labrador','nunavut',
             'quebec','manitoba','northwest_territories','ontario','saskatchewan'))

for(province in provinces) {

   prov <- switch(province,
               alberta='AB',new_brunswick='NB',
               northwest_territories='NT',nunavut='NU',
               prince_edward_island='PE',saskatchewan='SK',
               british_columbia='BC',manitoba='MB',
               newfoundland_and_labrador='NL',
               nova_scotia='NS',ontario='ON',
               quebec='QC',yukon='YK')
   
   cwec.2016.files <- list.files(path=paste0(epw.dir,province,'/'),pattern='_CWEC2016.epw')

   stats.matrix <- matrix(0,nrow=length(cwec.2016.files)+1,ncol=27)
   stats.matrix[1,] <- stats.names
   row.names <- rep('A',length(cwec.2016.files))
   for (i in seq_along(cwec.2016.files)) {
      file <- cwec.2016.files[i] 
      print(file)
      epw.coords <- get_epw_coordinates(paste0(epw.dir,province,'/'),file)
      names(epw.coords) <- c('lon','lat')
      epw.name <- strsplit(file,'_')[[1]][3]
      epw.stn.names <- strsplit(epw.name,'\\.')[[1]]
      new.location <- paste(epw.stn.names[-length(epw.stn.names)],collapse='_')  ##remove the stn id
      row.names[i] <- new.location
      cwec.stats <- calc_cwec_values(file,paste0(epw.dir,province,'/'))
      cwec.stats <- c(epw.coords,cwec.stats)
      
      stats.matrix[i+1,] <- cwec.stats     
      
   }
   rownames(stats.matrix) <- c('',row.names)
   write.file <- paste0(province,'_CWEC2016_statistics.csv')
   write.table(stats.matrix,file=paste0(epw.dir,write.file),sep=',',quote=F,row.name=T,col.name=F)
}