##Compare the time series of temperature data from the EPW file, 
##the available NOAA ISD data, and the ERA5-Land for the same location.

##----------------------------------------------------------------------

library(PCICt)
library(ncdf4)

source('/storage/home/ssobie/code/repos/epw/read.write.epw.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/epw.support.functions.r',chdir=T)

##------------------------------------------------------------
##Pull the time series
time_series_from_nc <- function(nc) { 

  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units

  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                     cal=time.calendar)
  time.values <- ncvar_get(nc,'time')
  time.series <- origin + time.values*3600

  time.min <- as.character(format(head(time.series,1),'%Y%m%d'))
  time.max <- as.character(format(tail(time.series,1),'%Y%m%d'))

  rv <- list(units=time.units,
             values=time.values,
             series=time.series,
             calendar=time.calendar,
             tmax=time.max,
             tmin=time.min)
  return(rv)
}


##----------------------------------------------------------------------
##Fetch EPW data

get_epw_series <- function(wmo.code,province,epw.var) {

   epw.dir <- paste0('/storage/data/projects/rci/weather_files/canada_cwec_files/canada_epw_files/',
                      province,'/')
   cwec.2016 <- list.files(epw.dir,pattern='CWEC2016') 
   epw.file <- cwec.2016[grep(wmo.code,cwec.2016)]
   print(epw.file)
   epw.coords <- get_epw_coordinates(epw.dir,epw.file)
   epw.name <- strsplit(epw.file,'_')[[1]][3]

   epw.info <- read.epw.file(epw.dir,epw.file)
   var.ix <- get_field_index(epw.var)

   epw.series <- epw.info$data[,var.ix]
   dates <- as.Date(paste(epw.info$data[,1],
                    sprintf('%02d',epw.info$data[,2]),
                    sprintf('%02d',epw.info$data[,3]),sep='-'))
   hours <- sprintf('%02d',epw.info$data[,4])

   rv <- list(date=dates,hour=hours,series=epw.series,coords=epw.coords,ename=epw.name)
   return(rv)    

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

get_era5_land_series <- function(epw.name,mon,yr) {

   save.dir <- "/storage/data/projects/rci/weather_files/era5_land_at_epw_sites/"
   
   load.file <- paste0(save.dir,epw.name,"_era5-land_tas_hourly_19810101-20191231.RData")
   load(load.file)

   mon.ix <- format(as.Date(era.series$date),'%Y-%m') == paste0(yr,'-',sprintf('%02d',mon))
   era.hours <- format(as.POSIXct(era.series$date),'%H')[mon.ix]

   era.tas <- era.series$series[mon.ix]   
   rv <- list(date=era.series$date[mon.ix],hour=era.hours,series=era.tas)
   return(rv)      
}


##----------------------------------------------------------------------
##**********************************************************************

##Set a site to compare

site <- 'Vancouver.Intl'
wmo.code <- '718920'
province <- 'british_columbia'
epw.var <- 'dry_bulb_temperature'
isd.var <- 'TAS'

epw.series <- get_epw_series(wmo.code,province,epw.var)
isd.series <- get_isd_series(wmo.code,province,isd.var)

##-----------------------------------------------------
##Plot one month at a time for now

par(mar=c(0,0,0,0),oma=c(2,2,2,2))
par(mfrow=c(4,3))

for (i in 1:12) {
month.ix <- i

ex <- format(epw.series$date,'%m') == sprintf('%02d',month.ix)
##Subtract 3600 to align the hour starting at 0 instead of 1
epw.time <- as.POSIXct(paste0(epw.series$date[ex],' ',epw.series$hour[ex],':00:00'))-3600

ey <- unique(format(epw.series$date,'%Y')[ex])
print('Year')
print(ey)
ix <- format(isd.series$date,'%Y-%m') == paste0(ey,'-',sprintf('%02d',month.ix))
##Offset by 8 hours as ISD seems to be UTC
isd.time <- as.POSIXct(paste0(isd.series$date[ix],' ',isd.series$hour[ix],':00:00'))-(8*3600)

era.series <- get_era5_land_series(epw.series$ename,month.ix,ey)
era.time <- as.POSIXct(era.series$date)-8*3600

xlim <- range(c(range(epw.time,na.rm=T),range(isd.time,na.rm=T)))

ylim <- range(pretty(c(epw.series$series[ex],isd.series$series[ix])))

plot(c(),xlim=xlim,ylim=ylim,xlab='',ylab='',axes=FALSE)
lines(epw.time,epw.series$series[ex],col='darkgreen',lwd=1.5)
points(epw.time,epw.series$series[ex],col='darkgreen',pch=19)
lines(isd.time,isd.series$series[ix],col='blue',lwd=1.5)
##points(isd.time,isd.series$series[ix],col='blue',pch=19)

lines(era.time,era.series$series,col='red',lwd=1.5)
##points(era.series$date,era.series$series,col='red',pch=19)
text(epw.time[round(length(epw.time)*0.95)],0.95*ylim[2],month.abb[i])

box(which='plot')

##Stats for the months
era.match <- era.time %in% epw.time
epw.era.match <- epw.time %in% era.time

isd.match <- isd.time %in% epw.time
epw.isd.match <- epw.time %in% isd.time

epw.subset <- epw.series$series[ex]
isd.subset <- isd.series$series[ix]
era.subset <- era.series$series

if (length(isd.subset) != sum(is.na(isd.subset))) {
print(paste0('Corr for ISD and EPW: ',round(cor(isd.subset[isd.match],epw.subset[epw.isd.match],use='complete.obs'),3)))
print(paste0('Corr for ERA5 and EPW: ',round(cor(era.subset[era.match],epw.subset[epw.era.match],use='complete.obs'),3)))

print(paste0('Bias for ISD and EPW: ',round(mean(isd.subset[isd.match]-epw.subset[epw.isd.match],na.rm=T),3)))
print(paste0('Bias for ERA5 and EPW: ',round(mean(era.subset[era.match]-epw.subset[epw.era.match],na.rm=T),3)))

}
}