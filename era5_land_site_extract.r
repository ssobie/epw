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
##Fetch ERA5-Land Data data

get_era5_land_series <- function(epw.coords,var.name,nc,lon,lat,time) {
  
   lon.ix <- which.min(abs(epw.coords[1]-lon))
   lat.ix <- which.min(abs(epw.coords[2]-lat))

   bf <- 3

   ##Check for ocean cells 
   era.one <- ncvar_get(nc,'tashour',start=c(lon.ix,lat.ix,1),count=c(1,1,1))
   if (is.na(era.one)) {
      print('EPW site is on an ocean cell')      
      print(paste0('Current coordinate index: ',lon.ix,', ',lat.ix))
      era.buffer <- ncvar_get(nc,'tashour',start=c(lon.ix-bf,lat.ix-bf,1),count=c(2*bf+1,2*bf+1,1))
      lon.mat <- matrix(lon[(lon.ix-bf):(lon.ix+bf)],nrow=2*bf+1,ncol=2*bf+1,byrow=TRUE)
      lat.mat <- matrix(rev(lat[(lat.ix-bf):(lat.ix+bf)]),nrow=2*bf+1,ncol=2*bf+1,byrow=F)

      cell.dist <- round(((lon.mat-lon[lon.ix])^2 + (lat.mat-lat[lat.ix])^2) * (era.buffer/era.buffer),4)
      min.dist <- min(cell.dist,na.rm=T)
      min.arr <- which(cell.dist==min.dist,arr.ind=TRUE)[1,]
      lon.ix <- ((lon.ix-bf):(lon.ix+bf))[min.arr[1]]
      lat.ix <- ((lat.ix-bf):(lat.ix+bf))[min.arr[2]]
      print(paste0('New coordinate index: ',lon.ix,', ',lat.ix))

      era.tas <- ncvar_get(nc,'tashour',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))     
   } else {
      era.tas <- ncvar_get(nc,'tashour',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
   }

   rv <- list(date=as.character(time$series),series=era.tas)
   return(rv)      
}

##----------------------------------------------------------------------
##**********************************************************************

##Set a site to compare

##Need rest of Canada ERA5-Land to extract the other provinces
##provinces <- c('alberta','new_brunswick','northwest_territories','nunavut','prince_edward_island','saskatchewan',
##               'british_columbia','manitoba','newfoundland_and_labrador','nova_scotia','ontario','quebec','yukon')

testing <- TRUE

if (testing) {
   province <- 'british_columbia'
   tmp.dir <- '/local_temp/ssobie/era5/'
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
   }
}

era.dir <- "/storage/data/climate/observations/reanalysis/ERA5-Land/"
era.file <- "tas_hour_ERA5-Land_BC_19810101-20191231.nc"

file.copy(from=paste0(era.dir,era.file),to=tmp.dir)
print('Done copying')
nc <- nc_open(paste0(tmp.dir,era.file))
lon <- ncvar_get(nc,'lon')
lat <- ncvar_get(nc,'lat')
time <- time_series_from_nc(nc)


epw.dir <- '/storage/data/projects/rci/weather_files/canada_cwec_files/canada_epw_files/'
save.dir <- "/storage/data/projects/rci/weather_files/era5_land_at_epw_sites/"

##-----------------------------------------------------

   prov <- switch(province,
               alberta='AB',new_brunswick='NB',
               northwest_territories='NT',nunavut='NU',
               prince_edward_island='PE',saskatchewan='SK',
               british_columbia='BC',manitoba='MB',
               newfoundland_and_labrador='NL',
               nova_scotia='NS',ontario='ON',
               quebec='QC',yukon='YK')

   cwec.2016.files <- list.files(path=paste0(epw.dir,province,'/'),pattern='_CWEC2016.epw')
   cwec.2016.files <- cwec.2016.files[68:77]

   for (file in cwec.2016.files) {
      epw.coords <- get_epw_coordinates(paste0(epw.dir,province,'/'),file)
      epw.name <- strsplit(file,'_')[[1]][3]
      epw.stn.names <- strsplit(epw.name,'\\.')[[1]]
      epw.wmo <-  gsub("[^0-9]","",epw.name)
      new.location <- paste(epw.stn.names[-length(epw.stn.names)],collapse='_')  ##remove the stn id
      print(epw.name)
      era.series <- get_era5_land_series(epw.coords,var.name='tashour',nc,lon,lat,time)
 
      save.file <- paste0(save.dir,epw.name,"_era5-land_tas_hourly_19810101-20191231.RData")      
      save(era.series,file=save.file)
   }

nc_close(nc)
