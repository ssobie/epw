##Script to plot the EPW file series

library(ncdf4)
library(PCICt)

source('/storage/home/ssobie/code/repos/crcm5/read.write.epw.r',chdir=T)
source('/storage/home/ssobie/code/repos/crcm5/epw.belcher.functions.r',chdir=T)


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
  if (grepl('HadGEM',gcm) & yrs[2]=='2100') {
    en <- length(years)
  }
  nc_close(nc)
  rv <- list(data=data[st:en],time=new.series[st:en])
  return(rv)
}


check_for_gcm_data <- function(lonc,latc,gcm.dir,gcm,scenario) {

  scen.files <- list.files(path=paste0(gcm.dir,gcm,'/'),pattern=scenario)
  tasmax.files <- scen.files[grep('tasmax_day_BCCAQ2',scen.files)]
  tasmax.file <- tasmax.files[grep('1951-2000',tasmax.files)]
  nc <- nc_open(paste0(gcm.dir,gcm,'/',tasmax.file))
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  lon.ix <- which.min(abs(lonc-lon))
  lat.ix <- which.min(abs(latc-lat))
  data.raw <- ncvar_get(nc,'tasmax',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))  
  rv <- sum(is.na(data.raw)) == length(data.raw)
  return(rv)
}

##------------------------------------------------------------------------------

##**************************************************************************************

scenario <- 'rcp85'

##lon <- -123.0684
##lat <- 49.3209

##method <- 'roll'
##rlen <- '21'

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

if (method!='roll') { 
  rlen <- ''
}

print(scenario)
print(method)
print(rlen)
print(infile)
print(lon)
print(lat)

epw.dir <- '/storage/data/projects/rci/weather_files/wx_files/offsets/'
write.dir <- '/storage/data/projects/rci/weather_files/wx_files/morphed_files/' ##tas_only/'

tas.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'

##infile ##'CAN_BC_1st_and_Clark_offset_from_VANCOUVER-INTL-A_1108395_CWEC.epw'

epw.files <- infile 


full.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')
sub.list <- c('ACCESS1-0','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MRI-CGCM3')

tmp.dir <- '/local_temp/ssobie/epw/'
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

for (gcm in full.list) {
  ##print(paste0('Copying: ',gcm))
  dir.create(paste0(tmp.dir,gcm,'/'))
  scen.files <- list.files(path=paste0(tas.dir,gcm,'/'),pattern=scenario)
  tasmax.files <- scen.files[grep('tasmax_day_BCCAQ2',scen.files)]
  file.copy(from=paste0(tas.dir,gcm,'/',tasmax.files),to=paste0(tmp.dir,gcm,'/'),overwrite=TRUE)
  tasmin.files <- scen.files[grep('tasmin_day_BCCAQ2',scen.files)]
  file.copy(from=paste0(tas.dir,gcm,'/',tasmin.files),to=paste0(tmp.dir,gcm,'/'),overwrite=TRUE)
}

##------------------------------------------------------------------------------
intervals <- c('2011-2040','2041-2070','2071-2100')
for (present.epw.file in epw.files) {
  for (interval in intervals) {
    print(interval)
    print(present.epw.file)  
    pef.split <- strsplit(present.epw.file,'_')[[1]]
    print(pef.split[3])
    future.epw.file <- paste0('MORPHED_',toupper(method),rlen,'_TAS_CAN_BC_',pef.split[3],'_',interval,'_CWEC.epw')
    epw.present <- read.epw.file(epw.dir,present.epw.file) 
    ##lat <- as.numeric(strsplit(epw.present$header[1],',')[[1]][7])
    ##lon <- as.numeric(strsplit(epw.present$header[1],',')[[1]][8])  
    print(lon)
    print(lat)
 
    flag <- check_for_gcm_data(lon,lat,tas.dir,'CanESM2',scenario)
    if (flag) {
      print('No data available for this location')
    } else {

    epw.morphed.tas <- morph_dry_bulb_temp(epw.present,
                              lon,lat,
                              gcm.list=full.list,
                              gcm.dir=tmp.dir, ##'/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/',##
                              scenario,interval=interval,
                              method=method,rlen=rlen)
##      browser()
    write.epw.file(epw.morphed.tas$data,epw.morphed.tas$header,write.dir,future.epw.file) 

    future.epw.file <- paste0('MORPHED_',toupper(method),rlen,'_TAS_RHS_CAN_BC_',pef.split[3],'_',interval,'_CWEC.epw')

    epw.morphed.rhs <- generate_stretched_series(epw.morphed.tas,'relative_humidity','rhs',
                            lon,lat,
                            gcm.list=sub.list,
                            gcm.dir="/storage/data/climate/downscale/CMIP5/building_code/",
                            scenario,interval=interval,
                            method=method,rlen=rlen)
##browser()
    write.epw.file(epw.morphed.rhs$data,epw.morphed.rhs$header,write.dir,future.epw.file)     

    } ##if data available
  }  ## for intervals
}  ##for epw files



##epw.morphed.dwpt <- morph_dew_point_temp(epw.present,lon,lat,gcm.list,
##                        gcm.dir="/storage/data/climate/downscale/CMIP5/building_code/",
##                        scenario,interval,method=method,rlen=rlen)


##epw.morphed.dnr <- generate_stretched_series(epw.present,'direct_normal_radiation','clt',
##                        lon,lat,
##                        gcm.list=sub.list,
##                        gcm.dir="/storage/data/climate/downscale/CMIP5/building_code/",
##                        scenario,interval=interval,
##                        method=method,rlen=rlen)



for (gcm in full.list) {
##  scen.files <- list.files(path=paste0(tmp.dir,gcm,'/'),pattern=scenario)
##  tasmax.files <- scen.files[grep('tasmax_day_BCCAQ2',scen.files)]
##  file.remove(paste0(tmp.dir,gcm,'/',tasmax.files))
##  tasmin.files <- scen.files[grep('tasmin_day_BCCAQ2',scen.files)]
##  file.remove(paste0(tmp.dir,gcm,'/',tasmin.files))
}
