##Script to plot the EPW file series

library(ncdf4)
library(PCICt)

source('/storage/home/ssobie/code/repos/crcm5/read.write.epw.r',chdir=T)
source('/storage/home/ssobie/code/repos/crcm5/epw.precomputed.belcher.functions.r',chdir=T)
source('/storage/home/ssobie/code/repos/crcm5/epw.prism.offset.r',chdir=T)


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

read_cell <- function(var.name,lonc,latc,input.file,read.dir) {

  ##print(input.file)              
  nc <- nc_open(paste(read.dir,input.file,sep=''))
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)
  time.series <- origin.pcict + time.values*86400

  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  
  lon.ix <- which.min(abs(lonc-lon))
  lat.ix <- which.min(abs(latc-lat))
  print(input.file)
  data <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  
  if (grepl('Zad',input.file)) {
     data.fix <- 1:365
     ix <- seq(40,365,70)
     flags <- data.fix %in% ix
     data.fix[!flags] <- data
     for (i in ix) {
       data.fix[i] <- (data.fix[i-1]+data.fix[i+1])/2
     }
     data <- data.fix
     time.series <- seq(from=as.Date('1995-01-01'),by='day',to=as.Date('1995-12-31'))
     browser()
  }
  nc_close(nc)

  rv <- list(data=data,time=time.series)
  return(rv)
}

##------------------------------------------------------------------------------

##**************************************************************************************

##---------------------------------------
##Info passed in from submission script
scenario <- 'rcp85'
new.location <- 'GoldenHinde'

lon <- -125.748631
lat <- 49.663237

method <- 'roll'
rlen <- '21'

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}

if (method!='roll') { 
  rlen <- ''
}

##--------------------------------------
##Base directories and info

prism.dir <- '/storage/data/climate/PRISM/dataportal/'
epw.dir <- '/storage/data/projects/rci/weather_files/wx_files/' 
offsets.dir <- '/storage/data/projects/rci/weather_files/wx_files/offsets/'
write.dir <- paste0('/storage/data/projects/rci/weather_files/wx_files/morphed_files/',new.location,'/')
if (!file.exists(write.dir)) {
  dir.create(write.dir,recursive=TRUE)
}

morph.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/epw_factors/'

gcm.list <- c('ACCESS1-0','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MRI-CGCM3')

tmp.dir <- '/local_temp/ssobie/epw/'
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
   dir.create(paste0(tmp.dir,'gcm_epws/'),recursive=TRUE)
}

intervals <- c('2011-2040','2041-2070','2071-2100')

##--------------------------------------
##First compute the offset file with 
##PRISM climatologies

epw.files <- generate_prism_offset(lon,lat,epw.dir,prism.dir,new.location)

##--------------------------------------
##Copy the morphing factors to temporary 
##directory
print('Copying EPW factors to tmp (31Gb)')
##file.copy(from=morph.dir,to=tmp.dir,recursive=TRUE,overwrite=TRUE)

##--------------------------------------
##For each interval, calculate the morphed
##epw files for each (10) GCMs for all
##variables possible

for (present.epw.file in epw.files) {
  for (interval in intervals) {
    print(interval)
    print(present.epw.file)  
    pef.split <- strsplit(present.epw.file$file,'_')[[1]]
    print(pef.split[3])
    future.epw.file <- paste0('MORPHED_',toupper(method),rlen,'_ALL_CAN_BC_',pef.split[3],'_',interval,'_CWEC.epw')
    epw.present <- read.epw.file(present.epw.file$dir,present.epw.file$file) 
    lat <- as.numeric(strsplit(epw.present$header[1],',')[[1]][7])
    lon <- as.numeric(strsplit(epw.present$header[1],',')[[1]][8])  
    print(lon)
    print(lat)

    epw.morphed.tas <- morph_dry_bulb_temp(epw.present,lon,lat,
                                           gcm=gcm.list,gcm.dir=paste0(tmp.dir,'epw_factors/'), 
                                           scenario,interval='2071-2100',
                                           method=method,rlen=rlen)
    ##write.epw.file(epw.morphed.tas$data,epw.morphed.tas$header,write.dir,future.epw.file) 
    ##future.epw.file <- paste0('MORPHED_',toupper(method),rlen,'_TAS_RHS_CAN_BC_',pef.split[3],'_',interval,'_CWEC.epw')

    epw.morphed.rhs <- generate_stretched_series(epw.present,'relative_humidity','rhs',
                            lon,lat,
                            gcm.list=gcm.list,
                            gcm.dir=paste0(tmp.dir,'epw_factors/'), 
                            scenario,interval=interval,
                            method=method,rlen=rlen)
    epw.morphed.dwpt <- morph_dew_point_temp(epw.present,lon,lat,gcm.list,
                            gcm.dir=paste0(tmp.dir,'epw_factors/'), 
                            scenario,interval,method=method,rlen=rlen)
    epw.morphed.rad <- generate_horizontal_radiation(epw.present,lon,lat,gcm.list,
                            gcm.dir=paste0(tmp.dir,'epw_factors/'), 
                            scenario,interval,method=method,rlen=rlen)
    epw.morphed.psl <- generate_stretched_series(epw.present,'atmospheric_station_pressure','psl',
                            lon,lat,
                            gcm.list=gcm.list,
                            gcm.dir=paste0(tmp.dir,'epw_factors/'), 
                            scenario,interval=interval,
                            method=method,rlen=rlen)
    epw.morphed.dnr <- generate_stretched_series(epw.present,'direct_normal_radiation','clt',
                            lon,lat,
                            gcm.list=gcm.list,
                            gcm.dir=paste0(tmp.dir,'epw_factors/'), 
                            scenario,interval=interval,
                            method=method,rlen=rlen)
    epw.morphed.wspd <- generate_stretched_series(epw.present,'wind_speed','wspd',
                            lon,lat,
                            gcm.list=gcm.list,
                            gcm.dir=paste0(tmp.dir,'epw_factors/'), 
                            scenario,interval=interval,
                            method=method,rlen=rlen)
    epw.morphed.tsc <- generate_stretched_series(epw.present,'total_sky_cover','clt',
                            lon,lat,
                            gcm.list=gcm.list,
                            gcm.dir=paste0(tmp.dir,'epw_factors/'), 
                            scenario,interval=interval,
                            method=method,rlen=rlen)
    epw.morphed.osc <- generate_stretched_series(epw.present,'opaque_sky_cover','clt',
                            lon,lat,
                            gcm.list=gcm.list,
                            gcm.dir=paste0(tmp.dir,'epw_factors/'), 
                            scenario,interval=interval,
                            method=method,rlen=rlen)

    for (g in seq_along(gcm.list)) {
       gcm <- gcm.list[g]
       epw.gcm.present <- epw.present 
       epw.gcm.present$data[,get_field_index('dry_bulb_temperature')] <- epw.morphed.tas[g,]
       epw.gcm.present$data[,get_field_index('dew_point_temperature')] <- epw.morphed.dwpt[g,]
       epw.gcm.present$data[,get_field_index('relative_humidity')] <- epw.morphed.rhs[g,]
       epw.gcm.present$data[,get_field_index('global_horizontal_radiation')] <- epw.morphed.rad$global[g,]
       epw.gcm.present$data[,get_field_index('diffuse_horizontal_radiation')] <- epw.morphed.rad$diffuse[g,]
       epw.gcm.present$data[,get_field_index('direct_normal_radiation')] <- epw.morphed.dnr[g,]
       epw.gcm.present$data[,get_field_index('atmospheric_station_pressure')] <- epw.morphed.psl[g,]     
       epw.gcm.present$data[,get_field_index('wind_speed')] <- epw.morphed.wspd[g,]      
       epw.gcm.present$data[,get_field_index('total_sky_cover')] <- epw.morphed.tsc[g,]      
       epw.gcm.present$data[,get_field_index('opaque_sky_cover')] <- epw.morphed.osc[g,]        
       epw.gcm.present$header[1] <- paste0(epw.gcm.present$header[1],',',gcm)
       epw.gcm.file <- paste0(gcm,'_',present.epw.file$file)
       write.epw.file(epw.gcm.present$data,epw.gcm.present$header,paste0(tmp.dir,'gcm_epws/'),epw.gcm.file)     

       ##Calculate the stats for each file
       ##cwec.site[g,] <- calc_cwec_values(epw.gcm.file,paste0(tmp.dir,'gcm_epws/'),var.names)
       
       
    }

    epw.ens.present <- epw.present 
    epw.ens.present$data[,get_field_index('dry_bulb_temperature')] <- round(apply(epw.morphed.tas,2,mean),1)
    epw.ens.present$data[,get_field_index('dew_point_temperature')] <- round(apply(epw.morphed.dwpt,2,mean),1)
    epw.ens.present$data[,get_field_index('relative_humidity')] <- round(apply(epw.morphed.rhs,2,mean),0)
    epw.ens.present$data[,get_field_index('atmospheric_station_pressure')] <- round(apply(epw.morphed.psl,2,mean),0)
    epw.ens.present$data[,get_field_index('global_horizontal_radiation')] <- round(apply(epw.morphed.rad$global,2,mean),0)
    epw.ens.present$data[,get_field_index('diffuse_horizontal_radiation')] <- round(apply(epw.morphed.rad$diffuse,2,mean),0)
    epw.ens.present$data[,get_field_index('direct_normal_radiation')] <- round(apply(epw.morphed.dnr,2,mean),0)
    epw.ens.present$data[,get_field_index('wind_speed')] <- round(apply(epw.morphed.wspd,2,mean),1)
    epw.ens.present$data[,get_field_index('total_sky_cover')] <- round(apply(epw.morphed.tsc,2,mean),0)
    epw.ens.present$data[,get_field_index('opaque_sky_cover')] <- round(apply(epw.morphed.osc,2,mean),0)
    epw.ens.present$header[1] <- paste0(epw.ens.present$header[1],',ENSEMBLE')

    write.epw.file(epw.ens.present$data,epw.ens.present$header,paste0(tmp.dir,'gcm_epws/'),epw.ens.file)
    write.epw.file(epw.ens.present$data,epw.ens.present$header,write.dir,future.epw.file)   



  }  ## for intervals
}  ##for epw files

