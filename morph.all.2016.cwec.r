##Script to plot the EPW file series

library(ncdf4)
library(PCICt)

source('/storage/home/ssobie/code/repos/epw/read.write.epw.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/epw.precomputed.belcher.functions.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/epw.prism.offset.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/epw.support.functions.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/epw.projected.stats.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/figures.for.summary.tables.r',chdir=T)

##------------------------------------------------------------------------------
##Read specified cell from the precomputed morphing values

read_cell <- function(var.name,lonc,latc,input.file,read.dir) {

  if (length(input.file) != 1) {
    print(input.file)              
    stop('More than one input file')
  }

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
  print(c(lon[lon.ix],lat[lat.ix]))
  data <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  
  if (grepl('Zad',input.file)) { ##Zad to skip for now this is to fix Hadley 360
     data.fix <- 1:365
     ix <- seq(40,365,70)
     flags <- data.fix %in% ix
     data.fix[!flags] <- data
     for (i in ix) {
       data.fix[i] <- (data.fix[i-1]+data.fix[i+1])/2
     }
     data <- data.fix
     time.series <- seq(from=as.Date('1995-01-01'),by='day',to=as.Date('1995-12-31'))
     print('Zad')
     browser()
  }
  nc_close(nc)

  rv <- list(data=data,time=time.series)
  return(rv)
}

##------------------------------------------------------------------------------
##Select the morphing function based on the variable

get_morphing_function <- function(epw.var) {

   fxn <- list(dry_bulb_temperature=generate_dry_bulb_temp,
               dew_point_temperature=generate_dew_point_temp,
               diffuse_horizontal_radiation=generate_horizontal_radiation,
               global_horizontal_radiation=generate_horizontal_radiation,
               direct_normal_radiation=generate_stretched_series,
               relative_humidity=generate_stretched_series,
               atmospheric_station_pressure=generate_stretched_series,
               wind_speed=generate_stretched_series,
               total_sky_cover=generate_stretched_series,
               opaque_sky_cover=generate_stretched_series,
               liquid_precip_quantity=generate_stretched_series)
    rv <- fxn[[epw.var]]
   return(rv)
}

##------------------------------------------------------------------------------
##Morph GCM list

morph_epw_by_each_gcm <- function(epw.file,variable.list,lon,lat,
                                  gcm.list,gcm.dir,
                                  scenario,
                                  past.int,proj.int,
                                  method,rlen) {

   pef.split <- strsplit(epw.file$file,'_')[[1]]
   epw.present <- read.epw.file(epw.file$dir,epw.file$file) 
 
   ##-------------------------------------------------------
   ##Calculate the morphed time series of hourly values from
   ##the EPW file for each possible variable
   morphed.gcm.list <- vector(mode='list',length=length(variable.list))
   for (i in seq_along(variable.list)) {
      var.names <- variable.list[[i]]
      print(var.names$epw)
      morphing_fxn <- get_morphing_function(var.names$epw)

      morphed.gcm.list[[i]] <- morphing_fxn(epw.present=epw.present,
                                       epw.var=var.names$epw,gcm.var=var.names$gcm,
                                       lon=lon,lat=lat,
                                       gcm=gcm.list,gcm.dir=gcm.dir, 
                                       scenario=scenario,
                                       past.int=past.int,
                                       proj.int=proj.int,
                                       method=method,rlen=rlen)
      names(morphed.gcm.list)[i] <- var.names$epw
   }
   return(morphed.gcm.list) 
}

##------------------------------------------------------------------------------
##Create individual EPW files for each GCM morphing and 
##save temporarily

create_gcm_morphed_epw_files <- function(morphed.gcm.list,variable.list,
                                         gcm.list,epw.file,interval,tmp.dir) {

   epw.present <- read.epw.file(epw.file$dir,epw.file$file) 
   for (g in seq_along(gcm.list)) {
      gcm <- gcm.list[g]
      epw.gcm.present <- epw.present 
      for (v in seq_along(variable.list)) {
         var.names <- variable.list[[v]]
         morphed.gcm <- morphed.gcm.list[[var.names$epw]]
         epw.gcm.present$data[,get_field_index(var.names$epw)] <- morphed.gcm[g,]
      }
      epw.gcm.present$header[1] <- paste0(epw.gcm.present$header[1],',',gcm)
      epw.gcm.file <- paste0(gcm,'_',interval,'_',epw.file$file)
      write.epw.file(epw.gcm.present$data,epw.gcm.present$header,paste0(tmp.dir,'gcm_epws/'),epw.gcm.file)           
   }
}

##------------------------------------------------------------------------------
##Create an ensemble average EPW file that will be written
##to file and be made available 

create_ensemble_average_morphed_epw <- function(epw.file,variable.list,morphed.gcm.list,file.version,
                                                write.dir,future.epw.file) {
    epw.present <- read.epw.file(epw.file$dir,epw.file$file) 
    epw.ens.present <- epw.present 
    for (v in seq_along(variable.list)) {
       var.names <- variable.list[[v]]
       rd <- 0
       if (grepl('(temperature|speed)',var.names$epw)) {
         rd <- 1
       }
       morphed.gcm <- morphed.gcm.list[[var.names$epw]]
       epw.ens.present$data[,get_field_index(var.names$epw)] <- round(apply(morphed.gcm,2,mean),rd)
    }
    short.names <- paste(unlist(lapply(variable.list,function(x){get_short_name(x$epw)})),collapse=',')
    epw.ens.present$header[1] <- paste0(epw.present$header[1],
                                        ' | Morphed:',short.names,
                                        ' | File Version: ',file.version,
                                        ' | Creation Date: ',format(Sys.time(),'%Y-%m-%d'))
 
##    epw.ens.present$header[1] <- paste0(epw.present$header[1],
##                                        ' Morphed:',short.names,
##                                        ' with an ENSEMBLE of 10 GCMs, Daily Morphing Factors with 21-Day Rolling Mean')
    
    write.epw.file(epw.ens.present$data,epw.ens.present$header,paste0(tmp.dir,'gcm_epws/'),future.epw.file)
    ##This copy below is passed to be used externally
    write.epw.file(epw.ens.present$data,epw.ens.present$header,write.dir,future.epw.file)   
}

##------------------------------------------------------------------------------
##Run the EPW File stats here 

calculate_epw_file_statistics <- function(epw.file,interval,gcm.list,stats.names) {

  glen <- length(gcm.list)
  epw.present <- read.epw.file(epw.file$dir,epw.file$file) 

  past.cwec <- calc_cwec_values(epw.file$file,epw.file$dir,stats.names)             
  proj.cwec <- matrix(NA,nrow=glen,ncol=length(stats.names))
  
  for (g in seq_along(gcm.list)) {
     gcm <- gcm.list[g]
     print(paste0(gcm,'_',interval,'_',epw.file$file))
     proj.cwec[g,] <- calc_cwec_values(paste0(gcm,'_',interval,'_',epw.file$file),paste0(tmp.dir,'gcm_epws/'),stats.names)              
  }       
  rv <- list(past=past.cwec,proj=proj.cwec)
  return(rv)
}

##------------------------------------------------------------------------------
##

create_cwec_table_sheets <- function(epw.file,past.int,intervals,lon,lat,
                                     gcm.list,variable.list,stats.list,
                                     tmp.dir,scenario,file.version,
                                     method,rlen,write.dir,fig.dir) {
   coords <- c(lon,lat)
   ###base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
   base.dir <- '/storage/data/climate/downscale/BCCAQ2/bccaqv2_climatologies/'
   gcm.site.tas  <- calc_gcm_tas_stats(coords,scenario,gcm.list,past.int,base.dir)

   pef.split <- strsplit(epw.file$file,'_')[[1]]
   epw.present <- read.epw.file(epw.file$dir,epw.file$file) 
   dates <- paste(sprintf('%04d',epw.present$data[,1]),
                  sprintf('%02d',epw.present$data[,2]),sep='-')
   tmy.years <- substr(unique(dates),1,4)

   stats.names <- unlist(lapply(stats.list,function(x){return(x$name)}))
   for (interval in intervals) {
      morphed.gcm.list <- morph_epw_by_each_gcm(epw.file=epw.file,
                              variable.list=variable.list,lon,lat,
                              gcm.list=gcm.list,gcm.dir=paste0(tmp.dir,'epw_factors/'),
                              scenario=scenario,
                              past.int=past.int,proj.int=interval,
                              method=method,rlen=rlen)

      create_gcm_morphed_epw_files(morphed.gcm.list=morphed.gcm.list,
                                 variable.list=variable.list,
                                 gcm.list=gcm.list,epw.file=epw.file,
                                 interval=interval,
                                 tmp.dir=tmp.dir)
      interval.name <- paste0(as.numeric(strsplit(interval,'-')[[1]][2]) - 20,'s')

      future.epw.file <- paste0(interval.name,'_CAN_',prov,'_',pef.split[3],'_CWEC2016.epw')
      create_ensemble_average_morphed_epw(epw.file=epw.file,
                                       variable.list=variable.list,
                                       morphed.gcm.list=morphed.gcm.list,
                                       file.version=file.version,
                                       write.dir=write.dir,future.epw.file=future.epw.file)
   }  

   cwec.2020s <- calculate_epw_file_statistics(epw.file=epw.file,interval='2011-2040',
                                               gcm.list=gcm.list,
                                               stats.names=stats.names)
   cwec.2050s <- calculate_epw_file_statistics(epw.file=epw.file,interval='2041-2070',
                                               gcm.list=gcm.list,
                                               stats.names=stats.names)
   cwec.2080s <- calculate_epw_file_statistics(epw.file=epw.file,interval='2071-2100',
                                               gcm.list=gcm.list,
                                               stats.names=stats.names)

   cwec.entries <- vector(mode='list',length=length(stats.names))
   glen <- length(gcm.list)
   for (v in seq_along(stats.names)) {
      vals <- cbind(rep(cwec.2020s$past[v],glen),
                    cwec.2020s$proj[,v],
                    cwec.2050s$proj[,v],
                    cwec.2080s$proj[,v])
      rv <- get.round.val(stats.names[v])
      cwec.entries[[v]] <- make_table_row(vals,rv)
   } 
   names(cwec.entries) <- stats.names

   ##Create plots here using the ensemble values

   figs <- make_summary_figures(cwec.2020s,cwec.2050s,cwec.2080s,pef.split[3],fig.dir)

   check.dir <- paste0(base.dir,'ACCESS1-0/rcp85/annual/climatologies/')
   rv <- cwec.entries
   rv <- list(cwec=cwec.entries,
              model=gcm.site.tas,
              years=tmy.years,
              figs=figs)

   return(rv)
}



##------------------------------------------------------------------------------

##**************************************************************************************
##---------------------------------------
##Info passed in from submission script
scenario <- 'rcp85'

method <- 'roll'
rlen <- '21'

if (method!='roll') { 
  rlen <- ''
}

##--------------------------------------
##Base directories and info

##prism.dir <- '/storage/data/climate/PRISM/dataportal/'
##offsets.dir <- '/storage/data/projects/rci/weather_files/wx_2016/offsets/'
##wx.morph.dir <- '/storage/data/projects/rci/weather_files/wx_2016/wx_2016_morphed_1998_2014/'

provinces <- c('alberta','new_brunswick','northwest_territories','nunavut','prince_edward_island','saskatchewan',
               'british_columbia','manitoba','newfoundland_and_labrador','nova_scotia','ontario','quebec','yukon')

epw.dir <- '/storage/data/projects/rci/weather_files/canada_cwec_files/canada_epw_files/' 
wx.morph.dir <- '/storage/data/projects/rci/weather_files/canada_cwec_files/canada_morphed_epw_files/' 

##morph.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/epw_factors/'
tas.morph.dir <- '/storage/data/climate/downscale/BCCAQ2/epw_factors/rolled_bccaq2_1998-2014/'
other.morph.dir <- '/storage/data/climate/downscale/CMIP5/epw_factors/rolled_gcm_1998-2014/'

fig.dir <- '/storage/data/projects/rci/weather_files/canada_cwec_files/summary_cwec_figures/'

gcm.list <- c('ACCESS1-0','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MRI-CGCM3')

##Primary variable list
## !!If changing this list, also change the summary file text
variable.list <- list(list(epw='dry_bulb_temperature',gcm='tas'),
                       list(epw='relative_humidity',gcm='rhs'),
                       list(epw='dew_point_temperature',gcm='dewpoint'),
                       list(epw='atmospheric_station_pressure',gcm='psl'))                       

##Secondary variable list
##                       list(epw='diffuse_horizontal_radiation',gcm='rsds'),
##                       list(epw='global_horizontal_radiation',gcm='rsds'),
##                       list(epw='direct_normal_radiation',gcm='clt'),
##                       list(epw='wind_speed',gcm='wspd'),
##                       list(epw='total_sky_cover',gcm='clt'),
##                       list(epw='opaque_sky_cover',gcm='clt'))

stats.list <- list(list(name='hdd',type='annual',title='HDD'),
                 list(name='tnnETCCDI',type='annual',title='TNN'),
                 list(name='tasmin.annual_quantile_010',type='annual',title='Heating 99.0%'),
                 list(name='tasmin.annual_quantile_025',type='annual',title='Heating 97.5%'),
                 list(name='wetbulb.annual_quantile_010',type='annual',title='Heating (Wetbulb) 99.0%'),
                 list(name='wetbulb.annual_quantile_025',type='annual',title='Heating (Wetbulb) 97.5%'),
                 list(name='cdd',type='annual',title='CDD'),
                 list(name='cdd_10',type='annual',title='CDD10'),
                 list(name='txxETCCDI',type='annual',title='TXX'),
                 list(name='tasmax.annual_quantile_975',type='annual',title='Cooling 2.5%'),
                 list(name='wetbulb.annual_quantile_975',type='annual',title='Cooling (Wetbulb) 2.5%'),
                 list(name='tas_jan',type='monthly',title='January Tmean'),
                 list(name='tas_feb',type='monthly',title='February Tmean'),
                 list(name='tas_mar',type='monthly',title='March Tmean'),
                 list(name='tas_apr',type='monthly',title='April Tmean'),
                 list(name='tas_may',type='monthly',title='May Tmean'),
                 list(name='tas_jun',type='monthly',title='June Tmean'),
                 list(name='tas_jul',type='monthly',title='July Tmean'),
                 list(name='tas_aug',type='monthly',title='August Tmean'),
                 list(name='tas_sep',type='monthly',title='September Tmean'),
                 list(name='tas_oct',type='monthly',title='October Tmean'),
                 list(name='tas_nov',type='monthly',title='November Tmean'),
                 list(name='tas_dec',type='monthly',title='December Tmean'),
                 list(name='tas_win',type='monthly',title='Winter Tmean'),
                 list(name='tas_spr',type='monthly',title='Spring Tmean'),
                 list(name='tas_sum',type='monthly',title='Summer Tmean'),
                 list(name='tas_fal',type='monthly',title='Fall Tmean'),
                 list(name='tas_ann',type='monthly',title='Annual Tmean'))


tmp.dir <- '/local_temp/ssobie/morph/'
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
   dir.create(paste0(tmp.dir,'gcm_epws/'),recursive=TRUE)
   dir.create(paste0(tmp.dir,'epw_factors/'),recursive=TRUE)
}

file.version <- '2.2'
past.int <- '1998-2014'
intervals <- c('2011-2040','2041-2070','2071-2100')


##--------------------------------------------------------------------------------

print('Copying EPW factors to tmp (31Gb)')
##if (!file.exists(tmp.dir)) {
  tas.files <- list.files(path=tas.morph.dir,pattern='tas',full.name=TRUE)
  tas.interval.files <- tas.files[grep(past.int,tas.files)]
##  file.copy(from=tas.interval.files,to=paste0(tmp.dir,'epw_factors/'),overwrite=T)

  other.files <- list.files(path=other.morph.dir,pattern='roll21',full.name=TRUE)
  other.interval.files <- other.files[grep(past.int,other.files)]
##  file.copy(from=other.interval.files,to=paste0(tmp.dir,'epw_factors/'),overwrite=T)
##}

provinces <- 'british_columbia'

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
   ###cwec.2016.files <- cwec.2016.files[65:75] ##cwec.2016.files[c(-5,-55)] ##Remove lighthouses without grid data

   cwec.2016.files <- cwec.2016.files[42]   

   for (file in cwec.2016.files) {
      epw.coords <- get_epw_coordinates(paste0(epw.dir,province,'/'),file)
      epw.name <- strsplit(file,'_')[[1]][3]
      epw.stn.names <- strsplit(epw.name,'\\.')[[1]]
      new.location <- paste(epw.stn.names[-length(epw.stn.names)],collapse='_')  ##remove the stn id

      ##tas.fig <- paste0(fig.dir,province,'/',new.location,'/',epw.name,'_mean_temperature_boxplots.png')
      ##wx.figures <- list(tas=tas.fig)

      ##Check for presence of temperature data (coastal sites may coincide with ocean cells)
      tas.check <- read_cell(var.name='tas',lonc=epw.coords[1],latc=epw.coords[2],
                             input.file=paste0('roll21_delta_tas_CanESM2_',past.int,'_2041-2070.nc'),
                             read.dir=paste0(tmp.dir,'epw_factors/'))

      if (all(is.na(tas.check$data))) {
         print(paste0('EPW location is an ocean cell in ANUSPLIN. Unable to morph',new.location,'.'))        
      } else {
         write.dir <- paste0(wx.morph.dir,province,'/',new.location,'/')
         wx.fig.dir <- paste0(fig.dir,province,'/',new.location,'/')
         if (!file.exists(write.dir)) {
           dir.create(write.dir,recursive=TRUE)
         }
         if (!file.exists(wx.fig.dir)) {
           dir.create(wx.fig.dir,recursive=TRUE)
         }

         n.lon <- epw.coords[1]
         n.lat <- epw.coords[2]

         ##--------------------------------------
         ##For each interval, calculate the morphed
         ##epw files for each (10) GCMs for all
         ##variables possible

         ##Apply to each epw file separately

         sheets.closest <- create_cwec_table_sheets(list(file=file,dir=paste0(epw.dir,province,'/')),
                                              past.int,intervals,n.lon,n.lat,
                                              gcm.list,variable.list,stats.list,
                                              tmp.dir,scenario,file.version,
                                              method,rlen,write.dir,wx.fig.dir)
         make_formated_stats_table(nearest=epw.name,site=new.location,
                                var.list=stats.list,
                                sheets.closest=sheets.closest,
                                sheets.offset=NULL,
                                file.version=file.version,
                                method=method,rlen=rlen,write.dir)
         browser()
      }

   }
}
clean.files <- list.files(path=tmp.dir,recursive=TRUE,full.name=T)
file.remove(clean.files)