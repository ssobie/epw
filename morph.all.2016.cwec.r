##Script to plot the EPW file series

source('/storage/home/ssobie/code/repos/epw/morph.cwec.main.functions.r',chdir=T)


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

###provinces <- c('alberta','new_brunswick','northwest_territories','nunavut','prince_edward_island','saskatchewan',
###               'british_columbia','manitoba','newfoundland_and_labrador','nova_scotia','ontario','quebec','yukon')

###'alberta','new_brunswick','northwest_territories','nunavut',
provinces <- c('prince_edward_island','saskatchewan',
               'manitoba','newfoundland_and_labrador','nova_scotia','ontario','quebec','yukon')

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

file.version <- '2.1'
past.int <- '1998-2014'
intervals <- c('2011-2040','2041-2070','2071-2100')


##--------------------------------------------------------------------------------

print('Copying EPW factors to tmp (31Gb)')
##if (!file.exists(tmp.dir)) {
  tas.files <- list.files(path=tas.morph.dir,pattern='tas',full.name=TRUE)
  tas.interval.files <- tas.files[grep(past.int,tas.files)]
  file.copy(from=tas.interval.files,to=paste0(tmp.dir,'epw_factors/'),overwrite=T)

  other.files <- list.files(path=other.morph.dir,pattern='roll21',full.name=TRUE)
  other.interval.files <- other.files[grep(past.int,other.files)]
  file.copy(from=other.interval.files,to=paste0(tmp.dir,'epw_factors/'),overwrite=T)
##}

###provinces <- 'ontario'

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
   ###cwec.2016.files <- cwec.2016.files[20] 

   ###cwec.2016.files <- cwec.2016.files[42]   

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
         ##browser()
      }

   }
}
clean.files <- list.files(path=tmp.dir,recursive=TRUE,full.name=T)
file.remove(clean.files)