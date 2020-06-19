##Script to plot the EPW file series

source('/storage/home/ssobie/code/repos/epw/morph.cwec.main.functions.r',chdir=T)

##------------------------------------------------------------------------------

##**************************************************************************************
##---------------------------------------
##Info passed in from submission script

new.location <- 'Royal_Columbian_Hospital'
new.name <- gsub('_','-',new.location)
province <- 'british_columbia'
lon <- -122.890853 ## -123.31 ###
lat <- 49.226618 ##48.5 ###

scenario <- 'rcp85'
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

prov <- switch(province,
               alberta='AB',new_brunswick='NB',
               northwest_territories='NT',nunavut='NU',
               prince_edward_island='PE',saskatchewan='SK',
               british_columbia='BC',manitoba='MB',
               newfoundland_and_labrador='NL',
               nova_scotia='NS',ontario='ON',
               quebec='QC',yukon='YK')

epw.dir <- paste0('/storage/data/projects/rci/weather_files/canada_cwec_files/canada_epw_files/',province,'/')

wx.morph.dir <- offsets.dir <- '/storage/data/projects/rci/weather_files/canada_cwec_files/canada_morphed_epw_files/'

tas.morph.dir <- '/storage/data/climate/downscale/BCCAQ2/epw_factors/rolled_bccaq2_1998-2014/'
other.morph.dir <- '/storage/data/climate/downscale/CMIP5/epw_factors/rolled_gcm_1998-2014/'
fig.dir <- '/storage/data/projects/rci/weather_files/canada_cwec_files/summary_cwec_figures/'
wx.fig.dir <- paste0(fig.dir,province,'/',new.location,'/')
if (!file.exists(wx.fig.dir)) {
   dir.create(wx.fig.dir,recursive=TRUE)
}

write.dir <- paste0(wx.morph.dir,province,'/',new.location,'/')
if (!file.exists(write.dir)) {
   dir.create(write.dir,recursive=TRUE)
}


gcm.list <- c('ACCESS1-0','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MRI-CGCM3')

##Primary variable list
## !!If changing this list, also change the summary file text
variable.list <- list(list(epw='dry_bulb_temperature',gcm='tas'),
                       list(epw='relative_humidity',gcm='rhs'),
                       list(epw='dew_point_temperature',gcm='dewpoint'),
                       list(epw='atmospheric_station_pressure',gcm='psl'))

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
##First compute the offset file with 
##PRISM climatologies

epw.files <- generate_prism_offset(lon,lat,epw.dir,prism.dir,write.dir,new.name)

nearest.epw <- strsplit(epw.files$closest$file,'_')[[1]][3]
nearest.epw.coords <- epw.files$closest$coords
n.lon <- nearest.epw.coords[1]
n.lat <- nearest.epw.coords[2]

##--------------------------------------

##Copy the morphing factors to temporary 
##directory
print('Copying EPW factors to tmp')
##if (!file.exists(tmp.dir)) {
  tas.files <- list.files(path=tas.morph.dir,pattern='tas',full.name=TRUE)
  tas.interval.files <- tas.files[grep(past.int,tas.files)]
##  file.copy(from=tas.interval.files,to=paste0(tmp.dir,'epw_factors/'),overwrite=T)

  other.files <- list.files(path=other.morph.dir,pattern='roll21',full.name=TRUE)
  other.interval.files <- other.files[grep(past.int,other.files)]
##  file.copy(from=other.interval.files,to=paste0(tmp.dir,'epw_factors/'),overwrite=T)
##}

epw.coords <- get_epw_coordinates(epw.files$closest$dir,epw.files$closest$file)

n.lon <- epw.coords[1]
n.lat <- epw.coords[2]

if (1==0) {

sheets.closest <- create_cwec_table_sheets(epw.files$closest,
                                           past.int,intervals,n.lon,n.lat,
                                           gcm.list,variable.list,stats.list,
                                           tmp.dir,scenario,file.version,
                                           method,rlen,write.dir,wx.fig.dir)
if (!is.null(epw.files$offset)) {
   sheets.offset <- create_cwec_table_sheets(epw.files$offset,
                                             past.int,intervals,lon=lon,lat=lat,
                                             gcm.list,variable.list,stats.list,
                                             tmp.dir,scenario,file.version,
                                             method,rlen,write.dir,wx.fig.dir)
} else {
   sheets.offset <- NULL
}


}

make_formated_stats_table(nearest=nearest.epw,site=new.location,
                          var.list=stats.list,
                          sheets.closest=sheets.closest,
                          sheets.offset=sheets.offset,
                          file.version=file.version,
                          method=method,rlen=rlen,write.dir)

##--------------------------------------
###clean.files <- list.files(path=tmp.dir,recursive=TRUE,full.name=T)
###file.remove(clean.files)