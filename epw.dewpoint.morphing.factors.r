##Script to generate shift and stretch factors for temperature
##from the BCCAQ2-TPS downscaled fields

library(ncdf4)
library(PCICt)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)

source('/storage/home/ssobie/code/repos/building_code/epw.morphing.factor.functions.r')
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##**************************************************************************************
##**************************************************************************************

scenario <- 'rcp85'
agg.fxn <- mean
gcm <- 'HadGEM2-ES'

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}

dwpt.dir <- '/storage/data/climate/downscale/CMIP5/building_code/'
epw.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/epw_factors/'

tmp.dir <- '/local_temp/ssobie/epw/'
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

##--------------------------------------------------------------------------------------------
##Copy temperature files to temporary directory

##dir.create(paste0(tmp.dir,'/',gcm,'/'))
scen.files <- list.files(path=paste0(dwpt.dir,gcm,'/'),pattern=scenario)
dewpoint.file <- scen.files[grep('dewpoint_day',scen.files)]
file.copy(from=paste0(dwpt.dir,gcm,'/',dewpoint.file),to=paste0(tmp.dir,'/'),overwrite=TRUE)

##------------------------------------------------------------------------------


dewpoint.nc <- nc_open(paste0(tmp.dir,'/',dewpoint.file))


##
dewpoint.1980s.avg <- daily_aggregate(dewpoint.nc,'dewpoint',
                        gcm=gcm,interval='1971-2000',
                        agg.fxn=mean)
dewpoint.1980s.sd <- daily_aggregate(dewpoint.nc,'dewpoint',
                        gcm=gcm,interval='1971-2000',
                        agg.fxn=sd)

intervals <- c('2011-2040','2041-2070','2071-2100')

for (interval in intervals) {

  dewpoint.proj.avg <- daily_aggregate(dewpoint.nc,'dewpoint',
                        gcm=gcm,interval=interval,
                        agg.fxn=mean)
  dewpoint.proj.sd <- daily_aggregate(dewpoint.nc,'dewpoint',
                        gcm=gcm,interval=interval,
                        agg.fxn=sd)

  delta.dewpoint <- dewpoint.proj.avg - dewpoint.1980s.avg

  delta.file <- paste0('delta_dewpoint_',gcm,'_1971-2000_',interval,'.nc')
  create_factor_file(dewpoint.nc,'dewpoint','degC',
                     delta.dewpoint,
                     tmp.dir,delta.file)
  file.copy(from=paste0(tmp.dir,delta.file),to=epw.dir,overwrite=TRUE)

  alpha.dewpoint <- dewpoint.proj.sd / dewpoint.1980s.sd
  alpha.file <- paste0('alpha_dewpoint_',gcm,'_1971-2000_',interval,'.nc')
  create_factor_file(dewpoint.nc,'alpha_dewpoint','ratio',
                   alpha.dewpoint,
                   tmp.dir,alpha.file)
  file.copy(from=paste0(tmp.dir,alpha.file),to=epw.dir,overwrite=TRUE)

  file.remove(paste0(tmp.dir,alpha.file))
  file.remove(paste0(tmp.dir,delta.file))

}

nc_close(dewpoint.nc)

file.remove(paste0(tmp.dir,dewpoint.file))

