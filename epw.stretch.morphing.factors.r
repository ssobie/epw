##Script to generate shift and stretch factors for temperature
##from the BCCAQ2-TPS downscaled fields

library(ncdf4)
library(PCICt)
library(doParallel)
registerDoParallel(cores=4)
library(foreach)

source('/storage/home/ssobie/code/repos/building_code/epw.morphing.factor.functions.r')
##source('/storage/home/ssobie/code/repos/crcm5/epw.belcher.functions.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##**************************************************************************************
##**************************************************************************************

scenario <- 'rcp85'

agg.fxn <- mean

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}

var.name <- 'dewpoint' ##varname
var.units <- 'ratio' ##units
gcm <- 'HadGEM2-ES'

gcm.dir <- '/storage/data/climate/downscale/CMIP5/building_code/'
epw.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/epw_factors/'

tmp.dir <- '/local_temp/ssobie/epw/'
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

##--------------------------------------------------------------------------------------------
##Copy gcm files to temporary directory

scen.files <- list.files(path=paste0(gcm.dir,gcm,'/'),pattern=scenario)
gcm.file <- scen.files[grep(paste0(var.name,'_day'),scen.files)]
file.copy(from=paste0(gcm.dir,gcm,'/',gcm.file),to=paste0(tmp.dir,'/'),overwrite=TRUE)

##------------------------------------------------------------------------------

gcm.nc <- nc_open(paste0(tmp.dir,'/',gcm.file))

##
gcm.2080s <- daily_aggregate(gcm.nc,var.name,
                        gcm=gcm,interval='2071-2100',
                        agg.fxn=agg.fxn)

gcm.1980s <- daily_aggregate(gcm.nc,var.name,
                        gcm=gcm,interval='1971-2000',
                        agg.fxn=agg.fxn)

gcm.2020s <- daily_aggregate(gcm.nc,var.name,
                        gcm=gcm,interval='2011-2040',
                        agg.fxn=agg.fxn)

gcm.2050s <- daily_aggregate(gcm.nc,var.name,
                        gcm=gcm,interval='2041-2070',
                        agg.fxn=agg.fxn)


alpha.2020s <- gcm.2020s / gcm.1980s
alpha.2050s <- gcm.2050s / gcm.1980s
alpha.2080s <- gcm.2080s / gcm.1980s

#alpha.2020s[alpha.2020s > 5] <- 1
#alpha.2050s[alpha.2050s > 5] <- 1
#alpha.2080s[alpha.2080s > 5] <- 1

alpha.2020s[is.na(alpha.2020s)] <- 1
alpha.2050s[is.na(alpha.2050s)] <- 1
alpha.2080s[is.na(alpha.2080s)] <- 1

alpha.2020s.file <- paste0('alpha_',var.name,'_',gcm,'_1971-2000_2011-2040.nc')
create_factor_file(gcm.nc,var.name,var.units,
                   alpha.2020s,
                   tmp.dir,alpha.2020s.file)
file.copy(from=paste0(tmp.dir,alpha.2020s.file),to=epw.dir,overwrite=TRUE)

alpha.2050s.file <- paste0('alpha_',var.name,'_',gcm,'_1971-2000_2041-2070.nc')
create_factor_file(gcm.nc,var.name,var.units,
                   alpha.2050s,
                   tmp.dir,alpha.2050s.file)
file.copy(from=paste0(tmp.dir,alpha.2050s.file),to=epw.dir,overwrite=TRUE)

alpha.2080s.file <- paste0('alpha_',var.name,'_',gcm,'_1971-2000_2071-2100.nc')
create_factor_file(gcm.nc,var.name,var.units,
                   alpha.2080s,
                   tmp.dir,alpha.2080s.file)
file.copy(from=paste0(tmp.dir,alpha.2080s.file),to=epw.dir,overwrite=TRUE)

nc_close(gcm.nc)

file.remove(paste0(tmp.dir,gcm.file))
file.remove(paste0(tmp.dir,alpha.2020s.file))
file.remove(paste0(tmp.dir,alpha.2050s.file))
file.remove(paste0(tmp.dir,alpha.2080s.file))

