##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)
library(zoo)


##-----------------------------------------------------------------
##Options for create ensemble average morphing factors

##1) Calculate 21-day rolling mean of daily series at each cell
##   Then calculate the ensemble average for each day at each cell
##2) Calculate the ensemble average for each day at each cell
##   Then apply the rolling mean to the ensemble average.

## Note for the GCM factors the models will need to be regridded
## to a common grid.


make_ensemble_average_netcdf_file <- function(file.name,ens.file,var.name,tmp.dir) {

    file.copy(from=paste0(tmp.dir,file.name),paste0(tmp.dir,ens.file),overwrite=T)

    nc <- nc_open(paste0(tmp.dir,ens.file),write=TRUE)
    ncatt_put(nc,varid=0,attname='history',attval='')
    ncatt_put(nc,varid=0,attname='time_averaging',
                      attval='Rolling 21-Day mean, ends extended')
    ncatt_put(nc,varid=0,attname='driving_experiment',
                      attval='PCIC 10 Ensemble Average, historical+rcp85')
    ncatt_put(nc,varid=0,attname='driving_model_id',
                      attval='PCIC 10 Ensemble Average')
    return(nc)
}

##--------------------------------------------------------------
##Create files with rolling mean applied to smooth the daily series

create_ensemble_average <- function(gcm.files,ens.file,var.name,tmp.dir) {
  
   gcm.one <- gcm.files[1]
   gcm.nc <- nc_open(paste0(tmp.dir,gcm.one))
   
   ens.nc <- make_ensemble_average_netcdf_file(gcm.one,ens.file,var.name,tmp.dir) 

   lon <- ncvar_get(gcm.nc,'lon')
   lat <- ncvar_get(gcm.nc,'lat')
   
   nc_close(gcm.nc)

   n.lon <- length(lon)
   n.lat <- length(lat)
   n.time <- 365

   gcm.ncs <- vector(mode='list',length=length(gcm.files))
   for (g in seq_along(gcm.files)) {
       gcm.ncs[[g]] <- nc_open(paste0(tmp.dir,gcm.files[g]))
   }

   for (j in 1:n.lat) {
      lat.ix <- j
      ltm <- proc.time()
      print(paste0('Latitude: ',j,' of ',n.lat))

      gcm.subsets <- array(NA,c(n.lon,n.time,length(gcm.files)))

      for (g in seq_along(gcm.files)) {
         gcm.subsets[,,g] <- ncvar_get(gcm.ncs[[g]],var.name,start=c(1,j,1),count=c(-1,1,-1))
      }
      ens.matrix <- apply(gcm.subsets,c(1,2),mean,na.rm=T)
      rm(gcm.subsets)
      ncvar_put(ens.nc,varid=var.name,vals=ens.matrix,
              start=c(1,lat.ix,1),count=c(-1,1,-1))
      rm(ens.matrix)
   }
   nc_close(ens.nc)
   for (g in seq_along(gcm.files)) {
     nc_close(gcm.ncs[[g]])
   }
}

##--------------------------------------------------------------


##*************************************************************************
##*************************************************************************

gcm.list <- c('ACCESS1-0','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MRI-CGCM3')

##morph.name <- 'alpha_tasmax_tasmin' ##To match the file name morphing factor
##var.name <- 'alpha_tas'

##morph.name <- 'delta_tas' ##To match the file name morphing factor
##var.name <- 'tas'

##morph.name <- 'alpha_dewpoint'
##var.name <- 'alpha_dewpoint'

past.int <- '2004-2018'

testing <- FALSE
if (testing) {
   type <- 'BCCAQ2'
   tmpdir <- '/local_temp/ssobie/ens'
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
   }
}

morph.name <- morph
proj.int <- proj ##'2071-2100'
var.name <- varname

##-------------------------------------------------------------------------

##For BCCAQ2 files
if (type=='BCCAQ2') {
   epw.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/epw_factors/rolled_bccaq2_',past.int,'/')
   write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/epw_factors/ens_avg_bccaq2_',past.int,'/')
}

##For GCM files
if(type=='GCM') {
   ##epw.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/epw_factors/rolled_gcm_',past.int,'/')
   ##write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/epw_factors/ens_avg_gcm_',past.int,'/')
   epw.dir <- paste0('/storage/data/climate/downscale/CMIP5/epw_factors/rolled_gcm_',past.int,'/')
   write.dir <- paste0('/storage/data/climate/downscale/CMIP5/epw_factors/ens_avg_gcm_',past.int,'/')
}

if (!file.exists(write.dir)) {
   dir.create(write.dir,recursive=TRUE)
}

tmp.dir <- paste0(tmpdir,'/',morph.name,'_',type,'_',proj.int,'/')
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

var.files <- list.files(path=epw.dir,pattern=morph.name)
int.files <- var.files[grep(paste0(past.int,'_',proj.int),var.files)]

gcm.files <- c()
for (gcm in gcm.list) {
   gcm.file <- int.files[grep(gcm,int.files)]
   gcm.files <- c(gcm.files,gcm.file)
}


ens.file <- gsub('ACCESS1-0','ENSEMBLE',gcm.files[1])

file.copy(from=paste0(epw.dir,gcm.files),to=tmp.dir,overwrite=T)

create_ensemble_average(gcm.files,ens.file,var.name,tmp.dir)

file.copy(from=paste0(tmp.dir,ens.file),to=write.dir,overwrite=TRUE)
