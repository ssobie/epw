##Script to add attributes to the aggregated downscaled output.

library(ncdf4)
library(PCICt)
library(zoo)

library(foreach)
library(doParallel)
registerDoParallel(cores=2) # or some other number that you're comfortable with.

##-----------------------------------------------------------------
##Options for create ensemble average morphing factors

##1) Calculate 21-day rolling mean of daily series at each cell
##   Then calculate the ensemble average for each day at each cell
##2) Calculate the ensemble average for each day at each cell
##   Then apply the rolling mean to the ensemble average.

## Note for the GCM factors the models will need to be regridded
## to a common grid.

##-----------------------------------------------------------------
##Remap GCM files to common grid
regrid_gcm_to_bccaq2_grid <- function(input.file,output.file,tmp.dir) {

   work <- paste0("cdo remapbil,/storage/home/ssobie/grid_files/gcm.epw.bccaq.grid.txt ",
                   tmp.dir,input.file," ",tmp.dir,output.file)
   system(work)
   Sys.sleep(2)
     
}

##-----------------------------------------------------------------

gcm_unit_conversion <- function(varname,gcm.subset,gcm.nc) {

   units <- switch(varname,
                   pr='kg m-2 d-1',
                   tasmax='degC',
                   tasmin='degC',
                   tas='degC',
                   alpha_tas='degC')
   var.units <- ncatt_get(gcm.nc,varname,'units')$value
   if (!(varname %in% c('pr','tasmax','tasmin','tas','alpha_tas'))) {
      print(paste0('Varname: ',varname,' Var.units: ',var.units))
      stop('Incorrect units for conversion')
   }


   if (var.units != units) {
      rv <- ud.convert(gcm.subset,var.units,units)
   } else {
      rv <- gcm.subset
   }
   return(rv)
}


make_rolling_mean_netcdf_file <- function(file.name,roll.name,var.name,tmp.dir) {

    file.copy(from=paste0(tmp.dir,file.name),paste0(tmp.dir,roll.name),overwrite=T)

    nc <- nc_open(paste0(tmp.dir,roll.name),write=TRUE)
    ncatt_put(nc,varid=0,attname='history',attval='')
    ncatt_put(nc,varid=0,attname='time_averaging',
                      attval='Rolling 21-Day mean, ends extended')
    ncatt_put(nc,varid=0,attname='driving_experiment',
                      attval='PCIC 10 Ensemble Average, historical+rcp85')
    ncatt_put(nc,varid=0,attname='driving_model_id',
                      attval='PCIC 10 Ensemble Average')
    return(nc)
}


#----------------------------------------------------------------------------------------------
rolling_mean_of_days <- function(var.name,var.ncs,lat.ix,n.lon,time.len,flag,
                                 day.list) {
   ##Variables
   flen <- sum(!flag)
   n.col <- time.len
   roll.matrix <- matrix(NA,nrow=n.lon,ncol=n.col)

   if (flen!=0) { ##Some Real Values
       sub.list <- day.list[!flag]
       roll.values <- foreach(
                            data=sub.list,
                            .export=c('rollmean')
                            ) %do% {
                                 rolled.values <- rollmean(data,k=21,fill='extend')
                            }
       sub.matrix <- matrix(unlist(roll.values),nrow=flen,ncol=n.col,byrow=TRUE)
       rm(roll.values)
       roll.matrix[!flag] <- sub.matrix
       rm(sub.matrix)
       rm(sub.list)
    } else {
         print('All NA values')
    }
    rm(day.list)

    ncvar_put(var.ncs,varid=var.name,vals=roll.matrix,
              start=c(1,lat.ix,1),count=c(-1,1,-1))
    rm(roll.matrix)
    gc()
}


##--------------------------------------------------------------
##Create files with rolling mean applied to smooth the daily series

create_rolling_means <- function(file.name,roll.name,var.name,tmp.dir) {

   gcm.nc <- nc_open(paste0(tmp.dir,file.name))
   
   roll.nc <- make_rolling_mean_netcdf_file(file.name,roll.name,var.name,tmp.dir) 

   lon <- ncvar_get(gcm.nc,'lon')
   lat <- ncvar_get(gcm.nc,'lat')
   n.lon <- length(lon)
   n.lat <- length(lat)

   for (j in 1:n.lat) {
      lat.ix <- j
      ltm <- proc.time()

      print(paste0('Latitude: ',j,' of ',n.lat))
      gcm.subset <- ncvar_get(gcm.nc,var.name,start=c(1,j,1),count=c(-1,1,-1))
      gcm.converted <- gcm.subset ##gcm_unit_conversion(var.name,gcm.subset,gcm.nc)
      rm(gcm.subset)

      flag <- is.na(gcm.converted[,1])
      n.time <- dim(gcm.converted)[2]  
      day.list <- vector(mode='list',length=n.lon)
      day.list <- lapply(seq_len(nrow(gcm.converted)), function(k) gcm.converted[k,])
      rm(gcm.converted)

      rolling_mean_of_days(var.name,roll.nc,lat.ix,n.lon,n.time,flag,
                           day.list)
      rm(day.list)
   }
   nc_close(roll.nc)
   nc_close(gcm.nc)
}

##--------------------------------------------------------------


##*************************************************************************
##*************************************************************************

##morph.name <- 'alpha_tasmax_tasmin' ##To match the file name morphing factor
##var.name <- 'alpha_tas'

##morph.name <- 'delta_tas' ##To match the file name morphing factor
##var.name <- 'tas'


past.int <- '1998-2014'

testing <- FALSE
if (testing) {
   gcm <- 'inmcm4'
   proj <- '2011-2040'
   type <- 'GCM'
   varname <- 'dewpoint'
   tmpdir <- '/local_temp/ssobie/morph'
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
   }
}

var.name <- varname ##'rhs'

morph.name <- switch(var.name,
                     rsds='alpha_rsds',
                     rhs='alpha_rhs',
                     clt='alpha_clt',
                     psl='alpha_psl',
                     wspd='alpha_wspd',
                     dewpoint='delta_dewpoint',
                     alpha_dewpoint='alpha_dewpoint')
##type <- 'GCM'
proj.int <- proj ##'2041-2070'

##-------------------------------------------------------------------------

##For BCCAQ2 files
if (type=='BCCAQ2') {
   epw.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/epw_factors/bccaq2_',past.int,'/')
   write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/epw_factors/rolled_bccaq2_',past.int,'/')
}

##For GCM files
if(type=='GCM') {
   ##epw.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/epw_factors/gcm_',past.int,'/')
   ##write.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/epw_factors/rolled_gcm_',past.int,'/')
   epw.dir <- paste0('/storage/data/climate/downscale/CMIP5/epw_factors/gcm_',past.int,'/')
   write.dir <- paste0('/storage/data/climate/downscale/CMIP5/epw_factors/rolled_gcm_',past.int,'/')
}

if (!file.exists(write.dir)) {
   dir.create(write.dir,recursive=TRUE)
}

tmp.dir <- paste0(tmpdir,'/',morph.name,'_',gcm,'_',proj.int,'/')
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

gcm.files <- list.files(path=epw.dir,pattern=gcm)
int.files <- gcm.files[grep(paste0(past.int,'_',proj.int),gcm.files)]
var.file <- int.files[grep(morph.name,int.files)]
if (type=='GCM') { roll.file <- gsub(gcm,paste0(gcm,'_GCM'),paste0('roll21_',var.file)) }
if (type=='BCCAQ2') { roll.file <- gsub(gcm,paste0(gcm,'_GCM'),paste0('roll21_',var.file))}

file.copy(from=paste0(epw.dir,var.file),to=tmp.dir,overwrite=T)

if (type=='GCM') {
   regrid.file <- gsub(gcm,paste0(gcm,'_at_BCCAQ2_Res_'),var.file)
   regrid_gcm_to_bccaq2_grid(var.file,regrid.file,tmp.dir)
   create_rolling_means(regrid.file,roll.file,var.name,tmp.dir)
}

if (type=='BCCAQ2') {
   create_rolling_means(var.file,roll.file,var.name,tmp.dir)
}

file.copy(from=paste0(tmp.dir,roll.file),to=write.dir,overwrite=TRUE)
