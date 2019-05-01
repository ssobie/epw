##Script to generate shift and stretch factors for temperature
##from the BCCAQ2-TPS downscaled fields

library(ncdf4)
library(PCICt)

source('/storage/home/ssobie/code/repos/crcm5/epw.belcher.functions.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##------------------------------------------------------------------------------
create_factor_file <- function(nc,var.name,var.units,
                               input,
                               tmp.dir,write.file) {

  time.calendar <- '365_day'
  time.units <- 'days since 1950-01-01 00:00:00'
  full.dates <- seq(from=as.Date('1950-01-01'),by='day',to=as.Date('1995-12-31'))
  fd.leap.flag <- grep('02-29',full.dates)
  noleap.dates <- full.dates[-fd.leap.flag]
  dates.95 <- which(grepl('1995',noleap.dates))-1

  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  lon.atts <- ncatt_get(nc,'lon')
  lat.atts <- ncatt_get(nc,'lat')
  n.lon <- length(lon)
  n.lat <- length(lat)

  global.atts <- ncatt_get(nc,varid=0)
  ##--------------------------------------------------------------
  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, dates.95,
                      unlim=FALSE, calendar=time.calendar)
  var.geog <- ncvar_def(var.name, units=var.units, dim=list(x.geog, y.geog, t.geog),
                        missval=-32768)
  file.nc <- nc_create(paste(tmp.dir,write.file,sep=''), var.geog)
  ##File Attributes
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(file.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  ##Time attributes
  ncatt_put(file.nc,varid='time',attname='units',attval=time.units)
  ncatt_put(file.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(file.nc,varid='time',attname='calendar',attval=time.calendar)

  lon.names <- names(lon.atts)
  for (j in 1:length(lon.atts))
    ncatt_put(file.nc,varid='lon',attname=lon.names[j],attval=lon.atts[[j]])

  lat.names <- names(lat.atts)
  for (j in 1:length(lat.atts))
    ncatt_put(file.nc,varid='lat',attname=lat.names[j],attval=lat.atts[[j]])

  ##Climdex Attributes
  ncatt_put(file.nc,varid=var.name,attname='units',attval=var.units)
  ncatt_put(file.nc,varid=var.name,attname='_FillValue',attval=-32768)
  ncatt_put(file.nc,varid=var.name,attname='standard_name',attval=var.name)
  ncatt_put(file.nc,varid=var.name,attname='long_name',attval=var.name)

  ncatt_put(file.nc,varid=0,attname='history',attval='')
  ncvar_put(file.nc,var.name,input)
  ncvar_put(file.nc,'lon',lon)
  ncvar_put(file.nc,'lat',lat)

  nc_close(file.nc)

}

##------------------------------------------------------------------------------
make_average_series <- function(series,fac,agg.fxn) {
                   
  agg.daily <- tapply(series,fac,agg.fxn)
  return(agg.daily)
}

##----------------------------------------------------------------------------------------------
##Compute the 365 day averaging factor
make_averaging_factor <- function(interval,time,cal) {
   yrs <- as.numeric(strsplit(interval,'-')[[1]])
   full.dates <- seq(from=as.Date(paste0(yrs[1],'-01-01')),by='day',to=as.Date(paste0(yrs[2],'-12-31')))
   fd.leap.flag <- grep('02-29',full.dates)
   noleap.dates <- full.dates[-fd.leap.flag]
   ymax <- max(format(time,'%Y'))
   if (ymax=='2099') 
      yrs[2] <- 2099
   ylen <- yrs[2]-yrs[1]+1
   
   fac <- as.factor(rep(sprintf('%03d',1:365),ylen))
   time <- as.Date(format(time,'%Y-%m-%d'))
   if (cal=='365') {
     flag <- noleap.dates %in% time
     fac <- fac[flag]
   }  
   if (grepl('(gregorian|standard)',cal)) {
     leap.flag <- grep('02-29',time)
     noleap.time <- time[-leap.flag]
     flag <- noleap.dates %in% noleap.time
     fac <- fac[flag]
   }  
   return(fac)
}
##----------------------------------------------------------------------------------------------
##Find year indices

find_year_indices <- function(time,interval) {
   yrs <- strsplit(interval,'-')[[1]]
   years <- format(time,'%Y')
   st <- head(grep(yrs[1],years),1)
   en <- tail(grep(yrs[2],years),1)
   if (grepl('HadGEM',gcm) & yrs[2]=='2100') {
     en <- length(years)
   }
   cnt <- en-st+1
  rv <- list(st=st,en=en,cnt=cnt)
  return(rv)
}


##----------------------------------------------------------------------------------------------
separate_365_into_list <- function(nc,var.name,j,time,time.bnds,interval) {

    data.subset <- ncvar_get(nc,var.name,start=c(1,j,time.bnds$st),count=c(-1,1,time.bnds$cnt))
    data.list <- lapply(seq_len(nrow(data.subset)), function(k) data.subset[k,])
    rm(data.subset)
    return(data.list)
}

##----------------------------------------------------------------------------------------------
separate_gregorian_into_list <- function(nc,var.name,j,time,time.bnds,interval) {

    data.raw <- ncvar_get(nc,var.name,start=c(1,j,time.bnds$st),count=c(-1,1,time.bnds$cnt))
    time.series <- netcdf.calendar(nc)[time.bnds$st:time.bnds$en]
    feb.flag <- grep('-02-29',time.series)
    data.subset <- data.raw[,-feb.flag]
    data.list <- lapply(seq_len(nrow(data.subset)), function(k) data.subset[k,])
    rm(data.subset)
    rm(data.raw)
    return(data.list)
}

##----------------------------------------------------------------------------------------------
separate_360_into_list <- function(nc,var.name,j,time,time.bnds,interval) {

    data.raw <- ncvar_get(nc,var.name,start=c(1,j,time.bnds$st),count=c(-1,1,time.bnds$cnt))
    yrs <- as.numeric(strsplit(interval,'-')[[1]])
    ymax <- max(format(time,'%Y'))
    if (ymax=='2099' & yrs[1]==2071) 
       yrs[2] <- 2099
    ylen <- yrs[2]-yrs[1]+1    
    inx <- seq(30,365*ylen,73)
    insrt <- 1:(ylen*365) %in% inx
    data.fix <- matrix(NA,nrow=nrow(data.raw),ncol=ylen*365)

    data.fix[,!insrt] <- data.raw

    for (i in seq_along(inx)) {
      ix <- inx[i]
      data.fix[,ix] <- (data.fix[,ix-1] + data.fix[,ix+1])/2
    }
    data.list <- lapply(seq_len(nrow(data.fix)), function(k) data.fix[k,])

    rm(data.raw)
    return(data.list)
}

##----------------------------------------------------------------------------------------------

##Calculate morphing factors

daily_aggregate <- function(nc,var.name,
                            gcm,interval,
                            agg.fxn) {
   print(paste0(gcm, ' for ',interval))
   time <- netcdf.calendar(nc)
   n.lat <- nc$dim$lat$len ##Latitude Length
   n.lon <- nc$dim$lon$len ##Longitude Length
   n.time <- length(time)
   time.bnds <- find_year_indices(time,interval)
   cal <- attr(time,'cal')
   sub.fac <- make_averaging_factor(interval,time[time.bnds$st:time.bnds$en],cal) 


   sep_fxn <- separate_365_into_list
   if (grepl('360',cal)) {
      sep_fxn <- separate_360_into_list
   }
   if (grepl('(gregorian|standard)',cal)) {
      sep_fxn <- separate_gregorian_into_list
   }

   agg.array <- array(NA,c(n.lon,n.lat,365))

   ##Iterate along latitude indices
   for (j in 1:n.lat) {
      ##print(paste0('Latitude: ',j,' of ',n.lat))

      data.list <- sep_fxn(nc,var.name,j,time,time.bnds,interval)

      na.flag <- unlist(lapply(data.list,function(x){any(is.na(x))}))
      data.sub.list <- data.list[!na.flag]

      data.test <- make_average_series(series=data.sub.list[[1]],
                                         fac=sub.fac,agg.fxn=agg.fxn)           

      data.result <- rep(list(rep(NA,length(levels(sub.fac)))),n.lon)

      ##Compute the aggregate values (default to daily)
      data.agg <- foreach(
                      data=data.sub.list,
                      .export=c('make_average_series','sub.fac','agg.fxn')
                    ) %dopar% {
                      objects <- make_average_series(series=data,
                                      fac=sub.fac,agg.fxn=agg.fxn)
                              }
      data.result[!na.flag] <- data.agg
      ncol <- length(data.result[[1]])
      data.matrix <- matrix(unlist(data.result),nrow=n.lon,ncol=ncol,byrow=TRUE)
      agg.array[,j,] <- data.matrix
   }
   return(agg.array)
}

##----------------------------------------------------------------------------------------------

