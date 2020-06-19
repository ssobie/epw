##Script to plot the EPW file series

library(ncdf4)
library(PCICt)

source('/storage/home/ssobie/code/repos/epw/read.write.epw.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/epw.precomputed.belcher.functions.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/epw.prism.offset.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/epw.support.functions.r',chdir=T)
##source('/storage/home/ssobie/code/repos/epw/epw.projected.stats.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/epw.with.about.tab.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/figures.for.summary.tables.r',chdir=T)

##------------------------------------------------------------------------------
##Read specified cell from the precomputed morphing values

read_cell <- function(var.name,lonc,latc,input.file,read.dir) {

  if (length(input.file) != 1) {
    print(read.dir)
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
   gcm.base.tas  <- calc_gcm_tas_stats(coords,scenario,gcm.list,'1971-2000',base.dir)

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
              base=gcm.base.tas,
              years=tmy.years,
              figs=figs)

   return(rv)
}



##------------------------------------------------------------------------------
