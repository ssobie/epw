##Script to produce table of statistics with the EPW files and GCM-PRISM projections

library(ncdf4)

source('/storage/home/ssobie/code/repos/crcm5/read.write.epw.r',chdir=T)
source('/storage/home/ssobie/code/repos/building_code/bc.building.code.fxns.r',chdir=T)
source('/storage/home/ssobie/code/repos/crcm5/epw.stats.formatted.table.functions.r',chdir=T)

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
##Return the lon/lat coordinates for the supplied EPW file

get_epw_coordinates <- function(epw.dir,epw.file) {

   epw <- read.epw.file(epw.dir,epw.file)
   epw.header <- epw$header
   epw.first <- strsplit(epw.header[1],',')[[1]]
   lon <- as.numeric(epw.first[8]) ##Fixed location
   lat <- as.numeric(epw.first[7])
   if (lon < -180 | lon > 180) {
      stop('Ill defined longitude coordinate')
   }
   if (lat < 40 | lat > 90) {
      stop('Ill defined latitude coordinate')
   }

   rv <- c(lon,lat)
   return(rv)
}

##----------------------------------------------------------------------- 

check_for_gcm_data <- function(lonc,latc,gcm.dir,gcm,scenario) {

  tasmax.files <- list.files(path=gcm.dir,pattern="tasmax_average_annual_climatology")
  tasmax.file <- tasmax.files[grep('1971-2000',tasmax.files)]
  nc <- nc_open(paste0(gcm.dir,tasmax.file))
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  lon.ix <- which.min(abs(lonc-lon))
  lat.ix <- which.min(abs(latc-lat))
  data.raw <- ncvar_get(nc,'tasmax',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  rv <- sum(is.na(data.raw)) == length(data.raw)
  nc_close(nc)
  return(rv)
}

##----------------------------------------------------------------------- 

##Degree Day Values
dd <- function(tmean,tbase) {
  g <- tmean - tbase
  days <- sum(g[g>0], na.rm=T)
  return(round(days))
}

##----------------------------------------------------------------------- 

specific_humidity <- function(dwpt,pas) {
  dwpt <- dwpt + 273
  vape.press <- sat.vape * exp( (lh.vape/R.vapor) * (1/T.zero - 1/dwpt))
  sp.hum <- vape.press * epsilon / pas
}

##----------------------------------------------------------------------- 

retrieve_closest_cell <- function(lonc,latc,var.name,file.dir,file.name) {
  if (grepl('\\.',var.name)) {
     var.name <- strsplit(var.name,'\\.')[[1]][1]
  }
  nc <- nc_open(paste0(file.dir,file.name))
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  lon.ix <- which.min(abs(lonc-lon))
  lat.ix <- which.min(abs(latc-lat))
  cell <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  nc_close(nc)
  return(cell)
}

##------------------------------------------------------------------------------ 

get_var_clims <- function(var.info,var.dir,coords) {
  
  intervals <- c('1971-2000','2011-2040','2041-2070','2071-2100')
  var.cell <- rep(NA,length(intervals))
  var.ref <- var.info$name
  if (grepl('\\.',var.info$name)) {
     var.split <- strsplit(var.info$name,'\\.')[[1]]
     var.ref <- paste0(var.split[1],'_',var.split[2])
  }

  for (i in seq_along(intervals)) {
    var.files <- list.files(path=var.dir,pattern=var.ref)
    var.file <- var.files[grep(intervals[i],var.files)]
    var.cell[i] <- retrieve_closest_cell(coords[1],coords[2],var.info$name,var.dir,var.file)
  }
  if (grepl('wetbulb',var.info$name)) {
     var.cell <- var.cell - 273 
  }
  return(var.cell)
}

##------------------------------------------------------------------------------

calc_cwec_values <- function(present.epw.file,epw.dir,var.names) {

  pef.split <- strsplit(present.epw.file,'_')[[1]]
  new.epw.file <- paste0('STATS_CAN_BC_',pef.split[3],'_CWEC.epw')
  epw.present <- read.epw.file(epw.dir,present.epw.file)

  tas.ix <- get_field_index('dry_bulb_temperature')
  epw.tas <- epw.present$data[,tas.ix]
  dates <- as.Date(paste('1999',sprintf('%02d',epw.present$data[,2]),sprintf('%02d',epw.present$data[,3]),sep='-'))

  fac <- as.factor(format(dates,'%Y-%m-%d'))
  epw.tas.daily <- tapply(epw.tas,fac,mean)

  epw.hdd <- dd(-1*epw.tas.daily,-18)
  epw.cdd <- dd(epw.tas.daily,18)

  epw.txx <- max(epw.tas)
  epw.tnn <- min(epw.tas)

  epw.975 <- quantile(epw.tas,0.975,names=F)
  epw.001 <- quantile(epw.tas,0.01,names=F)

  dwpt.ix <- get_field_index('dew_point_temperature')
  epw.dwpt <- epw.present$data[,dwpt.ix]
  pas.ix <- get_field_index('atmospheric_station_pressure')
  epw.pas <- epw.present$data[,pas.ix]

  epw.sph <- specific_humidity(epw.dwpt,epw.pas)
  epw.twb <- temp.wet.bulb(epw.tas+273,epw.dwpt,epw.pas,epw.sph) - 273

  twb.975 <- round(quantile(epw.twb,0.975,names=F),1)
  twb.001 <- round(quantile(epw.twb,0.01,names=F),1)

  cwec.col <- list(epw.hdd,
                   epw.cdd,
                   epw.txx,
                   epw.tnn,
                   epw.975,
                   epw.001,
                   twb.975,
                   twb.001)
  names(cwec.col) <- var.names 
  return(cwec.col)
}

##------------------------------------------------------------------------------
get_gcm_directory <- function(var.info,gcm,base.dir) {

  if (grepl('(hdd|cdd|fdd|gdd)',var.info$name))
    gcm.dir <- paste0(base.dir,gcm,'/rcp85/degree_days/climatologies/')

  if (grepl('ETCCDI',var.info$name))
    gcm.dir <- paste0(base.dir,gcm,'/rcp85/climdex/climatologies/')
 
  if (grepl('quantile',var.info$name))
    gcm.dir <- paste0(base.dir,gcm,'/rcp85/annual_quantiles/climatologies/')

  if (grepl('wetbulb',var.info$name))
    gcm.dir <- paste0('/storage/data/climate/downscale/CMIP5/building_code/',gcm,'/climatologies/')  
  return(gcm.dir)
}

##------------------------------------------------------------------------------
calc_gcm_stats <- function(var.info,coords,scenario,model.list,
                           check.dir,base.dir) {

  ##GCM Component
  flag <- check_for_gcm_data(coords[1],coords[2],check.dir,'ACCESS1-0',scenario)

  vals <- matrix(NA,nrow=length(model.list),ncol=4)

  for (g in seq_along(model.list)) {
    gcm <- model.list[g]
    ##print(gcm)
    gcm.dir <- get_gcm_directory(var.info,gcm,base.dir)

    if (grepl('wetbulb',var.info$name) & (gcm == 'CCSM4' | gcm == 'MPI-ESM-LR')) {
      vals[g,] <- vals[g,] <- rep(NA,4)
    } else {
      vals[g,] <- get_var_clims(var.info,gcm.dir,coords)
    }
  }
  rv <- get.round.val(var.info$name)
  val.row <- make_table_row(vals,rv)
  return(val.row)
}

##---==-=
old_calc_gcm_stats <- function(coords,scenario,model.list,
                           check.dir,base.dir) {

  ##GCM Component
  flag <- check_for_gcm_data(coords[1],coords[2],check.dir,'ACCESS1-0',scenario)

  hdd.vals <- matrix(NA,nrow=length(model.list),ncol=4)
  cdd.vals <- matrix(NA,nrow=length(model.list),ncol=4)
  txx.vals <- matrix(NA,nrow=length(model.list),ncol=4)
  tnn.vals <- matrix(NA,nrow=length(model.list),ncol=4)
  tx975.vals <- matrix(NA,nrow=length(model.list),ncol=4)
  tn010.vals <- matrix(NA,nrow=length(model.list),ncol=4)

  wb975.vals <- matrix(NA,nrow=length(model.list),ncol=4)
  wb010.vals <- matrix(NA,nrow=length(model.list),ncol=4)

  for (g in seq_along(model.list)) {
    gcm <- model.list[g]
    print(gcm)
    dd.dir <- paste0(base.dir,gcm,'/rcp85/degree_days/climatologies/')
    hdd.vals[g,] <- get_var_clims('hdd','annual',dd.dir,coords)
    cdd.vals[g,] <- get_var_clims('cdd','annual',dd.dir,coords)

    clim.dir <- paste0(base.dir,gcm,'/rcp85/climdex/climatologies/')
    txx.vals[g,] <- get_var_clims('txxETCCDI','annual',clim.dir,coords)
    tnn.vals[g,] <- get_var_clims('tnnETCCDI','annual',clim.dir,coords)

    quant.dir <- paste0(base.dir,gcm,'/rcp85/annual_quantiles/climatologies/')
    tx975.vals[g,] <- get_var_clims('tasmax','annual_quantile_975',quant.dir,coords)
    tn010.vals[g,] <- get_var_clims('tasmin','annual_quantile_010',quant.dir,coords)

    if (gcm == 'CCSM4' | gcm == 'MPI-ESM-LR') {
      wb975.vals[g,] <- wb010.vals[g,] <- rep(NA,4)
    } else {
      wbt.dir <- paste0('/storage/data/climate/downscale/CMIP5/building_code/',gcm,'/climatologies/')  
      wb975.vals[g,] <- get_var_clims('wetbulb','annual_quantile_975',wbt.dir,coords)-273
      wb010.vals[g,] <- get_var_clims('wetbulb','annual_quantile_010',wbt.dir,coords)-273
    }
  }
  hdd.row <- make_table_row(hdd.vals,0)
  cdd.row <- make_table_row(cdd.vals,0)
  txx.row <- make_table_row(txx.vals,1)
  tnn.row <- make_table_row(tnn.vals,1)
  tx975.row <- make_table_row(tx975.vals,1)
  tn010.row <- make_table_row(tn010.vals,1)
  wb975.row <- make_table_row(wb975.vals,1) 
  wb010.row <- make_table_row(wb010.vals,1)

  all.rows <- rbind(hdd.row,cdd.row,txx.row,tnn.row,
                    tx975.row,tn010.row,wb975.row,wb010.row)
  return(all.rows)
}

##------------------------------------------------------------------------------

make_table_row <- function(data.vals,rv) {
   
  val.past <- mean(data.vals[,1],na.rm=T)
  anoms.2020s <- data.vals[,2]-data.vals[,1]
  anoms.2050s <- data.vals[,3]-data.vals[,1]
  anoms.2080s <- data.vals[,4]-data.vals[,1]
  rv <- 0
  val.row <- c(round(val.past,rv),
               round(mean(anoms.2020s,na.rm=T),rv),
               paste0(' (',round(quantile(anoms.2020s,0.1,na.rm=T),rv),' to ',round(quantile(anoms.2020s,0.9,na.rm=T),rv),')'),
               round(mean(data.vals[,2],na.rm=T),rv),
               paste0(' (',round(quantile(data.vals[,2],0.1,na.rm=T),rv),' to ',round(quantile(data.vals[,2],0.9,na.rm=T),rv),')'),
               round(mean(anoms.2050s,na.rm=T),rv),
               paste0(' (',round(quantile(anoms.2050s,0.1,na.rm=T),rv),' to ',round(quantile(anoms.2050s,0.9,na.rm=T),rv),')'), 
               round(mean(data.vals[,3],na.rm=T),rv),
               paste0(' (',round(quantile(data.vals[,3],0.1,na.rm=T),rv),' to ',round(quantile(data.vals[,3],0.9,na.rm=T),rv),')'),
               round(mean(anoms.2080s,na.rm=T),rv),
               paste0(' (',round(quantile(anoms.2080s,0.1,na.rm=T),rv),' to ',round(quantile(anoms.2080s,0.9,na.rm=T),rv),')'),
               round(mean(data.vals[,4],na.rm=T),rv),
               paste0(' (',round(quantile(data.vals[,4],0.1,na.rm=T),rv),' to ',round(quantile(data.vals[,4],0.9,na.rm=T),rv),')'))
  return(val.row)
}

##------------------------------------------------------------------------------
  ##Formatted Table
make_formated_stats_table <- function(var.list,cwec.site,gcm.site.list,cwec.offset,gcm.offset.list) {

  sorted.vars <- filter.input.variables(var.list)
  row.locs <- find.row.locations(sorted.vars)
  
  wb <- createWorkbook()
  addWorksheet(wb, "Requested Site (Offset)")
  setColWidths(wb, sheet = 1, cols = 1:15, widths = 14) ##Set fixed width
  create.frozen.top.pane(wb,sheet=1)
  create.title.panes(wb,sheet=1,var.name='tas',start.row=row.locs$tas[[1]][1])
  write.variables(wb,sheet=1,sorted.vars$tas,row.locs$tas,'tas',cwec.offset,gcm.offset.list)
  freezePane(wb,sheet=1,firstActiveCol=2,firstActiveRow=3)

  addWorksheet(wb, "Weather File Site")
  setColWidths(wb, sheet = 2, cols = 1:15, widths = 14) ##Set fixed width
  create.frozen.top.pane(wb,sheet=2)
  create.title.panes(wb,sheet=2,var.name='tas',start.row=row.locs$tas[[1]][1])
  write.variables(wb,sheet=2,sorted.vars$tas,row.locs$tas,'tas',cwec.site,gcm.site.list)
  freezePane(wb,sheet=2,firstActiveCol=2,firstActiveRow=3)
  saveWorkbook(wb, paste0('/storage/home/ssobie/general/assessment_tables/downtown_victoria_variable_table_rcp85.xlsx'), overwrite = TRUE)
}



##----------------------------------------------------------------------- 

##***********************************************************************************
##------------------------------------------------------------------------------------


epw.dir <- '/storage/data/projects/rci/weather_files/wx_files/'
offset.dir <- '/storage/data/projects/rci/weather_files/wx_files/offsets/'
write.dir <- '/storage/data/projects/rci/weather_files/wx_files/stats/'

full.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')
sub.list <- c('ACCESS1-0','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MRI-CGCM3')

var.list <- list(list(name='hdd',type='annual',title='HDD'),     
                 list(name='cdd',type='annual',title='CDD'),
                 list(name='txxETCCDI',type='annual',title='TXX'),
                 list(name='tnnETCCDI',type='annual',title='TNN'),
                 list(name='tasmax.annual_quantile_975',type='annual',title='TX 97.5%'),
                 list(name='tasmin.annual_quantile_010',type='annual',title='TN 1.0%'),
                 list(name='wetbulb.annual_quantile_975',type='annual',title='WB 97.5%'),
                 list(name='wetbulb.annual_quantile_010',type='annual',title='WB 1.0%'))
var.names <- unlist(lapply(var.list,function(x){return(x$name)}))
 
base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
check.dir <- paste0(base.dir,'ACCESS1-0/rcp85/annual/climatologies/')

tmp.dir <- '/local_temp/ssobie/epw/'
if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=TRUE)
}

scenario <- 'rcp85'
offset.epw.file <- 'CAN_BC_DowntownVictoria_offset_from_VICTORIA-GONZALES-CS_1018611_CWEC.epw'
present.epw.file <- 'CAN_BC_VICTORIA-GONZALES-CS_1018611_CWEC.epw'

cwec.offset <- calc_cwec_values(offset.epw.file,offset.dir,var.names)
coords <- get_epw_coordinates(offset.dir,offset.epw.file)
gcm.offset.list <- lapply(var.list,calc_gcm_stats,coords,scenario,full.list,check.dir,base.dir)
names(gcm.offset.list) <- var.names

cwec.site <- calc_cwec_values(present.epw.file,epw.dir,var.names)
coords <- get_epw_coordinates(epw.dir,present.epw.file)
gcm.site.list <- lapply(var.list,calc_gcm_stats,coords,scenario,full.list,check.dir,base.dir)
names(gcm.site.list) <- var.names


##Send to xlsx means of writing
make_formated_stats_table(var.list,cwec.offset,gcm.offset.list,
                                   cwec.site,gcm.site.list)

