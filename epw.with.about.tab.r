##Script to produce table of statistics with the EPW files and GCM-PRISM projections

source('/storage/home/ssobie/code/repos/building_code/bc.building.code.fxns.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/epw.stats.formatted.table.functions.r',chdir=T)

##------------------------------------------------------------------------------
##Make sure GCM values are available for this location
check_for_gcm_data <- function(lonc,latc,gcm.dir,gcm,scenario,past.int) {

  tasmax.files <- list.files(path=gcm.dir,pattern="tasmax_average_annual_climatology")
  tasmax.file <- tasmax.files[grep(past.int,tasmax.files)]
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

##Midday Hours: Number of hours between 8am and 4pm where dry bulb
##is between 13C and 21C

midday <- function(t) {
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

get_var_clims <- function(var.info,var.dir,coords,past.int) {
  
  intervals <- c(past.int,'2011-2040','2041-2070','2071-2100')
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
  hours <- sprintf('%02d',epw.present$data[,4])
  fac <- as.factor(format(dates,'%Y-%m-%d'))
  day.dates <- as.Date(levels(fac))
  mon.fac <- as.factor(format(day.dates,'%m'))
  seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
  seas.fac <- factor(seasons[mon.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))

  ##Mid-day temperature index
  epw.tas.thresh <- epw.tas >= 13 & epw.tas <= 21
  mid.day.hours <- hours %in% 8:16 ##Hours between 8:00am and 4:00pm
  epw.tas.thresh[!mid.day.hours] <- FALSE
  epw.tas.midday <- tapply(epw.tas.thresh,fac,sum)

  epw.tas.daily <- tapply(epw.tas,fac,mean)

  epw.hdd <- dd(-1*epw.tas.daily,-18)
  epw.cdd <- dd(epw.tas.daily,18)
  epw.cdd.10 <- dd(epw.tas.daily,10)

  epw.txx <- max(epw.tas)
  epw.tnn <- min(epw.tas)

  epw.975 <- quantile(epw.tas,0.975,names=F)
  epw.010 <- quantile(epw.tas,0.01,names=F)
  epw.025 <- quantile(epw.tas,0.025,names=F)

  dwpt.ix <- get_field_index('dew_point_temperature')
  epw.dwpt <- epw.present$data[,dwpt.ix]
  pas.ix <- get_field_index('atmospheric_station_pressure')
  epw.pas <- epw.present$data[,pas.ix]

  epw.sph <- specific_humidity(epw.dwpt,epw.pas)
  epw.twb <- temp.wet.bulb(epw.tas+273,epw.dwpt,epw.pas,epw.sph) - 273

  if (any(is.na(epw.twb))) { browser()}
  twb.975 <- round(quantile(epw.twb,0.975,names=F),1)
  twb.010 <- round(quantile(epw.twb,0.01,names=F),1)
  twb.025 <- round(quantile(epw.twb,0.025,names=F),1)

  epw.tas.daily <- tapply(epw.tas,fac,mean)
  
  tas.monthly <- tapply(epw.tas.daily,mon.fac,mean)
  tas.seasonal <- tapply(epw.tas.daily,seas.fac,mean)
  tas.annual <- mean(epw.tas.daily)

  cwec.col <- c(epw.hdd,
                epw.cdd,
                epw.cdd.10,
                epw.txx,
                epw.tnn,
                epw.975,
                epw.010,
                epw.025,
                twb.975,
                twb.010,
                twb.025,
                tas.monthly,
                tas.seasonal,
                tas.annual)

  names(cwec.col) <- c("hdd",
                       "cdd",
                       "cdd_10", 
                       "txxETCCDI",
                       "tnnETCCDI",
                       "tasmax.annual_quantile_975",
                       "tasmin.annual_quantile_010",
                       "tasmin.annual_quantile_025",
                       "wetbulb.annual_quantile_975",
                       "wetbulb.annual_quantile_010",
                       "wetbulb.annual_quantile_025",
                       'tas_jan','tas_feb','tas_mar',
                       'tas_apr','tas_may','tas_jun',
                       'tas_jul','tas_aug','tas_sep',
                       'tas_oct','tas_nov','tas_dec',
                       'tas_win','tas_spr','tas_sum',
                       'tas_fal','tas_ann')                      
  re.ix <- sapply(var.names,function(x,y){which(y %in% x)},names(cwec.col))
  return(cwec.col[re.ix])
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

calc_gcm_stats <- function(var.info,coords,scenario,model.list,past.int,
                           check.dir,base.dir) {

  ##GCM Component
  flag <- check_for_gcm_data(coords[1],coords[2],check.dir,'ACCESS1-0',scenario,past.int)

  vals <- matrix(NA,nrow=length(model.list),ncol=4)

  for (g in seq_along(model.list)) {
    gcm <- model.list[g]
    ##print(gcm)
    gcm.dir <- get_gcm_directory(var.info,gcm,base.dir)

    if (grepl('wetbulb',var.info$name) & (gcm == 'CCSM4' | gcm == 'MPI-ESM-LR')) {
      vals[g,] <- vals[g,] <- rep(NA,4)
    } else {
      vals[g,] <- get_var_clims(var.info,gcm.dir,coords,past.int)
    }
  }
  rv <- get.round.val(var.info$name)
  val.row <- make_table_row(vals,rv)

  return(val.row)
}

##------------------------------------------------------------------------------

calc_gcm_tas_stats <- function(coords,scenario,model.list,interval,
                               base.dir) {

  intervals <- interval ##'1998-2014' ##c('1971-2000','1998-2014')

  ##GCM Component
  vals <- array(NA,c(length(model.list),length(intervals),17))
  
  for (g in seq_along(model.list)) {
     gcm <- model.list[g]
     print(gcm)
     clim.dir <- mon.dir <- ann.dir <- '/storage/data/climate/downscale/BCCAQ2/bccaqv2_climatologies/'
     gcm.files <- list.files(path=clim.dir,pattern=gcm)
     ##seas.dir <- paste0(base.dir,gcm,'/rcp85/seasonal/climatologies/')
     ##mon.dir <- paste0(base.dir,gcm,'/rcp85/monthly/climatologies/')
     ##ann.dir <- paste0(base.dir,gcm,'/rcp85/annual/climatologies/')

     ##print(gcm)
     for (j in seq_along(intervals)) {
        print(intervals[j])
        mon.files <- gcm.files[grep('tas_monthly',gcm.files)]
        mon.file <- mon.files[grep(intervals[j],mon.files)]
        vals[g,j,1:12] <- retrieve_closest_cell(coords[1],coords[2],'tas',clim.dir,mon.file)
        seas.files <- gcm.files[grep('tas_seasonal',gcm.files)]
        seas.file <- seas.files[grep(intervals[j],seas.files)]
        vals[g,j,13:16] <- retrieve_closest_cell(coords[1],coords[2],'tas',clim.dir,seas.file)
        ###ann.files <- list.files(path=ann.dir,pattern='tas_average_annual')
        ann.files <- gcm.files[grep('tas_annual',gcm.files)]
        ann.file <- ann.files[grep(intervals[j],ann.files)]
        vals[g,j,17] <- retrieve_closest_cell(coords[1],coords[2],'tas',clim.dir,ann.file)
     }
  }
  rv <- get.round.val('tas')
  ##val.row <- matrix('A',nrow=17,ncol=13)
  ##for (k in 1:17) {
  ##   val.row[k,] <- make_table_row(vals[,,k],rv)
  ##}
  vals.ens <- round(t(apply(vals,c(2,3),mean)),rv)

  return(vals.ens)
}


##------------------------------------------------------------------------------

make_table_row <- function(data.vals,rv) {
   
  val.past <- mean(data.vals[,1],na.rm=T)
  anoms.2020s <- data.vals[,2]-data.vals[,1]
  anoms.2050s <- data.vals[,3]-data.vals[,1]
  anoms.2080s <- data.vals[,4]-data.vals[,1]
  ##rv <- 0
  val.row <- c(round(val.past,rv),

               round(mean(data.vals[,2],na.rm=T),rv),
               paste0(' (',round(quantile(data.vals[,2],0.1,na.rm=T),rv),' to ',round(quantile(data.vals[,2],0.9,na.rm=T),rv),')'),
               round(mean(anoms.2020s,na.rm=T),rv),
               paste0(' (',round(quantile(anoms.2020s,0.1,na.rm=T),rv),' to ',round(quantile(anoms.2020s,0.9,na.rm=T),rv),')'),

               round(mean(data.vals[,3],na.rm=T),rv),
               paste0(' (',round(quantile(data.vals[,3],0.1,na.rm=T),rv),' to ',round(quantile(data.vals[,3],0.9,na.rm=T),rv),')'),
               round(mean(anoms.2050s,na.rm=T),rv),
               paste0(' (',round(quantile(anoms.2050s,0.1,na.rm=T),rv),' to ',round(quantile(anoms.2050s,0.9,na.rm=T),rv),')'), 

               round(mean(data.vals[,4],na.rm=T),rv),
               paste0(' (',round(quantile(data.vals[,4],0.1,na.rm=T),rv),' to ',round(quantile(data.vals[,4],0.9,na.rm=T),rv),')'),
               round(mean(anoms.2080s,na.rm=T),rv),
               paste0(' (',round(quantile(anoms.2080s,0.1,na.rm=T),rv),' to ',round(quantile(anoms.2080s,0.9,na.rm=T),rv),')'))
  return(val.row)
}

##------------------------------------------------------------------------------
  ##Formatted Table
make_formated_stats_table <- function(nearest,site,var.list,sheets.closest,sheets.offset,
                                      file.version,
                                      method,rlen,write.dir) {
               
  sorted.vars <- filter.input.variables(var.list)
  row.locs <- find.row.locations(sorted.vars)
  desc.start <- tail(unlist(row.locs$derived),1)+3
  wb <- createWorkbook()

  addWorksheet(wb, "Weather File Station")
  sheet <- 1
  setColWidths(wb, sheet = sheet, cols = 1:14, widths = c(30,rep(14,13))) ##Set fixed width
  create.frozen.top.pane(wb,sheet=sheet,site=nearest)
  ##create.title.panes(wb,sheet=sheet,var.name='tas',start.row=row.locs$tas[[1]][1])
  write_variables(wb,sheet=sheet,sorted.vars$derived,row.locs$derived,'derived',sheets.closest$cwec,
                     data.cols=2:14,highlights=c(5,6,9,10,13,14))
  textstyle <- createStyle(fgFill = 'white', halign = "LEFT",
                           fontColour = "black",bgFill='white')
  titlestyle <- createStyle(fgFill = 'white', halign = "LEFT",textDecoration = "Bold",
                           fontColour = "black",bgFill='white')
  linkstyle <- createStyle(fgFill = 'white', halign = "LEFT",textDecoration = "Underline",
                           fontColour = "blue",bgFill='white')

  ##File Description
  file.desc <- c("This tab summarizes the data from the ",
                 paste0(nearest," CWEC2016 Weather File"),
                 "A description of each variable can be seen by hovering over the variable name in ",
                 "Column A. The 'Past (TMY)' column summarizes data from the CWEC2016 weather file",
                 "while the remaining columns summarize future-shifted values. More details can be ",
                 "found in the 'About This File' tab.")
  writeData(wb, sheet=sheet, 'Data Description', startRow = desc.start-1, startCol = 2, headerStyle = titlestyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,titlestyle,rows=tail(unlist(row.locs$derived),1)+2,cols=2:7,gridExpand=FALSE,stack=FALSE)
  writeData(wb, sheet=sheet, file.desc, startRow = desc.start, startCol = 2, headerStyle = textstyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,textstyle,rows=seq(desc.start,by=1,length.out=length(file.desc)),cols=2:7,gridExpand=TRUE,stack=TRUE)
  ##-----------

  freezePane(wb,sheet=sheet,firstActiveCol=2,firstActiveRow=3)

  ##Adjusted Weather file Station Table
  if (!is.null(sheets.offset)) {

     offset.description <- c("This tab summarizes past and future versions of the same weather file as the",
                             "previous tab, but with temperatures adjusted using the difference between",
                             paste0("the historical climatology of ",site),
                             paste0("and ",nearest," according to average PRISM temperatures. More details"),
                             "can be found in the 'About This File' tab.")

     addWorksheet(wb, "Adjusted Weather File Station")
     sheet <- 2
     setColWidths(wb, sheet = sheet, cols = 1:14, widths = c(30,rep(14,13))) ##Set fixed width
     create.frozen.top.pane(wb,sheet=sheet,site=site)
     ##create.title.panes(wb,sheet=sheet,var.name='tas',start.row=row.locs$tas[[1]][1])
     write_variables(wb,sheet=sheet,sorted.vars$derived,row.locs$derived,'derived',sheets.offset$cwec,
                        data.cols=2:14,highlights=c(5,6,9,10,13,14))
     writeData(wb, sheet=sheet, 'Data Description', startRow = tail(unlist(row.locs$derived),1)+2, startCol = 2, headerStyle = titlestyle,
            colNames=FALSE)
     addStyle(wb,sheet=sheet,titlestyle,rows=tail(unlist(row.locs$derived),1)+2,cols=2:7,gridExpand=FALSE,stack=FALSE)
     writeData(wb, sheet=sheet, offset.description, startRow = tail(unlist(row.locs$derived),1)+3, startCol = 2, headerStyle = textstyle,
            colNames=FALSE)
     addStyle(wb,sheet=sheet,textstyle,rows=seq(desc.start,by=1,length.out=length(offset.description)),cols=2:7,gridExpand=TRUE,stack=TRUE)
     freezePane(wb,sheet=sheet,firstActiveCol=2,firstActiveRow=3)

     sheet <- 3
     sheets.tas <- sheets.offset
  } else {
     sheet <- 2
     sheets.tas <- sheets.closest
  }

  ##TAS Climatologies Tab
  ##Insert the TMY Years and model climatologies into the TAS rows
  tas.flag <- unlist(lapply(var.list,function(x){grepl('tas_',x$name)}))  
  tas.names <- var.list[tas.flag]
  years <- c(sheets.tas$years,rep('NA',5))
  for (t in 1:sum(tas.flag)) {
     var.name <- tas.names[[t]]$name

     line <- sheets.tas$cwec[[var.name]]
     sheets.tas$cwec[[var.name]] <- append(line,values=c(years[t],sheets.tas$base[t,],sheets.tas$model[t,]),after=1)
  }

  addWorksheet(wb, "Temperature Climatologies") ##"Requested Site (CWEC)")
  setColWidths(wb, sheet = sheet, cols = 1:17, widths = c(30,rep(14,16))) ##Set fixed width
  create_temp_climatologies_top_pane(wb,sheet=sheet,site=site)
  ##create.title.panes(wb,sheet=sheet,var.name='tas',start.row=row.locs$tas[[1]][1])
  write_variables(wb,sheet=sheet,sorted.vars$tas,row.locs$tas,'tas',sheets.tas$cwec,
                              data.cols=2:17,highlights=c(8,9,12,13,16,17))

  ##Add image of temperature boxplots below table values
  tas.fig.description <- c("The plot to the left displays temperature climatologies",
                           "from the above table. The blue 'Past' lines correspond to",
                           "the 'Past (TMY)' column. The boxplots are computed from",
                           "the ensemble of future projections for the 2020s, 2050s,",
                           "and 2080s that is summarized in the 'Future' columns.") 
  image.row <- as.numeric(tail(row.locs$tas,1)) + 2
  insertImage(wb,sheet=sheet,file=sheets.closest$figs$tas,width=5.25,height=3,units='in',dpi=600,
                startRow=image.row,startCol=2)
  writeData(wb, sheet=sheet, 'Plot Description', 
                startRow = image.row, startCol = 7, headerStyle = titlestyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,titlestyle,
              rows=image.row,cols=7,gridExpand=FALSE,stack=FALSE)

  writeData(wb, sheet=sheet, tas.fig.description, 
                startRow = image.row+1, startCol = 7, headerStyle = textstyle,
                colNames=FALSE)
  addStyle(wb,sheet=sheet,textstyle,rows=seq(image.row+1,by=1,length.out=length(tas.fig.description)),
           cols=7:11,gridExpand=TRUE,stack=TRUE)
  
  freezePane(wb,sheet=sheet,firstActiveCol=2,firstActiveRow=3)

  ##------------------------
  ##About this file sheet
  if (!is.null(sheets.offset)) {
     sheet <- 4
  } else {
     sheet <- 3
  }
  addWorksheet(wb, "About This File") ##"Requested Site (CWEC)")
  start.row <- 1
  tab.cols <- 1:11
  writeData(wb, sheet=sheet, 'About This File', startRow = start.row, startCol = 1, headerStyle = titlestyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,titlestyle,rows=1,cols=tab.cols,gridExpand=FALSE,stack=FALSE)
  addStyle(wb,sheet=sheet,titlestyle,rows=2,cols=tab.cols,gridExpand=FALSE,stack=FALSE)
  start.row <- start.row + 2
  writeData(wb, sheet=sheet, 'File Description', startRow = start.row, startCol = 1, headerStyle = titlestyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,titlestyle,rows=start.row,cols=tab.cols,gridExpand=FALSE,stack=FALSE)

  start.row <- start.row + 1

  file.description1 <- c("This file contains summary statistics for the ",
                        paste0(nearest,"CWEC2016 Weather File"),
                       "Further about the variables in the rows can be seen by hovering over the row headings",
                       "(Column A) of each tab. This summary file for the past and future-shifted weather",
                       paste0("files (Version ",file.version,") was downloaded from:"))
  writeData(wb, sheet=sheet, file.description1, startRow = start.row, startCol = 1, headerStyle = textstyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,textstyle,rows=seq(start.row,by=1,length.out=length(file.description1)),cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  start.row <- start.row + length(file.description1)
  epw.file.link <- "https://pacificclimate.org/data/weather-files"
  names(epw.file.link) <- "https://pacificclimate.org/data/weather-files"
  class(epw.file.link) <- "hyperlink"
  writeData(wb, sheet=sheet, x=epw.file.link, startRow = start.row, startCol = 1, headerStyle = titlestyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,linkstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  start.row <- start.row + 1
  file.description2 <- c("More information about each tab is provided in the descriptions below.")
  writeData(wb, sheet=sheet, file.description2, startRow = start.row, startCol = 1, headerStyle = textstyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,textstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  addStyle(wb,sheet=sheet,textstyle,rows=start.row+1,cols=tab.cols,gridExpand=TRUE,stack=TRUE)

  start.row <- start.row + 2
  ##Methods and more information
  writeData(wb, sheet=sheet, 'Methods and more information', startRow = start.row, startCol = 1, headerStyle = titlestyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,titlestyle,rows=start.row,cols=tab.cols,gridExpand=FALSE,stack=FALSE)

  start.row <- start.row + 1

  meth.description1 <- c("These future-shifted weather files were created by applying morphing factors (Belcher,",
                         "et al., 2005) from an ensemble of 10 global climate models selected for Western Canada",
                         "(Cannon, 2015). Daily morphing factors (smoothed with a 21-day rolling mean) are",
                         "applied to the following variables:",
                         "dry bulb temperature, dewpoint temperature, relative humidity, surface pressure",
                         "Dry bulb temperature morphing factors are obtained from downscaled (BCCAQv2) climate",
                         "models (Cannon et al., 2015). Morphing factors for other variables are obtained from",
                         "GCMs interpolated to the BCCAQv2 resolution (about 10 km).",
                         "For a description of the TMY methodology see:")
  writeData(wb, sheet=sheet, meth.description1, startRow = start.row, startCol = 1, headerStyle = textstyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,textstyle,rows=seq(start.row,by=1,length.out=length(meth.description1)),cols=tab.cols,gridExpand=TRUE,stack=TRUE)

  start.row <- start.row + length(meth.description1)
  tmy.file.link <- "https://climate.onebuilding.org/papers"
  names(tmy.file.link) <- "https://climate.onebuilding.org/papers"
  class(tmy.file.link) <- "hyperlink"
  writeData(wb, sheet=sheet, x=tmy.file.link, startRow = start.row, startCol = 1, headerStyle = titlestyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,linkstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  start.row <- start.row + 1

  meth.description2 <- c("For the complete set of CWEC2016 files see:")
  writeData(wb, sheet=sheet, meth.description2, startRow = start.row, startCol = 1, headerStyle = textstyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,textstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)

  start.row <- start.row + length(meth.description2)
  cwec.file.link <- "https://climate.weather.gc.ca/prods_servs/engineering_e.html"
  names(cwec.file.link) <- "https://climate.weather.gc.ca/prods_servs/engineering_e.html"
  class(cwec.file.link) <- "hyperlink"
  writeData(wb, sheet=sheet, x=cwec.file.link, startRow = start.row, startCol = 1, headerStyle = titlestyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,linkstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  start.row <- start.row + 1

  meth.description3 <- c("For the more information on the 'morphing' procedure see:")
  writeData(wb, sheet=sheet, meth.description3, startRow = start.row, startCol = 1, headerStyle = textstyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,textstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)

  start.row <- start.row + length(meth.description3)
  morph.file.link <- "www.pacificclimate.org/sites/default/files/Eketal_2018_Proceedings_22_Feb_2019.pdf"
  names(morph.file.link) <- "www.pacificclimate.org/sites/default/files/Eketal_2018_Proceedings_22_Feb_2019.pdf"
  class(morph.file.link) <- "hyperlink"
  writeData(wb, sheet=sheet, x=morph.file.link, startRow = start.row, startCol = 1, headerStyle = titlestyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,linkstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  start.row <- start.row + 1

  meth.description4 <- c("For a detailed description of BCCAQv2 see:")
  writeData(wb, sheet=sheet, meth.description4, startRow = start.row, startCol = 1, headerStyle = textstyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,textstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)

  start.row <- start.row + length(meth.description4)
  bccaqv2.file.link <- "https://pacificclimate.org/data/statistically-downscaled-climate-scenarios"
  names(bccaqv2.file.link) <- "https://pacificclimate.org/data/statistically-downscaled-climate-scenarios"
  class(bccaqv2.file.link) <- "hyperlink"
  writeData(wb, sheet=sheet, x=bccaqv2.file.link, startRow = start.row, startCol = 1, headerStyle = titlestyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,linkstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  start.row <- start.row + 1
  
  ##Weather File Station Tab
  addStyle(wb,sheet=sheet,textstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  start.row <- start.row + 1
  ##Methods and more information
  writeData(wb, sheet=sheet, 'Weather File Station Tab', startRow = start.row, startCol = 1, headerStyle = titlestyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,titlestyle,rows=start.row,cols=tab.cols,gridExpand=FALSE,stack=FALSE)

  start.row <- start.row + 1

  tab.description1 <- c("The 'Past (TMY)' column summarizes the data from ",
                        paste0("CAN_BC_",nearest,"_CWEC2016.epw"),
                        "while the remaining columns summarize future-shifted values (i.e. data from files:",
                        paste0("2020s_CAN_BC_",nearest,"_CWEC2016.epw"),
                        paste0("2050s_CAN_BC_",nearest,"_CWEC2016.epw"),
                        paste0("2080s_CAN_BC_",nearest,"_CWEC2016.epw)."),
                        "The 2020s, 2050s, and 2080s columns summarize information averaged over the 30-year",
                        "period, based on the change according to the average, 10th percentile, and 90th",
                        "percentile from the 10-model ensemble. The percentiles are provided to encapsulate",
                        "the range of  plausible future change under the RCP8.5 (high emissions) scenario.")

  writeData(wb, sheet=sheet, tab.description1, startRow = start.row, startCol = 1, headerStyle = textstyle,
            colNames=FALSE)
  addStyle(wb,sheet=sheet,textstyle,rows=seq(start.row,by=1,length.out=length(tab.description1)),cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  start.row <- start.row + length(tab.description1)

  if (!is.null(sheets.offset)) {
     ##Adjusted Weather File Station Tab
     addStyle(wb,sheet=sheet,textstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)
     start.row <- start.row + 1
     writeData(wb, sheet=sheet, 'Adjusted Weather File Station Tab', startRow = start.row, startCol = 1, headerStyle = titlestyle,
               colNames=FALSE)
     addStyle(wb,sheet=sheet,titlestyle,rows=start.row,cols=tab.cols,gridExpand=FALSE,stack=FALSE)

     start.row <- start.row + 1

     adjust.description <- c("The adjusted tab summarizes past and future versions of",
                              paste0("CAN_BC_",site,"_offset-from_",nearest,"_CWEC2016.epw"),
                              "following the same format as the 'Weather File Station Tab'. Temperatures in this file",
                              "are adjusted using the difference between the monthly historical climatology of",
                              paste0(nearest," and ",site," according to PRISM"),
                              "temperature climatologies. For more information about the high-resolution PRISM",
                              "climatologies see:")
      writeData(wb, sheet=sheet, adjust.description, startRow = start.row, startCol = 1, headerStyle = textstyle,
                colNames=FALSE)
      addStyle(wb,sheet=sheet,textstyle,rows=seq(start.row,by=1,length.out=length(adjust.description)),cols=tab.cols,gridExpand=TRUE,stack=TRUE)

      start.row <- start.row + length(adjust.description)
      prism.file.link <-  "https://www.pacificclimate.org/data/prism-climatology-and-monthly-timeseries-portal"
      names(prism.file.link) <- "https://www.pacificclimate.org/data/prism-climatology-and-monthly-timeseries-portal"
      class(prism.file.link) <- "hyperlink"
      writeData(wb, sheet=sheet, x=prism.file.link, startRow = start.row, startCol = 1, headerStyle = titlestyle,
                    colNames=FALSE)
      addStyle(wb,sheet=sheet,linkstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)
      start.row <- start.row + 1
  }

  ##Temperature Climatologies
  addStyle(wb,sheet=sheet,textstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  start.row <- start.row + 1
  writeData(wb, sheet=sheet, 'Temperature Climatologies Tab', startRow = start.row, startCol = 1, headerStyle = titlestyle,
  colNames=FALSE)
  addStyle(wb,sheet=sheet,titlestyle,rows=start.row,cols=tab.cols,gridExpand=FALSE,stack=FALSE)

  start.row <- start.row + 1


  temp.description <- c("The temperature climatologies tab shows monthly, seasonal, and annual temperature",
                        "averages. Temperatures for the past are shown from both the TMY file and from the",
                        "gridded historical past (BCCAQv2) averaged over the years 1998-2014. This is the period",
                        "from which TMY months are chosen for the CWEC2016 files. Future changes in monthly,",
                        "seasonal and annual average temperatures from the ensemble of 10 BCCAQv2 downscaled",
                        "scenarios are also shown.")
  writeData(wb, sheet=sheet, temp.description, startRow = start.row, startCol = 1, headerStyle = textstyle,
                colNames=FALSE)
  addStyle(wb,sheet=sheet,textstyle,rows=seq(start.row,by=1,length.out=length(temp.description)),cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  start.row <- start.row + length(temp.description)
  ##--------------------------
  ##References
  addStyle(wb,sheet=sheet,textstyle,rows=start.row,cols=tab.cols,gridExpand=TRUE,stack=TRUE)
  start.row <- start.row + 1
  writeData(wb, sheet=sheet, 'References', startRow = start.row, startCol = 1, headerStyle = titlestyle,
  colNames=FALSE)
  addStyle(wb,sheet=sheet,titlestyle,rows=start.row,cols=tab.cols,gridExpand=FALSE,stack=FALSE)
  start.row <- start.row + 1

  ref.description <- c("Belcher, S.E., Hacker, J. and Powell, D.S. (2005) Constructing design",
                       "weather data for future climates. Building Services Engineering Research and",
                       "Technology. 26. 10.1191/0143624405bt112oa.",
                        "Cannon, A.J.,Sobie, S.R., and Murdock, T.Q. (2015) Bias correction of GCM ",
                        "precipitation by quantile mapping: How well do methods preserve changes in",
                        "quantiles and extremes? J. Clim., 28, pp. 6938-6959",
                        "Cannon,A.J. (2015) Selecting GCM scenarios that span the range of changes",
                        "in a multimodel ensemble: application to CMIP5 climate extremes indices",
                        "J. Clim. 28, pp. 1260-1267")

  writeData(wb, sheet=sheet, ref.description, startRow = start.row, startCol = 1, headerStyle = textstyle,
                colNames=FALSE)
  addStyle(wb,sheet=sheet,textstyle,rows=seq(start.row,by=1,length.out=length(ref.description)),cols=tab.cols,gridExpand=TRUE,stack=TRUE)

   
  saveWorkbook(wb, paste0(write.dir,site,'_Summary_rcp85.xlsx'), overwrite = TRUE)
  print(paste0(write.dir,site,'_Summary_rcp85.xlsx'))

}

