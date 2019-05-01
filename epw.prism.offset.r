##Script to adjust EPW file temperature using PRISM

library(ncdf4)

source('/storage/home/ssobie/code/repos/crcm5/read.write.epw.r',chdir=T)


##----------------------------------------------------------

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

##----------------------------------------------------------

##Create a list of coordinates for all available EPW files

list_of_epw_coordinates <- function(epw.dir,epw.files) {
   epw.coords <- matrix(NA,nrow=length(epw.files),ncol=2)
   for (i in seq_along(epw.files)) {
      epw.coords[i,] <- get_epw_coordinates(epw.dir,epw.files[i])
   }
   return(epw.coords)
}

##----------------------------------------------------------

##

##----------------------------------------------------------

##Search the available EPW files and find the nearest to the
##supplied coordinates

find_closest_epw_file <- function(coords,
                                  epw.dir='/storage/data/projects/rci/weather_files/wx_files/') 
                         {
  epw.files <- list.files(path=epw.dir,pattern='CWEC.epw')
  epw.coords <- list_of_epw_coordinates(epw.dir,epw.files)
  wx.ix <- which.min(apply((epw.coords - matrix(coords,nrow=nrow(epw.coords),ncol=2,byrow=T))^2,1,sum))
  wx.selected <- epw.files[wx.ix]
  ###wx.selected <- 'CAN_BC_EDWARD-MILNE_9999999_CWEC.epw'
  return(list(dir=epw.dir,file=wx.selected))
}

##----------------------------------------------------------

##Open prism climatologies file

prism_nc <- function(var.name,prism.dir) {

   prism.filename <- paste0(prism.dir,var.name,'_monClim_PRISM_historical_run1_198101-201012.nc')
   nc <- nc_open(prism.filename)  
   return(nc)
}

##----------------------------------------------------------

##Find PRISM cell containing coordinates

get_prism_indices <- function(nc,coords) {

  lon <- ncvar_get(nc,'lon')
  lon.ix <- which.min(abs(coords[1]-lon))
  lat <- ncvar_get(nc,'lat')
  lat.ix <- which.min(abs(coords[2]-lat))
  rv <- c(lon.ix,lat.ix)
  return(rv)
}

##----------------------------------------------------------

##Calculate the TAS climatologies for the selected PRISM 
##cell

prism_tas <- function(cell,prism.dir) {

  tx.nc <- prism_nc('tmax',prism.dir)
  tasmax <- ncvar_get(tx.nc,'tmax',start=c(cell,1),count=c(1,1,-1))
  tn.nc <- prism_nc('tmin',prism.dir)
  tasmin <- ncvar_get(tn.nc,'tmin',start=c(cell,1),count=c(1,1,-1))
  tas <- (tasmax+tasmin)/2
  return(tas)  
}

##----------------------------------------------------------

##Adjust the epw temperature series with the PRISM 
##climatologies - col 7 for dry bulb temperature

adjust_epw_with_prism <- function(epw.data,prism.diff) {
   new.epw <- epw.data
   mons <- 1:12
   mon.dates <- epw.data[,2]
   for (mn in mons) {
     mn.ix <- mon.dates==mn
     new.epw[mn.ix,7] <- round(epw.data[mn.ix,7] + prism.diff[mn],1)
   }     
   return(new.epw)
}

##----------------------------------------------------------

##Given coordinates find the nearest weather file and adjust
##the temperature series based on the PRISM climatologies

generate_prism_offset <- function(lon,lat,epw.dir,prism.dir,new.location) {
   coords <- c(lon,lat)
   epw.closest <- find_closest_epw_file(coords,epw.dir)
   epw.closest.coords <- get_epw_coordinates(epw.closest$dir,
                                             epw.closest$file)

   prism.nc <- prism_nc('tmax',prism.dir)
   prism.cell <- get_prism_indices(prism.nc,coords)
   print(paste0('PRISM Cells: ',prism.cell))
   epw.cell <- get_prism_indices(prism.nc,epw.closest.coords)
   print(paste0('EPW Cells: ',prism.cell))
   prism.loc.tas <- prism_tas(prism.cell,prism.dir)
   prism.epw.tas <- prism_tas(epw.cell,prism.dir)
   prism.diff <- prism.loc.tas - prism.epw.tas
   print(prism.diff)
   epw <- read.epw.file(epw.closest$dir,
                        epw.closest$file)   
   epw.offset <- adjust_epw_with_prism(epw$data,prism.diff)
   epw.split <- strsplit(epw.closest$file,'_')[[1]]
   write.epw.file <- paste(epw.split[1],epw.split[2],
                           new.location,'offset_from',
                           epw.split[3],epw.split[4],
                           epw.split[5],sep='_')
   line.split <- strsplit(epw$header[1],',')[[1]]
   line.split[7] <- round(lat,3)
   line.split[8] <- round(lon,3)
   epw$header[1] <- paste0(line.split,collapse=',')
   write.epw.file(epw.offset,epw$header,
                  paste0(epw.dir,'offsets/'),write.epw.file)
   print(write.epw.file)
   epw.offset <- list(file=write.epw.file,dir=paste0(epw.dir,'offsets/'))
   rv <- list(closest=epw.closest,offset=epw.offset)
   return(rv)
}

##----------------------------------------------------------
##epw.dir <- '/storage/data/projects/rci/weather_files/wx_files/'
##prism.dir <- '/storage/data/climate/PRISM/dataportal/'
##new.location <- 'UBCWestbrook'

##test <- generate_prism_offset(-123.234678,49.253745,epw.dir,prism.dir,new.location) ##UBCWesbrook

##test <- generate_prism_offset(-123.094656,49.314414,epw.dir,prism.dir,new.location) ##Harbourside (BC Housing)
##test <- generate_prism_offset(-123.69667,48.38444,epw.dir,prism.dir,new.location) ##Sooke Reserve
##test <- generate_prism_offset(-123.722775,48.786147,epw.dir,prism.dir,new.location) ##Cowichan Hospital
##test <- generate_prism_offset(-122.891355,49.226239,epw.dir,prism.dir,new.location) ##RoyalColumbian
##test <- generate_prism_offset(-123.117856,49.218684,epw.dir,prism.dir,new.location) ##Dogwood
##test <- generate_prism_offset(-123.97049,49.18498,epw.dir,prism.dir,new.location) ##NRGH
##test <- generate_prism_offset(-123.1469,49.1687,epw.dir,prism.dir,new.location) ##Richmond
##test <- generate_prism_offset(-123.0684,49.3209,epw.dir,prism.dir,new.location) ##Lions gate
##test <- generate_prism_offset(-123.307119,48.464223,epw.dir,prism.dir,new.location)
##test <- generate_prism_offset(-123.0774,49.2697,epw.dir,prism.dir,new.location) ##1st and Clark

