##Plot a projected map of Canada with the locations of all Canadian CWEC files
##identified on the map
##Have an option to represent sites as coloured circles so quantities
##can be added to the map.

source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)

library(graticule)
library(raster)
library(rgdal)

add_graticules <- function(lons,lats,crs) {

  xl <-  range(lons)
  yl <- range(lats)
  grat <- graticule(lons, lats, proj = CRS(crs),xlim=c(-180,0),ylim=c(30,89))
  rv <- grat
  return(rv)
}

##X-Axis Ticks
get.proj.xaxis <- function(lons,crs,plot.window.ylim) {

  y <- seq(0,80,0.1)
  xm <- sapply(lons,rep,length(y))
  S <- apply(xm,2,function(x) {SpatialPoints(cbind(x,y), proj4string = CRS("+proj=longlat +datum=WGS84"))})
  S2<- lapply(S,spTransform, crs)
  indices <- lapply(S2,function(x){which.min(abs(x@coords[,'y']-plot.window.ylim[1]))})
  xticks <- mapply(FUN=function(s,i){s@coords[,'x'][i]},S2,indices)
  return(xticks)
}


 ##Y-Axis Ticks
get.proj.yaxis <- function(lats,crs,plot.window.xlim) {

  x <- seq(-180,-80,0.1)
  ym <- sapply(lats,rep,length(x))
  S <- apply(ym,2,function(y) {SpatialPoints(cbind(x,y), proj4string = CRS("+proj=longlat +datum=WGS84"))})
  S2<- lapply(S,spTransform, crs)
  indices <- lapply(S2,function(x){which.min(abs(x@coords[,'x']-plot.window.xlim[1]))})
  yticks <- mapply(FUN=function(s,i){s@coords[,'y'][i]},S2,indices)
  return(yticks)
}

convert_to_proj_coords <- function(lon,lat,proj.crs="+init=epsg:3005") {

  d <- data.frame(lon=lon,lat=lat)
  coordinates(d) <- c('lon','lat')
  proj4string(d) <- CRS("+init=epsg:4326")
  d.albers <- spTransform(d,CRS(proj.crs))
  rv <- d.albers@coords
  return(rv)
}

##-------------------------------------------------------------------

make_canada_plot <- function(morph.name,var.name,plot.file,plot.title,
                             morph.proj,epw.proj.coords,can.crs,
                             class.breaks=NULL) {

   map.range <- range(as.matrix(morph.proj),na.rm=T)
   print(paste0('Map range: ',map.range))
   if (is.null(class.breaks)) {
      class.breaks <- get.class.breaks(var.name,type='anomaly',map.range,manual.breaks='')
   }
   
   colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=map.range,
                                       my.bp=0,class.breaks=class.breaks,
                                       type='anomaly')
   map.class.breaks.labels <- get.class.break.labels(class.breaks)


   plot.window.xlim <- c(-2500000,3100000)
   plot.window.ylim <- c(100000,4800000)

   lons <- seq(-180,0,by=10)
   lats <- seq(30,90,by=5)
   grats <- add_graticules(lons,lats,can.crs)
  
   xtks <- get.proj.xaxis(lons,can.crs,plot.window.ylim)
   ytks <- get.proj.yaxis(lats,can.crs,plot.window.xlim)

   shape.dir <- '/storage/data/gis/basedata/base_layers'

   provinces.region <- 'canada_provinces'
   can.shp <- readOGR(shape.dir, provinces.region, stringsAsFactors=F, verbose=F)

   na.region <- 'north_america_state_provincial_boundaries'
   na.shp <- readOGR(shape.dir, na.region, stringsAsFactors=F, verbose=F)
   na.proj.shp <- spTransform(na.shp,CRS(can.crs))
   ocean.shp <- readOGR(shape.dir, 'ocean_sym_diff_continents', stringsAsFactors=F, verbose=F)
   us.shp <- readOGR(shape.dir, 'united_states', stringsAsFactors=F, verbose=F)
   greenland.shp <- readOGR('/storage/home/ssobie/general/', 'greenland', stringsAsFactors=F, verbose=F)
   lakes.shp <- readOGR(shape.dir, 'lakes', stringsAsFactors=F, verbose=F)
   rivers.shp <- readOGR(shape.dir, 'rivers', stringsAsFactors=F, verbose=F)

   ocean.transformed <- spTransform(ocean.shp, CRS(can.crs))
   us.transformed <- spTransform(us.shp, CRS(can.crs))
   can.transformed <- spTransform(can.shp, CRS(can.crs))
   greenland.transformed <- spTransform(greenland.shp, CRS(can.crs))
   lakes.transformed <- spTransform(lakes.shp, CRS(can.crs))
   rivers.transformed <- spTransform(rivers.shp, CRS(can.crs))

   plot.dir <- '/storage/data/projects/rci/weather_files/plots/'
   write.file <- paste0(plot.dir,plot.file)
  
   cx <- 2.5

   png(file=write.file,width=6,height=5,units='in',res=600,pointsize=6,bg='white')
   par(mar=c(5.1,5.6,4.1,2.1))
   plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
   bg='white',# 'gray94',
   xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main=plot.title,
   cex.axis=cx,cex.lab=cx,cex.main=1.5,mgp=c(3.5,2,0),axes=F)
   rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightgray')
  
   image(morph.proj,col=colour.ramp,breaks=class.breaks,axes=FALSE,
         xlim=plot.window.xlim,ylim=plot.window.ylim,bg='white',
        main='',xlab='',ylab='',add=T)

   lw <- 0.4


   plot(us.transformed,col='gray',add=TRUE,lwd=lw)
   plot(lakes.transformed,col='lightblue',border='lightblue',add=TRUE,lwd=lw)
   plot(can.transformed,add=TRUE,lwd=lw)

   plot(greenland.transformed,col='lightgray',add=TRUE,lwd=lw)
   plot(ocean.transformed,col='lightblue',add=TRUE,lwd=lw)
   plot(grats,add=TRUE,lty=5,col='gray',lwd=lw)

   points(epw.proj.coords,pch=18,cex=0.75,col='black')
   points(epw.proj.coords,pch=18,cex=0.5,col='green')

   ##plot(rivers.transformed,col='lightblue',add=TRUE,lwd=lw)
   
   axis(1,at=xtks,label=rep('',length(lons)),cex.axis=cx-1)
   axis(1,at=xtks,label=lons,cex.axis=cx-1)
   axis(2,at=ytks,label=rep('',length(lats)),cex.axis=cx-1)
   axis(2,at=ytks,label=lats,cex.axis=cx-1)

   leg.title <- switch(morph.name,
                       delta_tas='\u00B0C',
                       alpha_tasmax_tasmin='Ratio')

   legend('topright', col = "black", legend=rev(map.class.breaks.labels), 
           pch=22, pt.bg = rev(colour.ramp),
           pt.cex=4.0, y.intersp=0.8, title.adj=0.2, title=leg.title, xjust=0, cex=1.25)
   box(which='plot',lwd=2)

   dev.off()
}


##---------------------------------------------------------------------
##*********************************************************************


##Canadian Projection - Lambert Conformal Conic
can.crs <- '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'

##-------------------------------------------------------------------
##EPW locations
epw.metadata.file <- paste0('/storage/data/projects/rci/weather_files/canada_cwec_files/',
                            'canada_cwec2016_files_metadata.csv')

epw.info <- read.csv(epw.metadata.file,header=T)

epw.lons <- epw.info$lon
epw.lats <- epw.info$lat
epw.proj.coords <- matrix(NA,nrow=length(epw.lons),ncol=2)
for (i in seq_along(epw.lons)) {
   epw.proj.coords[i,] <- convert_to_proj_coords(epw.lons[i],epw.lats[i],
                          proj.crs=can.crs)
}




##-------------------------------------------------------------------
##Morphing Factors for illustration

var.name <- 'tas' ##'alpha_tas' ##
morph.title <- 'Delta Tas' ##'Alpha Tasmax/Tasmin' ##
morph.name <- 'delta_tas' ##'alpha_tasmax_tasmin'  ##
past.int <- '1998-2014'
proj.int <- '2041-2070'


##day <- 335 ##Dec 1
##day <- 32 ##Feb 1
##day <- 244 ##Sept 1
day <- 152 ##June 1 Day of the year 

year.date <- as.Date(paste0('1995-',sprintf('%03d',day)),format='%Y-%j')
text.date <- format(year.date,'%B %d')

shared.breaks <- round(seq(0.4,4,0.4),1)
diff.breaks <- round(seq(-0.6,0.6,by=0.2),1)

type <- 'GCM'
gcm.plot.title <- paste0(morph.title,' ',type,' Ensemble Average\n',text.date)
gcm.plot.file <- paste0(tolower(type),'_ens_avg_',tolower(morph.name),'_',
                        past.int,'_',proj.int,'_',gsub(' ','_',tolower(text.date)),'.png')
gcm.morph.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/epw_factors/ens_avg_',tolower(type),'_',past.int,'/')
gcm.morph.file <- paste0('roll21_',morph.name,'_ENSEMBLE_',type,'_',past.int,'_',proj.int,'.nc')
gcm.morph.brick <- brick(paste0(gcm.morph.dir,gcm.morph.file))
gcm.morph.day <- subset(gcm.morph.brick,day)
gcm.morph.proj <- projectRaster(gcm.morph.day,crs=CRS(can.crs))
if (1==1) {
make_canada_plot(morph.name,var.name,
                gcm.plot.file,gcm.plot.title,
                gcm.morph.proj,
                epw.proj.coords,can.crs,
                class.breaks=shared.breaks)
}

type <- 'BCCAQ2'
bccaq2.plot.title <- paste0(morph.title,' ',type,' Ensemble Average\n',text.date)
bccaq2.plot.file <- paste0(tolower(type),'_ens_avg_',tolower(morph.name),'_',
                           past.int,'_',proj.int,'_',gsub(' ','_',tolower(text.date)),'.png')
bccaq2.morph.dir <- paste0('/storage/data/climate/downscale/BCCAQ2/epw_factors/ens_avg_',tolower(type),'_',past.int,'/')
bccaq2.morph.file <- paste0('roll21_',morph.name,'_ENSEMBLE_',type,'_',past.int,'_',proj.int,'.nc')
bccaq2.morph.brick <- brick(paste0(bccaq2.morph.dir,bccaq2.morph.file))
bccaq2.morph.day <- subset(bccaq2.morph.brick,day)
bccaq2.morph.proj <- projectRaster(bccaq2.morph.day,crs=CRS(can.crs))

if (1==1) {
make_canada_plot(morph.name,var.name,
                 bccaq2.plot.file,bccaq2.plot.title,
                 bccaq2.morph.proj,
                 epw.proj.coords,can.crs,
                 class.breaks=shared.breaks)
}

##Difference
diff.plot.title <- paste0(morph.title,' Difference (BCCAQ2-GCM) in Ensemble Average\n',text.date)
diff.plot.file <- paste0('Difference_ens_avg_',tolower(morph.name),'_',
                         past.int,'_',proj.int,'_',gsub(' ','_',tolower(text.date)),'.png')
diff.morph.file <- paste0('roll21_',morph.name,'_ENSEMBLE_DIFFERENCE_BCCAQ2-GCM_',past.int,'_',proj.int,'.nc')
diff.morph.day <- bccaq2.morph.day - gcm.morph.day
diff.morph.proj <- projectRaster(diff.morph.day,crs=CRS(can.crs))

make_canada_plot(morph.name,var.name,
                 diff.plot.file,diff.plot.title,
                 diff.morph.proj,
                 epw.proj.coords,can.crs,
                 class.breaks=diff.breaks)



