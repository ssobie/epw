##Plot a projected map of Canada with the locations of all Canadian CWEC files
##identified on the map
##Have an option to represent sites as coloured circles so quantities
##can be added to the map.

library(ncdf4)
library(maps)
source('/storage/home/ssobie/code/repos/assessments/resource.region.map.support.r',chdir=T)
source('/storage/home/ssobie/code/repos/LHASA/bc.albers.map.support.r')

##-------------------------------------------------------------------

make_epw_bc_albers_plot <- function(var.name,plot.type,seas,plot.title,
                                    plot.data,epw.proj.coords=NULL,
                                    shared.range=NULL) {

   col.var <- var.name

   if (!is.null(shared.range)) {
      map.range <- shared.range
   } else {
      map.range <- range(as.matrix(plot.data),na.rm=T)
   }

   class.breaks <- get.class.breaks(col.var,plot.type,map.range,manual.breaks='')
   print('Default Breaks')
   print(get.class.breaks(col.var,plot.type,range(as.matrix(plot.data),na.rm=T),manual.breaks=''))

   colour.ramp <- get.legend.colourbar(var.name=col.var,map.range=map.range,
                                       my.bp=0,class.breaks=class.breaks,
                                       type=plot.type)
   map.class.breaks.labels <- get.class.break.labels(class.breaks,type=plot.type)

   alb.crs <- "+init=epsg:3005"
   lons <- c(-145.0, -140.0,-135.0,-130.0,-125.0,-120.0,-115.0,-110.0)
   lats <- c(  48.0,  50.0,  52.0,  54.0,  56.0,  58.0,  60.0,  62.0)
   grats <- add.graticules(lons,lats,alb.crs)

   shape.dir <- '/storage/data/projects/rci/data/assessments/shapefiles/bc'
   bc.shp <- get.region.shape('bc',shape.dir)

   map.extent <- extent(c(-136.0,-114.0,48.0,60.5))
   crop.extent <- extent(c(-145,-100,45.0,65))

   map.extent <- extent(c(-136.0,-114.0,48.0,60.5))
   crop.extent <- extent(c(-145,-100,45.0,65))

   plot.data.albers <- projectRaster(plot.data,crs=CRS(alb.crs))
   map.bounds <- extent(plot.data.albers)
   plot.bounds <- make.plot.window(map.bounds,spTransform(bc.shp,CRS(alb.crs)))

   plot.window.xlim <- plot.bounds$xlim ##c(map.bounds@xmin,map.bounds@xmax)
   plot.window.ylim <- plot.bounds$ylim ##c(map.bounds@ymin,map.bounds@ymax)

   par(mar=c(5,5,2,2))

   plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
   bg='white',col.main='gray22',col.lab='gray22', # 'gray94',
   xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main=plot.title,
   cex.axis=2,cex.lab=1.75,cex.main=2,mgp=c(2.5,2,0),axes=F)
   rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightgray')

   xtks <- get.proj.xaxis(lons,alb.crs,plot.window.ylim)
   ytks <- get.proj.yaxis(lats,alb.crs,plot.window.xlim)

   axis(1,at=xtks,label=rep('',length(lons)),cex.axis=1.75,col='gray22')
   axis(1,at=xtks,label=lons,cex.axis=1.75,col='gray22',col.axis='gray22')
   axis(2,at=ytks,label=rep('',length(lats)),cex.axis=1.75,col='gray22')
   axis(2,at=ytks,label=lats,cex.axis=1.75,col='gray22',col.axis='gray22')

   bc.overlay <- 'north_america_state_provincial_boundaries'
   borders.shp <- readOGR(shape.dir, bc.overlay, stringsAsFactors=F, verbose=F)
   wco.shp <- readOGR('/storage/data/projects/rci/data/assessments/shapefiles/bc_common',
                           'ocean_mask', stringsAsFactors=F, verbose=F)

##   can.shp <- readOGR('/storage/data/gis/basedata/base_layers',
##                      'canada_provinces', stringsAsFactors=F, verbose=F)
   na.shp <- readOGR('/storage/data/projects/rci/data/assessments/shapefiles/bc_common',
                      'bc_surroundings', stringsAsFactors=F, verbose=F)

   image(plot.data.albers,add=T,col=colour.ramp,breaks=class.breaks)

   ##plot(spTransform(can.shp,CRS(alb.crs)),add=TRUE,border='lightgray',lwd=0.25)
   plot(spTransform(wco.shp,CRS(alb.crs)),add=TRUE,col='aliceblue',lwd=0.25)
   plot(spTransform(na.shp,CRS(alb.crs)),add=TRUE,col='lightgray',lwd=0.25)
##   plot(spTransform(borders.shp,CRS(alb.crs)),add=TRUE,border='lightgray',lwd=0.25)
   if (!is.null(epw.proj.coords)) {
##      points(epw.proj.coords$lon,epw.proj.coords$lat,
##             pch=21,cex=1.5,col='gray5',lwd=0.5,
##             bg=epw.proj.coords$col)
      points(epw.proj.coords,pch=18,cex=2,col='royalblue')
   }

   plot(grats,add=TRUE,lty=3,col='gray',lwd=0.7)



##   text.albers <- convert.to.alb.coords(lon=-132,lat=49,crs=alb.crs)
##   text(text.albers[1],text.albers[2],plot.title,cex=2.5)

   leg.title <- 'Percent'
   if ((var.name=='pr' | var.name=='rx5dayETCCDI') & plot.type=='past') {
     leg.title <- 'mm'
   }
   if ((var.name=='tasmax' | var.name=='tasmin') & plot.type=='past') {
     leg.title <- '\u00B0C'
   }

   legend('topright', col = "black", legend=rev(map.class.breaks.labels), pch=22, pt.bg = rev(colour.ramp),
         pt.cex=2.95, y.intersp=0.8, title.adj=0.2, title='\u00B0C', xjust=0, cex=1.55) ##leg.title
   box(which='plot',lwd=2,col='gray22')
   rv <- list(labels=rev(map.class.breaks.labels),cols=rev(colour.ramp))
   return(rv)
}
##-------------------------------------------------------------------
convert_to_proj_coords <- function(lon,lat,proj.crs="+init=epsg:3005") {

  d <- data.frame(lon=lon,lat=lat)
  coordinates(d) <- c('lon','lat')
  proj4string(d) <- CRS("+init=epsg:4326")
  d.albers <- spTransform(d,CRS(proj.crs))
  rv <- d.albers@coords
  return(rv)
}

##-------------------------------------------------------------------

read_climatologies <- function(var.name,clim,ex,
                               clim.dir,clim.file) {

   clim.brick <- brick(paste0(clim.dir,clim.file))
   clim.mean <- calc(clim.brick,mean)
   clim.data <- subset(clim.brick,1)
   clim.crop <- crop(clim.data,ex)
   return(clim.crop)
}

##-------------------------------------------------------------------
##*********************************************************************
obs.dir <- '/storage/data/climate/observations/gridded/'
plot.dir <- '/storage/data/projects/rci/data/cas/wx_files/'

##BC Extent from PRISM extent
ex <- extent(x=c(-141.0,-105.0),y=c(47,63))

plot.file <- paste0(plot.dir,'tasmaxx_975.epw.locations.for.bc.png')
var.name <- 'tasmax'
plot.type <- 'past'
plot.title <- 'EPW Locations'


##-----------------------------------------------------------
##PNWNAmet Climatologies
pnw.dir <- paste0(obs.dir,'PNWNAmet/Derived/annual_quantiles/')

pnw.file <- "tasmax_annual_quantile_975_PNWNAmet_observations_1945-2012.nc"

pnw.data <- read_climatologies(var.name,clim,ex,
                               pnw.dir,pnw.file)
##-----------------------------------------------------------

##EPW locations
epw.metadata.file <- paste0('/storage/data/projects/rci/weather_files/canada_cwec_files/',
                            'bc_cwec2016_files_metadata.csv')

epw.info <- read.csv(epw.metadata.file,header=T)

epw.lons <- epw.info$lon
epw.lats <- epw.info$lat
epw.proj.coords <- matrix(NA,nrow=length(epw.lons),ncol=2)
for (i in seq_along(epw.lons)) {
   epw.proj.coords[i,] <- convert_to_proj_coords(epw.lons[i],epw.lats[i],
                          proj.crs="+init=epsg:3005")
}


##-----------------------------------------------------------

png(file=plot.file,width=5,height=5,units='in',res=600,pointsize=6,bg='white')
##par(mar=c(0,0,0,0),oma=c(6,6,4,8))
##par(mgp=c(4,1.5,0))

rv <- make_epw_bc_albers_plot(var.name,plot.type,seas,plot.title,
                              pnw.data,epw.proj.coords=epw.proj.coords,
                              shared.range=NULL)
leg.title <- '\u00B0C'
par(xpd=NA)
legend('topright',inset=c(-0.32,0), col = "black",
                                    legend=rv$labels, pch=22, pt.bg = rv$cols,
        pt.cex=4, y.intersp=0.8, title.adj=0.2, title=leg.title, xjust=0, cex=2.2,box.lwd=2)
par(xpd=FALSE)

dev.off()



