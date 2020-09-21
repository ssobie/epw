##Plot a projected map of Canada with the locations of all Canadian CWEC files
##identified on the map
##Have an option to represent sites as coloured circles so quantities
##can be added to the map.

source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/canada.epw.site.plot.r',chdir=T)


library(graticule)
library(raster)
library(rgdal)
library(scales)

convert_to_proj_coords <- function(lon,lat,proj.crs="+init=epsg:3005") {

  d <- data.frame(lon=lon,lat=lat)
  coordinates(d) <- c('lon','lat')
  proj4string(d) <- CRS("+init=epsg:4326")
  d.albers <- spTransform(d,CRS(proj.crs))
  rv <- d.albers@coords
  return(rv)
}

##---------------------------------------------------------------------
##*********************************************************************
##Canadian Projection - Lambert Conformal Conic
can.crs <- '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'

var.name <- 'cdd'
col.var <- 'tasmax'
class.breaks <- NULL


##-------------------------------------------------------------------
##EPW locations
provinces <- sort(c('alberta','new_brunswick','nova_scotia','prince_edward_island',
             'yukon','british_columbia','newfoundland_and_labrador','nunavut',
             'quebec','manitoba','northwest_territories','ontario','saskatchewan'))
##provinces <- 'british_columbia'
stats.proj.coords <- c(NA,NA,NA)

for (province in provinces) {
   print(province)
   stats.file <- paste0('/storage/data/projects/rci/weather_files/canada_cwec_files/canada_epw_files/',
                        province,'_CWEC2016_statistics.csv')

   stats.info <- read.csv(stats.file,header=T)

   stats.lons <- stats.info$lon
   stats.lats <- stats.info$lat
   prov.proj.coords <- matrix(NA,nrow=length(stats.lons),ncol=2)
   for (i in seq_along(stats.lons)) {
      prov.proj.coords[i,] <- convert_to_proj_coords(stats.lons[i],stats.lats[i],
                             proj.crs=can.crs)
   }
   prov.proj.coords <- cbind(prov.proj.coords,stats.info[[var.name]]) ##,as.character(stats.info[['X']]))
   stats.proj.coords <- rbind(stats.proj.coords,prov.proj.coords)  

}


##-------------------------------------------------------------------


stats.vals <- stats.proj.coords[,3]
map.range <- range(stats.vals,na.rm=T)
print(paste0('Map range: ',map.range))
if (is.null(class.breaks)) {
   class.breaks <- get.class.breaks(var.name,type='past',map.range,manual.breaks='')
}
map.class.breaks.labels <- get.class.break.labels(class.breaks)
col.ix <- unlist(sapply(stats.vals,function(x,y){tail(which(x>=y),1)},class.breaks))

bp <- 0
colour.ramp <- get.legend.colourbar(var.name=col.var,map.range=map.range,
                                    my.bp=bp,class.breaks=class.breaks,
                                    type='past')
stats.colours <- colour.ramp[col.ix]
##          delta_tas='\u00B0C',
leg.title <- 'Degree Days'

plot.points <- list(lon=stats.proj.coords[-1,1],
                    lat=stats.proj.coords[-1,2],
                    col=stats.colours)



plot.title <- 'Cooling Degree Days'
plot.dir <- '/storage/data/projects/rci/weather_files/plots/'
plot.file <- 'cdd.test.png'

make_canada_plot(var.name,col.name,plot.file,plot.title,
                 morph.proj=NULL,epw.proj.coords=plot.points,
                 class.breaks=class.breaks,colour.ramp=colour.ramp,
                 map.class.breaks.labels=map.class.breaks.labels,
                 leg.title='Degree Days')







