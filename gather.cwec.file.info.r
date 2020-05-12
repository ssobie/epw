##Go through the EPW files for Canada and create a csv file 
##with the location, elevation, and station info from each file.

source('/storage/home/ssobie/code/repos/epw/read.write.epw.r',chdir=T)


##------------------------------------------------------------------------------
##Return the information for the supplied EPW file

get_epw_info <- function(epw.dir,epw.file) {

   epw <- read.epw.file(epw.dir,epw.file)
   epw.header <- epw$header
   epw.first <- strsplit(epw.header[1],',')[[1]]
   lon <- as.numeric(epw.first[8]) ##Fixed location
   lat <- as.numeric(epw.first[7])
   stn.name <- gsub(' ','_',epw.first[2])
   stn.id <- epw.first[6]
   elev <- epw.first[10]

   if (lon < -180 | lon > 180) {
      stop('Ill defined longitude coordinate')
   }
   if (lat < 40 | lat > 90) {
      stop('Ill defined latitude coordinate')
   }

   rv <- list(lon=lon,lat=lat,name=stn.name,id=stn.id,elev=elev)
   return(rv)
}

##-----------------------------------------------------------------------
##
epw.dir <- '/storage/data/projects/rci/weather_files/canada_cwec_files/canada_epw_files/'

regions <- sort(c('alberta','new_brunswick','nova_scotia','prince_edward_island',
             'yukon','british_columbia','newfoundland_and_labrador','nunavut',
             'quebec','manitoba','northwest_territories','ontario','saskatchewan'))

file.type <- 'CWEC2016'

##----------------------------------------------------------------------
##Format of output file
## Station Name, Region, StationID,Lon,Lat,Elevation 
titles <- c('station_name','region','station_id','lon','lat','elev')

epw.metadata <- titles

for (region in regions) {
   epw.files <- sort(list.files(path=paste0(epw.dir,region),pattern=file.type))
   for (i in seq_along(epw.files)) {
      epw.file <- epw.files[i]
      epw.info <- get_epw_info(paste0(epw.dir,region,'/'),epw.file)
      print(epw.info$name)
      epw.line <- c(epw.info$name,region,epw.info$id,
                    epw.info$lon,epw.info$lat,epw.info$elev)
      epw.metadata <- rbind(epw.metadata,epw.line)                          
   }

}

write.dir <- '/storage/data/projects/rci/weather_files/canada_cwec_files/'
write.file <- 'canada_cwec2016_files_metadata.csv'
write.table(epw.metadata,file=paste0(write.dir,write.file),
            sep=',',quote=FALSE,row.name=F,col.name=F)
