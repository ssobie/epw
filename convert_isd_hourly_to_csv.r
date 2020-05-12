##Read the downloaded ISD gz files and create single csv files for each site

##--------------------------------------------------------------------------

metadata.file <- paste0('/storage/data/projects/rci/weather_files/canada_cwec_files/',
                        'canada_cwec2016_files_metadata.csv')
region.metadata <- read.csv(metadata.file,header=T)


regions <- sort(c('alberta','new_brunswick','nova_scotia','prince_edward_island',
             'yukon','british_columbia','newfoundland_and_labrador','nunavut',
             'quebec','manitoba','northwest_territories','ontario','saskatchewan'))

download.dir <- "/storage/data/climate/observations/station/noaa_isd/downloads/"

csv.title <- c('Year','Month','Day','Hour','TAS','Dewpoint','SLP',
               'WindDir','Windspeed','SkyCover','Precip','Precip_6_Hour')
csv.template <- rep(NA,length(csv.title))


for (region in regions) {

   region.dir <- paste0(download.dir,region,'/')
   write.dir <- paste0("/storage/data/climate/observations/station/noaa_isd/",region,"/")
   if (!file.exists(write.dir)) {
      dir.create(write.dir,recursive=T)
   }
   sites <- list.files(path=region.dir)
   
   for (site in sites) {
      site.ix <- region.metadata$station_name %in% site
      site.info <- region.metadata[site.ix,]

      site.dir <- paste0(region.dir,site,'/')
      year.files <- sort(list.files(path=site.dir,pattern='.gz'))
print('Check whether the years are sequential')
browser()
      for (i in seq_along(year.files)) {
         file <- year.files[i]
         site.data <- read.table(gzfile(paste0(site.dir,file)),header=FALSE)
         site.data[site.data==-9999] <- NA
         site.data[,5] <- site.data[,5] / 10 ##TAS Scaling
         site.data[,6] <- site.data[,6] / 10 ##Dewpoint Scaling
         site.data[,7] <- site.data[,7] / 10 ##SLP Scaling
         site.data[,9] <- site.data[,9] / 10 ##Windspeed Scaling
         site.data[,11] <- site.data[,11] / 10 ##Hourly Precip Scaling
         site.data[,12] <- site.data[,12] / 10 ##6-Hour Precip Scaling       
         csv.template <- rbind(csv.template,site.data)
      }
      
      write.data <- rbind(csv.title,csv.template[-1,])
      year.range <- range(csv.template[-1,1])
      month.range <- range(csv.template[-1,2])
      day.range <- range(csv.template[-1,3])

      write.site.file <- paste0(write.dir,
                           site.info$station_name,'_',
                           site.info$station_id,'_',
                           region,'_ISD_hourly_',
                           paste0(paste0(year.range[1],month.range[1],day.range[1]),
                                  paste0(year.range[2],month.range[2],day.range[2]),
                           collapse='-'),'.csv')
                         
      write.table(write.data,file=write.site.file,quote=F,row.name=F,col.name=F,sep=',')
      browser()
   }


}






