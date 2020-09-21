##Read the downloaded ISD gz files and create single csv files for each site

##--------------------------------------------------------------------------

metadata.file <- paste0('/storage/data/projects/rci/weather_files/canada_cwec_files/',
                        'canada_cwec2016_files_metadata.csv')
region.metadata <- read.csv(metadata.file,header=T)


regions <- sort(c('alberta','new_brunswick','nova_scotia','prince_edward_island',
             'yukon','british_columbia','newfoundland_and_labrador','nunavut',
             'quebec','manitoba','northwest_territories','ontario','saskatchewan'))
regions <- 'british_columbia'

download.dir <- "/storage/data/climate/observations/station/noaa_isd/downloads/"

csv.title <- c('Year','Month','Day','Hour','TAS','Dewpoint','SLP',
               'WindDir','Windspeed','SkyCover','Precip','Precip_6_Hour')

for (region in regions) {
   print(regions)
   region.dir <- paste0(download.dir,region,'/')
   write.dir <- paste0("/storage/data/climate/observations/station/noaa_isd/",region,"/")
   if (!file.exists(write.dir)) {
      dir.create(write.dir,recursive=T)
   }
   sites <- list.files(path=region.dir)
   sites <- sites[66] ##YVR
   
   for (site in sites) {
      csv.template <- rep(NA,length(csv.title))
      print(site)
      print(which(sites %in% site))
      site.ix <- region.metadata$station_name %in% site
      site.info <- region.metadata[site.ix,]

      site.dir <- paste0(region.dir,site,'/')
      year.files <- sort(list.files(path=site.dir,pattern='.gz')) 
      year.files <- year.files[33:45]

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
      year.range <- range(csv.template[-1,1])
      month.range <- range(csv.template[-1,2])
      day.range <- range(csv.template[-1,3])
      dlen <- dim(csv.template)[1]
      start.date <- paste0(paste(csv.template[2,1:3],collapse='-'),' ',csv.template[2,4],":00:00")
      end.date <- paste0(paste(csv.template[dlen,1:3],collapse='-'),' ',csv.template[dlen,4],":00:00")
      ##Fill in missing years

      hourly.series <- seq(from=as.POSIXct(start.date,tz='UTC'),
                           by='hour',
                           to=as.POSIXct(end.date,tz='UTC'))
      hourly.compare <- format(hourly.series,'%Y-%m-%d-%H') 
      series.years <- format(hourly.series,'%Y') 
      series.months <- format(hourly.series,'%m') 
      series.days <- format(hourly.series,'%d') 
      series.hours <- format(hourly.series,'%H') 

      site.hours <- apply(cbind(csv.template[-1,1],
                                sprintf('%02d',csv.template[-1,2]),
                                sprintf('%02d',csv.template[-1,3]),
                                sprintf('%02d',csv.template[-1,4])),1,paste,collapse='-')
      hours.missing <- (hourly.compare %in% site.hours)
                    
      full.matrix <- matrix(NA,nrow=length(hourly.compare),ncol=12)
      full.matrix[hours.missing,] <- as.matrix(csv.template[-1,])
      full.matrix[,1] <- series.years
      full.matrix[,2] <- series.months
      full.matrix[,3] <- series.days
      full.matrix[,4] <- series.hours

      write.data <- rbind(csv.title,full.matrix)

      write.site.file <- paste0(write.dir,
                           site.info$station_name,'_',
                           site.info$station_id,'_',
                           region,'_ISD_hourly_',
                           format(as.Date(start.date),'%Y%m%d'),'-',
                           format(as.Date(end.date),'%Y%m%d'),'.csv')

      write.table(write.data,file=write.site.file,quote=F,row.name=F,col.name=F,sep=',')
browser()
   }


}






