##Download NOAA Integrated Surface Data for EPW CWEC Canada sites


download_isd_file <- function(stn.id,year,stn.name,region,write.dir) {
  isd.path <- paste0("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-lite/",year,"/")
  isd.file <- paste0(stn.id,'-99999-',year,'.gz')

  dest.file <- paste0(stn.name,'_',stn.id,'_',region,'_ISD_hourly_',year,'.gz')
  dest <- paste0(write.dir,dest.file)
  url <- paste0(isd.path,isd.file)
  tryCatch(download.file(url,destfile=dest,quiet=TRUE),
     error = function(e) print("No file available"))
}


##---------------------------------------------------------

##*********************************************************

years <- 1945:2019

download.dir <- '/storage/data/climate/observations/station/noaa_isd/downloads/'

metadata.file <- paste0('/storage/data/projects/rci/weather_files/canada_cwec_files/',
                        'canada_cwec2016_files_metadata.csv')
region.metadata <- read.csv(metadata.file,header=T)

##'alberta',
regions <- sort(c('new_brunswick','nova_scotia','prince_edward_island',
             'yukon','british_columbia','newfoundland_and_labrador','nunavut',
             'quebec','manitoba','northwest_territories','ontario','saskatchewan'))
regions <- 'alberta'
##Turn of the warnings from failing to find files
oldw <- getOption('warn')
options(warn = -1)

for (region in regions) {
    print(region)
    region.ix <- region.metadata$region %in% region
    sites <- region.metadata[region.ix,]   

    for (i in seq_along(sites[,1])) {
       site <- sites[i,]
       print(site)
       site.dir <- paste0(download.dir,region,'/',site$station_name,'/')
       if (!file.exists(site.dir)) {
          dir.create(site.dir,recursive=T)
       }
       dl.status <- rep(NA,length(years))
       for (y in seq_along(years)) {
          print(years[y])
          dl.status[y] <- download_isd_file(site$station_id,years[y],site$station_name,region,site.dir)
       }
    }
}

options(warn=oldw)
