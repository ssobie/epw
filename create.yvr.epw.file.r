##Script to generate an hourly temperature series by combining BCCAQ2-PRISM 
##with the hourly temperatures from the YVR EPW File

library(MASS)
library(ncdf4)
source('/storage/home/ssobie/code/repos/epw/epw.support.functions.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/read.write.epw.r',chdir=T)

##-----------------------------------------------------------------
##YVR Coordinates
lonc <- -123.1815
latc <- 49.1967
 
##-----------------------------------------------------------------
##Subset dates for pdf

subset_dates <- function(input,interval) {

   bnds <- strsplit(interval,'-')[[1]]
   gcm.data <- input$data
   gcm.dates <- input$dates

   yst <- head(grep(bnds[1],gcm.dates),1)
   yen <- tail(grep(bnds[2],gcm.dates),1)

   data.sub <- gcm.data[yst:yen]
   dates.sub <- gcm.dates[yst:yen]
   return(list(data=data.sub,dates=dates.sub))
}

##-----------------------------------------------------------------
##Read EPW file for YVR
epw.file <- 'CAN_BC_Vancouver.Intl.AP.718920_2015.epw'
epw.dir <- '/storage/data/projects/rci/weather_files/wx_2016/'
epw.data <- read.epw.file(epw.dir,epw.file)
epw.tas <- epw.data$data[,7]
epw.years <- epw.data$data[,1]
epw.mons <- epw.data$data[,2]
epw.days <- epw.data$data[,3]
epw.hours <- epw.data$data[,4]
epw.dates <- paste0(paste(epw.years,
                          sprintf('%02d',epw.mons),
                          sprintf('%02d',epw.days),sep='-'),
                          ' ',sprintf('%02d',epw.hours-1),':00')
new.epw <- epw.data
new.epw.data <- epw.data$data

##-----------------------------------------------------------------
##Read hourly data for YVR (2017)

yvr.dir <- '/storage/data/projects/rci/weather_files/YVR/yvr_hourly_2017/'
yvr.files <- list.files(path=yvr.dir,pattern='eng',full.name=T)

yvr.matrix <- c()
for (yvr.file in yvr.files) {
   yvr.data <- read.csv(yvr.file,skip=15,header=T,as.is=T)
   yvr.matrix <- rbind(yvr.matrix,yvr.data)
}

yvr.matrix[['Stn.Press..kPa.']] <- yvr.matrix[['Stn.Press..kPa.']]*1000
yvr.matrix[['Wind.Dir..10s.deg.']] <- yvr.matrix[['Wind.Dir..10s.deg.']]*10

new.years <- rep(2017,length(epw.tas)) ##

yvr.vars <- c('Temp...C.','Dew.Point.Temp...C.','Rel.Hum....',
              'Wind.Dir..10s.deg.','Wind.Spd..km.h.',
              'Stn.Press..kPa.')
epw.vars <- c('dry_bulb_temperature','dew_point_temperature','relative_humidity',       
              'wind_direction','wind_speed','atmospheric_station_pressure')

##-----------------------------------------------------------------
##
date.ix <- substr(yvr.matrix$Date.Time,6,16) %in% substr(epw.dates,6,16)

for (v in seq_along(yvr.vars)) {
    new.series <- yvr.matrix[[yvr.vars[v]]]
    na.ix <- is.na(new.series)
    if (sum(na.ix) !=0) {
       for (n in which(na.ix)) {
          new.series[n] <- (new.series[n-1] + new.series[n+1])/2
       }
    }
    epw.ix <- get_field_index(epw.vars[v])
    new.epw.data[,epw.ix] <- new.series
}

new.epw.data[,1] <- new.years

new.epw.file <- 'CAN_BC_Vancouver.Intl.AP.718920_2017.epw'
write.epw.file(new.epw.data,new.epw$header,epw.dir,new.epw.file)


