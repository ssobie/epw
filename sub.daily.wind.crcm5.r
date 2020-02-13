##Script to generate an hourly temperature series by combining BCCAQ2-PRISM 
##with the hourly temperatures from the YVR EPW File

library(MASS)
library(ncdf4)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##-----------------------------------------------------------------
##Gonzales Coordinates
lonc <- -123.325463
latc <- 48.413693

##-----------------------------------------------------------------------

make_cartoon <- function(series.past,series.future,
                         xlim,ylim,xlab,ylab='Density',brk,title,
                         plot.file) {
#   norm.past <- fitdistr(unlist(series.past),'weibull')
#   norm.future <- fitdistr(unlist(series.future),'weibull')
#   para.past <- norm.past$estimate
#   para.future <- norm.future$estimate

   png(file=plot.file,width=4,height=2.5,units='in',res=600,pointsize=6,bg='white')
   plot(c(),xaxs='i',yaxs='i',xlim=xlim,ylim=ylim,main=title,
           xlab=xlab,
           ylab=ylab)

   breaks <- seq(0,ceiling(max(c(unlist(series.past),unlist(series.future)),na.rm=T)/10)*10,by=brk)

   hist(breaks=breaks,unlist(series.past),border='blue',add=T,prob=TRUE)
   hist(breaks=breaks,unlist(series.future),border='red',add=T,prob=TRUE)

##   curve(dweibull(x, para.past[1], para.past[2]), from=0,col ='blue', add = TRUE,lwd=2)
##   curve(dweibull(x, para.future[1], para.future[2]),from=0, col ='red', add = TRUE,lwd=2)

   legend('topright',leg=c('1981-2010','2021-2050'),col=c('blue','red'),pch=16)
   box(which='plot')

   dev.off()

}

##-----------------------------------------------------------------

retrieve_derived_series <- function(var.name,scenario,
                                    type,var.file,
                                    lonc,latc,
                                    read.dir) {

  nc <- nc_open(paste0(read.dir,var.file))
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  lon.ix <- 28 ##which.min(abs(lonc-lon))
  lat.ix <- 7 ##which.min(abs(latc-lat))

  var.data <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  var.dates <- netcdf.calendar(nc)
  nc_close(nc)
  rv <- list(data=var.data,dates=var.dates)
  return(rv)
}

 
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
##Subset dates for pdf

subset_season <- function(tas,tasmax,season) {

   dates <- tasmax$dates

   seas.ix <- grepl(season,dates)
   data.sub <- tas[seas.ix]
   dates.sub <- dates[seas.ix]
   return(data.sub) ###list(data=data.sub,dates=dates.sub))
}



##-----------------------------------------------------------------
make_annual <- function(tas,tasmax) {

   yfac <- as.factor(format(tasmax$dates,'%Y'))
   tas.ann <- tapply(tas,yfac,mean)
   return(tas.ann)
}

##-----------------------------------------------------------------
make_seasonal <- function(tas,tasmax) {
   yfac <- as.factor(format(tasmax$dates,'%Y'))           
   fac <- as.factor(format(tasmax$dates,'%m'))
   seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
   fac <- factor(seasons[fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))
   tas.seas <- tapply(tas,list(yfac,fac),mean)
   return(tas.seas)
}


##-----------------------------------------------------------------
##Read Gonzales cell and save as RData (for 3 RCMs)

create_rcm_series <- function(var.name,type) {
   crcm5.dir <- '/storage/data/climate/downscale/RCM/CRCM5/reconfig/hourly/'
   save.dir <- '/storage/data/projects/rci/weather_files/Gonzales/'

   rcm.list <- c('ERA-Interim','CanESM2','MPI')
   var.files <- c('wspd_hour_WC011_ERA-Interim+CRCM5_historical_19800101-20141231.nc',
                  'wspd_hour_WC011_CanESM2+CRCM5_historical+rcp85_19800101-20501231.nc',
                  'wspd_hour_WC011_MPI+CRCM5_historical+rcp85_19800101-20501231.nc')
   rcm.file <- paste0(save.dir,var.name,'_',type,'_Gonzales_CRCM5.RData')
   if (!file.exists(rcm.file)) {
   rcm.series <- vector(mode='list',length=length(rcm.list))
      for (g in seq_along(rcm.list)) {
         gcm <- rcm.list[g]
         var.file <- var.files[g]
         rcm.series[[g]] <- retrieve_derived_series(
                                    var.name=var.name,scenario='rcp85',
                                    type=type,var.file=var.file,
                                    lonc=lonc,latc=latc,
                                    read.dir=crcm5.dir)     
      }
      save(rcm.series,file=rcm.file)
   } else {
      load(rcm.file)
   }
   return(rcm.series)
}

test <- create_rcm_series(var.name='wspd',type='RCM')

for (i in 2:3) {
   rcm.list <- c('ERA-Interim','CanESM2','MPI')
   slct <- test[[i]]
   series <- round(slct$data,2)
   dates <- slct$dates
   years <- format(dates,'%Y')
   months <- format(dates,'%m')
   days <- format(dates,'%d')
   hours <- format(dates,'%H')
   csv.version <- rbind(c('Year','Month','Day','Hour','Windspeed (m/s)'),
                         cbind(years,months,days,hours,series))
   csv.file <- paste0('/storage/data/projects/rci/weather_files/',
                      'Gonzales/windspeed_hourly_Gonzales_',rcm.list[i],'+CRCM5_historical+rcp85_19800101-20501231.csv')
   write.table(csv.version,file=csv.file,quote=F,row.name=F,col.name=F,sep=',')


wspd.past <- subset_dates(test[[i]],'1981-2010')
wspd.future <- subset_dates(test[[i]],'2021-2050')

plot.file <- paste0('/storage/data/projects/rci/weather_files/Gonzales/windspeed.',rcm.list[i],'.gonzales.pdfs.png')
xlab <- '3-Hour Windspeed (km/h)'

make_cartoon(wspd.past$data*3.6,wspd.future$data*3.6,
             xlim=c(0,45),ylim=c(0,0.11),xlab,ylab='Density',brk=0.5,title=paste0(rcm.list[i],'+CRCM5 3-Hourly Windspeed at Gonzales'),
             plot.file)

}


