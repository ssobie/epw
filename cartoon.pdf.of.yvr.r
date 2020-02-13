##Script to generate an hourly temperature series by combining BCCAQ2-PRISM 
##with the hourly temperatures from the YVR EPW File

library(MASS)
library(ncdf4)
source('/storage/home/ssobie/code/repos/assessments/single.site.bccaq.to.prism.r',chdir=T)
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
##Read YVR cell and save as RData (for 12 GCMs)

create_bccaq2_prism_series <- function(var.name) {
   bccaq2.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
   save.dir <- '/storage/data/projects/rci/weather_files/YVR/'

   gcm.list <- c('ACCESS1-0','CCSM4','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')

   bccaq2.file <- paste0(save.dir,var.name,'_YVR_BCCAQ2-PRISM_daily_series_1951-2100.RData')
   if (!file.exists(bccaq2.file)) {

      bccaq2.series <- retrieve_bccaq2_prism_series(
                                var.name=var.name,scenario='rcp85',
                                gcm.list=gcm.list,
                                lonc=lonc,latc=latc,
                                base.dir=bccaq2.dir)     

      save(bccaq2.series,file=bccaq2.file)
   } else {
      load(bccaq2.file)
   }
   return(bccaq2.series)
}

##-----------------------------------------------------------------
##Read EPW hourly series for YVR
epw.file <- 'CAN_BC_Vancouver.Intl.AP.718920_2015.epw'
epw.dir <- '/storage/data/projects/rci/weather_files/wx_2016/'
epw.data <- read.epw.file(epw.dir,epw.file)
epw.tas <- epw.data$data[,7]
epw.years <- rep(1995,length(epw.tas)) ##epw.data$data[,1]
epw.mons <- epw.data$data[,2]
epw.days <- epw.data$data[,3]

##-----------------------------------------------------------------
##

tasmax.series <- create_bccaq2_prism_series(var.name='tasmax') 
tasmin.series <- create_bccaq2_prism_series(var.name='tasmin') 



##-----------------------------------------------------------------
##

epw.dates <- as.Date(paste(epw.years,epw.mons,epw.days,sep='-'))
day.fac <- as.factor(format(epw.dates,'%Y-%m-%d'))
epw.tas.daily <- tapply(epw.tas,day.fac,mean)

mon.fac <- as.factor(format(epw.dates,'%m'))
seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
seas.fac <- factor(seasons[mon.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))

epw.tas.seas <- tapply(epw.tas,seas.fac,mean)


##-----------------------------------------------------------------
##

tasmax.sub.1971.2000 <- lapply(tasmax.series,subset_dates,'1971-2000')
tasmin.sub.1971.2000 <- lapply(tasmin.series,subset_dates,'1971-2000')
tas.sub.1971.2000 <- mapply(FUN=function(x,y){(x$data+y$data)/2},x=tasmax.sub.1971.2000,y=tasmin.sub.1971.2000)
tas.summer.1971.2000 <- mapply(FUN=subset_season,tas.sub.1971.2000,tasmax.sub.1971.2000,'(-06-|-07-|-08-)')

tas.ann.1971.2000 <- mapply(FUN=make_annual,tas.sub.1971.2000,tasmax.sub.1971.2000)
tas.seas.1971.2000 <- mapply(FUN=make_seasonal,tas.sub.1971.2000,tasmax.sub.1971.2000)
jja.ix <- seq(3,120,4)
tas.sum.1971.2000 <- tas.seas.1971.2000[jja.ix,]

tasmax.sub.1986.2015 <- lapply(tasmax.series,subset_dates,'1986-2015')
tasmin.sub.1986.2015 <- lapply(tasmin.series,subset_dates,'1986-2015')
tas.sub.1986.2015 <- mapply(FUN=function(x,y){(x$data+y$data)/2},x=tasmax.sub.1986.2015,y=tasmin.sub.1986.2015)
tas.summer.1986.2015 <- mapply(FUN=subset_season,tas.sub.1986.2015,tasmax.sub.1986.2015,'(-06-|-07-|-08-)')

tas.seas.1986.2015 <- mapply(FUN=make_seasonal,tas.sub.1986.2015,tasmax.sub.1986.2015)
tas.sum.1986.2015 <- tas.seas.1986.2015[jja.ix,]

norm.1971.2000 <- fitdistr(unlist(tas.sum.1971.2000),'normal')
norm.1986.2015 <- fitdistr(unlist(tas.sum.1986.2015),'normal')
para.1971.2000 <- norm.1971.2000$estimate
para.1986.2015 <- norm.1986.2015$estimate


plot.file <- '/storage/data/projects/rci/weather_files/cartoon.yvr.2015.png'

##png(file=plot.file,width=5,height=4,units='in',res=600,pointsize=6,bg='white')
plot(c(),xlim=c(-8,27),ylim=c(0,0.20),xaxs='i',yaxs='i',
        xlab='Average Summer Temperature (\u00B0C)',
        ylab='Density')
hist(unlist(tas.sum.1971.2000),border='blue',add=T,prob=TRUE)
hist(unlist(tas.sum.1986.2015),border='red',add=T,prob=TRUE)

curve(dnorm(x, para.1971.2000[1], para.1971.2000[2]), col ='blue', add = TRUE,lwd=2)
curve(dnorm(x, para.1986.2015[1], para.1986.2015[2]), col ='red', add = TRUE,lwd=2)
abline(v=mean(epw.tas.seas[3]))
text(x=21,y=0.08,'YVR 2015')
legend('topleft',leg=c('1971-2000','1986-2015'),col=c('blue','red'),pch=16)
box(which='plot')

##dev.off()