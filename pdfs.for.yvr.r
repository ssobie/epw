##Script to generate an hourly temperature series by combining BCCAQ2-PRISM 
##with the hourly temperatures from the YVR EPW File

library(MASS)
library(ncdf4)
###source('/storage/home/ssobie/code/repos/assessments/single.site.bccaq.to.prism.r',chdir=T)
source('/storage/home/ssobie/code/repos/epw/read.write.epw.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##-----------------------------------------------------------------
##YVR Coordinates
lonc <- -123.1815
latc <- 49.1967

##-----------------------------------------------------------------------

##Degree Day Values
dd <- function(tmean,tbase) {
  g <- tmean - tbase
  days <- sum(g[g>0], na.rm=T)
  return(round(days))
}

##-----------------------------------------------------------------------
##MAke the cartoon

make_cartoon <- function(series.past,series.future,epw.lines,epw.line.hourly=NULL,
                         xlim,ylim,xlab,ylab='Density',brk,title,
                         plot.file) {
   norm.past <- fitdistr(unlist(series.past),'normal')
   norm.future <- fitdistr(unlist(series.future),'normal')
   para.past <- norm.past$estimate
   para.future <- norm.future$estimate

   png(file=plot.file,width=4,height=2.5,units='in',res=600,pointsize=6,bg='white')
   plot(c(),xaxs='i',yaxs='i',xlim=xlim,ylim=ylim,main=title,
           xlab=xlab,
           ylab=ylab)
   breaks <- seq(0,ceiling(max(c(unlist(series.past),unlist(series.future)),na.rm=T)/10)*10,by=brk)

   hist(breaks=breaks,unlist(series.past),border='blue',add=T,prob=TRUE)
   hist(breaks=breaks,unlist(series.future),border='red',add=T,prob=TRUE)

   curve(dnorm(x, para.past[1], para.past[2]), from=0,col ='blue', add = TRUE,lwd=2)
   curve(dnorm(x, para.future[1], para.future[2]),from=0, col ='red', add = TRUE,lwd=2)
   abline(v=epw.lines)
   text(x=epw.lines[1],y=0.75*ylim[2],'YVR 2015')
   text(x=epw.lines[2],y=0.8*ylim[2],'YVR 2017')
   if (!is.null(epw.line.hourly)) {
      abline(v=epw.line.hourly,lty=2)
      text(x=epw.line.hourly[1],y=0.75*ylim[2],'YVR 2015 Hourly')
      text(x=epw.line.hourly[2],y=0.85*ylim[2],'YVR 2017 Hourly')
   }
   legend('topleft',leg=c('1971-2000','1986-2015'),col=c('blue','red'),pch=16)
   box(which='plot')

   dev.off()

}

##-----------------------------------------------------------------

retrieve_derived_series <- function(var.name,scenario,
                                    type,gcm,
                                    lonc,latc,
                                    base.dir) {

  read.dir <- paste0(base.dir,gcm,'/',scenario,'/',type,'/')
  var.file <- list.files(path=read.dir,pattern=var.name)
  nc <- nc_open(paste0(read.dir,var.file))
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  lon.ix <- which.min(abs(lonc-lon))
  lat.ix <- which.min(abs(latc-lat))

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
##Read YVR cell and save as RData (for 12 GCMs)

create_bccaq2_derived_series <- function(var.name,type) {
   bccaq2.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/high_res_downscaling/bccaq_gcm_bc_subset/'
   save.dir <- '/storage/data/projects/rci/weather_files/YVR/'

   gcm.list <- c('ACCESS1-0','CCSM4','CanESM2','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')

   bccaq2.file <- paste0(save.dir,var.name,'_',type,'_YVR_BCCAQ2-PRISM_1951-2100.RData')
   if (!file.exists(bccaq2.file)) {
   bccaq2.series <- vector(mode='list',length=length(gcm.list))
      for (g in seq_along(gcm.list)) {
         gcm <- gcm.list[g]
         bccaq2.series[[g]] <- retrieve_derived_series(
                                    var.name=var.name,scenario='rcp85',
                                    type=type,
                                    gcm=gcm,
                                    lonc=lonc,latc=latc,
                                    base.dir=bccaq2.dir)     
      }
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
epw.hours <- epw.data$data[,3]
epw.dates <- as.Date(paste(epw.years,epw.mons,epw.days,sep='-'))
epw.hour.dates <- as.Date(paste0(epw.years,'-',epw.mons,'-',epw.days,' ',epw.hours,':00:00'))
day.fac <- as.factor(format(epw.dates,'%Y-%m-%d'))
epw.tas.daily <- tapply(epw.tas,day.fac,mean)
epw.tas.max <- tapply(epw.tas,day.fac,max)

mon.fac <- as.factor(format(epw.dates,'%m'))
seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
seas.fac <- factor(seasons[mon.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))

epw.tas.seas <- tapply(epw.tas,seas.fac,mean)


epw.file.2017 <- 'CAN_BC_Vancouver.Intl.AP.718920_2017.epw'
epw.data.2017 <- read.epw.file(epw.dir,epw.file.2017)
epw.tas.2017 <- epw.data.2017$data[,7]
epw.tas.2017.daily <- tapply(epw.tas.2017,day.fac,mean)
epw.tas.2017.max <- tapply(epw.tas.2017,day.fac,max)

##-----------------------------------------------------------------
##Cooling Degree Days

cdd.series <- create_bccaq2_derived_series(var.name='cdd',type='degree_days') 
epw.cdd <- dd(epw.tas.daily,18)
epw.cdd.hourly <- dd(epw.tas,18)/24

epw.cdd.2017 <- dd(epw.tas.2017.daily,18)
epw.cdd.2017.hourly <- dd(epw.tas.2017,18)/24

cdd.1971.2000 <- lapply(cdd.series,subset_dates,'1971-2000')
cdd.1986.2015 <- lapply(cdd.series,subset_dates,'1986-2015')

cdd.1971.2000.data <- lapply(cdd.1971.2000,function(x){x$data})
cdd.1986.2015.data <- lapply(cdd.1986.2015,function(x){x$data})

plot.file <- '/storage/data/projects/rci/weather_files/cdd.cartoon.yvr.2015.png'
xlab <- 'Cooling Degree Days (Degree Days)'

make_cartoon(cdd.1971.2000.data,cdd.1986.2015.data,c(epw.cdd,epw.cdd.2017),c(epw.cdd.hourly,epw.cdd.2017.hourly),
             xlim=c(0,300),ylim=c(0,0.035),xlab,ylab='Density',brk=5,title='Cooling Degree Days',
             plot.file)


series.file <- '/storage/data/projects/rci/weather_files/cdd.series.yvr.2015.png'

png(file=series.file,width=5,height=2.5,units='in',res=600,pointsize=6,bg='white')
plot(x=1971:2015,y=rep(0,45),col='white',yaxs='i',xlim=c(1971,2015),ylim=c(0,310),
         xlab='Years',main='Cooling Degree Days at YVR (Modelled)',
         ylab='CDD')

lapply(cdd.1986.2015,function(x,y){points(y,x$data,col='red',pch=18)},1986:2015)
lapply(cdd.1971.2000,function(x,y){points(y,x$data,cex=0.8,pch=18)},1971:2000)

abline(h=epw.cdd,col='gray')
##abline(h=epw.cdd.hourly,col='gray',lty=2)
text(x=2010,y=epw.cdd,'2015')
abline(h=epw.cdd.2017,col='gray')
text(x=2010,y=epw.cdd.2017,'2017')

#abline(h=epw.cdd.hourly,col='gray',lty=2)
#abline(h=epw.cdd.2017.hourly,col='gray',lty=2)



dev.off()



##-----------------------------------------------------------------
##SU25

su.series <- create_bccaq2_derived_series(var.name='suETCCDI',type='climdex') 
epw.su <- sum(epw.tas.max >= 25)
epw.su.2017 <- sum(epw.tas.2017.max >= 25)

su.1971.2000 <- lapply(su.series,subset_dates,'1971-2000')
su.1986.2015 <- lapply(su.series,subset_dates,'1986-2015')

su.1971.2000.data <- lapply(su.1971.2000,function(x){x$data})
su.1986.2015.data <- lapply(su.1986.2015,function(x){x$data})

plot.file <- '/storage/data/projects/rci/weather_files/su.cartoon.yvr.2015.png'
xlab <- 'Days above 25\u00B0C'

make_cartoon(su.1971.2000.data,su.1986.2015.data,c(epw.su,epw.su.2017),epw.line.hourly=NULL,
             xlim=c(0,60),ylim=c(0,0.09),xlab,ylab='Density',brk=1,title='No. Days above 25 degC',
             plot.file)

series.file <- '/storage/data/projects/rci/weather_files/su.series.yvr.2015.png'

png(file=series.file,width=4,height=2.5,units='in',res=600,pointsize=6,bg='white')
plot(x=1971:2015,y=rep(0,45),col='white',yaxs='i',xlim=c(1971,2015),ylim=c(0,65),
         xlab='Years',main='No. Days above 25degC at YVR (Modelled)',
         ylab='SU')

lapply(su.1986.2015,function(x,y){points(y,x$data,col='red',pch=18)},1986:2015)
lapply(su.1971.2000,function(x,y){points(y,x$data,cex=0.8,pch=18)},1971:2000)

abline(h=epw.su,col='gray')
text(x=2010,y=epw.su-3,'2015')
abline(h=epw.su.2017,col='gray')
text(x=2010,y=epw.su.2017+3,'2017')

dev.off()

tas.file <- '/storage/data/projects/rci/weather_files/tas.yvr.2015.png'
png(file=tas.file,width=4,height=2.5,units='in',res=600,pointsize=6,bg='white')
plot(epw.hour.dates,epw.tas,col='black',yaxs='i',ylim=c(-6,30),pch=18,
         xlab='Date',cex=0.5,
         ylab='TAS (degC)')
abline(h=25,lwd=0.75)
dev.off()


##-----------------------------------------------------------------
##

epw.file <- 'CAN_BC_Vancouver.Intl.AP.718920_2017.epw'
epw.dir <- '/storage/data/projects/rci/weather_files/wx_2016/'
epw.data <- read.epw.file(epw.dir,epw.file)
epw.tas <- epw.data$data[,7]
epw.years <- rep(1995,length(epw.tas)) ##epw.data$data[,1]
epw.mons <- epw.data$data[,2]
epw.days <- epw.data$data[,3]
epw.hours <- epw.data$data[,3]

epw.dates <- as.Date(paste(epw.years,epw.mons,epw.days,sep='-'))

epw.hour.dates <- as.Date(paste0(epw.years,'-',epw.mons,'-',epw.days,' ',epw.hours,':00:00'))

tas.file <- '/storage/data/projects/rci/weather_files/tas.yvr.2017.png'
png(file=tas.file,width=5,height=2.5,units='in',res=600,pointsize=6,bg='white')
plot(epw.hour.dates,epw.tas,col='black',yaxs='i',ylim=c(-10,30),pch=18,
         xlab='Date',cex=0.5,
         ylab='TAS (degC)')
abline(h=25,lwd=0.75)
dev.off()
