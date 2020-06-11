## Create figures to embed in the EPW summary excel table

library(zoo)

##-----------------------------------------------------------

new_box <- function(at,data,add,boxwex,axes,col) { #     new.box(at=i+0.2,seas.2080s[[i]],add=T,boxwex=0.25,axes=F)##col='red',
  bb <- boxplot(data,plot=FALSE)
  bb$stats[c(1,5),] <- quantile(data,c(0.1,0.9),na.rm=TRUE)
  bxp(z=bb,at=at,add=add,boxwex=boxwex,axes=axes,boxfill=col)
}
  
mean_temperature_boxplot <- function(cwec.2020s,cwec.2050s,cwec.2080s,
                                     site,fig.dir) {

   months <- 1:12
   tas.ix <- grep('tas_',names(cwec.2020s$past))   
   cwec.tas <- c(cwec.2020s$past[tas.ix],cwec.2020s$proj[,tas.ix],
                 cwec.2050s$past[tas.ix],cwec.2050s$proj[,tas.ix],
                 cwec.2080s$past[tas.ix],cwec.2080s$proj[,tas.ix])
   yvals <- pretty(cwec.tas) ##c(-5,40,5)
   
   plot.file <- paste0(fig.dir,site,'_mean_temperature_boxplots.png')    
   plot.title <- 'Monthly Average Temperatures'
   y.label <- 'Temperature (\u00B0C)'


   png(file=plot.file,width=7,height=4,units='in',res=600,pointsize=6,bg='white')
   par(mar=c(5,5,5,5))
   plot(c(),xlim=c(1,12),ylim=c(min(yvals),max(yvals)),main=plot.title,
       xlab='',ylab=y.label,cex.axis=2,cex.lab=2,cex.main=2.5,axes=F,yaxs='i')

   axis(1,at=1:12,label=month.abb,cex.axis=2)
   axis(2,at=yvals,yvals,cex.axis=2)

   for (i in seq_along(months)) {

      var.name <- paste0('tas_',tolower(month.abb[i]))
      ix <- grep(var.name,names(cwec.2020s$past))
      past.tas <- cwec.2020s$past[ix]
      lines(x=c(i-0.5,i+0.5),y=c(past.tas,past.tas),col='blue',lwd=2,lend=1)

      new_box(at=i-0.2,cwec.2020s$proj[,ix],add=T,boxwex=0.25,axes=F,col='goldenrod')
      new_box(at=i,cwec.2050s$proj[,ix],add=T,boxwex=0.25,axes=F,col='orange')
      new_box(at=i+0.2,cwec.2080s$proj[,ix],add=T,boxwex=0.25,axes=F,col='red')
    }

    abline(v=seq(0.5,12.5,by=1))
    ##grid(ny=NULL,nx=NA,col='gray',lwd=1)
    box(which='plot')
    legend('topright',legend=c('Past','2020s','2050s','2080s'),col=c('blue','goldenrod','orange','red'),cex=2,pch=15)

    dev.off()  

    rv <- plot.file
    return(rv)
}

##-----------------------------------------------------------------------------------
##**********************************************************************************


make_summary_figures <- function(cwec.2020s,cwec.2050s,cwec.2080s,
                                 site,fig.dir) {

   ##Mean temperature
   plot.file <- mean_temperature_boxplot(cwec.2020s,cwec.2050s,cwec.2080s,
                                         site,fig.dir)
   rv <- list(tas=plot.file)
   return(rv)

}