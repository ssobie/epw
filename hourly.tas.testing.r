##Script to test methods of simulating hourly temperature

library(ncdf4)
library(fBasics)
library(MASS)
library(pracma)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

rgbeta <- function(n, mean, var, min = 0, max = 1)
{
  dmin <- mean - min
  dmax <- max - mean

  if (dmin <= 0 || dmax <= 0)
  {
    stop(paste("mean must be between min =", min, "and max =", max)) 
  }

  if (var >= dmin * dmax)
  {
    stop(paste("var must be less than (mean - min) * (max - mean) =", dmin * dmax))
  }

  # mean and variance of the standard beta distributed variable
  mx <- (mean - min) / (max - min)
  vx <- var / (max - min)^2

  # find the corresponding alpha-beta parameterization
  a <- ((1 - mx) / vx - 1 / mx) * mx^2
  b <- a * (1 / mx - 1)

  # generate standard beta observations and transform
  x <- rbeta(n, a, b)
  y <- (max - min) * x + min

  return(y)
}

##-----------------------------------------------------------------------

daily_start_and_end_values <- function(tas.mean,time) {
  tlen <- length(time)
  tas.start <- tas.end <- rep(0,tlen)
  tas.start[1] <- tas.mean[1]
  tas.end[tlen] <- tas.mean[tlen]
  x <- 1:24
  for (i in 2:tlen) {
     tas.slope <- (tas.mean[i] - tas.mean[i-1])/24
     tas.line <- tas.slope*x + tas.mean[i-1]
     tas.start[i] <- tas.line[13]
     tas.end[i-1] <- tas.line[12]     
  }

  rv <- list(start=tas.start,end=tas.end)
  return(rv)
}

##-----------------------------------------------------------------------


nc <- nc_open('/storage/data/climate/observations/reanalysis/ERA5/tas_hour_ERA5_BC_19800101-20181231.nc')
tas <- ncvar_get(nc,'tas',start=c(50,25,1),count=c(1,1,-1))
time <- netcdf.calendar(nc)
days <- 31
st <- seq(1,by=24,length.out=days) ##seq(1,145,24)
en <- seq(24,by=24,length.out=days) ##seq(24,168,24)

daily.fac <- as.factor(format(time,'%Y-%m-%d'))
time.day <- levels(daily.fac)

tas.mean <- tapply(tas,daily.fac,mean)
tas.max <-  tapply(tas,daily.fac,max)
tas.min <-  tapply(tas,daily.fac,min)
tas.var <-  tapply(tas,daily.fac,var)
tas.range <- tas.max-tas.min
tas.skew <- (tas.mean-tas.min)/tas.range

midday <- time[grep('12:00:00',time)]

plot.dir <- '/storage/data/projects/rci/weather_files/hourly_testing/'
plot.file <- paste0(plot.dir,'hourly_disaggregation_testing.png')
###png(file=plot.file,width=10,height=6,units='in',res=600,pointsize=6,bg='white')

###plot(c(time[1],time[days*24]),c(0,0),ylim=c(-35,0),col='white',cex.axis=2,cex.lab=2,
###        xlab='Date',ylab='Temperature')

for (d in 1:days) {
    day <- time[grep(time.day[d],time)]
###    lines(day,rep(tas.mean[d],length(day)))
###    lines(day,rep(tas.max[d],length(day)),col='red')
###    lines(day,rep(tas.min[d],length(day)),col='blue')
}
###lines(time[1:(days*24)],tas[1:(days*24)],col='green')
###lines(midday[1:days],tas.mean[1:days],col='orange')

mids <- daily_start_and_end_values(tas.mean[1:days],time[1:days]) 
###points(time[st],mids$start,pch=24)
###points(time[en],mids$end,pch=25)

##Try a loess curve to the mean temperatures and mid values

mid.c <- c(midday[1:days],time[st],time[en])
mid.ord <- rank(mid.c)
mid.vals <- mid.ord*0
mid.vals[mid.ord] <- c(tas.mean[1:days],mids$start,mids$end)

##loess.data <- list(x=1:length(mid.c),y=mid.vals)
##model <- lm(y~poly(x,10),data=loess.data)
##predicted <- predict(model,newdata=data.frame(x=seq(1,21,length.out=168)))
##lines(time[1:168],predicted)

##Add autocorrelation by simulating mvrnorm with acf lags from
##historical days and reranking the sorted beta dist values 
##using the ranks from the mvrnorm series.
##Still need to ensure continuity at the day boundaries and 
##Replace the max and min from the beta with the actual max
##and minimum values

##Use a fitted exponential distribution relating temperature variance
##to temperature range to produce the required variance for the beta
##distribution

##A potential way to ensure continuity is to fix the 1st and 24th values
##for each day (using the linear fit between daily averages) and generate
##22 values from the lagged tas and Beta distribution. Then fix the max
##and min, and adjust the remaining 20 values to match the mean required.



tas.anom <- tas[1:168]*0
len <- 23
tas.acf <- matrix(0,nrow=days,ncol=len)
for (i in 1:days) {
  tas.anom[st[i]:en[i]] <- tas[st[i]:en[i]] - mean(tas[st[i]:en[i]])
  a <- acf(tas[1:744],plot=FALSE,lag.max=23) ##[st[i]:en[i]],plot=FALSE,lag.max=23)
  tas.acf[i,] <- a$acf[1:len]
}

##avg.acf <- apply(tas.acf,2,mean)
tas.lag <- matrix(0,nrow=days,ncol=24)
for (l in 1:days) {
   acm <- diag(rep(1,len))
   for (m in 1:len) {
      acm[m:len,m] <- tas.acf[l,1:(len-m+1)]
      acm[m,m:len] <- tas.acf[l,1:(len-m+1)]
   }
   tas.lag[l,] <- c(tas.mean[l],mvrnorm(n=1, mu=rep(tas.mean[l],len), Sigma=acm))
   day <- time[grep(time.day[l],time)]


}

set.seed(1)
n <- 10
beta.list <- vector(mode='list',length=days)
beta.dist <- matrix(0,nrow=10,ncol=days)
hours <- rep(1:31,each=24)
days <- 31
for (i in 1:days) {
   beta.matrix <- matrix(0,nrow=n,ncol=24)
   for (j in 1:n) {
      beta.matrix[j,] <- rgbeta(24, mean = tas.mean[i], 
                            var = tas.var[i], 
                            min = tas.min[i], 
                            max = tas.max[i])
   }
   beta.list[[i]] <- beta.matrix

   beta.max <- apply(beta.matrix,1,max)
   beta.min <- apply(beta.matrix,1,min)
   beta.mean <- apply(beta.matrix,1,mean)
 
   beta.dist[,i] <- (beta.max-tas.max[i])^2 + 
               (beta.min-tas.min[i])^2 +
               (beta.mean-tas.mean[i])^2
   beta.ix <- which.min(beta.dist[,i])               
   beta.slct <- beta.matrix[beta.ix,]           
   lag.rank <- rank(tas.lag[i,])   
   beta.sort <- sort(beta.slct)
   ##beta.sort[1] <- tas.min[i]
   ##beta.sort[24] <- tas.max[i]
   beta.ord <- beta.sort[lag.rank]
   beta.past <- beta.ord
   beta.ord[1] <- mids$start[i]
   beta.ord[24] <- mids$end[i]

   beta.ord[which.max(beta.ord)] <- tas.max[i]
   beta.ord[which.min(beta.ord)] <- tas.min[i]

   print(paste0('Max diff ',max(beta.ord)-max(beta.past)))
   print(paste0('Mean diff ',mean(beta.ord)-mean(beta.past)))
   print(paste0('Min diff ',min(beta.ord)-min(beta.past)))

   day <- time[grep(time.day[i],time)]
   ###lines(day,beta.ord,col='goldenrod',lwd=2)
   lines(which(hours==i),beta.ord,col='goldenrod',lwd=2)
   print('Skew compare')
   print(paste0('TAS Skew: ',tas.skew[i]))
   print(paste0('Beta Skew: ',(mean(beta.ord)-min(beta.ord))/(max(beta.ord)-min(beta.ord))))

}


##PCHIP method mimmicking the MetSim python package

tas.min.time <- paste0(names(tas.min[1:days]),' 05:00:00')
tas.max.time <- paste0(names(tas.max[1:days]),' 17:00:00')

tas.ix <- order(c(tas.min.time,tas.max.time))
tas.min.max <- c(tas.min[1:days],tas.max[1:days])[tas.ix]
int.time <- time[1:(days*24)]
tas.time <- c(tas.min.time,tas.max.time)[tas.ix]

tas.xi <- seq_along(tas.time)
for (i in seq_along(tas.time)) {
    tas.xi[i] <- which(time %in% as.PCICt(tas.time[i],cal="proleptic_gregorian"))
}

metsim <- pchip(xi=tas.xi,yi=tas.min.max,
                x=seq_along(int.time))

###lines(time[1:(days*24)],metsim,col='purple')

###box(which='plot')
###dev.off()        
##-------------------------------------------
##Add other methods of disaggregation for comparison

##Adding diurnal cycle from EPW file
##MetSim using Hermite Polynomials
##Sinusoidal fitting






##-------------------------------------------




##January 1 only data 
if (1==0) {
ix <- grep('[0-9]{4}-01-01 *',time)
tas.jan <- tas[ix]
st <- seq(1,913,24)
en <- seq(24,936,24)

jan.var <- jan.max <- jan.min <- jan.mean <- jan.skew <- rep(0,39)
jan.anom <- tas.jan*0
len <- 23
jan.acf <- matrix(0,nrow=39,ncol=len)

for (i in 1:39) {
  jan.max[i] <- max(tas.jan[st[i]:en[i]])
  jan.min[i] <- min(tas.jan[st[i]:en[i]])
  jan.var[i] <- var(tas.jan[st[i]:en[i]])
  jan.mean[i] <- mean(tas.jan[st[i]:en[i]])
  jan.anom[st[i]:en[i]] <- tas.jan[st[i]:en[i]] - mean(tas.jan[st[i]:en[i]])
  a <- acf(tas.jan[st[i]:en[i]],plot=FALSE,lag.max=23)
  jan.acf[i,] <- a$acf[1:len]
}

jan.skew <- (jan.mean-jan.min)/(jan.max-jan.min)
avg.acf <- apply(jan.acf,2,mean)
acm <- diag(rep(1,len))

for (m in 1:len) {
  acm[m:len,m] <- avg.acf[1:(len-m+1)]
  acm[m,m:len] <- avg.acf[1:(len-m+1)]
}

##x <- mvrnorm(jan.mean[1], rep(0,len), acm)
x <- mvrnorm(n=1, mu=rep(jan.mean[1],len), Sigma=acm)



set.seed(1)
n <- 10
beta.matrix <- matrix(0,nrow=n,ncol=24)

for (j in 1:n) {
 beta.matrix[j,] <- rgbeta(24, mean = jan.mean[1], 
                            var = jan.var[1], 
                            min = jan.min[1], 
                            max = jan.max[1])
}

beta.max <- apply(beta.matrix,1,max)
beta.min <- apply(beta.matrix,1,min)
beta.mean <- apply(beta.matrix,1,mean)

beta.dist <- (beta.max-jan.max[1])^2 + 
             (beta.min-jan.min[1])^2 +
             (beta.mean-jan.mean[1])^2

}

nc_close(nc)