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

get_beta_series <- function(epw.tas,epw.dates,tasmax.data,tasmin.data) {

tas <- epw.tas
time <- epw.dates

days <- 31
st <- seq(1,by=24,length.out=days) 
en <- seq(24,by=24,length.out=days) 

time.day <- levels(day.fac)

tas.mean <- (tasmax.data+tasmin.data)/2 ##tapply(epw.tas,day.fac,mean)
tas.max <-  tasmax.data ##tapply(epw.tas,day.fac,max)
tas.min <-  tasmin.data ##tapply(epw.tas,day.fac,min)
tas.var <-  tapply(epw.tas,day.fac,var)[1:31]

##Tas Var must be less than (tas.mean -tas.min) * (tas.max - tas.mean) for beta dist
tas.check <- (tas.mean -tas.min) * (tas.max - tas.mean)
flag <- tas.check < as.numeric(tas.var)

tas.var[flag] <- tas.check[flag] * 0.25

tas.range <- tas.max-tas.min
tas.skew <- (tas.mean-tas.min)/tas.range

midday <- time[grep('12:00:00',time)]

mids <- daily_start_and_end_values(tas.mean[1:days],time[1:days]) 

mid.c <- c(midday[1:days],time[st],time[en])
mid.ord <- rank(mid.c)
mid.vals <- mid.ord*0
mid.vals[mid.ord] <- c(tas.mean[1:days],mids$start,mids$end)

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

tas.anom <- epw.tas*0
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
beta.series <- hours*0

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

   ##print(paste0('Max diff ',max(beta.ord)-max(beta.past)))
   ##print(paste0('Mean diff ',mean(beta.ord)-mean(beta.past)))
   ##print(paste0('Min diff ',min(beta.ord)-min(beta.past)))

   day <- time[grep(time.day[i],time)]
   ###lines(day,beta.ord,col='goldenrod',lwd=2)
   ###lines(which(hours==i),beta.ord,col='goldenrod',lwd=2)
   beta.series[which(hours==i)] <- beta.ord
}

return(beta.series)

}


##-------------------------------------------



