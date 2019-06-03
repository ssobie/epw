##Means of testing read and write of epw data


read.epw.file <- function(epw.dir,epw.file) {
   
   ##epw.dir <- '/storage/home/ssobie/code/repos/bc-projected-weather/bcweather/tests/data/'
   ##epw.file <- 'CAN_BC_ABBOTSFORD-A_1100031_CWEC.epw'
   epw.filename <- paste0(epw.dir,epw.file)

   epw.data <- read.csv(epw.filename,skip=8,header=F,as.is=T)
   epw.header <- readLines(epw.filename,n=8)
   rv <- list(data=epw.data,header=epw.header)
}

write.table_with_header <- function(x, file, header, ...){
  cat(header[1], '\n',  file = file)                          
  for (l in 2:length(header)) {                         
    cat(header[l], '\n',  file = file,append=T)
  }
  write.table(x, file, append = T, col.names=F,row.names=F,quote=F,...)
}

write.epw.file <- function(epw.data,epw.header,epw.dir,epw.file) {
   
   epw.filename <- paste0(epw.dir,epw.file)

   ##Format the different columns as needed
   epw.data[,7] <- format(epw.data[,7],nsmall=1,trim=TRUE)
   epw.data[,8] <- format(epw.data[,8],nsmall=1,trim=TRUE)
   epw.data[,22] <- format(epw.data[,22],nsmall=1,trim=TRUE)
   epw.data[,25] <- format(epw.data[,25],nsmall=1,trim=TRUE)
   epw.data[,28] <- sprintf('%08d',epw.data[,28])
   epw.data[,30] <- format(epw.data[,30],nsmall=4,trim=TRUE)
   epw.data[,33] <- format(epw.data[,33],nsmall=3,trim=TRUE)
   epw.data[,34] <- format(epw.data[,34],nsmall=1,trim=TRUE)
   epw.data[,35] <- format(epw.data[,35],nsmall=1,trim=TRUE)

   write.table_with_header(epw.data,epw.filename,epw.header,sep=',')
}
