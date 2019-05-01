##Script to sort out the data requested for Okanagan

source('/storage/home/ssobie/code/repos/assessments/summary.table.comments.r',chdir=T)

library(openxlsx)

##------------------------------------------------------------- 
##------------------------------------------------------------- 

get.units.pane <- function(var.name) {
  leg.label <- NA
  if (grepl("(tas|txx|tnn|txn|tnx|tmax|tmin|dtr)", var.name))
    leg.label <- rep('degC',14)
  if (grepl("(pr|rx|r10|r20|r9|RP|sdii|prcptot)", var.name))
    leg.label <- rep('mm',14)
  if (grepl("(dd)", var.name))
    leg.label <- rep('degree days',14)
  if (grepl("(fdE|cddE|cdd90|cddmax|cwd|su|gsl|id|trE|su30|r95daysE|r99daysE)", var.name))
    leg.label <- rep('days',14)
  if (grepl("(dtr)", var.name))
    leg.label <- rep('degC',14)
  if (grepl("(r95sep|r95dist|r95days|r99days)", var.name))
    leg.label <- rep('days',14)
  return(leg.label)
} 

##------------------------------------------------------------- 
get.round.val <- function(var.name) {
  rd <- 0
  if (grepl("(dd)", var.name))
    rd <- 0    
  if (grepl("(tas|txx|tnn|tnx|txn|tmax|tmin)", var.name))
    rd <- 1
  if (grepl("(pr|rx|r9|RP|rp|tx90|tn10|trE|cddE|cdd90|cddmax|cwdE|idE|dtrE|wsdiE|csdiE|r95sep)", var.name))
    rd <- 0
  if (grepl("(pas|snowdepth)", var.name))
    rd <- 0
  return(rd)
} 

##------------------------------------------------------------- 
##Separate variables into precip and temperature
## 'E' denotes a climdex index
##Set the order that the variables should be arranged here.

filter.input.variables <- function(table.vars) {
  all.vars <- unlist(lapply(table.vars,function(x){return(x[1])}))
  pr.vars <- c('pr',
               'prcptotETCCDI','sdiiETCCDI','r10mmETCCDI','r20mmETCCDI',
               'rx1dayETCCDI','rx2dayETCCDI','rx5dayETCCDI',
               'r95pETCCDI','r95daysETCCDI','r99pETCCDI','r99daysETCCDI',
               'pr.maximum','pr.minimum','pr.standard_deviation',
               'pr_rp5','pr_rp20','pr_rp50',
               'cddETCCDI','cdd90ETCCDI','cddmaxETCCDI','cwdETCCDI') 
  pr.ix <- match(pr.vars,all.vars)
  pr.selected <- table.vars[pr.ix[!is.na(pr.ix)]]

  tas.vars <- c('tasmax','tas','tasmin',
                'txxETCCDI','tnnETCCDI','txnETCCDI','tnxETCCDI',
                'dtrETCCDI','suETCCDI','su30ETCCDI','trETCCDI',
                'idETCCDI','fdETCCDI','gslETCCDI',
                'cdd','gdd','hdd','fdd',
                'tasmax.annual_quantile_975','tasmax.annual_quantile_990','tasmax.annual_quantile_996',                        
                'tasmin.annual_quantile_004','tasmin.annual_quantile_010','tasmin.annual_quantile_025',
                'wetbulb.annual_quantile_975','wetbulb.annual_quantile_990','wetbulb.annual_quantile_996',
                'wetbulb.annual_quantile_004','wetbulb.annual_quantile_010','wetbulb.annual_quantile_025',
                'tasmax_rp5','tasmax_rp20',
                'tasmin_rp5','tasmin_rp20')                          
  tas.ix <- match(tas.vars,all.vars)
  tas.selected <- table.vars[tas.ix[!is.na(tas.ix)]]
  rv <- list(pr=pr.selected,tas=tas.selected)
  return(rv)
}

##------------------------------------------------------------- 

find.row.locations <- function(sorted.vars) {

   tas.vars <- sorted.vars$tas
   tas.rows <- vector(length=length(tas.vars)+1,mode='list')

   tas.rows[[1]] <- 3:5 ##Header Rows
   for (i in 1:length(tas.vars)) {
      row.len <- 2 ##switch(tas.vars[[i]][[2]],annual=2,seasonal=6)
      tas.rows[[i+1]] <- seq(tail(tas.rows[[i]],1)+1,length.out=row.len)
   }

   if (length(sorted.vars$pr!=0)) {
     pr.start <- tail(tas.rows[[length(tas.vars)+1]],1)
     pr.vars <- sorted.vars$pr
     pr.rows <- vector(length=length(pr.vars)+1,mode='list')
     pr.rows[[1]] <- (1:3)+pr.start
     for (i in 1:length(pr.vars)) {
        row.len <- switch(pr.vars[[i]][[2]],annual=2,seasonal=6)
        pr.rows[[i+1]] <- seq(tail(pr.rows[[i]],1)+1,length.out=row.len)
     }
   } else {
     pr.rows <- NULL
   }
    
   rv <- list(pr=pr.rows,tas=tas.rows)
   return(rv)

}
##---------------------------------------------

seasonal.table <- function(var.name,scenario='rcp85',rp=NULL) {

  no.percent <- '(wb|tas|tasmax|tasmin|txxETCCDI|tnnETCCD|tnxETCCDI|txnETCCDI|trETCCDI|suETCCDI|su30ETCCDI)'
  result <- get.seasonal.data(var.name,scenario)
  if (grepl(no.percent,var.name))
    result[,10:15] <- 'NA'
  return(as.data.frame(result))
}


annual.table <- function(var.name,scenario='rcp85',rp=NULL) {
  
  no.percent <- '(wb|tas|tasmax|tasmin|txxETCCDI|tnnETCCD|trETCCDI|r95sep|r99days|r95days)'
  result <- get.annual.data(var.name,scenario,rp) 

  if (grepl(no.percent,var.name))
    result[10:15] <- 'NA'
  return(as.data.frame(t(result)))
}


##----------------------------------------------------------------------------------------
write.variables <- function(wb,sheet,sorted.vars,row.locs,type,cwec.col,gcm.stats.list) {
  len <- length(sorted.vars)
  for (i in 1:len) {               
    var.info <- sorted.vars[[i]]
    var.name <- var.info$name
    print(var.name)
    season <- grepl('season',var.info$type)
    var.title <- var.info$title
    print(var.title)
    current.row <- row.locs[[i+1]]
    var.entry <- c(var.info$type,cwec.col[[var.name]],gcm.stats.list[[var.name]])
    
    s1 <- createStyle(fontSize = 12, fontColour = "black", textDecoration = c("BOLD"))
    s2 <- createStyle(fontSize = 12, fontColour = "black")
    c1 <- createComment(comment = variable.comment(var.name),style=c(s1,s2),visible=FALSE)
    writeComment(wb, sheet, col = 1, row = current.row[1], comment = c1)

    pane.colour <- switch(type,pr='lightblue',tas='tan1')             

    hdstyle <- createStyle(fgFill = pane.colour, halign = "CENTER", textDecoration = "Bold",
                             border = "TopBottomLeftRight", fontColour = "black")  
    units.pane <- c(var.title,get.units.pane(var.name))                             
    writeData(wb, sheet=sheet, as.data.frame(t(units.pane)), startRow = current.row[1], startCol = 1, headerStyle = hdstyle,
              borders = "rows", borderStyle = "medium",colNames=FALSE)
    addStyle(wb,sheet=sheet,hdstyle,rows=current.row[1],cols=1:15,gridExpand=FALSE,stack=FALSE)

    datastyle <- createStyle(fgFill = 'white', halign = "RIGHT",
                             border = "TopBottomLeftRight", fontColour = "black")                              
    writeData(wb, sheet=sheet, as.data.frame(t(var.entry)), startRow = current.row[2], startCol =1, headerStyle = hdstyle,
              borders = "rows", borderStyle = "medium",colNames=FALSE)
    addStyle(wb,sheet=sheet,datastyle,rows=current.row[-1],cols=2:15,gridExpand=TRUE,stack=FALSE)

    highlight <- createStyle(fgFill = 'lightyellow', halign = "CENTER", 
                             border = "TopBottomLeftRight", fontColour = "black")  
    addStyle(wb,sheet=sheet,highlight,rows=current.row[-1],cols=3,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row[-1],cols=6,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row[-1],cols=7,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row[-1],cols=10,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row[-1],cols=11,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row[-1],cols=14,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row[-1],cols=15,gridExpand=FALSE,stack=FALSE)

  }
     
}
##----------------------------------------------------------------------------------------
      ##Top Frozen Pane

create.frozen.top.pane <- function(wb,sheet) {

      pane.titles <- list(' ','CWEC',
                              'Model Past',
                              '2020s Change',' ',
                              '2020s Future',' ',
                              '2050s Change',' ',
                              '2050s Future',' ',
                              '2080s Change',' ',
                              '2080s Future',' ')

      fz1 <- createStyle(fgFill = "gray94", halign = "CENTER", textDecoration = "Bold",
                         border = "TopBottomLeftRight", fontColour = "black")
      writeData(wb, sheet=sheet, pane.titles, startRow = 1, startCol = 1, headerStyle = fz1,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=sheet,fz1,rows=1,cols=1:15,gridExpand=FALSE,stack=FALSE)
      mergeCells(wb,sheet=sheet,cols=c(4,5),rows=1)
      mergeCells(wb,sheet=sheet,cols=c(6,7),rows=1)
      mergeCells(wb,sheet=sheet,cols=c(8,9),rows=1)
      mergeCells(wb,sheet=sheet,cols=c(10,11),rows=1)
      mergeCells(wb,sheet=sheet,cols=c(12,13),rows=1)
      mergeCells(wb,sheet=sheet,cols=c(14,15),rows=1)
      ##freezePane(wb,sheet=sheet,firstRow=TRUE)

      prct.header <- list(' ',' ',' ','Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                                  'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%')          
      writeData(wb, sheet=sheet, prct.header, startRow = 2, startCol = 1, headerStyle = fz1,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=sheet,fz1,rows=2,cols=1:15,gridExpand=FALSE,stack=FALSE)
      ##freezePane(wb,sheet=sheet,firstActiveRow=3)
}

##-----------------------------------------------------------------------------------------
      ##Precipitation Header Rows
create.title.panes <- function(wb,var.name,start.row,sheet) {

      pane.titles <- list(' ','CWEC',
                              'Past',
                              '2020s Change',' ',
                              '2050s Change',' ',
                              '2080s Change',' ',
                              '2020s Percent Change',' ',
                              '2050s Percent Change',' ',
                              '2080s Percent Change',' ')

      pane.colour <- switch(var.name,pr='lightblue',tas='tan1')             
      hdstyle <- createStyle(fgFill = pane.colour, halign = "CENTER", textDecoration = "Bold",
                             border = "TopBottomLeftRight", fontColour = "black") 
      writeData(wb, sheet=sheet, pane.titles, startRow = start.row, startCol = 1, headerStyle = hdstyle,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=sheet,hdstyle,rows=start.row,cols=1:15,gridExpand=FALSE,stack=FALSE)
      mergeCells(wb,sheet=sheet,cols=c(4,5),rows=start.row)
      mergeCells(wb,sheet=sheet,cols=c(6,7),rows=start.row)
      mergeCells(wb,sheet=sheet,cols=c(8,9),rows=start.row)
      mergeCells(wb,sheet=sheet,cols=c(10,11),rows=start.row)
      mergeCells(wb,sheet=sheet,cols=c(12,13),rows=start.row)
      mergeCells(wb,sheet=sheet,cols=c(14,15),rows=start.row)
      prct.header <- list(' ',' ',' ','Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                                  'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%')          
      writeData(wb, sheet=sheet, prct.header, startRow = start.row+1, startCol = 1, headerStyle = hdstyle,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=sheet,hdstyle,rows=start.row+1,cols=1:15,gridExpand=FALSE,stack=FALSE)

      mergeCells(wb,sheet=sheet,cols=1:15,rows=start.row+2)
      pane.header <- list(switch(var.name,pr='Precipitation',tas='Temperature'))
      writeData(wb, sheet=sheet, pane.header, startRow = start.row+2, startCol = 1, headerStyle = hdstyle,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=sheet,hdstyle,rows=start.row+2,cols=1:15,gridExpand=FALSE,stack=FALSE)
}

##-----------------------------------------------------------------------------------------