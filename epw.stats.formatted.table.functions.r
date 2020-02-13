##Script to sort out the data requested for Okanagan

source('/storage/home/ssobie/code/repos/assessments/summary.table.comments.r',chdir=T)

library(openxlsx)

##------------------------------------------------------------- 
##------------------------------------------------------------- 

get.units.pane <- function(var.name) {
  leg.label <- NA
  if (grepl("(tas|txx|tnn|txn|tnx|tmax|tmin|dtr|wetbulb)", var.name))
    leg.label <- '\u00B0C' ##rep('degC',13)
  if (grepl("(pr|rx|r10|r20|r9|RP|sdii|prcptot)", var.name))
    leg.label <- 'mm' ##rep('mm',13)
  if (grepl("(dd)", var.name))
    leg.label <- 'Degree Days' ##rep('degree days',13)
  if (grepl("(fdE|cddE|cdd90|cddmax|cwd|su|gsl|id|trE|su30|r95daysE|r99daysE)", var.name))
    leg.label <- 'Days' ##rep('days',13)
  if (grepl("(dtr)", var.name))
    leg.label <- '\u00B0C' ##rep('degC',13)
  if (grepl("(r95sep|r95dist|r95days|r99days)", var.name))
    leg.label <- 'Days' ##rep('days',13)
  return(leg.label)
} 

##------------------------------------------------------------- 
get.round.val <- function(var.name) {
  rd <- 0
  if (grepl("(dd)", var.name))
    rd <- 0    
  if (grepl("(tas|txx|tnn|tnx|txn|tmax|tmin|wetbulb)", var.name))
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
  derived.vars <- c('hdd','tnnETCCDI','tasmin.annual_quantile_010','wetbulb.annual_quantile_010','cdd',
                    'txxETCCDI','tasmax.annual_quantile_975','wetbulb.annual_quantile_975')
  derived.ix <- match(derived.vars,all.vars)
  derived.selected <- table.vars[derived.ix[!is.na(derived.ix)]]

  tas.vars <- c(paste0('tas_',tolower(month.abb)),
                       'tas_win','tas_spr','tas_sum','tas_fal',
                       'tas_ann')
  tas.ix <- match(tas.vars,all.vars)
  tas.selected <- table.vars[tas.ix[!is.na(tas.ix)]]
  rv <- list(derived=derived.selected,tas=tas.selected)
  return(rv)
}

##------------------------------------------------------------- 

find.row.locations <- function(sorted.vars) {
 
   tas.vars <- sorted.vars$tas
   tas.rows <- vector(length=length(tas.vars)+1,mode='list')

   tas.rows[[1]] <- 1:2 ##Header Rows
   for (i in 1:length(tas.vars)) {
      row.len <- 1 ##One for the single row ##switch(tas.vars[[i]][[2]],annual=2,seasonal=6)
      tas.rows[[i+1]] <- seq(tail(tas.rows[[i]],1)+1,length.out=row.len)
   }

   derived.vars <- sorted.vars$derived
   derived.rows <- vector(length=length(derived.vars)+1,mode='list')

   derived.rows[[1]] <- 1:2 ##Header Rows
   for (i in 1:length(derived.vars)) {
      row.len <- 1 ##One for the single row ##switch(derived.vars[[i]][[2]],annual=2,seasonal=6)
      derived.rows[[i+1]] <- seq(tail(derived.rows[[i]],1)+1,length.out=row.len)
   }
   
   rv <- list(derived=derived.rows,tas=tas.rows)
   return(rv)

}
##---------------------------------------------

seasonal.table <- function(var.name,scenario='rcp85',rp=NULL) {

  no.percent <- '(wb|tas|tasmax|tasmin|txxETCCDI|tnnETCCD|tnxETCCDI|txnETCCDI|trETCCDI|suETCCDI|su30ETCCDI)'
  result <- get.seasonal.data(var.name,scenario)
  if (grepl(no.percent,var.name))
    result[,19:14] <- 'NA'
  return(as.data.frame(result))
}


annual.table <- function(var.name,scenario='rcp85',rp=NULL) {
  
  no.percent <- '(wb|tas|tasmax|tasmin|txxETCCDI|tnnETCCD|trETCCDI|r95sep|r99days|r95days)'
  result <- get.annual.data(var.name,scenario,rp) 

  if (grepl(no.percent,var.name))
    result[9:14] <- 'NA'
  return(as.data.frame(t(result)))
}


##----------------------------------------------------------------------------------------
write.variables <- function(wb,sheet,sorted.vars,row.locs,type,sheet.data) {
  len <- length(sorted.vars)
  for (i in 1:len) {               
    var.info <- sorted.vars[[i]]
    var.name <- var.info$name
    print(var.name)
    season <- grepl('season',var.info$type)
    var.title <- var.info$title
    print(var.title)
    current.row <- row.locs[[i+1]]
    title.cell <- paste0(var.title,' (',get.units.pane(var.name),')')                             
    var.entry <- c(title.cell,sheet.data[[var.name]])

    s1 <- createStyle(fontSize = 12, fontColour = "black", textDecoration = c("BOLD"))
    s2 <- createStyle(fontSize = 12, fontColour = "black")
    c1 <- createComment(comment = variable.comment(var.name),style=c(s1,s2),visible=FALSE)
    writeComment(wb, sheet, col = 1, row = current.row[1], comment = c1)

    hdstyle <- createStyle(fgFill = 'white', halign = "CENTER", textDecoration = "Bold",
                             border = "TopBottomLeftRight", fontColour = "black")
    datastyle <- createStyle(fgFill = 'white', halign = "RIGHT",
                             border = "TopBottomLeftRight", fontColour = "black")                              
    writeData(wb, sheet=sheet, as.data.frame(t(var.entry)), startRow = current.row, startCol =1, headerStyle = datastyle,
              borders = "rows", borderStyle = "medium",colNames=FALSE)
    addStyle(wb,sheet=sheet,datastyle,rows=current.row,cols=2:14,gridExpand=TRUE,stack=FALSE)
    addStyle(wb,sheet=sheet,hdstyle,rows=current.row,cols=1,gridExpand=TRUE,stack=FALSE)

    highlight <- createStyle(fgFill = 'lightyellow', halign = "CENTER", 
                             border = "TopBottomLeftRight", fontColour = "black")  
    ##addStyle(wb,sheet=sheet,highlight,rows=current.row,cols=2,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row,cols=5,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row,cols=6,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row,cols=9,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row,cols=10,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row,cols=13,gridExpand=FALSE,stack=FALSE)
    addStyle(wb,sheet=sheet,highlight,rows=current.row,cols=14,gridExpand=FALSE,stack=FALSE)

  }
     
}
##----------------------------------------------------------------------------------------
      ##Top Frozen Pane

create.frozen.top.pane <- function(wb,sheet,site) {

      pane.titles <- list(site,'Past (TMY)',
                              '2020s Future',' ',
                              '2020s Change',' ',
                              '2050s Future',' ',
                              '2050s Change',' ',
                              '2080s Future',' ',
                              '2080s Change',' ')

      fz1 <- createStyle(fgFill = "gray94", halign = "CENTER", textDecoration = "Bold",
                         border = "TopBottomLeftRight", fontColour = "black")
      writeData(wb, sheet=sheet, pane.titles, startRow = 1, startCol = 1, headerStyle = fz1,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=sheet,fz1,rows=1,cols=1:14,gridExpand=FALSE,stack=FALSE)
      mergeCells(wb,sheet=sheet,cols=c(3,4),rows=1)
      mergeCells(wb,sheet=sheet,cols=c(5,6),rows=1)
      mergeCells(wb,sheet=sheet,cols=c(7,8),rows=1)
      mergeCells(wb,sheet=sheet,cols=c(9,10),rows=1)
      mergeCells(wb,sheet=sheet,cols=c(11,12),rows=1)
      mergeCells(wb,sheet=sheet,cols=c(13,14),rows=1)
      ##freezePane(wb,sheet=sheet,firstRow=TRUE)

      prct.header <- list(' ',' ','Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                                  'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%')          
      writeData(wb, sheet=sheet, prct.header, startRow = 2, startCol = 1, headerStyle = fz1,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=sheet,fz1,rows=2,cols=1:14,gridExpand=FALSE,stack=FALSE)
      ##freezePane(wb,sheet=sheet,firstActiveRow=3)
}

##-----------------------------------------------------------------------------------------
##Header Rows
create.title.panes <- function(wb,var.name,start.row,sheet) {

      pane.titles <- list(' ','Past',
                              '2020s Future',' ',
                              '2020s Change',' ',
                              '2050s Future',' ',
                              '2050s Change',' ',
                              '2080s Future',' ',
                              '2080s Change',' ')

      pane.colour <- switch(var.name,pr='lightblue',tas='tan1')             
      hdstyle <- createStyle(fgFill = pane.colour, halign = "CENTER", textDecoration = "Bold",
                             border = "TopBottomLeftRight", fontColour = "black") 
      writeData(wb, sheet=sheet, pane.titles, startRow = start.row, startCol = 1, headerStyle = hdstyle,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=sheet,hdstyle,rows=start.row,cols=1:14,gridExpand=FALSE,stack=FALSE)
      mergeCells(wb,sheet=sheet,cols=c(3,4),rows=start.row)
      mergeCells(wb,sheet=sheet,cols=c(5,6),rows=start.row)
      mergeCells(wb,sheet=sheet,cols=c(7,8),rows=start.row)
      mergeCells(wb,sheet=sheet,cols=c(9,10),rows=start.row)
      mergeCells(wb,sheet=sheet,cols=c(11,12),rows=start.row)
      mergeCells(wb,sheet=sheet,cols=c(13,14),rows=start.row)
      prct.header <- list(' ',' ','Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                                  'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%')          
      writeData(wb, sheet=sheet, prct.header, startRow = start.row+1, startCol = 1, headerStyle = hdstyle,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=sheet,hdstyle,rows=start.row+1,cols=1:14,gridExpand=FALSE,stack=FALSE)

      mergeCells(wb,sheet=sheet,cols=1:14,rows=start.row+2)
      ##pane.header <- list(switch(var.name,pr='Precipitation',tas='Temperature'))
      ##writeData(wb, sheet=sheet, pane.header, startRow = start.row+2, startCol = 1, headerStyle = hdstyle,
      ##          borders = "rows", borderStyle = "medium")
      ##addStyle(wb,sheet=sheet,hdstyle,rows=start.row+2,cols=1:14,gridExpand=FALSE,stack=FALSE)
}

##-----------------------------------------------------------------------------------------