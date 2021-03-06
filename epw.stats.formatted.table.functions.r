##Script to sort out the data requested for Okanagan

source('/storage/home/ssobie/code/repos/assessments/summary.table.comments.r',chdir=T)

library(openxlsx)

##------------------------------------------------------------- 
##------------------------------------------------------------- 

get.units.pane <- function(var.name) {
  leg.label <- NA
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
  if (grepl("(tas|txx|tnn|txn|tnx|tmax|tmin|dtr|wetbulb)", var.name))
    leg.label <- '\u00B0C' ##rep('degC',13)
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
  derived.vars <- c('hdd','tnnETCCDI','tasmin.annual_quantile_010','tasmin.annual_quantile_025',
                    'wetbulb.annual_quantile_010','wetbulb.annual_quantile_025','cdd','cdd_10',
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
write_variables <- function(wb,sheet,sorted.vars,row.locs,type,sheet.data,data.cols,highlights) {
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
    print(title.cell)
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
    addStyle(wb,sheet=sheet,datastyle,rows=current.row,cols=data.cols,gridExpand=TRUE,stack=FALSE)
    addStyle(wb,sheet=sheet,hdstyle,rows=current.row,cols=1,gridExpand=TRUE,stack=FALSE)

    highlight <- createStyle(fgFill = 'lightyellow', halign = "CENTER", 
                             border = "TopBottomLeftRight", fontColour = "black")  
    ##addStyle(wb,sheet=sheet,highlight,rows=current.row,cols=2,gridExpand=FALSE,stack=FALSE)
    for (h in highlights) { ##Add yellow shade for change columns
       addStyle(wb,sheet=sheet,highlight,rows=current.row,cols=h,gridExpand=FALSE,stack=FALSE)
    }

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
      s1 <- createStyle(fontSize = 12, fontColour = "black")
      c1 <- createComment(comment = 'Variables calculated from the CWEC2016 TMY file',style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 2, row = 1, comment = c1)
      mergeCells(wb,sheet=sheet,cols=c(3,4),rows=1)
      com <- createComment(comment = 'Variables calculated from the morphed TMY file for the 2020s',style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 3, row = 1, comment = com)
      mergeCells(wb,sheet=sheet,cols=c(5,6),rows=1)
      com <- createComment(comment = 'Differences in the variables between 2020s morphed TMY file and the CWEC2016 TMY File',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 5, row = 1, comment = com)
      mergeCells(wb,sheet=sheet,cols=c(7,8),rows=1)
      com <- createComment(comment = 'Variables calculated from the morphed TMY file for the 2050s',style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 7, row = 1, comment = com)
      mergeCells(wb,sheet=sheet,cols=c(9,10),rows=1)
      com <- createComment(comment = 'Differences in the variables between 2050s morphed TMY file and the CWEC2016 TMY File',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 9, row = 1, comment = com)
      mergeCells(wb,sheet=sheet,cols=c(11,12),rows=1)
      com <- createComment(comment = 'Variables calculated from the morphed TMY file for the 2080s',style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 11, row = 1, comment = com)
      mergeCells(wb,sheet=sheet,cols=c(13,14),rows=1)
      com <- createComment(comment = 'Differences in the variables between 2080s morphed TMY file and the CWEC2016 TMY File',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 13, row = 1, comment = com)
      ##freezePane(wb,sheet=sheet,firstRow=TRUE)

      prct.header <- list(' ',' ','Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                                  'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%')          
      writeData(wb, sheet=sheet, prct.header, startRow = 2, startCol = 1, headerStyle = fz1,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=sheet,fz1,rows=2,cols=1:14,gridExpand=FALSE,stack=FALSE)

      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2020s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 3, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2020s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 4, row = 2, comment = com)

      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2020s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 5, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2020s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 6, row = 2, comment = com)

      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2050s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 7, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2050s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 8, row = 2, comment = com)
 
      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2050s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 9, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2050s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 10, row = 2, comment = com)

      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2080s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 11, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2080s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 12, row = 2, comment = com)

      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2080s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 13, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2080s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 14, row = 2, comment = com)




      ##freezePane(wb,sheet=sheet,firstActiveRow=3)
}

##-----------------------------------------------------------------------------------------
create_temp_climatologies_top_pane <- function(wb,sheet,site) {

##                              'Model Past',
      pane.titles <- list(site,'Past (TMY)',
                              'TMY Year',
                              'BCCAQv2 Past',
                              'BCCAQv2 Past',
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
      addStyle(wb,sheet=sheet,fz1,rows=1,cols=1:17,gridExpand=FALSE,stack=FALSE)
      s1 <- createStyle(fontSize = 12, fontColour = "black")
      c1 <- createComment(comment = 'Mean temperatures calculated from the CWEC2016 TMY file',style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 2, row = 1, comment = c1)
      com <- createComment(comment = 'Representative year selected for each month in the CWEC2016 TMY file',style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 3, row = 1, comment = c1)
      com <- createComment(comment = 'Mean temperatures calculated from the model scenario past',style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 4, row = 1, comment = com)
      writeComment(wb, sheet=sheet, col = 5, row = 1, comment = com)
      mergeCells(wb,sheet=sheet,cols=c(6,7),rows=1)
      com <- createComment(comment = 'Mean temperatures calculated from the morphed TMY file for the 2020s',style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 6, row = 1, comment = com)
      mergeCells(wb,sheet=sheet,cols=c(8,9),rows=1)
      com <- createComment(comment = 'Differences in mean temperatures between 2020s morphed TMY file and the CWEC2016 TMY File',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 8, row = 1, comment = com)
      mergeCells(wb,sheet=sheet,cols=c(10,11),rows=1)
      com <- createComment(comment = 'Mean temperatures calculated from the morphed TMY file for the 2050s',style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 10, row = 1, comment = com)
      mergeCells(wb,sheet=sheet,cols=c(12,13),rows=1)
      com <- createComment(comment = 'Differences in mean temperatures between 2050s morphed TMY file and the CWEC2016 TMY File',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 12, row = 1, comment = com)
      mergeCells(wb,sheet=sheet,cols=c(14,15),rows=1)
      com <- createComment(comment = 'Mean temperatures calculated from the morphed TMY file for the 2080s',style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 14, row = 1, comment = com)
      mergeCells(wb,sheet=sheet,cols=c(16,17),rows=1)
      com <- createComment(comment = 'Differences in mean temperatures between 2080s morphed TMY file and the CWEC2016 TMY File',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 16, row = 1, comment = com)
      ##freezePane(wb,sheet=sheet,firstRow=TRUE)

      ##'(1971-2000)',
      prct.header <- list(' ',' ',' ','(1971-2000)','(1998-2014) ',
                          'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%',
                          'Average','10th%-90th%','Average','10th%-90th%','Average','10th%-90th%')          
      writeData(wb, sheet=sheet, prct.header, startRow = 2, startCol = 1, headerStyle = fz1,
                borders = "rows", borderStyle = "medium")
      addStyle(wb,sheet=sheet,fz1,rows=2,cols=1:17,gridExpand=FALSE,stack=FALSE)
      ##freezePane(wb,sheet=sheet,firstActiveRow=3)

      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2020s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 6, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2020s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 7, row = 2, comment = com)

      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2020s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 8, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2020s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 9, row = 2, comment = com)

      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2050s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 10, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2050s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col =11, row = 2, comment = com)
 
      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2050s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 12, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2050s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 13, row = 2, comment = com)

      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2080s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 14, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2080s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 15, row = 2, comment = com)

      com <- createComment(comment = 'Average of the 10 model ensemble of future scenarios for the 2080s',
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 16, row = 2, comment = com)
      com <- createComment(comment = paste0('Range of the 10 model ensemble of future scenarios for the 2080s',
                                            ' (10th and 90th percentiles of the ensemble are shown)'),
                           style=s1,visible=FALSE)
      writeComment(wb, sheet=sheet, col = 17, row = 2, comment = com)

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