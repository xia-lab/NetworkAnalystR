##################################################
## R script for NetworkAnalyst
## Description: functions only for list data analysis
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

GetNumOfLists <- function(){
  return(numOfLists)
}

# parse a list file
ReadListFile <- function(fileName) {
  dat1 <- data.table::fread(fileName, header=FALSE, check.names=FALSE, data.table=FALSE);
  dataSet$name <- fileName
  rowNms <- dat1[,1]
  if(length(dat1) == 1){
    dat1[,1] <- 0
  }else{
    dat1[,1] <- dat1[,2]
    dat1 <- dat1[,-2];
  }
  dataSet$prot.mat <- as.matrix(dat1)
  rownames(dataSet$prot.mat) <- rowNms;
  qs::qsave(dataSet, file=fileName); # keep original copy, not in mem
  return(1)
}
