##################################################
## R script for NetworkAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

Set.Config <-function(anal.mode="web"){
  globalConfig <- list();
  globalConfig$anal.mode <- anal.mode
  globalConfig <<- globalConfig;
}

# init resources for analysis

#'Initialize resources for analysis
#'@description call this function before performing any analysis
#'@param path path pointing to different built-in resources
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
Init.Data<-function(onWeb=T, path="../../"){
  node.infoU <<- data.frame()

  dataSet <- list(annotated=FALSE);
  analSet <- list(annotated=FALSE);
  paramSet <- list(annotated=FALSE);
  msgSet <- list(annotated=FALSE);
  cmdSet <- list(annotated=FALSE);
  
  if(onWeb){
  Set.Config("web");
  }else{
  Set.Config("local");
  }

  dataSet$params <- list();
  partialToBeSaved <<- c("Rload.RData", "Rhistory.R")
  Sys.setenv("OMP_NUM_THREADS" = 2); # need to control parallel computing for some packages
  init.lib <<- "kegg";
  regids <<- ""
  selectedFactorInx <<- 1;
  data.idType <<- "NA";
  chord_count <<- 0;
  netUploadU <<- 0;
  net.stats <<- as.data.frame(matrix(0, ncol = 3, nrow = 1));
  enr.mat <<- NULL;
  rankOptGlobal <<- "pval";
  data.org <<- "hsa";
  keggpw.count <<- 0;
  pvalu <<- 0.05;
  reg.count <<- 0; 


  dataSet$name <- "";

  if(file.exists("/home/glassfish/sqlite/")){
    sqlite.path <- "/home/glassfish/sqlite/";  #public server
  }else if(file.exists("/Users/xia/Dropbox/sqlite/")){
    sqlite.path <- "/Users/xia/Dropbox/sqlite/"; #xia local
  }else if(file.exists("/Users/jeffxia/Dropbox/sqlite/")){
    sqlite.path <- "/Users/jeffxia/Dropbox/sqlite/"; #xia local2
  }else if(file.exists("/media/zzggyy/disk/sqlite/")){  
    sqlite.path <-"/media/zzggyy/disk/sqlite/"; #zgy local
  }else if(file.exists("/home/le/sqlite/networkanalystdatabase/")){
    sqlite.path <- "/home/le/sqlite/networkanalystdatabase/"; #le local
  }else if(file.exists("/home/zgy/sqlite/")){
    sqlite.path <-"/home/zgy/sqlite/"; #zgy local
  }else if(file.exists("/Users/jessicaewald/sqlite/")){ # ewald local
    sqlite.path <- "/Users/jessicaewald/sqlite/";
  }else{
    sqlite.path <-"";
  }
  
  data.org <<- NULL;
  module.count <<- 0;
  msg.vec <<- vector(mode="character");
  current.msg <<- "";
  
  paramSet$partialToBeSaved <- c("Rload.RData", "Rhistory.R", "paramSet.qs", "msgSet.qs", "analSet.qs", "cmdSet.qs");

  Sys.setenv("OMP_NUM_THREADS" = 2); # need to control parallel computing for some packages
  paramSet$init.lib <- "kegg";
  paramSet$selectedFactorInx <- 1; #in multi comparison (i.e pairwise, time-series) which contrast is used
  analSet$net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = 1));
  msgSet$summaryVec <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "NA"); 
  analSet$enr.mat <- NULL;
  paramSet$numOfLists <- 1;
  paramSet$data.idType <- "";
  paramSet$pvalu <- 0.05;
  paramSet$selDataNm <- "meta_default";
  paramSet$mdata.all <- list();
  paramSet$anal.type <- "onedata";
  paramSet$api.bool <- F;
  paramSet$jsonNms <- list();

  paramSet$sqlite.path <- sqlite.path;
  paramSet$lib.path <- paste0(path, "data/");

  if(!onWeb) {
    paramSet$sqlite.path <- paste0(getwd(), "/");
    paramSet$lib.path <- "https://www.expressanalyst.ca/ExpressAnalyst/resources/data/";
    .on.public.web <<- T;
  } 

  paramSet$on.public.web <- onWeb;
  paramSet$data.org <- "hsa";
  paramSet$module.count <- 0;
  msgSet$current.msg <- vector(mode="character");
  msgSet$msg.list <- list(); #numbered list, each element: function name, line number, time stamp, severity

  # preload some general package
  require('Cairo');
  CairoFonts("Arial:style=Regular","Arial:style=Bold","Arial:style=Italic","Helvetica","Symbol")
  require('igraph');
  print("called networkanalyst init!");

  dataSets <<- list();
  saveSet(paramSet, "paramSet");
  saveSet(msgSet, "msgSet");
  saveSet(analSet, "analSet");
  saveSet(cmdSet, "cmdSet");

  SetAnalType("genelist");

  return(RegisterData(dataSet,1));
}

.set.mSet <- function(dataSetObj=NA){
  dataSet <<- dataSetObj;
  if(.on.public.web){
    return (1);
  }else{
    return(dataSetObj);
  }
}

.get.mSet <- function(dataSetObj=NA){
  if(.on.public.web){
    return(dataSet)
  }else{
    return(dataSetObj);
  }
}

# genelist, onedata, metadata
# also set up or clear the other global objects
SetAnalType <- function(analType){
  paramSet <- readSet(paramSet, "paramSet");
  paramSet$anal.type <- analType;
  paramSet$mdata.all <- list();
  paramSet$meta.selected <- TRUE;
  paramSet$meta.upload <- FALSE; # when upload merged data from meta-analysis b4
  if(analType == "metadata"){
    paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, "inmex_meta.qs")
  }
  saveSet(paramSet, "paramSet");
  return(paste0("Set to ",analType));
}

SetNetType <- function(netType){
  net.type <<- netType;
}

# When multiple genelists/datasets/results, record their name and save the data as .RDS file
# a) Current dataSet object
# Note, the memory will only contain one dataSet object. By default, the last one will be the current dataSet object;
# Users can switch this (from the interface) to specify which data is load into memory (dataSet object)
# b) Include for certain analysis
# For chord and heatmap analysis, users can do multiple selection (include)
# All datasets are selected by default (1 for selected, 0 for unselected)

# note, dataSet need to have "name" property
# note, dataSet need to have "name" property
RegisterData <- function(dataSet, output=1){
  dataName <- dataSet$name;

  paramSet <- readSet(paramSet, "paramSet");
  if(length(dataName)>0){
  mdata.all <- paramSet$mdata.all;
  mdata.all[[dataName]] <- 1;
  paramSet$mdata.all <- mdata.all;
  saveSet(paramSet, "paramSet");
  }

  if(globalConfig$anal.mode == "web"){
    dataSets[[dataName]] <- dataSet;
    dataSets <<- dataSets;
    return(output);
  }else{
    if( globalConfig$anal.mode == "api"){
        qs::qsave(dataSet, file=dataName);
        return(output);
    }else{
        dataSets[[dataName]] <- dataSet;
        dataSets <<- dataSets;
        return(dataSets);
    }
  }
} 

# only for switching single expression data results
SetCurrentData <- function(nm){
  if(dataSet$name != nm){
    dataSet <- qs::qread(nm);
  }
  return(.set.mSet(dataSet));
}

# remove data object, the current dataSet will be the last one by default
RemoveData <- function(dataName){
  paramSet <- readSet(paramSet, "paramSet");
  if(!is.null(mdata.all[[dataName]])){
    mdata.all[[dataName]] <- NULL;
  }
  paramSet$mdata.all <- mdata.all;
  saveSet(paramSet, "paramSet");

}

# users can select one or more data for analysis
# note, we use 1 to indicate this is selected
# and by default is all selected.
SelectDataSet <- function(){
  if(!exists('nm.vec')){
    msgSet <- readSet(msgSet, "msgSet");
    msgSet$current.msg <-"No dataset is selected for analysis!";
    saveSet(msgSet, "msgSet");
    return(0);
  }

  paramSet <- readSet(paramSet, "paramSet");
  mdata.all <- paramSet$mdata.all;

  all.nms <- names(mdata.all);
  for(nm in all.nms){
    if(nm %in% nm.vec){
      mdata.all[[nm]] <- 1;
    }else{
      mdata.all[[nm]] <- 0;
    }
  }
  paramSet$mdata.all <- mdata.all;
  
  rm('nm.vec', envir = .GlobalEnv);
  return(1);
}


GetAllDataNames <- function(){
  names(mdata.all);
}

SetOrganism <- function(org){
  paramSet <- readSet(paramSet, "paramSet");
  paramSet$data.org <- org;
  saveSet(paramSet, "paramSet");
}

SetSelectedFactorInx <- function(inx){
  selectedFactorInx <<- inx;
}

SetSelNetDataset <- function(type){
  paramSet <- readSet(paramSet, "paramSet");
  paramSet$selectedNetDataset <- type;
  saveSet(paramSet, "paramSet");
}

SetSelMultiNet <- function(type){
  selectedMultiNet <<- type;
}

SetRankingMetric <- function(opt){
  rankOptGlobal <<- opt;
}


SetListNms <- function(){
  newDat <- list();
  tot.count <- 0;
  listSizes <- list();
  
  # convert to entrez

  en.ids <- rownames(dataSet$resTable)
  nm <- "dataSet"
  
  names(en.ids) <- doEntrez2SymbolMapping(en.ids)
  
  listSizes[[1]] <- list(
    name = nm,
    label = nm,
    size = length(en.ids)
  );
  
  list.genes <<- en.ids;
  listSizes <<- listSizes;
}

GetDataListNames <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(names(paramSet$mdata.all));
}

GetCurrentJson <-function(type){
  return(paramSet$jsonNms[[type]]);
}

GetFilesToBeSaved <-function(naviString){
  return(unique(partialToBeSaved));
}

PrepareJsonFromR <- function(fileNm, type, jsonString, dataSetString){
    # rjson bug use RJSONIO
    dataSet <- RJSONIO::fromJSON(dataSetString);
    dataSet <<- dataSet
    sink(fileNm);
    cat(jsonString);
    sink();
    return(1)
}



#'Record R Commands
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param cmd Commands 
#'@export
RecordRCommand <- function(cmd){
  cmdSet <- readSet(cmdSet, "cmdSet"); 
  cmdSet$cmdVec <- c(cmdSet$cmdVec, cmd);
  saveSet(cmdSet, "cmdSet");
  return(1);
}

SaveRCommands <- function(){
  cmdSet <- readSet(cmdSet, "cmdSet"); 
  cmds <- paste(cmdSet$cmdVec, collapse="\n");
  pid.info <- paste0("# PID of current job: ", Sys.getpid());
  cmds <- c(pid.info, cmds);
  write(cmds, file = "Rhistory.R", append = FALSE);
}

#'Export R Command History
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@export
GetRCommandHistory <- function(){
  cmdSet <- readSet(cmdSet, "cmdSet"); 
  if(length(cmdSet$cmdVec) == 0){
    return("No commands found");
  }
  return(cmdSet$cmdVec);
}
