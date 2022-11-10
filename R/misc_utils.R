##################################################
## R scripts for NetworkAnalyst 
## Various utility methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# given a data with duplicates, dups is the one with duplicates
.removeDuplicates <- function(data, lvlOpt, quiet=T){
  
  all.nms <- rownames(data);
  colnms <- colnames(data);
  dup.inx <- duplicated(all.nms);
  dim.orig  <- dim(data);
  data <- apply(data, 2, as.numeric); # force to be all numeric
  dim(data) <- dim.orig; # keep dimension (will lost when only one item) 
  rownames(data) <- all.nms;
  colnames(data) <- colnms;
  if(sum(dup.inx) > 0){
    uniq.nms <- all.nms[!dup.inx];
    uniq.data <- data[!dup.inx,,drop=F];
    
    dup.nms <- all.nms[dup.inx];
    uniq.dupnms <- unique(dup.nms);
    uniq.duplen <- length(uniq.dupnms);
    
    for(i in 1:uniq.duplen){
      nm <- uniq.dupnms[i];
      hit.inx.all <- which(all.nms == nm);
      hit.inx.uniq <- which(uniq.nms == nm);
      
      # average the whole sub matrix 
      if(lvlOpt == "mean"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, mean, na.rm=T);
      }else if(lvlOpt == "median"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, median, na.rm=T);
      }else if(lvlOpt == "max"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, max, na.rm=T);
      }else{ # sum
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, sum, na.rm=T);
      }
    }
    if(!quiet){
      if(numOfLists == 1){
        current.msg <<- paste(current.msg, paste("A total of ", sum(dup.inx), " of duplicates were replaced by their ", lvlOpt, ".", sep=""), collapse="\n");
      }else{
        current.msg <<- paste(current.msg, paste0("<b>", listInxU, "</b> : ", length(data), " genes;"), collapse="\n");
      }
    }
    return(uniq.data);
  }else{
    if(!quiet){
      if(numOfLists == 1){
        current.msg <<- paste(current.msg, "All IDs are unique.", collapse="\n");
      }else{
        current.msg <<- paste(current.msg, paste0("<b>", listInxU, "</b> : ", length(data), " genes;"), collapse="\n");
      }
    }
    return(data);
  }
} 

# need to obtain the full path to convert (from imagemagik) for cropping images
GetBashFullPath<-function(){
  path <- system("which bash", intern=TRUE);
  if((length(path) == 0) && (typeof(path) == "character")){
    print("Could not find bash in the PATH!");
    return("NA");
  }
  return(path);
}


cleanMem <- function(n=8) { for (i in 1:n) gc() }

###########
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(format(utils::object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  print(lapply(dataSet, object.size));
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
ShowMemoryUse <- function(..., n=30) {
  library(pryr);
  sink(); # make sure print to screen
  print(mem_used());
  print(sessionInfo());
  print(.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n));
  print(warnings());
}

GetListEnrGeneNumber <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  mdata.all <- paramSet$mdata.all;
  anal.type <- paramSet$anal.type;
  all.enIDs <- NULL;
  listSizes <- list();
  if(anal.type == "genelist"){
    if(paramSet$numOfLists > 1){
      newDat <- list();
      tot.count <- 0;
      all.nms <- paramSet$listNms;
      for(i in 1:length(all.nms)){
        dataNm <- all.nms[i];
        dataSet <- readDataset(dataNm);
        gene.mat <- dataSet$prot.mat;
        
        # convert to entrez
        expr.val <- gene.mat[,1];
        en.ids <- rownames(gene.mat);
        
        names(expr.val) <- en.ids;
        newDat[[dataNm]] <- expr.val;
        names(en.ids) <- doEntrez2SymbolMapping(en.ids, paramSet$data.org, paramSet$data.idType)
        all.enIDs <- c(all.enIDs, en.ids);
        listSizes[[i]] <- list(
          name = dataNm,
          label = dataNm,
          size = length(en.ids)
          #val = de.prct[i]
        )
      }
      
    }else{
      
      all.enIDs <- rownames(dataSet$prot.mat);
      names(all.enIDs ) <- doEntrez2SymbolMapping(all.enIDs, paramSet$data.org, paramSet$data.idType)
      listSizes[[1]] <- list(
        name = "datalist1",
        label = "datalist1",
        size = length(all.enIDs)
        #val = de.prct[i]
      )
    }
  }else if(anal.type == "onedata"){
    all.enIDs <- rownames(dataSet$sig.mat);
    names(all.enIDs) <- doEntrez2SymbolMapping(all.enIDs, paramSet$data.org, paramSet$data.idType)
    listSizes[[1]] <- list(
      name = "dataSet1",
      label = "dataSet1",
      size = length(all.enIDs)
      #val = de.prct[i]
    )
  }else{
    newDat <- list();
    tot.count <- 0;
    listSizes <- list();
    all.nms <- names(mdata.all);
    for(i in 1:length(all.nms)){
      dataNm <- all.nms[i];
      dataSet <- readDataset(dataNm);
      gene.mat <- dataSet$sig.mat;
      
      # convert to entrez
      expr.val <- gene.mat[,1];
      en.ids <- rownames(gene.mat);
      
      names(expr.val) <- en.ids;
      newDat[[dataNm]] <- expr.val;
      names(en.ids) <- doEntrez2SymbolMapping(en.ids, paramSet$data.org)
      all.enIDs <- c(all.enIDs, en.ids);
      listSizes[[i]] <- list(
        name = dataNm,
        label = dataNm,
        size = length(en.ids)
      )
    }
  }
  analSet$list.genes <- all.enIDs;
  analSet$listSizes <- listSizes;
  saveSet(analSet, "analSet");
}

color_scale <- function(c1="grey", c2="red") {
  pal <- colorRampPalette(c(c1, c2))
  colors <- pal(100)
  return(colors)
}

# for project saving
PrepareSignatureOfNetworkAnalyst <- function(){
  
  if(anal.type == "genelist"){
    signature.gene <- dataSet$sig.mat;
    signature.gene.org <- data.org;
    save(signature.gene, signature.gene.org, file="RShare_networkanalyst.RData");  
  }else if(anal.type == "onedata"){
    if(!file.exists("express.res.t.qs")){
      return("-1");
    }
    resT <- qs::qread("express.res.t.qs");
    if(exists("P.Value", where=resT)){
      signature.gene <- as.matrix(resT$P.Value);
    }else if(exists("PValue", where=resT)){
      signature.gene <- as.matrix(resT$PValue);
    }
    rownames(signature.gene) <- rownames(resT);
    
    signature.gene.org <- data.org;
    save(signature.gene, signature.gene.org, file="RShare_networkanalyst.RData");  
  }
  
  return(.set.mSet(dataSet));
}

# in public web, this is done by microservice
.perform.computing <- function(){
  dat.in <- qs::qread("dat.in.qs"); 
  dat.in$my.res <- dat.in$my.fun();
  qs::qsave(dat.in, file="dat.in.qs");    
}

fast.write <- function(dat, file, row.names=TRUE){
    tryCatch(
        {
           if(is.data.frame(dat)){
                # there is a rare bug in data.table (R 3.6) which kill the R process in some cases 
                data.table::fwrite(dat, file, row.names=row.names);
           }else{
                write.csv(dat, file, row.names=row.names);  
           }
        }, error=function(e){
            print(e);
            fast.write.csv(dat, file, row.names=row.names);   
        }, warning=function(w){
            print(w);
            fast.write.csv(dat, file, row.names=row.names); 
        });
}

checkfac = function(fac) {

  if(is.numeric(fac)) {
    nrgrp = as.integer(max(fac, na.rm=TRUE)+1)
    fac   = as.integer(fac)
  }
  ## this must precede the factor test
  if(is.character(fac))
    fac = factor(fac)

  if (is.factor(fac)) {
    nrgrp = nlevels(fac)
    fac   = as.integer(as.integer(fac)-1)
  } 
  if(!is.integer(fac))
    stop("'fac' must be factor, character, numeric, or integer.")
  
  if(any(fac<0, na.rm=TRUE))
    stop("'fac' must not be negative.")
    
  return(list(fac=fac, nrgrp=nrgrp))
}

.getDynLoadPath <- function() {
    path = "../../rscripts/networkanalystr/src/NetworkAnalyst.so";
    return(path)
}

LoadRObjects <- function(path="", imgName, jsonName, fileNms){
    fileNms.vec <- unlist(strsplit(fileNms, ";"));
    for(i in 1:length(fileNms.vec)){
        link <- paste0(path, "/", fileNms.vec[i])
        download.file(link, fileNms.vec[i], quiet=T)
    }

    dataSet <- qs:::qread(imgName);
    data.org <<- dataSet$org;
    listSizes <<- dataSet$listSizes;
    ppi.comps <<- dataSet$ppi.comps;
    current.net.nm <<- dataSet$current.net.nm;
    if(dataSet$anal.type == "onedata"){
       rownames(dataSet$resTable) <- dataSet$resTableRowNames;
    }else if(dataSet$anal.type == "genelist"){
       dataSet$all.prot.mat[,1] <- as.numeric(dataSet$all.prot.mat[,1]);
    }else{
       meta.mat.all <<- dataSet$meta.mat.all
       performedDE <<- T;
    }
    anal.type <<- anal.type;
    all.prot.mat <<- dataSet$all.prot.mat;
    .set.mSet(dataSet);
}

ReadList <- function(dataSetObj=NA, fullPath, fileNm){
    fullUrl <- url(paste0(fullPath,"/", fileNm))
    all_str <- paste0(readLines(fullUrl),collapse="\n");
    return(all_str);
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


  GetGeneList <- function(dataSetObj=NA,fileNm, type){
    paramSet <- readSet(paramSet, "paramSet");
    mdata.all <- paramSet$mdata.all;
    all_str <- "";
    if(type == "genelist"){
      if(paramSet$numOfLists > 1){
        dataSet <- list();
        my.vec <- names(mdata.all);
        for(i in 1:length(my.vec)){
          datSet <- readDataset(my.vec[i]);
          if(i == 1){
            all_str = datSet$orig
          }else{
            all_str = paste0(all_str, "\n//", datSet$orig)
          }
        }
      }else{
        all_str = dataSet$orig;
      }
    }else{
      require(readr);
      my.vec <- names(mdata.all);
      for(i in 1:length(my.vec)){
        dataSet <- readDataset(my.vec[i]);
        sig.ids <- rownames(dataSet$sig.mat);
        stat.fc <- dataSet$sig.mat[,1];
        df <- data.frame(ids=sig.ids, fc=stat.fc);
        df_str <- readr:::format_tsv(df)
        df_str <- paste0("#", df_str);
        if(i == 1){
          all_str <- df_str;
        }else{
          all_str <- paste0(all_str, "\n//", df_str)
        }
      }
    }
    writeLines(all_str, fileNm)
    return(all_str);
  }


# given a data with duplicates, dups is the one with duplicates
RemoveDuplicates <- function(data, lvlOpt, quiet=T, paramSet, msgSet, listInx=1){
  paramSet <- readSet(paramSet, "paramSet");
  msgSet <- readSet(msgSet, "msgSet");

  all.nms <- rownames(data);
  colnms <- colnames(data);
  dup.inx <- duplicated(all.nms);
  dim.orig  <- dim(data);
  data <- apply(data, 2, as.numeric); # force to be all numeric
  dim(data) <- dim.orig; # keep dimension (will lost when only one item) 
  rownames(data) <- all.nms;
  colnames(data) <- colnms;
  if(sum(dup.inx) > 0){
    uniq.nms <- all.nms[!dup.inx];
    uniq.data <- data[!dup.inx,,drop=F];
    
    dup.nms <- all.nms[dup.inx];
    uniq.dupnms <- unique(dup.nms);
    uniq.duplen <- length(uniq.dupnms);
    
    for(i in 1:uniq.duplen){
      nm <- uniq.dupnms[i];
      hit.inx.all <- which(all.nms == nm);
      hit.inx.uniq <- which(uniq.nms == nm);
      
      # average the whole sub matrix 
      if(lvlOpt == "mean"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, mean, na.rm=T);
      }else if(lvlOpt == "median"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, median, na.rm=T);
      }else if(lvlOpt == "max"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, max, na.rm=T);
      }else{ # sum
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, sum, na.rm=T);
      }
    }
    if(!quiet){
      if(paramSet$numOfLists == 1){
        msgSet$current.msg <- paste(msgSet$current.msg, paste("A total of ", sum(dup.inx), " of duplicates were replaced by their ", lvlOpt, ".", sep=""), collapse="\n");
      }else{
        msgSet$current.msg <- paste(msgSet$current.msg, paste0("<b>", listInx, "</b> : ", length(data), " genes;"), collapse="\n");
      }
    }
    saveSet(msgSet, "msgSet");
    return(list(uniq.data,msgSet));
  }else{
    if(!quiet){
      if(paramSet$numOfLists == 1){
        msgSet$current.msg <- paste(msgSet$current.msg, "All IDs are unique.", collapse="\n");
      }else{
        msgSet$current.msg <- paste(msgSet$current.msg, paste0("<b>", listInx, "</b> : ", length(data), " genes;"), collapse="\n");
      }
    }
    saveSet(msgSet, "msgSet");
    return(list(data,msgSet));
  }
} 

###Gene list
GetNumOfLists <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(paramSet$numOfLists)
}
