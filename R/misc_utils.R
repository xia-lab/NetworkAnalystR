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

PrepareEnrichNet<-function(netNm, type, overlapType){
    if(!exists("my.prepareEnrichNet")){ 
        compiler::loadcmp("../../rscripts/networkanalystr/_utils_prepareEnrichNet.Rc");    
    }
    return(my.prepareEnrichNet(netNm, type, overlapType));
}

GetListEnrGeneNumber <- function(){
  all.enIDs <- NULL;
  listSizes <- list();
  if(anal.type == "genelist"){
    if(numOfLists > 1){
      newDat <- list();
      tot.count <- 0;
      all.nms <- listNms;
      for(i in 1:length(all.nms)){
        dataNm <- all.nms[i];
        dataSet <- qs::qread(dataNm);
        gene.mat <- dataSet$prot.mat;
        
        # convert to entrez
        expr.val <- gene.mat[,1];
        en.ids <- rownames(gene.mat);
        
        names(expr.val) <- en.ids;
        newDat[[dataNm]] <- expr.val;
        names(en.ids) <- doEntrez2SymbolMapping(en.ids)
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
      names(all.enIDs ) <- doEntrez2SymbolMapping(all.enIDs)
      listSizes[[1]] <- list(
        name = "datalist1",
        label = "datalist1",
        size = length(all.enIDs)
        #val = de.prct[i]
      )
    }
  }else if(anal.type == "onedata"){
    all.enIDs <- rownames(dataSet$sig.mat);
    names(all.enIDs) <- doEntrez2SymbolMapping(all.enIDs)
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
      dataSet <- qs::qread(dataNm);
      gene.mat <- dataSet$sig.mat;
      
      # convert to entrez
      expr.val <- gene.mat[,1];
      en.ids <- rownames(gene.mat);
      
      names(expr.val) <- en.ids;
      newDat[[dataNm]] <- expr.val;
      names(en.ids) <- doEntrez2SymbolMapping(en.ids)
      all.enIDs <- c(all.enIDs, en.ids);
      listSizes[[i]] <- list(
        name = dataNm,
        label = dataNm,
        size = length(en.ids)
      )
    }
  }
  list.genes <<- all.enIDs;
  listSizes <<- listSizes;
}

InitListEnrichment <- function(dataSet=NA, type){
  dataSet <- .get.mSet();
  GetListEnrGeneNumber();
  res <- .performEnrichAnalysis(dataSet, paste0("enrichment_", type), type, list.genes);
  if(res){
    res <- PrepareEnrichNet(paste0('enrichNet_', type), 'list', "mixed");
  }
  return(res)
}

PerformListEnrichmentView <- function(dataSetObj=NA, file.nm, fun.type, netNm, IDs){
  dataSet <- .get.mSet(dataSetObj);
  gene.vec <- unlist(strsplit(IDs, "; "));
  gene.vec <- unique(gene.vec);
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  names(gene.vec) <- sym.vec;
  list.genes <<- gene.vec
  res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, list.genes);
  if(res){
    res <- PrepareEnrichNet(netNm, "list", "mixed");
  }
  if(res == 0){
    return(0);
  }else{
    return(.set.mSet(dataSet));
  }
}

overlap_ratio <- function(x, y, type) {
  x <- unlist(x)
  y <- unlist(y)
  if(type == "mixed"){
    res <- 0.5 * length(intersect(x, y))/length(unique(y)) + 0.5 * length(intersect(x, y))/length(unique(c(x,y)))
  }else if(type == "overlap"){
    if(length(x)>length(y)){
      res=length(intersect(x, y))/length(unique(y))
    }else{
      res=length(intersect(x, y))/length(unique(x))
    }
  }else{
    res=length(intersect(x, y))/length(unique(c(x,y)))
  }
  return(res)
}

color_scale <- function(c1="grey", c2="red") {
  pal <- colorRampPalette(c(c1, c2))
  colors <- pal(100)
  return(colors)
}


CalculateDEgeneSetEnr <- function(nms, operation, refNm, filenm){
  nms <- strsplit(nms, ";")[[1]];
  if(anal.type == "metadata" || anal.type == "onedata"){
    com.smbls <- PerformSetOperation_DataEnr(nms, operation, refNm);
  }else{
    com.smbls <- PerformSetOperation_ListEnr(nms, operation, refNm);
  }
  
  sink(filenm);
  cat(rjson::toJSON(com.smbls));
  sink();
}

PerformSetOperation_ListEnr <- function(nms, operation, refNm){
  all.nms <- names(mdata.all);
  include.inx <- all.nms %in% nms;
  my.vec <- all.nms[include.inx];
  if(anal.type == "onedata"){
    my.vec <- c("1");
  }
  if(operation == "diff"){ # make sure reference is the first
    inx <- which(my.vec == refNm);
    my.vec <- my.vec[-inx];
  }
  com.ids <- NULL;
  ids.list <- list()
  for(i in 1:length(my.vec)){
    if(anal.type != "onedata"){
      dataSet <- qs::qread(my.vec[i]);
    }
    if(operation == "diff"){
      ids.list[[i]]=dataSet$GeneAnotDB[,"gene_id"]
      #com.ids <- setdiff(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
    }else if(is.null(com.ids)){
      com.ids <- dataSet$GeneAnotDB[,"gene_id"];
    }else{
      if(operation == "intersect"){
        com.ids <- intersect(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
      }else if(operation == "union"){
        com.ids <- union(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
      }
    }
  }
  if(operation == "diff"){
    dataSet <- qs::qread(refNm);
    ids <- unique(unlist(ids.list));
    com.ids <-setdiff(dataSet$GeneAnotDB[,"gene_id"], ids);
  }
  
  com.ids <- unique(as.character(com.ids[!is.na(com.ids)])); # make sure it is interpreted as name not index
  com.symbols <- doEntrez2SymbolMapping(com.ids);
  names(com.symbols) <- com.ids
  
  com.symbols<-com.symbols[!is.null(com.symbols)];
  venn.genes <<- com.ids;
  return(com.symbols);
}

PerformSetOperation_DataEnr <- function(nms, operation, refNm){
  
  my.vec <- nms
  if(operation == "diff"){ # make sure reference is the first
    inx <- which(my.vec == refNm);
    my.vec <- my.vec[-inx];
  }
  com.ids <- NULL;
  ids.list <- list()
  if(anal.type == "onedata"){
    my.vec <- "dat"
  }
  for(nm in my.vec){
    if(anal.type != "onedata"){
      dataSet <- qs::qread(nm);
    }
    if(operation == "diff"){
      ids.list[[nm]]=rownames(dataSet$sig.mat);
      #com.ids <- setdiff(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
    }else if(is.null(com.ids)){
      com.ids <- rownames(dataSet$sig.mat);
    }else{
      if(operation == "intersect"){
        com.ids <- intersect(com.ids, rownames(dataSet$sig.mat));
      }else if(operation=="union"){
        com.ids <- union(com.ids, rownames(dataSet$sig.mat));
      }
    }
  }
  if(operation == "diff"){
    dataSet <- qs::qread(refNm);
    ids <- unique(unlist(ids.list));
    com.ids <-setdiff(rownames(dataSet$sig.mat), ids);
  } 
  com.symbols <- doEntrez2SymbolMapping(com.ids);
  names(com.symbols) <- com.ids;
  venn.genes <<- com.ids;
  return(com.symbols);
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

rowcolFt =  function(x, fac, var.equal, which = 1L) {
  
  if(!(which %in% c(1L, 2L)))
    stop(sQuote("which"), " must be 1L or 2L.")
  
  if(which==2L)
    x = t(x)

  if (typeof(x) == "integer")
      x[] <- as.numeric(x)

  sqr = function(x) x*x
  
  stopifnot(length(fac)==ncol(x), is.factor(fac), is.matrix(x))
  x   <- x[,!is.na(fac), drop=FALSE]
  fac <- fac[!is.na(fac)]

  ## Number of levels (groups)
  k <- nlevels(fac)

  ## xm: a nrow(x) x nlevels(fac) matrix with the means of each factor
  ## level
  xm <- matrix(
     sapply(levels(fac), function(fl) rowMeans(x[,which(fac==fl), drop=FALSE])),
     nrow = nrow(x),
     ncol = nlevels(fac))

  ## x1: a matrix of group means, with as many rows as x, columns correspond to groups 
  x1 <- xm[,fac, drop=FALSE]

  ## degree of freedom 1
  dff    <- k - 1

  if(var.equal){
    ## x0: a matrix of same size as x with overall means
    x0 <- matrix(rowMeans(x), ncol=ncol(x), nrow=nrow(x))
  
    ## degree of freedom 2
    dfr    <- ncol(x) - dff - 1

    ## mean sum of squares
    mssf   <- rowSums(sqr(x1 - x0)) / dff
    mssr   <- rowSums(sqr( x - x1)) / dfr

    ## F statistic
    fstat  <- mssf/mssr

  } else{

    ## a nrow(x) x nlevels(fac) matrix with the group size  of each factor
    ## level
    ni <- t(matrix(tapply(fac,fac,length),ncol=nrow(x),nrow=k))

    ## wi: a nrow(x) x nlevels(fac) matrix with the variance * group size of each factor
    ## level
    sss <- sqr(x-x1)
    x5 <- matrix(
       sapply(levels(fac), function(fl) rowSums(sss[,which(fac==fl), drop=FALSE])),
       nrow = nrow(sss),
       ncol = nlevels(fac))          
    wi <- ni*(ni-1) /x5

    ## u : Sum of wi
    u  <- rowSums(wi)

    ## F statistic
    MR <- rowSums(sqr((1 - wi/u)) * 1/(ni-1))*1/(sqr(k)-1)
    fsno <- 1/dff * rowSums(sqr(xm - rowSums(wi*xm)/u) * wi)
    fsdeno <- 1+ 2* (k-2)*MR
    fstat <- fsno/fsdeno

    ## degree of freedom 2: Vector with length nrow(x)
    dfr <- 1/(3 * MR)
  
  }
  
  res = data.frame(statistic = fstat,
                   p.value   = pf(fstat, dff, dfr, lower.tail=FALSE),
                   row.names = rownames(x))

  attr(res, "df") = c(dff=dff, dfr=dfr)
  return(res)
}

rowcoltt =  function(x, fac, tstatOnly, which, na.rm) {
    
  if(.on.public.web){
    dyn.load(.getDynLoadPath());
  }

  if (!missing(tstatOnly) && (!is.logical(tstatOnly) || is.na(tstatOnly)))
      stop(sQuote("tstatOnly"), " must be TRUE or FALSE.")
  
  f = checkfac(fac)
  if ((f$nrgrp > 2) || (f$nrgrp <= 0))
    stop("Number of groups is ", f$nrgrp, ", but must be >0 and <=2 for 'rowttests'.")

  if (typeof(x) == "integer")
      x[] <- as.numeric(x)

  cc = .Call("rowcolttests", x, f$fac, f$nrgrp, which-1L, na.rm)
    
  res = data.frame(statistic = cc$statistic,
                   dm        = cc$dm,
                   row.names = dimnames(x)[[which]])

  if (!tstatOnly)
    res = cbind(res, p.value = 2*pt(abs(res$statistic), cc$df, lower.tail=FALSE))

  attr(res, "df") = cc$df    
  return(res)
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

GetQEA.pathNames<-function(){
  current.geneset <- qs::qread("current_geneset.qs")
  hit.inx <- match(rownames(analSet$qea.mat),names(current.geneset));
  return(names(current.geneset)[hit.inx]);
}

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


  GetGeneList <- function(dataSetObj=NA,fileNm, type){
    all_str <- "";
      if(numOfLists > 1){
        dataSet <- list();
        my.vec <- names(mdata.all);
        for(i in 1:length(my.vec)){
          datSet <- qs::qread(my.vec[i]);
          if(i == 1){
            all_str = datSet$orig
          }else{
            all_str = paste0(all_str, "\n//", datSet$orig)
          }
        }
      }else{
        all_str = dataSet$orig;
      }
    writeLines(all_str, fileNm)
    return(all_str);
  }

GetDatasetNamesString <- function(){
    inmex.meta <- qs::qread("inmex_meta.qs");
    paste(unique(inmex.meta$data.lbl), collapse="||");
}

CalculateGsNet <- function(name, netNm, type, mType, db){
  res <- PerformMetaGseaNet(name, netNm, "pval", db, mType,0.05)
  return(.set.mSet(dataSet))
}

#For internal GSEA of meta-analysis of gene sets

PerformMetaGseaNet<- function(name, netNm, method, lib, mType, BHth=0.05){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  library(fgsea)
  curr.geneset <- .loadEnrichLib(lib);
  inmex.method <<- "effectsize";
  meta.stat <- "null";
  
  allMeta.mat <- qs::qread("allMeta.mat.qs");
  #allMeta.mat[allMeta.mat[,2]==0] <- 1e-20
  rankedVec <- as.vector(allMeta.mat[,1])*sign(allMeta.mat[,1]);
  
  names(rankedVec) <- rownames(allMeta.mat);
  rankedVec <- sort(rankedVec)
  rankedVec <- rankedVec[unique(names(rankedVec))]
  fgseaRes <- fgsea(pathways = curr.geneset, 
                    stats = rankedVec,
                    minSize=1,
                    maxSize=1000,
                    nperm=5000)
  fgseaRes <- fgseaRes[!duplicated(fgseaRes$pathway),]
  fgseaRes <- fgseaRes[,c("size","ES", "NES","padj", "pathway", "pval")]
  fgseaRes <- fgseaRes[order(-abs(fgseaRes$ES)),]
  fgseaRes <- fgseaRes[order(fgseaRes$pval),] 
  
  fgseaRe <- data.frame(fgseaRes)
  rownames(fgseaRes) <- fgseaRes$pathway
  es.mat <- as.matrix(fgseaRes[,c("ES","padj")]);
  colnames(es.mat) <- c("EnrichmentScore","Pvalue")
  rownames(es.mat) <- fgseaRes$pathway
  sig.inx <- which(es.mat[, "Pvalue"]<=BHth);
  if(length(sig.inx)<10){
    sig.inx <- c(1:10)
  }
  metaset.mat <<- es.mat[sig.inx,];
  metaset.mat.all <<- fgseaRes[sig.inx,]
  ii <- SetupMetaGSEAStats(name, netNm, BHth, mType, curr.geneset,lib);
  if(ii == 1){
    return(length(sig.inx));
  }else{
    return(0);
  }
}


SetupMetaGSEAStats <- function(name, netNm, BHth, mType, curr.geneset, lib){
  
  allmat <- qs::qread("allMeta.mat.qs");
  allmat.vec <- rownames(allmat);
  meta.mat <- metaset.mat
  metade.genes <- rownames(meta.mat);
  pval.mat <- fc.mat <- matrix(nrow=nrow(meta.mat), ncol=sum(mdata.all==1));
  fc.list <- split(rep(" ", length(allmat.vec)), allmat.vec);
  
  
  current.geneset <- curr.geneset[!duplicated(names(curr.geneset))]
  inx <- names(current.geneset) %in% rownames(meta.mat)  ;
  
  resTable <- meta.mat
  current.mset <- current.geneset[inx];
  
  inmex.meta <- qs::qread("inmex_meta.qs");

  ora.vec <- rownames(inmex.meta$data)
  ora.nms <- doEntrez2SymbolMapping(rownames(inmex.meta$data))
  
  hits.query <- lapply(current.mset, 
                       function(x) {
                         unique(ora.nms[ora.vec%in%unlist(x)]);
                       });  
  qs::qsave(hits.query, "hits_query.qs");
  
  set.num <- unlist(lapply(current.mset, function(x){length(unique(x))}), use.names=TRUE);
  names(hits.query) <- names(current.mset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=TRUE);
  
  vote.bool <- "false"
  meta.matcolinx <- 2;
  enr.score <- "NA"   
  
  padj <- p.adjust(as.vector(meta.mat[,meta.matcolinx]),method="BH");
  if(mType == "network"){
    json.res <- list(
      fun.anot = hits.query,
      fun.ids = as.vector(rownames(meta.mat)),
      fun.pval = as.vector(meta.mat[,meta.matcolinx]),
      fun.padj = padj,
      hit.num = hit.num,
      total= set.num
    );
  }else{
    json.res <- list(
      hits = hit.num,
      total= set.num,
      enr.pval= as.vector(meta.mat[,meta.matcolinx]),
      enr.padj= padj,
      enr.names= as.vector(rownames(meta.mat)),
      cls.lbl=inmex.meta$cls.lbl,
      smps.lbl=smps.vec,
      data.lbl = inmex.meta$data.lbl,
      path.lbl = rownames(meta.mat),
      enr.score = as.vector(meta.mat[,1]),
      isVote = vote.bool
    );
  }
  
  json.mat <- rjson::toJSON(json.res);
  json.nm <- paste0(name, ".json");
  
  sink(json.nm);
  cat(json.mat);
  sink();

  if(mType == "network"){
    res.mat<-matrix(0, nrow=length(names(current.mset)), ncol=4);
    rownames(res.mat)<-names(current.mset);
    colnames(res.mat)<-c("Total","Hits", "P.Value", "FDR");
    res.mat[,"Total"] <- set.num
    res.mat[,"Hits"] <- hit.num;
    res.mat[,"P.Value"] <- meta.mat[,meta.matcolinx];
    res.mat[,"FDR"] <- padj;
    enr.mat <<- res.mat
    res <- data.frame(Name=as.vector(rownames(meta.mat)), Total=set.num, Hits= hit.num, EnrichmentScore=as.vector(meta.mat[,1]), Pval=as.vector(meta.mat[,meta.matcolinx]), Padj = padj);
    list.genes <<- allmat.vec
    SetListNms();
    netnm <- paste0(netNm, ".json");
    PrepareEnrichNet(netNm, "meta", "mixed");
    fast.write(res.mat, file=paste0(name, ".csv"), row.names=F);

  }else{
    res.mat<-matrix(0, nrow=length(names(current.mset)), ncol=5);
    colnames(res.mat)<-c("Name", "Total","Hits", "P.Value", "FDR");
    res.mat[,"Name"] <- names(current.mset);
    res.mat[,"Total"] <- set.num
    res.mat[,"Hits"] <- hit.num;
    res.mat[,"P.Value"] <- meta.mat[,meta.matcolinx];
    res.mat[,"FDR"] <- padj;
    fast.write(res.mat, file=paste("meta_sig_genesets_", lib, ".csv", sep=""), row.names=F);
  }
  return(.set.mSet(dataSet))
}
