##################################################
## R script for NetworkAnalyst
## Description: Gene/Probe/Protein ID Annotation
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

.loadEnrichLib <- function(fun.type){
  
  
  my.path <- paste(lib.path, data.org, "/", fun.type, ".rds", sep="");
  
  my.lib <- readRDS(my.path);
  
  if(substr(fun.type, 0, 2)=="go"){
    if(is.null(names(my.lib))){ # some go lib does not give names
      names(my.lib) <- c("link", "term", "sets");
    }
  }
  
  current.geneset <- my.lib$sets;
  set.ids<- names(current.geneset); 
  names(set.ids) <- names(current.geneset) <- my.lib$term;
  
  if(substr(fun.type, 0, 2)=="go"){
    names(current.geneset) = firstup(names(current.geneset))
    names(current.geneset) = gsub("-", "_", names(current.geneset))
    names(set.ids) = firstup(names(set.ids));
    names(set.ids) = gsub("-", "_", names(set.ids));
  }
  qs::qsave(current.geneset, "current_geneset.qs");
  
  current.setlink <<- my.lib$link;
  current.setids <<- set.ids;
  return(current.geneset);
}

# geneIDs is text one string, need to make to vector
performGene2ProteinMapping <- function(listNm, geneIDs, org, type){
  
  dataSet <- list();
  dataSet$orig <- geneIDs;
  current.msg <<- NULL;
  data.org <<- org;
  listNms = vector();
  dataList <- .parseListInput(geneIDs);
  all.prot.mat <- list(); 
  for(i in 1:length(dataList)){
    dataSet$name = paste0("datalist", i);
    listNms[i] = dataSet$name;
    gene.mat <- prot.mat <- dataList[[i]];
    GeneAnotDB <-convertIdToEntrez(rownames(gene.mat), type);
    
    na.inx <- is.na(GeneAnotDB[,1]) | is.na(GeneAnotDB[,2]);
    if(sum(!na.inx) < 2){
      current.msg <<- "Less than two hits found in uniprot database. ";
      print(current.msg);
    }
    rownames(prot.mat) <- GeneAnotDB[,2];
    prot.mat <- prot.mat[!na.inx, , drop=F];
    
    # now merge duplicates
    prot.mat <- .removeDuplicates(prot.mat, "mean", quiet=T); 
    
    if(is.null(dim(prot.mat))){
      prot.mat <- matrix(prot.mat);
    }
    
    seed.proteins <- rownames(prot.mat);
    dataSet$GeneAnotDB <- GeneAnotDB;
    dataSet$sig.mat <- gene.mat;
    dataSet$prot.mat <- prot.mat;
    dataSet$seeds.proteins <- seed.proteins;
    if(i == 1){
      all.prot.mat <- prot.mat;
      totalseed.proteins <- seed.proteins;
    }else{
      totalseed.proteins <- c(totalseed.proteins, seed.proteins);
      all.prot.mat <- rbind(all.prot.mat, prot.mat)
    }
    RegisterData(dataSet); 
  }
  all.ent.mat <<- all.prot.mat
  rownames(all.prot.mat) = doEntrez2SymbolMapping(rownames(all.prot.mat))
  all.prot.mat <- data.frame(as.numeric(all.prot.mat[,1]), rownames(all.prot.mat));
  all.prot.mat <<- all.prot.mat
  listNms <<- listNms
  mdata.all <- list();
  for(i in 1:length(listNms)){
    mdata.all[i] <- 1;
  }
  names(mdata.all) <- listNms;
  partialToBeSaved <<- c(partialToBeSaved, listNms);
  mdata.all <<- mdata.all
  
  if(.on.public.web){
    return(.set.mSet(dataSet));
  }else{
    dataSet <<- dataSet;
    return(totalseed.proteins)
  }
}

doAnnotation <- function(id.vec, idType){
  if(idType %in% c("entrez", "symbol", "refseq", "gb", "embl_gene","embl_protein", "embl_transcript", "orf", "tair", "wormbase", "ko", "custom")){
    anot.id <- doGeneIDMapping(id.vec, idType);
  }else{
    anot.id <- doProbeMapping(id.vec, idType);
    names(anot.id) <- id.vec;
  }
  return(anot.id);        
}

# from probe ID to entrez ID 
doProbeMapping <- function(probe.vec, platform){
  platform.path <- paste(lib.path,  data.org, "/", platform, ".rds", sep="");
  probe.map <- readRDS(platform.path);
  if(is.null(probe.vec)){
    entrez <- probe.map[, "entrez"];
  }else{
    hit.inx <- match(probe.vec, probe.map[, "probe"]);
    entrez <- probe.map[hit.inx, "entrez"];
  }
  return(entrez);
}


queryGeneDB <- function(table.nm, data.org){
  if(table.nm == "custom" || data.org == "custom"){
    db.map <- qs::qread("anot_table.qs");
  }else{
    require('RSQLite');
    conv.db <- dbConnect(SQLite(), paste(sqlite.path, data.org, "_genes.sqlite", sep="")); 
    db.map <- dbReadTable(conv.db, table.nm);
    dbDisconnect(conv.db); cleanMem();
  }
  
  return(db.map)
}

# mapping between genebank, refseq and entrez
doGeneIDMapping <- function(q.vec, type){
  set.seed(1)
  if(is.null(q.vec)){
    db.map <-  queryGeneDB("entrez", data.org);
    q.vec <- db.map[, "gene_id"];
    type <-"entrez";
  }
  if(type == "symbol"){
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "symbol"]);
  }else if(type == "entrez"){
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else if(type == "ko"){ # only for ko
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else if(type == "custom"){ # only for ko
    db.map <-  queryGeneDB("custom", data.org);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else{
    # note, some ID can have version number which is not in the database
    # need to strip it off NM_001402.5 => NM_001402
    q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
    q.vec <- q.mat[,1];
    if(type == "gb"){
      db.map <-  queryGeneDB("entrez_gb", data.org);
    }else if(type == "refseq"){
      db.map <-  queryGeneDB("entrez_refseq", data.org);
    }else if(type == "embl_gene"){
      db.map <-  queryGeneDB("entrez_embl_gene", data.org);
    }else if(type == "embl_transcript"){
      db.map <-  queryGeneDB("entrez_embl_transcript", data.org);
    }else if(type == "embl_protein"){
      db.map <-  queryGeneDB("entrez_embl_protein", data.org);
    }else if(type == "orf"){ # only for yeast
      db.map <-  queryGeneDB("entrez_orf", data.org);
      db.path <- paste(lib.path, data.org, "/entrez_orf.rds", sep="");
    }else if(type == "tair"){ # only for ath
      db.map <-  queryGeneDB("tair", data.org);
    }else if(type == "wormbase"){ # only for cel
      db.map <-  queryGeneDB("entrez_wormbase", data.org);
    }else{
      print("Unknown data type");
      return(0);
    }
    hit.inx <- match(q.vec, db.map[, "accession"]);
  }
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  rm(db.map, q.vec); gc();
  #print(head(entrezs));
  return(entrezs);
}

doIdMapping <- function(q.vec, type){
  if(type %in% c("mir_acc", "mir_id", "mirnet")){
    require('RSQLite');
    path <- paste0(sqlite.path, "mir2gene.sqlite")
    mir.db <- dbConnect(SQLite(), path);
    if(type == "mir_id"){
      q.vec <- tolower(q.vec)
    }
    query <- paste (shQuote(q.vec),collapse=",");
    table.n <- data.org
    statement <- paste("SELECT * FROM ", data.org, " WHERE ((",type," IN (",query,")) OR (mir_id IN (", query, ")))", sep="");
    mirtable <- dbSendQuery(mir.db, statement);       
    mir.dic <- fetch(mirtable, n=-1);
    entrezs <- mir.dic[c("mir_acc", type)];
    paste("mirLength:", nrow(entrezs));
    entrezs <- entrezs[!duplicated(entrezs[,"mir_acc"]),]
    rownames(entrezs) = seq.int(nrow(entrezs));
    entrezs <- data.frame(lapply(entrezs, as.character), stringsAsFactors=FALSE)
    colnames(entrezs) = c("gene_id", "accession");
    
    return(entrezs[,1]);
  }
}

doEntrez2SymbolMapping <- function(entrez.vec){
  gene.map <-  queryGeneDB("entrez", data.org);
  gene.map[] <- lapply(gene.map, as.character)
  
  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  symbols <- gene.map[hit.inx, "symbol"];
  
  # if not gene symbol, use id by itself
  na.inx <- is.na(symbols);
  symbols[na.inx] <- entrez.vec[na.inx];
  return(symbols);
}

checkEntrezMatches <- function(entrez.vec){
  gene.map <-  queryGeneDB("entrez", data.org);
  gene.map[] <- lapply(gene.map, as.character)
  
  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  
  return(length(hit.inx));
}

doSymbol2EntrezMapping <- function(symbol.vec){
  db.map <-  queryGeneDB("entrez", data.org);
  db.map[] <- lapply(db.map, as.character)
  
  hit.inx <- match(symbol.vec, db.map[, "symbol"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  
  na.inx <- is.na(entrezs);
  entrezs[na.inx] <- symbol.vec[na.inx];
  
  mode(entrezs) <- "character";
  return(entrezs);
}


# note, entrez.vec could contain NA/null, cannot use rownames
doEntrezIDAnot <- function(entrez.vec){
  gene.map <-  queryGeneDB("entrez", data.org);
  
  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  anot.mat <- gene.map[hit.inx, c("gene_id", "symbol", "name")];
  
  na.inx <- is.na(hit.inx);
  anot.mat[na.inx, "symbol"] <- entrez.vec[na.inx];
  anot.mat[na.inx, "name"] <- 'NA';
  return(anot.mat);
}

convertIdToEntrez <- function(q.vec, type){ #convert user input ids to entrez
  if(type == "entrez"){
    # need to get only our data
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
    
  }else if(type == "symbol"){
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "symbol"]);
    
  }else if(type == "ko"){ # only for ko
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else if(type == "custom"){ # only for ko
    db.map <-  queryGeneDB("custom", data.org);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else{
    if(type == "gb"){
      # note, some ID can have version number which is not in the database
      # need to strip it off NM_001402.5 => NM_001402
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.map <-  queryGeneDB("entrez_gb", data.org);
    }else if(type == "refseq"){
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.map <-  queryGeneDB("entrez_refseq", data.org);
    }else if(type == "embl_gene"){
      db.map <-  queryGeneDB("entrez_embl_gene", data.org);
    }else if(type == "tair"){
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.map <-  queryGeneDB("tair", data.org);
    }else if(type == "embl_transcript"){
      db.map <-  queryGeneDB("entrez_embl_transcript", data.org);
    }else if(type == "embl_protein"){
      db.map <-  queryGeneDB("entrez_embl_protein", data.org);
    }else if(type == "orf"){ # only for yeast
      db.map <-  queryGeneDB("entrez_orf", data.org);
    }else if(type == "flybase"){
      db.map <-  queryGeneDB("entrez_flybase", data.org);
    }else if(type == "string"){ 
      db.map <-  queryGeneDB("entrez_string", data.org);
    }else if(type == "ecogene"){ # only for ecoli
      db.map <-  queryGeneDB("entrez_ecogene", data.org);
    }else if(type == "uniprot"){
      db.map <-  queryGeneDB("entrez_uniprot", data.org);
    }else if(type == "paelocus"){
      db.map <-  queryGeneDB("entrez", data.org);
    }else{
      print("Unknown data type");
      return(0);
    }
    hit.inx <- match(q.vec, db.map[, "accession"]);
    entrezs <- db.map[hit.inx, ];
  }
  entrezs <- db.map[hit.inx, ]; 
  if(type == "entrez"){
    entrezs = entrezs[,c(1,1)];
  }else{
    entrezs = entrezs[,c(2,1)];
  }
  na.inx <- is.na(entrezs[,1]);
  entrezs[,1][na.inx] <- q.vec[na.inx];
  na.inx <- is.na(entrezs[,2]);
  entrezs[,2][na.inx] <- q.vec[na.inx];
  colnames(entrezs) <- c("accession", "gene_id")
  #print(head(entrezs));
  #print("abc");
  return(entrezs);
}

PerformListAnnot <- function(listNm, org, geneIDs, type){
  dataSet <- list();
  dataSet$orig <- "";
  current.msg <<- NULL;
  data.org <<- org;
  listNms <- multiFileNamesU;
  numOfLists <<-length(multiFileNamesU);
  notOk = 0
  for(i in 1:length(listNms)){
    dataSet = qs::qread(listNms[i])
    dataSet$name = listNms[i];
    gene.mat <- prot.mat <- dataSet$prot.mat;
    GeneAnotDB <-convertIdToEntrez(rownames(gene.mat), type);
    
    na.inx <- is.na(GeneAnotDB[,1]) | is.na(GeneAnotDB[,2]);
    if(sum(!na.inx) < 2){
      current.msg <<- paste0("Less than two hits found in database for ", listNms[i]);
      print(current.msg);
      return(0);
    }
    rownames(prot.mat) <- GeneAnotDB[,2];
    prot.mat <- prot.mat[!na.inx, , drop=F];
    
    # now merge duplicates
    listInxU <<- listNms[i];
    prot.mat <- .removeDuplicates(prot.mat, "mean", quiet=F); 
    
    if(is.null(dim(prot.mat))){
      prot.mat <- matrix(prot.mat);
    }
    
    seed.proteins <- rownames(prot.mat);
    dataSet$GeneAnotDB <- GeneAnotDB;
    dataSet$sig.mat <- gene.mat;
    dataSet$prot.mat <- prot.mat;
    dataSet$seeds.proteins <- seed.proteins;
    if(i == 1){
      all.prot.mat <- prot.mat;
      totalseed.proteins = seed.proteins
    }else{
      totalseed.proteins  = c(totalseed.proteins, seed.proteins);
      all.prot.mat <- rbind(all.prot.mat, prot.mat)
    }
    RegisterData(dataSet); 
  }
  all.ent.mat <<- all.prot.mat
  rownames(all.prot.mat) = doEntrez2SymbolMapping(rownames(all.prot.mat))
  all.prot.mat <- data.frame(as.numeric(all.prot.mat[,1]), rownames(all.prot.mat));
  all.prot.mat <<- all.prot.mat
  listNms <<- listNms
  mdata.all <- list();
  for(i in 1:length(listNms)){
    mdata.all[i] <- 1;
  }
  names(mdata.all) <- listNms;
  mdata.all <<- mdata.all
  return(totalseed.proteins)
}

#########################################
##########################################
############# private utility methods #### 
##########################################
##########################################
# given a text from input area: one or two-col entries (geneID and logFC)
# parse out return the a matrix containing the logFc, with rows named by gene ID
# if no 2nd col (logFC), the value will be padded by 0s
.parseListInput <- function(geneIDs){
  
  spl <- unlist(strsplit(geneIDs, "\\//")[1]);
  spl <- spl[unlist(lapply(spl,function(x){!x %in% ""}))]
  spl <- lapply(spl,function(x){gsub("\\/", "",x)})
  numOfLists <<- length(spl)
  dataList <- list();
  inxU <- 0;
  for (i in 1:length(spl)){
    lines <- unlist(strsplit(spl[[i]], "\r|\n|\r\n")[1]);
    # remove the beginning & trailing space 
    lines <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", lines, perl=TRUE);
    if(substring(lines[1],1,1)=="#"){
      lines <- lines[-1];
    }
    gene.lists <- strsplit(lines, "\\s+");
    gene.mat <- do.call(rbind, gene.lists);
    
    if(dim(gene.mat)[2] == 1){ # add 0
      gene.mat <- cbind(gene.mat, rep(0, nrow(gene.mat)));
      current.msg <- "Only one column found in the list. Abundance will be all zeros. ";
    }else if(dim(gene.mat)[2] > 2){
      gene.mat <- gene.mat[,1:2];
      current.msg <- "More than two columns found in the list. Only first two columns will be used. ";
    }
    print(current.msg);
    
    rownames(gene.mat) <- gene.mat[,1];
    gene.mat <- gene.mat[,-1, drop=F];
    inxU <- inxU + 1;
    listInxU <<- paste0("datalist", inxU);
    gene.mat <- .removeDuplicates(gene.mat, "mean", quiet=F); 
    good.inx <- !is.na(gene.mat[,1]);
    gene.mat <- gene.mat[good.inx, , drop=F];
    dataList[[i]] <- gene.mat  
  }
  return(dataList)
}

##########################################
############# private utility methods #### 
##########################################

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
