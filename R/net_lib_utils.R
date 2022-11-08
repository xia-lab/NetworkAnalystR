##################################################
## R script for NetworkAnalyst
## Description: Functions to load various libraries for functional enrichment analysis during network visualization
## Authors:
## G. Zhou (guangyan.zhou@mail.mcgill.ca)
## J. Xia, jeff.xia@mcgill.ca
###################################################

# note, last two par only for STRING database
QueryPpiSQLite <- function(sqlite.path, table.nm, q.vec, requireExp, min.score){
  require('RSQLite')
  ppi.db <- dbConnect(SQLite(), paste(sqlite.path, "ppi.sqlite", sep="")); 
  query <- paste(shQuote(q.vec),collapse=",");
  
  if(grepl("string$", table.nm)){
    if(requireExp){
      statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND combined_score >=", min.score, " AND experimental > 0", sep="");
    }else{
      statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND combined_score >=", min.score, sep="");
    }
  }else{
    statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))", sep="");
  }
  ppi.res <- .query.sqlite(ppi.db, statement);
  
  # remove dupliated edges
  ppi.res <- ppi.res[!duplicated(ppi.res$row_id),]
  return(ppi.res);  
}

#for signaling pathway as well
QueryOtherPpiSQLite <- function(sqlite.path, table.nm, q.vec, data.org){
  require('RSQLite')
    if(grepl("signal$", table.nm)){
        ppi.db <- dbConnect(SQLite(), paste(sqlite.path, "signor.sqlite", sep="")); 
        table.nm <- data.org
    }else{
        ppi.db <- dbConnect(SQLite(), paste(sqlite.path, "cross_species_ppi.sqlite", sep="")); 
    }
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste('SELECT * FROM "',table.nm, '" WHERE ((id1 IN (', query, ')) OR (id2 IN (', query, ')))', sep='');
  ppi.res <- .query.sqlite(ppi.db, statement);
  return(ppi.res);  
}

# table name is org code, id.type is column name
QueryMirSQLite <- function(sqlite.path, org, id.type, q.vec, db.nm){
  require('RSQLite');
  mir.db <- dbConnect(SQLite(), paste(sqlite.path, "mir2gene.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  #statement <- paste("SELECT * FROM ", org, " WHERE ",id.type," IN (",query,")", " AND ", db.nm," == 1 ", sep="");
  statement <- paste("SELECT * FROM ", org, " WHERE ",id.type," IN (",query,")", " AND ", db.nm," == 1 ", sep="");
  return(.query.sqlite(mir.db, statement));
}

# table name is org code, id.type is column name
QueryDrugSQLite <- function(sqlite.path, q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "drug.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM human WHERE upid IN (",query,")", sep="");
  return(.query.sqlite(drug.db, statement));
}

QueryDiseaseSQLite <- function(sqlite.path, q.vec){
  require('RSQLite');
  dis.db <- dbConnect(SQLite(), paste(sqlite.path, "disease.sqlite", sep=""));
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM human WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(dis.db, statement));
}

QueryTfmirSQLite <- function(sqlite.path, q.vec, data.org){
  require('RSQLite');
  paramSet <- readSet(paramSet, "paramSet");
  sqlite.path <- paramSet$sqlite.path;
  tf.db <- dbConnect(SQLite(), paste(sqlite.path, "tfmirgene.sqlite", sep=""));
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", data.org, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  return(.query.sqlite(tf.db, statement));
}

QueryDiffNetSQLite <- function(sqlite.path, q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "tissuePPI.sqlite", sep=""));
  table.nm <- diffNetName;
  query <- paste (shQuote(q.vec),collapse=",");
  topPct <- 1-diffPct;
  botstatement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND rank <=", diffPct ,sep="");
  topstatement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND rank >=", topPct ,sep="");
  
  drug.dic1 <- .query.sqlite(drug.db, botstatement, FALSE);# no close db connection
  drug.dic2 <- .query.sqlite(drug.db, topstatement);
  drug.dic <- rbind(drug.dic1, drug.dic2);
  return(drug.dic);
}

QueryCellCoexSQLite <- function(sqlite.path, q.vec, data.org){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, data.org,"_immune.sqlite", sep=""));
  tblNm <- paste0(data.org,"_",cellCoexNumber);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", tblNm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  return(.query.sqlite(drug.db,statement));
}

QueryTissueCoexSQLite <- function(sqlite.path, q.vec, data.org){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "tissueCoex.sqlite", sep=""));
  tblNm <- paste0(data.org,"_",cellCoexNumber);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", tblNm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  return(.query.sqlite(drug.db, statement));
}

QueryChemSQLite<- function(sqlite.path, org, q.vec){
  require('RSQLite');
  chem.db <- dbConnect(SQLite(), paste(sqlite.path, "chem.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", org, " WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(chem.db, statement));
}

QueryTFSQLite<- function(sqlite.path, table.nm, q.vec){
  require('RSQLite');
  chem.db <- dbConnect(SQLite(), paste(sqlite.path, "tf2gene.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(chem.db, statement));
}

doPpiIDMapping <- function(sqlite.path, q.vec, data.org="entrez_swissprot"){
  paramSet <- readSet(paramSet, "paramSet");
  data.org <- paramSet$data.org;
  invertAcc <-F;
  idType <- ""
  if(data.org == "sce"){
    if(net.type == "string"){ # only for yeast
      idType <- "entrez";
    }else{
      idType <- "entrez_uniprot";
    }
  }else if(net.type %in% c("irefinx", "rolland")){
    idType <- "entrez_uniprot";
  }else{
    if(data.org == "pae" & net.type == "interactome"){
      idType <- "entrez_string";
    }else if(data.org %in% c("tbr", "cel", "pfa")){
      idType <- "entrez_string";
    }else{
      idType <- "entrez";
      invertAcc<- T;
    }
  }
  db.map <-  queryGeneDB(idType, data.org);
  hit.inx <- match(q.vec, db.map[, "gene_id"]);
  ppi.mat <- db.map[hit.inx, ];

if(idType == "entrez"){
  return(ppi.mat$gene_id);
  }else{
  return(ppi.mat$accession);
}
}

doUniprot2EntrezMapping<-function(uniprot.vec){
  db.map <-  queryGeneDB("entrez_uniprot", data.org);
  hit.inx <- match(uniprot.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  na.inx <- is.na(entrezs);
  entrezs[na.inx] <- uniprot.vec[na.inx];
  return(entrezs);
}

doEntrez2UniprotMapping<-function(entrez.vec){
  db.map <-  queryGeneDB("entrez_uniprot", data.org);
  hit.inx <- match(entrez.vec, db.map[, "gene_id"]);
  entrezs <- db.map[hit.inx, "accession"];
  mode(entrezs) <- "character";
  na.inx <- is.na(entrezs);
  entrezs[na.inx] <- entrez.vec[na.inx];
  return(entrezs);
}

doEntrez2UniprotMapping <- function(entrez.vec){
  db.map <-  queryGeneDB("entrez_uniprot", data.org);
  hit.inx <- match(entrez.vec, db.map[, "gene_id"]);
  unips <- db.map[hit.inx, "accession"];
  mode(unips) <- "character";
  return(unips);
}

doString2EntrezMapping <- function(string.vec){

  db.map <-  queryGeneDB("entrez_string", data.org);
  hit.inx <- match(string.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEmblGene2EntrezMapping <- function(emblgene.vec){
  db.map <-  queryGeneDB("entrez_embl_gene", data.org);
  hit.inx <- match(emblgene.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEmblProtein2EntrezMapping <- function(emblprotein.vec){
  db.map <-  queryGeneDB("entrez_embl_protein", data.org);
  hit.inx <- match(emblprotein.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

QueryMetSQLiteNet <- function(sqlite.path, table.nm, q.vec, inv){
    require('RSQLite');
    path <- paste0(sqlite.path, "met.sqlite")
    lnc.db <- dbConnect(SQLite(), path);
    query <- paste (shQuote(q.vec),collapse=",");
    
    if(inv == "inverse"){
        statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", sep="");
    }else{
        statement <- paste("SELECT * FROM ", table.nm, " WHERE kegg IN (",query,")", sep="");
    } 

    lnctable <- dbSendQuery(lnc.db, statement);
    lnc.dic <- fetch(lnctable, n=-1); # get all records
    lnc.dic <<- lnc.dic
    dbDisconnect(lnc.db);
    return(lnc.dic);
}


##########################################
############# private utility methods #### 
##########################################

.query.sqlite <- function(db.con, statement, offline=TRUE){
  rs <- dbSendQuery(db.con, statement);
  res <- fetch(rs, n=-1); # get all records
  dbClearResult(rs);
  if(offline){
    dbDisconnect(db.con);
  }
  cleanMem();
  return(res);
}