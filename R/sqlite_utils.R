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

.connect.sqlite <- function(db.path){
  if(!PrepareSqliteDB(db.path, paramSet$on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }
  return(dbConnect(SQLite(), db.path)); 
}

PrepareSqliteDB <- function(sqlite_Path, onweb = TRUE) {
  if(onweb) {return(TRUE)};
  if(file.exists(sqlite_Path)) {return(TRUE)};

  dbNM <- basename(sqlite_Path);
  DonwloadLink <- paste0("https://www.xialab.ca/resources/sqlite/", dbNM);
  download.file(DonwloadLink, sqlite_Path);
  return(TRUE)
}

queryGeneDB <- function(table.nm, data.org){
  paramSet <- readSet(paramSet, "paramSet");  
  if(length(table.nm) == 0){
    table.nm <- "";
  }

  if(table.nm == "custom" || data.org == "custom"){
    db.map <- qs::qread("anot_table.qs");
  }else{
    require('RSQLite');
    db.path <- paste(paramSet$sqlite.path, data.org, "_genes.sqlite", sep="")

    if(!PrepareSqliteDB(db.path, paramSet$on.public.web)){
      stop("Sqlite database is missing, please check your internet connection!");
    }
    conv.db <- dbConnect(SQLite(), db.path); 
    tbls <- dbListTables(conv.db)
    if(!table.nm %in% tbls){
        return(0);
    }
    db.map <- dbReadTable(conv.db, table.nm);
    dbDisconnect(conv.db); cleanMem();
  }
  return(db.map)
}