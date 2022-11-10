##################################################
## R scripts for NetworkAnalyst
## Description: biological network analysis methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# table.nm is the org code used for sqlite table (ppi)
# for chem type, table.nm is drugbank or ctd
# note, last two param only for STRING database

BuildMultiNet <- function(tblNm1="NA", tblNm2="NA", tblNm3="NA", tblNm4="NA", order=1){

  opts.vec <- c(tblNm1, tblNm2, tblNm3, tblNm4);
  opts.vec <- strsplit(opts.vec, "__");
  tbl.vec <- vector();
  dbNm.vec <- vector();

  for(i in 1:length(opts.vec)){
    dbNm.vec[i] <- opts.vec[[i]][1];
    tbl.vec[i] <- opts.vec[[i]][2];
  }

  node.resu <<- "empty";
  edge.resu <<- "empty";
  res.list <- list();
  for ( i in 1:length(tbl.vec)){
    if(tbl.vec[i] != "NA"){
      net.type <<- tbl.vec[i];
      if(dbNm.vec[i] == "ppi"){
        tblNm <- paste0(data.org, "_", tbl.vec[i]); 
      }else{
        tblNm <- tbl.vec[i];
      }
      res = SearchNetDB(dbNm.vec[i], tblNm, require.exp=TRUE, min.score = 900, order=1);
      res.list[[i]] = res;
    }
  }
  return(res);
}

#' Query database using input list or previously computed network
#'
#' @param db.type Input the name of the created dataSetObj (see Init.Data)
#' @param type Network type (gene, met, mir, tf, mic, peak, m2m, snp)
#' @param dbType Database name (i.e innatedb)
#' @param inputType Omics type of input features (gene, protein, met, tf, mic, peak, snp)
#'
#' @export
#'
SearchNetDB <- function(db.type, table.nm, require.exp=TRUE, min.score = 900, order=1){

    regids<<- "";
    db.typeu <<- db.type;
    result.list <- .prepareSigProteinJSON();

    paramSet <- readSet(paramSet, "paramSet");
    data.org <- paramSet$data.org;
    sqlite.path <- paramSet$sqlite.path;
    lib.path <- paramSet$lib.path;

    net.type <- gsub(paste0(data.org,"_"),'',table.nm);
    SetNetType(net.type);

    protein.vec <- result.list$protein.vec; # this actually is entrez IDs?
    seed.proteins <<- protein.vec;
    # now do the database search
    if(db.type == "ppi"){
        protein.vec <- doPpiIDMapping(sqlite.path, protein.vec, data.org);
        res <- QueryPpiSQLite(sqlite.path, table.nm, protein.vec, require.exp, min.score);
        # no hits
      if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,1],Target=res[,2]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);
        node.ids <- c(res[,1], res[,2])
        node.nms <- c(res[,3], res[,4]);
        node.types <- c(rep("Protein", nrow(res)), rep("Protein", nrow(res)))
    }else if(db.type == "tf"){
        table.nm <- paste(data.org,table.nm, sep="_");
        res <- QueryTFSQLite(sqlite.path, table.nm, protein.vec);
        # no hits
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,"entrez"],Target=res[,"tfid"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);
        node.ids <- c(res[,"entrez"], res[,"tfid"])
        node.nms <- c(res[,"symbol"], res[,"tfname"]);
        node.types <- c(rep("Protein", nrow(res)), rep("TF", nrow(res)))
    }else if(db.type == "crossppi"){
        res <- QueryOtherPpiSQLite(sqlite.path, table.nm, protein.vec, data.org);
        # no hits
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,"id1"],Target=res[,"id2"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);
        node.ids <- c(res[,"id1"], res[,"id2"])
        node.nms <- c(res[,"name1"], res[,"name2"]);
        node.types <- c(rep("Protein", nrow(res)), rep("Foreign_Protein", nrow(res)))
        if(table.nm == "microbiolink1"){
            edge.extra.info <- data.frame(domain=c(res[,"PFAM_domain1"], res[,"PFAM_domain2"]))
            node.extra.info <- data.frame(id=node.ids, species=c(rep("Human", nrow(res)), res[,"species2"]), domain=c(res[,"PFAM_domain1"], res[,"PFAM_domain2"]))
            node.infoU <<- node.extra.info[!duplicated(node.extra.info$id),];
        }else if(table.nm == "microbiolink2"){
            edge.extra.info <- data.frame(domain=c(res[,"linear_motif1"], res[,"PFAM_domain2"]))
            node.extra.info <- data.frame(id=node.ids, species=c(rep("Human", nrow(res)), res[,"species2"]), domain=c(res[,"linear_motif1"], res[,"PFAM_domain2"]))
            node.infoU <<- node.extra.info[!duplicated(node.extra.info$id),];
        }
    }else if(db.type == "signal"){
        res <- QueryOtherPpiSQLite(sqlite.path,table.nm, protein.vec, data.org);
        # no hits
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,"id1"],Target=res[,"id2"], Direction=res[, "EFFECT"], Mechanism=res[,"MECHANISM"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);
        node.ids <- c(res[,"id1"], res[,"id2"])
        node.nms <- c(res[,"name1"], res[,"name2"]);
        edge.extra.info <- res[,-c(1:4)];
        node.extra.info <- c(res[,"TYPEA"], res[,"TYPEB"])
        node.infoU <<- node.extra.info
        node.types <- c(res[,"TYPEA"], res[,"TYPEB"])

    }else if(db.type == "mir"){ # in miRNA, table name is org code, colname is id type
        db.nm <- table.nm;
        res <- QueryMirSQLite(sqlite.path, data.org, "entrez", protein.vec, db.nm);
        # no hits
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,"entrez"],Target=res[,"mir_acc"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);
        node.ids <- c(res[,"entrez"], res[,"mir_acc"])
        node.nms <- c(res[,"symbol"], res[,"mir_id"]);
        node.types <- c(rep("Protein", nrow(res)), rep("miRNA", nrow(res)))
    }else if(db.type == "drug"){
        # note, all drug data is on human,
        protein.vec <- doEntrez2UniprotMapping(protein.vec);
        protein.vec <- protein.vec[!is.na(protein.vec)];
        res <- QueryDrugSQLite(sqlite.path, protein.vec);
        # no hits
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=doUniprot2EntrezMapping(res[,"upid"]),Target=res[,"dbid"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);

        node.ids <- c(doUniprot2EntrezMapping(res[,"upid"]), res[,"dbid"]);
        node.nms <- c(res[,"symbol"], res[,"dbname"]);
        node.types <- c(rep("Protein", nrow(res)), rep("Drug", nrow(res)));
        protein.vec <- doUniprot2EntrezMapping(protein.vec);
    }else if(db.type == "disease"){
        res <- QueryDiseaseSQLite(sqlite.path, protein.vec);
        # no hits
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,"entrez"],Target=res[,"diseaseId"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);

        node.ids <- c(res[,"entrez"], res[,"diseaseId"])
        node.nms <- c(res[,"symbol"], res[,"diseaseName"]);
        node.types <- c(rep("Protein", nrow(res)), rep("Disease", nrow(res)))

    }else if(db.type == "tfmir"){
        res <- QueryTfmirSQLite(sqlite.path, protein.vec, data.org);
        # no hits
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,"id1"],Target=res[,"id2"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);

        node.ids <- c(res[,"id1"], res[,"id2"])
        node.nms <- c(res[,"name1"], res[,"name2"]);

        node.types <- c(rep("Protein", 2*nrow(res)))
        db.path <- paste(lib.path, data.org, "/mirlist.rds", sep="");
        db.map <-  readRDS(db.path);
        tf.inx <- node.ids %in% edge.res[,"Source"];
        mir.inx <- node.ids %in% db.map[,1];
        node.types[tf.inx] <- "TF";
        node.types[mir.inx] <- "miRNA";

    }else if(db.type == "cellcoex"){
        res <- QueryCellCoexSQLite(sqlite.path, protein.vec, data.org);
        # no hits
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,"id1"],Target=res[,"id2"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);

        node.ids <- c(res[,"id1"], res[,"id2"])
        node.nms <- c(res[,"name1"], res[,"name2"]);
        node.types <- c(rep("Protein", nrow(res)), rep("Protein", nrow(res)))

    }else if(db.type == "tissueppi"){
        res <- QueryDiffNetSQLite(sqlite.path, protein.vec);
        # no hits
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,"id1"],Target=res[,"id2"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);

        node.ids <- c(res[,"id1"], res[,"id2"])
        node.nms <- c(res[,"name1"], res[,"name2"]);
        node.types <- c(rep("Protein", nrow(res)), rep("Protein", nrow(res)))
    }else if(db.type == "tissuecoex"){
        # note, all drug data is on human,
        res <- QueryTissueCoexSQLite(sqlite.path, protein.vec, data.org);
        # no hits
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,"id1"],Target=res[,"id2"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);

        node.ids <- c(res[,"id1"], res[,"id2"])
        node.nms <- c(res[,"name1"], res[,"name2"]);
        node.types <- c(rep("Protein", nrow(res)), rep("Protein", nrow(res)))
    }else if(db.type == "chem"){
        res <- QueryChemSQLite(sqlite.path, data.org, protein.vec);
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,"entrez"],Target=res[,"ctdid"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);

        node.ids <- c(res[,"entrez"], res[,"ctdid"])
        node.nms <- c(res[,"symbol"], res[,"name"]);
        node.types <- c(rep("Protein", nrow(res)), rep("Chemical", nrow(res)))
    }else if(db.type == "met"){
        table.nm <- paste(data.org, met.type, sep="_");
        res <- QueryMetSQLiteNet(sqlite.path, table.nm, protein.vec, netInv);
        if(nrow(res)==0){ return(c(0,0)); }
        edge.res <- data.frame(Source=res[,"entrez"],Target=res[,"kegg"]);
        row.names(edge.res) <- 1:nrow(res);
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);

        node.ids <- c(res[,"kegg"], res[,"entrez"])
        node.nms <- c(res[,"met"], res[,"symbol"]);
        node.types <- c(rep("Metabolite", nrow(res)), rep("Protein", nrow(res)))
    }
    if(!is.null(node.ids)){
        seed.inx <- node.ids %in% unique(seed.proteins);
        node.types[seed.inx] <- "Seed"
    }

    if(order == 0 && db.type == "ppi"){
        keep.inx <- node.ids %in% protein.vec;
        node.ids <- node.ids[keep.inx];
        node.nms <- node.nms[keep.inx];
        node.types <- node.types[keep.inx];
        keep.inx <- edge.res[,1] %in% protein.vec;
        edge.res <- edge.res[keep.inx,];
        keep.inx <- edge.res[,2] %in% protein.vec;
        edge.res <- edge.res[keep.inx,];
    }

    node.res <- data.frame(Id=node.ids, Label=node.nms, Types=node.types);
    node.res <- node.res[!duplicated(node.res$Id),];

    if( grepl("microbiolink", table.nm, fixed = TRUE) || db.type == "signal"){
        edge.res <- cbind(edge.res, edge.extra.info);
        edge.infoU <<- edge.res;
    }

    nodeListu <<- node.res;
    table.nmu <<- table.nm;

    fast.write(node.res, file="orig_node_list.csv", row.names=FALSE);

    ppi.net <<- list(
        db.type=db.type,
        order=1,
        seeds=protein.vec,
        table.nm=table.nm,
        node.data = node.res,
        edge.data = edge.res,
        require.exp = require.exp,
        min.score = min.score
    );

    paramSet$nets[[db.type]] <- ppi.net
    saveSet(paramSet, "paramSet");
    analSet$ppi.net <- ppi.net;
    output <- c(nrow(node.res), nrow(res))
    return(saveSet(analSet, "analSet", output))
}

.prepareSigProteinJSON <- function(){
    paramSet <- readSet(paramSet, "paramSet");
    anal.type <- paramSet$anal.type;
    if(anal.type == "genelist"){
        result.list <- .prepareListSeeds();
    }else{ # single expression data or meta.mat
        result.list <- .prepareExpressSeeds();
    }
    return(result.list);
}

.prepareListSeeds <- function(){
    protein.list <- list();
    gene.list <- list();
    paramSet <- readSet(paramSet, "paramSet");
    selectedNetDataset <- paramSet$selectedNetDataset;
    msgSet <- readSet(msgSet, "msgSet");

    data.org <- paramSet$data.org;

    mdata.all <- paramSet$mdata.all;

    if(paramSet$numOfLists > 1){
        if(selectedNetDataset %in% c("intersect","union")){
            dataSet <- list();
            dataSet$name <- selectedNetDataset
            my.vec <- names(mdata.all);
            com.ids <- NULL;
            list.vec <- list()
            for(i in 1:length(my.vec)){
                datSet <- readDataset(my.vec[i]);
                if(is.null(com.ids)){
                  com.ids <- datSet$GeneAnotDB[,"gene_id"];
                  prot.mat <- datSet$prot.mat
                  list.vec[[i]] <- com.ids
                }else{
                  if(selectedNetDataset == "intersect"){
                    com.ids <- datSet$GeneAnotDB[,"gene_id"];
                    list.vec[[i]] <- com.ids
                  }else{
                    com.ids <- union(com.ids, datSet$GeneAnotDB[,"gene_id"]);
                  }
                    prot.mat <- rbind(prot.mat, as.matrix(datSet$prot.mat[rownames(datSet$prot.mat) %in% com.ids,]))
                }
           }
            if(selectedNetDataset == "intersect"){
            com.ids <- Reduce(intersect, list.vec)
            prot.mat <- as.matrix(datSet$prot.mat[rownames(datSet$prot.mat) %in% com.ids,])
            }else{
            com.ids <- unique(as.character(com.ids[!is.na(com.ids)])); # make sure it is interpreted as name not index
            }

            com.symbols <- doEntrez2SymbolMapping(com.ids, data.org, paramSet$data.idType);
            dataSet$GeneAnotDB <- data.frame(gene_id=com.ids, accession=com.symbols);
            dataSet$prot.mat <- prot.mat;
            dataSet$sig.mat <- prot.mat
            dataSet$seeds.proteins <- com.ids
        }else{
           my.vec <- names(mdata.all);
           # make sure reference is the first
           inx <- which(my.vec == selectedNetDataset);
           my.vec <- my.vec[-inx];
           com.ids <- NULL;
           ids.list <- list()
           for(i in 1:length(my.vec)){
                dataSet <- readDataset(my.vec[i]);
                ids.list[[i]]=dataSet$GeneAnotDB[,"gene_id"]
           }
            dataSet <- readDataset(selectedNetDataset);
            ids <- unique(unlist(ids.list));
            com.ids <-setdiff(dataSet$GeneAnotDB[,"gene_id"], ids);
            prot.mat <- as.matrix(dataSet$prot.mat[which(rownames(dataSet$prot.mat) %in% com.ids),])
            com.symbols <- doEntrez2SymbolMapping(com.ids, data.org, paramSet$data.idType);
            dataSet$GeneAnotDB <- data.frame(gene_id=com.ids, accession=com.symbols);
            dataSet$prot.mat <- prot.mat;
            dataSet$sig.mat <- prot.mat
            dataSet$seeds.proteins <- com.ids
        }
    }else{
        dataSet <- readDataset("datalist1");
    }
    

    # return a json array object
    # each object for a single dataset its sig proteins
    meta.vec <- meta.gene.vec <- meta.seed.expr <- NULL;
    file.create("seed_proteins.txt");
    GeneAnotDB <- NULL;

    gene.mat <- dataSet$sig.mat;
    prot.mat <- dataSet$prot.mat;
    write(paste("#DataSet:", dataSet$name),file="sig_genes.txt",append=TRUE);
    write.table(dataSet$sig.mat, file="sig_genes.txt", append=TRUE);

    meta.gene.vec <- c(meta.gene.vec, rownames(gene.mat));
    gene.list[[dataSet$name]] <- list(gene=rownames(gene.mat),logFC=unname(gene.mat[,1]));
    GeneAnotDB <- rbind(GeneAnotDB, dataSet$GeneAnotDB);
    meta.seed.expr <- c(meta.seed.expr, prot.mat[,1]);
    write(paste("#DataSet:", dataSet$name),file="seed_proteins.txt",append=TRUE);
    write.table(cbind(Emblprotein=rownames(prot.mat), Expression=prot.mat[,1]), file="seed_proteins.txt", row.names=F, quote=F,append=TRUE);
    protein.vec <- prot.mat[,1];
    meta.vec <- c(meta.vec, names(protein.vec));
    if(length(protein.vec) == 1){
        protein.vec <- as.matrix(protein.vec)
    }
    protein.list[[dataSet$name]] <- signif(protein.vec, 3);

    gene.list$name <- dataSet$name;
    paramSet$seed.genes <- unique(meta.gene.vec);

    meta.seed.df <- as.matrix(meta.seed.expr);
    rownames(meta.seed.df) <- names(meta.seed.expr);

    res <- RemoveDuplicates(meta.seed.df, "max", quiet=F, paramSet, msgSet);
    seed.expr <- res[[1]];
    msgSet <- res[[2]];
    paramSet$seed.expr <- seed.expr[,1];
    protein.vec <- unique(meta.vec);

    result <- list(
        gene.list = gene.list,
        protein.list = protein.list,
        protein.vec = protein.vec
    );

    saveSet(paramSet, "paramSet");
    return(result)
}


# 2nd order network, only for PPI
ExpandNetworkSearch <- function(){
    paramSet <- readSet(paramSet, "paramSet");
    data.org <- paramSet$data.org;
    sqlite.path <- paramSet$sqlite.path
    if(ppi.net$order == 0){
        SearchNetDB(ppi.net$db.type, ppi.net$table.nm);
        CreateGraph();
    }
    protein.vec <- V(overall.graph)$name;
    table.nm <- ppi.net$table.nm;
    if(db.typeu == "ppi"){
        res <- QueryPpiSQLite(sqlite.path, table.nm, protein.vec, ppi.net$require.exp, ppi.net$min.score);
    }else if (db.typeu == "tissuecoex"){
        res <- QueryTissueCoexSQLite(sqlite.path, protein.vec, data.org);
    }else if (db.typeu == "tissueppi"){
        res <- QueryDiffNetSQLite(sqlite.path, protein.vec);
    }else if (db.typeu == "cellcoex"){
        res <- QueryCellCoexSQLite(sqlite.path, protein.vec, data.org);
    }
    edge.res <- data.frame(Source=res[,1],Target=res[,2]);
    #row.names(edge.res) <- res[,5];

    node.ids <- c(res[,1], res[,2])
    node.nms <- c(res[,3], res[,4]);
    node.res <- data.frame(Id=node.ids, Label=node.nms);
    node.res <- node.res[!duplicated(node.res$Id),];
    if(nrow(node.res) < 10000){ # only overwrite if within the range
        fast.write(edge.res, file="orig_edge_list.csv",row.names=FALSE);
        fast.write(node.res, file="orig_node_list.csv", row.names=FALSE);
        ppi.net <<- list(db.type=db.typeu, order=2, seeds=protein.vec, table.nm=table.nm, node.data = node.res, edge.data = edge.res);
    }
    analSet$ppi.net <- ppi.net;
    output <- nrow(node.res);
    return(saveSet(analSet, "analSet", output));
}

# zero-order network - create ppi nets from only input (seeds)

#' Create network from only input (seeds)
#'
#' @export
#'
BuildSeedProteinNet <- function(){

    nodes <- V(overall.graph)$name;
    hit.inx <-  nodes %in% seed.proteins;
    nodes2rm <- nodes[!hit.inx];
    g <- simplify(delete.vertices(overall.graph, nodes2rm));

    nodeList <- get.data.frame(g, "vertices");
    nodeList <- nodeList[,1:2];
    colnames(nodeList) <- c("Id", "Label");

    fast.write(nodeList, file="orig_node_list.csv");
    nd.inx <- ppi.net$node.data[,1] %in% nodeList[,1];

    edgeList <- get.data.frame(g, "edges");
    edgeList <- edgeList[,1:2];
    colnames(edgeList) <- c("Source", "Target");
    fast.write(edgeList, file="orig_edge_list.csv");

    # update ppi.net
    ppi.net$order <- 0;
    ppi.net$node.data <- nodeList;
    ppi.net$edge.data <- edgeList;
    ppi.net <<- ppi.net;
    analSet$ppi.net <- ppi.net;
    return(saveSet(analSet, "analSet"));
}

# create igraph from the edgelist saved from graph DB
# and decompose into subnets

#' Create igraph object from edgelist created from database selection and decompose into connected subnetworks
#'
#' @export
#'
CreateGraph <- function(){
    require('igraph');
    paramSet <- readSet(paramSet, "paramSet");
    data.org <- paramSet$data.org;
    seed.expr <- paramSet$seed.expr;
    seed.genes <- paramSet$seed.genes;

    node.list <- ppi.net$node.data;
    edge.list <- ppi.net$edge.data;

    seed.proteins <- ppi.net$seeds;
    overall.graph <- simplify(graph.data.frame(edge.list, directed=FALSE, vertices=node.list));

    if("domain" %in% colnames(edge.list)){
        overall.graph <- simplify(graph.data.frame(edge.list, directed=F, vertices=node.list));
        overall.graph <- set.vertex.attribute(overall.graph, "species", index = V(overall.graph), value = node.infoU$species[which(V(overall.graph)$name %in% node.infoU$id)]);
        overall.graph <- set.vertex.attribute(overall.graph, "domain", index = V(overall.graph), value = node.infoU$domain[which(V(overall.graph)$name %in% node.infoU$id)]);
    }else if("Direction" %in% colnames(edge.list)){
        overall.graph <- simplify(graph.data.frame(edge.list, directed=F, vertices=node.list));
        overall.graph <-set_edge_attr(overall.graph, "direction",index = E(overall.graph), value = as.character(edge.list$Direction))
        overall.graph <-set_edge_attr(overall.graph, "mechanism",index = E(overall.graph), value = as.character(edge.infoU$Mechanism))
    }else if("TYPEA" %in% colnames(edge.list)){
        overall.graph <- simplify(graph.data.frame(edge.list, directed=F, vertices=node.list));
        overall.graph <-set_edge_attr(overall.graph, "direction",index = E(overall.graph), value = as.character(edge.infoU$EFFECT))
        overall.graph <-set_edge_attr(overall.graph, "mechanism",index = E(overall.graph), value = as.character(edge.list$Mechanism))

        overall.graph <- set.vertex.attribute(overall.graph, "type", index = V(overall.graph), value = node.infoU);
    }else {
        overall.graph <- simplify(graph.data.frame(edge.list, directed=FALSE, vertices=node.list));
        overall.graph <-set_edge_attr(overall.graph, "direction",index = E(overall.graph), value = "NA")
    }

    # add node expression value
    if(ppi.net$db.type == "ppi"){
        ## Temp fix seed.expr all in entrez?
        ## all converted to new IDs based on PPI from the given organism
        newIDs <- doPpiIDMapping(sqlite.path, names(seed.expr), data.org)
    }else{# all entrez in mirNet
        newIDs <- names(seed.expr);
    }
    match.index <- match(V(overall.graph)$name, newIDs);
    expr.vals <- seed.expr[match.index];
 #  expr.vec <- abs(expr.vals) # logFC can be negative!!
    expr.vec <- expr.vals;
    names(expr.vec)<- V(overall.graph)$name;
    expr.vec <<- expr.vec[!is.na(expr.vec)]
    overall.graph <- set.vertex.attribute(overall.graph, "abundance", index = V(overall.graph), value = expr.vals);

    hit.inx <- seed.proteins %in% node.list[,1];
    seed.proteins <<- seed.proteins[hit.inx];

    analSet <- DecomposeGraph(overall.graph,analSet);
    substats <- analSet$substats;
    overall.graph <<- overall.graph;
    analSet$overall.graph <- overall.graph;
    if(!is.null(substats)){
        output <- c(length(seed.genes), length(seed.proteins), nrow(node.list), nrow(edge.list), length(ppi.comps), substats);
    }else{
        output <- 0;
    }
       return( saveSet(analSet, "analSet", output) );

}

#' Filter bipartite network by degree and/or betweenness
#'
#' @export
#'
FilterBipartiNet <- function(nd.type, min.dgr, min.btw){
    paramSet <- readSet(paramSet, "paramSet");
    seed.genes <- paramSet$seed.genes;

    all.nms <- V(overall.graph)$name;
    edge.mat <- get.edgelist(overall.graph);
    dgrs <- degree(overall.graph);
    nodes2rm.dgr <- nodes2rm.btw <- NULL;

    if(nd.type == "gene"){
        hit.inx <- all.nms %in% edge.mat[,1];
    }else if(nd.type=="other"){
        hit.inx <- all.nms %in% edge.mat[,2];
    }else{ # all
        hit.inx <- rep(TRUE, length(all.nms));
    }

    if(min.dgr > 0){
        rm.inx <- dgrs <= min.dgr & hit.inx;
        nodes2rm.dgr <- V(overall.graph)$name[rm.inx];
    }
    if(min.btw > 0){
        btws <- betweenness(overall.graph);
        rm.inx <- btws <= min.btw & hit.inx;
        nodes2rm.btw <- V(overall.graph)$name[rm.inx];
    }

    nodes2rm <- unique(c(nodes2rm.dgr, nodes2rm.btw));
    overall.graph <- simplify(delete.vertices(overall.graph, nodes2rm));
    current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
    analSet <- DecomposeGraph(overall.graph,analSet);
    substats <- analSet$substats;
    if(!is.null(substats)){
        overall.graph <<- overall.graph;
        output <- c(length(seed.genes),length(seed.proteins), vcount(overall.graph), ecount(overall.graph), length(ppi.comps), substats);
    }else{
        output <- 0;
    }

    analSet$overall.graph <- overall.graph;
    return(saveSet(analSet, "analSet", output));
}

#' Prepare the json file for network visualization
#'
#' @param net.nm Name of subnetwork created from network building (i.e subnetwork1)
#' @param json.nm Name of json file to be exported for network visualization
#'
#' @export
#'
PrepareNetwork <- function(net.nm, json.nm){
   analSet <- readSet(analSet, "analSet");
   my.ppi <- analSet$ppi.comps[[net.nm]];
   nd.nms <- V(my.ppi)$name;


   analSet <- convertIgraph2JSON(net.nm, json.nm);
   current.net.nm <<- net.nm;
   return(saveSet(analSet, "analSet", 1));

}

#' Prepare the json file for network visualization from user uploaded graph file
#'
#' @param net.nm Name of subnetwork created from network building (i.e subnetwork1)
#' @param json.nm Name of json file to be exported for network visualization
#'
#' @export
#'
PrepareNetworkUpload <- function(net.nm, json.nm){

   my.ppi <- ppi.comps[[net.nm]];
   nd.nms <- V(my.ppi)$name;
   expr.vec <<- rep(0, length(nd.nms))
   table.nmu <<- "NA"
   convertIgraph2JSONFromFile(net.nm, json.nm, data.idType);
   current.net.nm <<- net.nm;
   return(saveSet(analSet, "analSet", 1));

}

GetNodeIDs <- function(){
    V(overall.graph)$name;
}

GetNodeNames <- function(){
    V(overall.graph)$Label;
}

GetNodeDegrees <- function(){
    degree(overall.graph);
}

GetNodeBetweenness <- function(){
    round(betweenness(overall.graph, directed=F, normalized=F), 2);
}

#' Compute minimum connect subnetwork based on input nodes using prize-collecting steiner forest approach
#'
#' @return check overall.graph object for result
#' @export
#'
ComputePCSFNet <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  seed.genes <- paramSet$seed.genes;
  edg <- as.data.frame(get.edgelist(overall.graph));
  edg$V3 <- rep(1, nrow(edg));
  colnames(edg) <- c("from", "to", "cost");

  node_names <- unique(c(as.character(edg[,1]),as.character(edg[,2])))
  ppi <- graph.data.frame(edg[,1:2],vertices=node_names,directed=F)
  E(ppi)$weight <- as.numeric(edg[,3])
  ppi <- simplify(ppi)

  if(sum(expr.vec) == 0){ # make sure weights are not 0?!
     expr.vec = expr.vec +1
  }
  g <- Compute.SteinerForest(ppi, expr.vec, w = 5, b = 100, mu = 0.0005);

  nodeList <- get.data.frame(g, "vertices");
  colnames(nodeList) <- c("Id", "Label");
  write.csv(nodeList, file="orig_node_list.csv", row.names=F, quote=F);
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  write.csv(edgeList, file="orig_edge_list.csv", row.names=F, quote=F);
  
  path.list <- NULL;
  analSet <- DecomposeGraph(g,analSet);
  substats <- analSet$substats;
  if(!is.null(substats)){
    overall.graph <<- g;
    current.msg<<- "Steiner Forest was completed successfully";
    output <- c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats);
  }else{
    output <- 0;
  }
  analSet$overall.graph <- overall.graph;
  return(saveSet(analSet, "analSet", output));
}

# Adapted from PCSF
# https://github.com/IOR-Bioinformatics/PCSF
Compute.SteinerForest <- function(ppi, terminals, w = 2, b = 1, mu = 0.0005, dummies){

  # Gather the terminal genes to be analyzed, and their scores
  terminal_names <- names(terminals)
  terminal_values <- as.numeric(terminals)

  # Incorporate the node prizes
  node_names <- V(ppi)$name
  node_prz <- vector(mode = "numeric", length = length(node_names))
  index <- match(terminal_names, node_names)
  percent <- signif((length(index) - sum(is.na(index)))/length(index)*100, 4)
  if (percent < 5){
    print("Less than 1% of your terminal nodes are matched in the interactome!");
    return(NULL);
  }
  paste0("  ", percent, "% of your terminal nodes are included in the interactome\n");
  terminal_names <- terminal_names[!is.na(index)]
  terminal_values <- terminal_values[!is.na(index)]
  index <- index[!is.na(index)]
  node_prz[index] <- terminal_values

  if(missing(dummies)||is.null(dummies)||is.na(dummies)){
    dummies <- terminal_names #re-assign this to allow for input
  }

  ## Prepare input file for MST-PCSF implementation in C++

  # Calculate the hub penalization scores
  node_degrees <- igraph::degree(ppi)
  hub_penalization <- - mu*node_degrees

  # Update the node prizes
  node_prizes <- b*node_prz
  index <- which(node_prizes==0)
  node_prizes[index] <- hub_penalization[index]

  # Construct the list of edges
  edges <- ends(ppi,es = E(ppi))
  from <- c(rep("DUMMY", length(dummies)), edges[,1])
  to <- c(dummies, edges[,2])

  cost <- c(rep(w, length(dummies)), E(ppi)$weight)

  #PCSF will faill if there are NAs in weights, this will check and fail gracefully
  if(any(is.na(E(ppi)$weight))){
    print("NAs found in the weight vector!");
    return (NULL);
  }

  ## Feed the input into the PCSF algorithm
  output <- XiaLabCppLib::call_sr(from,to,cost,node_names,node_prizes)

  # Check the size of output subnetwork and print a warning if it is 0
  if(length(output[[1]]) != 0){

    # Contruct an igraph object from the MST-PCSF output
    e <- data.frame(output[[1]], output[[2]], output[[3]])
    e <- e[which(e[,2]!="DUMMY"), ]
    names(e) <- c("from", "to", "weight")

    # Differentiate the type of nodes
    type <- rep("Steiner", length(output[[4]]))
    index <- match(terminal_names, output[[4]])
    index <- index[!is.na(index)]
    type[index] <- "Terminal"

    v <- data.frame(output[[4]], output[[5]], type)
    names(v) <- c("terminals", "prize", "type")
    subnet <- graph.data.frame(e,vertices=v,directed=F)
    E(subnet)$weight <- as.numeric(output[[3]])
    subnet <- delete_vertices(subnet, "DUMMY")
    subnet <- delete_vertices(subnet, names(which(degree(subnet)==0)));
    return(subnet)

  } else{
    print("Subnetwork can not be identified for a given parameter set")
    return(NULL);
  }
}


# for a given graph, obtain the smallest subgraphs that contain
# all the seed nodes. This is acheived by iteratively remove
# the marginal nodes (degree = 1) that are not in the seeds

#' Compute minimum connected network composed of seed nodes using shortest path based approach
#'
#' @param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#' @param max.len Maximum number of seeds, if more than this number, take top 100 seed nodes based on degrees
#'
#' @return check overall.graph object for result
#' @export
#'
GetMinConnectedGraphs <- function(max.len = 200){
    set.seed(8574);
    paramSet <- readSet(paramSet, "paramSet");
    seed.genes <- paramSet$seed.genes;
    # first get shortest paths for all pair-wise seeds
    my.seeds <- seed.proteins;
    sd.len <- length(my.seeds);
    paths.list <-list();

    # first trim overall.graph to remove no-seed nodes of degree 1
    dgrs <- degree(overall.graph);
    keep.inx <- dgrs > 1 | (names(dgrs) %in% my.seeds);
    nodes2rm <- V(overall.graph)$name[!keep.inx];
    overall.graph <-  simplify(delete.vertices(overall.graph, nodes2rm));

    # need to restrict the operation b/c get.shortest.paths is very time consuming
    # for top max.len highest degrees
    if(sd.len > max.len){
        hit.inx <- names(dgrs) %in% my.seeds;
        sd.dgrs <- dgrs[hit.inx];
        sd.dgrs <- rev(sort(sd.dgrs));
        # need to synchronize all (seed.proteins) and top seeds (my.seeds)
        seed.proteins <- names(sd.dgrs);
        if(max.len>table(hit.inx)[["TRUE"]]){
            sd.len <-  table(hit.inx)[["TRUE"]];
        }else{
            sd.len <-  max.len;
        }
        my.seeds <- seed.proteins[1:sd.len];
        current.msg <<- paste("The minimum connected network was computed using the top", sd.len, "seed proteins in the network based on their degrees.");
    }else{
        current.msg <<- paste("The minimum connected network was computed using all seed proteins in the network.");
    }
    # now calculate the shortest paths for
    # each seed vs. all other seeds (note, to remove pairs already calculated previously)
    for(pos in 1:sd.len){
            paths.list[[pos]] <- get.shortest.paths(overall.graph, my.seeds[pos], seed.proteins[-(1:pos)])$vpath;
    }
    nds.inxs <- unique(unlist(paths.list));
    nodes2rm <- V(overall.graph)$name[-nds.inxs];
    g <- simplify(delete.vertices(overall.graph, nodes2rm));

    nodeList <- get.data.frame(g, "vertices");
    colnames(nodeList) <- c("Id", "Label");
    fast.write(nodeList, file="orig_node_list.csv", row.names=F);

    edgeList <- get.data.frame(g, "edges");
    edgeList <- cbind(rownames(edgeList), edgeList);
    colnames(edgeList) <- c("Id", "Source", "Target");
    fast.write(edgeList, file="orig_edge_list.csv", row.names=F);

    path.list <- NULL;
    analSet <- DecomposeGraph(g,analSet);
    substats <- analSet$substats;
    net.stats<<-net.stats
    if(!is.null(substats)){
        overall.graph <<- g;
        output <- c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats);
    }else{
        output <- 0;
    }
    analSet$overall.graph <- overall.graph;
    return(saveSet(analSet, "analSet", output));
}

UpdateSubnetStats <- function(){
    old.nms <- names(ppi.comps);
    net.stats <- ComputeSubnetStats(ppi.comps);
    ord.inx <- order(net.stats[,1], decreasing=TRUE);
    net.stats <- net.stats[ord.inx,];
    rownames(net.stats) <- old.nms[ord.inx];
    net.stats <<- net.stats;
}

# exclude nodes in current.net (networkview)
ExcludeNodes <- function(nodeids, filenm){

    nodes2rm <- strsplit(nodeids, ";")[[1]];
    current.net <- ppi.comps[[current.net.nm]];
    current.net <- delete.vertices(current.net, nodes2rm);

    # need to remove all orphan nodes
    bad.vs<-V(current.net)$name[degree(current.net) == 0];
    current.net <- delete.vertices(current.net, bad.vs);

    # return all those nodes that are removed
    nds2rm <- paste(c(bad.vs, nodes2rm), collapse="||");

    # update topo measures
    node.btw <- as.numeric(betweenness(current.net));
    node.dgr <- as.numeric(degree(current.net));
    node.exp <- as.numeric(get.vertex.attribute(current.net, name="abundance", index = V(current.net)));
    nms <- V(current.net)$name;
    hit.inx <- match(nms, ppi.net$node.data[,1]);
    lbls <- ppi.net$node.data[hit.inx,2];

    nodes <- vector(mode="list");
    for(i in 1:length(nms)){
        nodes[[i]] <- list(
                  id=nms[i],
                  label=lbls[i],
                  degree=node.dgr[i],
                  between=node.btw[i],
                  expr = node.exp[i]
                );
    }
    # now only save the node pos to json
    netData <- list(deletes=nds2rm,nodes=nodes);
    sink(filenm);
    cat(RJSONIO::toJSON(netData));
    sink();

    ppi.comps[[current.net.nm]] <<- current.net;
    UpdateSubnetStats();

    # remember to forget the cached layout, and restart caching, as this is now different object (with the same name)
    #forget(PerformLayOut_mem);
    return(filenm);
}

# exclude nodes in overall net (network builder)
ExcludeNodesOverall <- function(nodeids, id.type, vismode){
    paramSet <- readSet(paramSet, "paramSet");
    seed.genes <- paramSet$seed.genes;
    data.org <- paramSet$data.org;
    # all convert to uniprot ID
    lines <- strsplit(nodeids, "\r|\n|\r\n")[[1]];
    lines<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", lines, perl=TRUE);
    if(vismode != "network"){
        prot.anots <- .doGeneIDMapping(lines, id.type, data.org, "matrix");
        nodes2rm <- unique(prot.anots$accession);
    }else{
        prot.anots <- lines
        nodes2rm <- unique(lines);
    }

    # now find the overlap
    nodes2rm <- nodes2rm[nodes2rm %in% V(overall.graph)$name];
    g <- delete.vertices(overall.graph, nodes2rm);

    nodeList <- get.data.frame(g, "vertices");
    colnames(nodeList) <- c("Id", "Label");

    fast.write(nodeList, file="orig_node_list.csv", row.names=F);

    edgeList <- get.data.frame(g, "edges");
    edgeList <- cbind(rownames(edgeList), edgeList);
    colnames(edgeList) <- c("Id", "Source", "Target");
    fast.write(edgeList, file="orig_edge_list.csv", row.names=F);

    analSet <- DecomposeGraph(g,analSet);
    substats <- analSet$substats;
    if(!is.null(substats)){
        overall.graph <<- g;
        output <- c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats);
    }else{
        output <- 0;
    }
    analSet$overall.graph <- g;
    return(saveSet(analSet,"analSet", output));
}

PrepareSubnetDownloads <- function(nm){
  g <- ppi.comps[[nm]];
  # need to update graph so that id is compound names rather than ID
  V(g)$name <- as.character(doID2LabelMapping(V(g)$name));
  saveNetworkInSIF(g, nm);
}

# adapted from BioNet
saveNetworkInSIF <- function(network, name){
    edges <- .graph.sif(network=network, file=name);
    sif.nm <- paste(name, ".sif", sep="");
    if(length(list.edge.attributes(network))!=0){
	edge.nms <- .graph.eda(network=network, file=name, edgelist.names=edges);
        sif.nm <- c(sif.nm, edge.nms);

    }
    if(length(list.vertex.attributes(network))!=0){
	node.nms <- .graph.noa(network=network, file=name);
        sif.nm <- c(sif.nm, node.nms);
    }
    # need to save all sif and associated attribute files into a zip file for download
    zip(paste(name,"_sif",".zip", sep=""), sif.nm);
}

# internal function to write cytoscape .sif file
.graph.sif <- function(network, file){
    edgelist.names <- igraph::get.edgelist(network, names=TRUE)
    edgelist.names <- cbind(edgelist.names[,1], rep("pp", length(E(network))), edgelist.names[,2]);
    write.table(edgelist.names, row.names=FALSE, col.names=FALSE, file=paste(file, ".sif", sep=""), sep="\t", quote=FALSE)
    return(edgelist.names)
}

doID2LabelMapping <- function(entrez.vec){
    if(exists("nodeListu")){
        hit.inx <- match(entrez.vec, nodeListu[, "Id"]);
        symbols <- nodeListu[hit.inx, "Label"];

        # if not gene symbol, use id by itself
        na.inx <- is.na(symbols);
        symbols[na.inx] <- entrez.vec[na.inx];
        return(symbols);
    }else{ # network upload
        return(entrez.vec);
    }
}

# internal method to write cytoscape node attribute files
.graph.noa <- function(network, file){
  all.nms <- c();
  attrib <- list.vertex.attributes(network)
  for(i in 1:length(attrib)){
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    noa <- cbind(V(network)$name, rep("=", length(V(network))), get.vertex.attribute(network, attrib[i]))
    first.line <- paste(attrib[i], " (class=java.lang.", type, ")", sep="")
    file.nm <- paste(file, "_", attrib[i], ".NA", sep="");
    write(first.line, file=file.nm, ncolumns = 1, append=FALSE, sep=" ")
    write.table(noa, row.names = FALSE, col.names = FALSE, file=file.nm, sep=" ", append=TRUE, quote=FALSE);
    all.nms <- c(all.nms, file.nm);
  }
  return(all.nms);
}

# internal method to write cytoscape edge attribute files
.graph.eda <- function(network, file, edgelist.names){
  all.nms <- c();
  attrib <- list.edge.attributes(network)
  for(i in 1:length(attrib)){
    if(is(get.edge.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(is(get.edge.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(is(get.edge.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    eda <- cbind(cbind(edgelist.names[,1], rep("(pp)", length(E(network))), edgelist.names[,3]), rep("=", length(E(network))), get.edge.attribute(network, attrib[i]))
    first.line <- paste(attrib[i], " (class=java.lang.", type, ")", sep="");
    file.nm <- paste(file, "_", attrib[i], ".EA", sep="");
    write(first.line, file=file.nm, ncolumns=1, append=FALSE, sep =" ")
    write.table(eda, row.names = FALSE, col.names = FALSE, file=file.nm, sep=" ", append=TRUE, quote=FALSE);
    all.nms <- c(all.nms, file.nm);
  }
  return(all.nms);
}

SetCellCoexNumber <-function(num){
    cellCoexNumber <<- num;
}

SetDiffNetName <-function(nm){
    diffNetName <<- nm;
}

SetDiffFilter <-function(pct){
    diffPct <<- pct/10;
}

# note: hit.query, resTable must synchronize

#' Perform gene enrichment analysis or identify gene regulatory targets
#'
#' @param file.nm File name of result table to be exported in csv format, do not include file extension
#' @param fun.type Enrichment database type
#' @param IDs String of ids to be tested separated by "; "
#'
#' @export
#'
PerformNetEnrichment <- function(file.nm, fun.type, IDs){
    paramSet <- readSet(paramSet, "paramSet");
    data.org <- paramSet$data.org;
    # prepare query
    ora.vec <- NULL;
    if(ppi.net$db.type == 'ppi'){
        if(data.org == "ath"){
                idtype <- "entrez"
        }else if(data.org %in% c("bsu", "tbr", "cel", "dme", "eco", "pfa") & net.type == "string"){
                idtype <- "string"
        }else if(data.org %in% c("bta","dre","rno","gga") & net.type == "string"){
                idtype <- "entrez";
        }else if(data.org %in% c("hsa","mmu", "pae") & net.type %in% c("string","innate", "huri", "interactome")){
                idtype <- "entrez";
        }else if(data.org %in% c("hsa","mmu", "cel", "dme","sce") & net.type %in% c("irefinx", "rolland")){
                idtype <- "uniprot";
        }else if(data.org == "sce" & net.type == "string"){ # only for yeast
                idtype <- "embl_gene";
        }
        if(idtype=="uniprot"){
            uniprot.vec <- unlist(strsplit(IDs, "; "));
            ora.vec <- doUniprot2EntrezMapping(uniprot.vec);
            names(ora.vec) <- uniprot.vec;
        }else if(idtype=="embl_protein"){
            emblprotein.vec <- unlist(strsplit(IDs, "; "))
            ora.vec <- doEmblProtein2EntrezMapping(emblprotein.vec);
            names(ora.vec) <- emblprotein.vec;
        }else if(idtype=="string"){
            string.vec <- unlist(strsplit(IDs, "; "))
            ora.vec <- doString2EntrezMapping(string.vec);
            names(ora.vec) <- string.vec;
        }else if(idtype=="embl_gene"){
            emblgene.vec <- unlist(strsplit(IDs, "; "))
            ora.vec <- doEmblGene2EntrezMapping(emblgene.vec);
            names(ora.vec) <- emblgene.vec;
        }else{
            ora.vec <- unlist(strsplit(IDs, "; "));
            names(ora.vec) <- ora.vec;
        }

        #if nothing matches, try original ids;
        if(all(is.na(ora.vec))){
            ora.vec <- unlist(strsplit(IDs, "; "));
            names(ora.vec) <- ora.vec;
        }

    }else{ # net is tf/mir/drug, they already in entrez
        ora.vec <- unlist(strsplit(IDs, "; "));
        names(ora.vec) <- as.character(ora.vec);
    }

    if(fun.type %in% c("trrust", "encode", "jaspar", "mirnet", "met", "drugbank", "disease")){
        res <- PerformRegEnrichAnalysis(file.nm, fun.type, ora.vec, "inverse");
    }else{
        res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, ora.vec);
    }

    if(checkEntrezMatches(unname(ora.vec)) > 0 && res == 0){
        res=2;
    }
    return(res);
}


PlotCorHistogram <- function(imgNm){
    Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
    library(ggplot2)
    cor <- E(overall.graph)$correlation
    p<-ggplot(data.frame(Correlation=abs(cor)), aes(x=Correlation)) + geom_histogram() + ggtitle("Correlation distribution")+
    theme(plot.title = element_text(hjust = 0.5));
    print(p);
    dev.off();
}

#' Filter interactions by tissue coexpression data
#'
#' @param type Name of database (encode, gtex, hpa)
#' @param tissue Selected tissue to filter (refer to web server interface for options, replace space with "." (i.e. Adipose.Tissue))
#'
#' @export
#'
FilterByTissue <- function(type, tissue){
    paramSet <- readSet(paramSet, "paramSet");
    seed.genes <- paramSet$seed.genes;
    all.nms <- V(overall.graph)$name;
    gene.nms <- all.nms[which(V(overall.graph)$Types %in% c("Gene", "Protein"))]

    tissue.mat <- queryFilterDB(type, data.org);
    hit.inx <- tissue.mat$Tissue %in% c(tissue, "All");

    tissue.genes <- tissue.mat[,1][!hit.inx];
    nodes2rm <- gene.nms[which(gene.nms %in% tissue.genes)]

    g <- simplify(delete.vertices(overall.graph, nodes2rm));

    nodeList <- get.data.frame(g, "vertices");
    nodeList <- nodeList[,1:2];
    colnames(nodeList) <- c("Id", "Label");

    fast.write(nodeList, file="orig_node_list.csv");
    nd.inx <- ppi.net$node.data[,1] %in% nodeList[,1];

    edgeList <- get.data.frame(g, "edges");
    edgeList <- edgeList[,1:2];
    colnames(edgeList) <- c("Source", "Target");
    fast.write(edgeList, file="orig_edge_list.csv");

    # update ppi.net
    ppi.net$order <- 1;
    ppi.net$node.data <- nodeList;
    ppi.net$edge.data <- edgeList;
    ppi.net <<- ppi.net;
    overall.graph <- g
    analSet <- DecomposeGraph(overall.graph,analSet);
    substats <- analSet$substats;
    if(!is.null(substats)){
        overall.graph <<- overall.graph;
        output <- c(length(seed.genes),length(seed.proteins), vcount(overall.graph), ecount(overall.graph), length(ppi.comps), substats);
    }else{
        output <- 0;
    }
    analSet$overall.graph <- overall.graph;
    return(saveSet(analSet, "analSet", output));
}

#' Filter interactions by correlation (only applicable to coexpression network)
#'
#' @param min.cor Correlation threshold (0.5 - 1)
#'
#' @export
#'
FilterNetByCor <- function( min.cor){
    paramSet <- readSet(paramSet, "paramSet");
    seed.genes <- paramSet$seed.genes;
    all.nms <- V(overall.graph)$name;
    edge.mat <- get.edgelist(overall.graph);
    cor <- E(overall.graph)$correlation;
    edges2rm.inx <-NULL;

    if(min.cor > 0.5){
        rm.inx <- cor <= min.cor;
    }

    overall.graph <- simplify(delete.edges(overall.graph, which(rm.inx)));
    current.msg <<- paste("A total of", length(which(rm.inx)) , "edges was reduced.");
    analSet <- DecomposeGraph(overall.graph,analSet);
    substats <- analSet$substats;
    if(!is.null(substats)){
        overall.graph <<- overall.graph;
        return(c(length(seed.genes),length(seed.proteins), vcount(overall.graph), ecount(overall.graph), length(ppi.comps), substats));
    }else{
        return(0);
    }
}

regEnrichment <- function(file.nm, fun.type, IDs, netInv){
  regBool <<- "false";
  res <- PerformNetEnrichment(file.nm, fun.type, IDs);
}

PerformRegEnrichAnalysis <- function(file.nm, fun.type, ora.vec, netInv){
    if(!exists("my.reg.enrich")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/networkanalystr/_utils_regenrich.Rc");    
    }
    return(my.reg.enrich(file.nm, fun.type, ora.vec, netInv));
}

PrepareTheraNet <- function(net.nm, json.nm, theraids, regids, regnms){
  require("igraph")
  net.nm = gsub(".json", "", net.nm);
  
  if(is.null(net.nm)){
    net.nm <- names(ppi.comps)[1];
  }
  
  reg.count <- reg.count + 1;
  reg.nm <- paste("addedReg", reg.count, sep="");
  
  my.ppi <- ppi.comps[[net.nm]];
  net.nm <<- net.nm;
  json.nm <<- json.nm;
  theraids <<- theraids
  regids<<-regids;

  theraids <- strsplit(theraids, ";");
  theraids <- theraids[[1]];
  theraids <- strsplit(theraids, ",");
  
  regids <- strsplit(regids, ",");
  regids <- regids[[1]]
  regnms <- strsplit(regnms, ",");
  regnms <- regnms[[1]]
  regnms<<-regnms;
  regids<<-regids;
  my.ppi <- ppi.comps[[net.nm]];
               
  nms <- V(my.ppi)$name                           
  remove.inx <-regids %in% nms
  nodesToAdd <- regids[!remove.inx];
    if(length(nodesToAdd)>0){
      my.ppi <- my.ppi + vertices(nodesToAdd);
    }
  for(i in 1:length(regids)){
    for(j in 1:length(theraids[[i]])){
      my.ppi <- my.ppi + edge(regids[[i]], theraids[[i]][j]);
    }
  }
  if(!is.null(V(my.ppi)$moltype)){
    moltype.vec <- V(my.ppi)$moltype
  }else{
    moltype.vec <- rep("NA", length(V(my.ppi)))
  }
    nms <- V(my.ppi)$name
    reg.inx <- nms %in% regids
    moltype.vec[reg.inx] <- "regulator"
    my.ppi <- set_vertex_attr(my.ppi, "moltype", value = moltype.vec)

  edge.data <- data.frame(regids, regnms, rep("regulator", length(regids)) ,stringsAsFactors=FALSE);
  
  colnames(edge.data) <- c("Id", "Label", "Types");
  colnames(ppi.net$node.data) <- c("Id", "Label", "Types");
  ppi.net$node.data <<-rbind(ppi.net$node.data, edge.data);
  filenm <- paste(reg.nm, ".json", sep="");
  reg.count <<- reg.count

  my.ppi <- simplify(my.ppi);
  ppi.comps[[reg.nm]] <<- my.ppi;
  current.net.nm <<- reg.nm;
  ppi.comps <<- ppi.comps;
  UpdateSubnetStats();
  #net.info$reg.ids <<- regids;
  convertIgraph2JSON(reg.nm, filenm);

  return(filenm);
  
}

GetSelectedNetworks <- function(){
    paramSet <- readSet(paramSet, "paramSet");
    return(names(paramSet$nets));
}

#' Merge or identify shared portions of multiple networks
#'
#' @param type Type of integrative network, either intersect, union or unique
#'
#' @export
#'
ComputeIntegNet <- function(type){
    gObjs <- list();
    paramSet <- readSet(paramSet, "paramSet");
    data.org <- paramSet$data.org;
    seed.expr <- paramSet$seed.expr;
    seed.genes <- paramSet$seed.genes;

    for(i in 1:length(paramSet$nets)){
        ppi.net <<- paramSet$nets[[i]];
        CreateGraph();
        gObj <- overall.graph;
        gObjs[[i]] <- gObj;
    }

    if(type == "intersect"){
       g <- do.call(igraph::intersection, gObjs)
    }else if(type == "union"){
       g <- do.call(igraph::union, gObjs)
    }else{
        gObjs.union <- list()
        gObj.unique <- "";
        for(i in 1:length(paramSet$nets)){
            gObj <- gObjs[[i]];
            nm <- names(paramSet$nets)[i]
            if(nm != type){
                gObjs.union[[nm]] <- gObj
            }else{
                gObj.unique <- gObj
            }
        }
        union.graph <- do.call(igraph::union, gObjs.union)
        g <- igraph::difference(gObj.unique, union.graph)  
    }

    nms = V(g)$name
    labels = rep("NA", length(V(g)$name))
    types = rep("NA", length(V(g)$name))

    for(i in 1:length(paramSet$nets)){
        ppi.net <- paramSet$nets[[i]];
        node.ids <- ppi.net$node.data[,1]
        node.nms <- ppi.net$node.data[,2]
        node.types <- ppi.net$node.data[,3]
        hit.inx <- nms %in% node.ids;
        match.inx <- which(node.ids %in% nms)
        labels[hit.inx] <- node.nms[match.inx]
        types[hit.inx] <- node.types[match.inx]
    }
    g <- set.vertex.attribute(g, "Label", index = V(g), value = labels);
    g <- set.vertex.attribute(g, "Types", index = V(g), value = types);

    node.res <- data.frame(Id=V(g)$name, Label=V(g)$Label, Types=V(g)$Types)
    edge.mat <- as_edgelist(g);
    edge.res <- data.frame(Source=edge.mat[,1], Target=edge.mat[,2])
    
    ppi.net$db.type <- "multi"
    ppi.net$node.data <- node.res;
    ppi.net$edge.data <- edge.res;

    node.list <- ppi.net$node.data;
    edge.list <- ppi.net$edge.data;
        
    overall.graph <- g

 # add node expression value
    if(ppi.net$db.type == "ppi"){
        ## Temp fix seed.expr all in entrez?
        ## all converted to new IDs based on PPI from the given organism
        newIDs <- doPpiIDMapping(sqlite.path, names(seed.expr), data.org);
    }else{# all entrez in mirNet
        newIDs <- names(seed.expr);
    }
    match.index <- match(V(overall.graph)$name, newIDs);
    expr.vals <- seed.expr[match.index];
    expr.vec <- expr.vals;
    names(expr.vec)<- V(overall.graph)$name;
    expr.vec <<- expr.vec[!is.na(expr.vec)]
    overall.graph <- set.vertex.attribute(overall.graph, "abundance", index = V(overall.graph), value = expr.vals);

    hit.inx <- seed.proteins %in% node.list[,1];
    seed.proteins <<- seed.proteins[hit.inx];

    analSet <- DecomposeGraph(overall.graph,analSet);
    substats <- analSet$substats;
    overall.graph <<- overall.graph;
    ppi.net <<- ppi.net
    if(!is.null(substats)){
        output <- c(length(seed.genes), length(seed.proteins), nrow(node.list), nrow(edge.list), length(ppi.comps), substats);
    }else{
        output <- 0;
    } 
    analSet$overall.graph <- overall.graph;
    analSet$ppi.net <- ppi.net;
    return(saveSet(analSet,"analSet", output));
}

ClearNetSelection <- function(type){
    paramSet <- readSet(paramSet, "paramSet");
    if(type %in% names(paramSet$nets)){
        inx = which( names(paramSet$nets) == type)
         paramSet$nets[[inx]] = NULL;
    }
    saveSet(paramSet, "paramSet");
    return(1)
}

queryFilterDB <- function(type, org){
  paramSet <- readSet(paramSet, "paramSet");
  sqlite.path <- paramSet$sqlite.path;
  require('RSQLite');
  conv.db <- dbConnect(SQLite(), paste(sqlite.path, "tissue_filter.sqlite", sep="")); 
  db.map <- dbReadTable(conv.db, paste0(data.org,"_",type));
  dbDisconnect(conv.db); cleanMem();

  return(db.map)
}

