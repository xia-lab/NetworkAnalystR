##################################################
## R script for NetworkAnalyst
## Description: Graph IO functions for network upload module
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################


#'Read graph file
#'@description read graph file and store into graph object
#'@param fileName file name of the graph file
#'@param fileType file type of the graph file (graphml, sif, txt, json)
#'@param org (optional, can be NA) organism
#'@param idType (optional, can be NA) identifier type
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'

ReadGraphFile <- function(fileName, fileType, org, idType) {

  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");

  json3d <- F;
  data.org<<- org;
  paramSet$data.org <- data.org;

  require("igraph");
  types_arr <<- "";
  ppi.comps <<- list();
  
  current.msg <<- NULL;
  
  if(fileType == "graphml"){
    graphX = tryCatch({
      read_graph(fileName, format = "graphml")
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {
      
    })
  }else if(fileType == "sif"){
    graphX = tryCatch({
      read.sif(fileName, format="igraph", directed = FALSE, sep="\t")
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {
      
    })
  }else if(fileType == "txt" || endsWith(fileName, "csv")){

    raw <- .readDataTable(fileName);
    node.types <- vector();
    if(ncol(raw)>2){
        df.list <- list();
        edge.list <- data.frame()
        for(i in 1:(ncol(raw)-1)){
            df.list[[i]] <- raw[,c(i,i+1)];
            colnames(df.list[[i]]) = c("Source", "Target");
        }
        for(i in 1:length(df.list)){
           edge.list <- rbind(edge.list, df.list[[i]])
        }
        df <- edge.list;
        df <- as.matrix(df)
    }else{
        df <- as.matrix(raw)
    }

    graphX = tryCatch({
      graph_from_edgelist(df)
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {
      
    })

    V(graphX)$moltype = V(graphX)$name
    if(colnames(raw)[1] != "Source"){
        nms = V(graphX)$name;
        for(i in 1:length(colnames(raw))){
            inx = which(nms %in% raw[,i])
            V(graphX)$moltype[inx] = colnames(raw[i])
        }
    }

  }else if(fileType == "json"){
    # rjson bug use RJSONIO
    dat = RJSONIO::fromJSON(fileName);
    dfn = unlist(dat$elements$nodes);
    conv = cbind_dif(list(id1=dfn[which(names(dfn)=='data.id')], name1=dfn[which(names(dfn)=='data.name')]));
    dfe = unlist(dat$elements$edges);
    dffe = data.frame(id1=dfe[which(names(dfe) == "data.source")], id2=dfe[which(names(dfe) == "data.target")]);
    dfint = merge(conv, dffe, by="id1");
    colnames(conv) = c("id2", "name2");
    df = merge(conv, dfint, by="id2");
    df = df[,c("id1", "id2")];
    df=as.matrix(df)
    
    
    graphX = tryCatch({
      graph_from_edgelist(df, directed=FALSE)
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {
      
    })
  }else if(fileType == "netjson"){
    
    # rjson bug use RJSONIO
    dat <- RJSONIO::fromJSON(fileName);
    
    org <- dat$org
    data.org <<- org
    if(!is.null(dat$threed)){
    json3d = T;
    }
    idType <<- dat$idType
    dfn <- unlist(dat$nodes);
    conv <- cbind_dif(list(id1=dfn[which(names(dfn)=='id')], name1=dfn[which(names(dfn)=='label')]));
    dfe <- unlist(dat$edges);
    dffe <- data.frame(id1=dfe[which(names(dfe) == "source")], id2=dfe[which(names(dfe) == "target")]);
    dfint <- merge(conv, dffe, by="id1");
    colnames(conv) <- c("id2", "name2");
    df <- merge(conv, dfint, by="id2");
    df <- df[,c("id1", "id2")];
    df <- as.matrix(df)
    
    graphX <- tryCatch({
      graph_from_edgelist(df, directed=FALSE)
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {
      
    })
    
    sink("networkanalyst_0.json");
    cat(RJSONIO::toJSON(dat));
    sink();
    
  }else{
    current.msg <<- "Unknown format, please make sure that the file is saved in the supported formats!";
    return(0)
  }
  fileTypeu <<- fileType;
  
  if(!is_igraph(graphX)){
    current.msg <<- "Failed to parse your file, please make sure that the file is formatted correctly";
    return(0)
  }
  current.msg <<- "Successfully parsed your graph file!";
  if(idType == "NA"){
    nms <- V(graphX)$name;
    if(length(nms)<1){
      nms <- V(graphX)$id;
      graphX = set_vertex_attr(graphX, "name", value=nms)
    }
    node.data = data.frame(nms, nms);
    seed.proteins <<- nms;
  }else{
    nms <- V(graphX)$name;
    if(length(nms)<1){
      nms <- V(graphX)$id;
      graphX = set_vertex_attr(graphX, "name", value=nms)
    }
    entrezs = .doGeneIDMapping(nms, idType, data.org, "matrix")
    nms = entrezs[,"gene_id"]
    symbols = doEntrez2SymbolMapping(nms)
    node.data = data.frame(nms, symbols);
    seed.proteins <<- entrezs;
    V(graphX)$name = nms;
  }
  seed.genes <<- seed.proteins;
  e=get.edgelist(graphX)
  edge.data= data.frame(Source=e[,1], Target=e[,2])
  
  paramSet$seed.expr <- rep(0, length(node.data));
  node.data <- cbind(node.data, rep("NA", nrow(node.data)) )
  colnames(node.data) = c("Id", "Label", "Types");
  analSet <- DecomposeGraph(graphX,analSet, 3, fileType);
  substats <- analSet$substats;

  net.nm <- names(ppi.comps)[1];
  net.nmu <<- net.nm;
  current.net.nm <<- net.nm;
  g <- ppi.comps[[net.nm]];

  ppi.net <<- list(db.type="uploaded", 
                   order=1, 
                   seeds=nms, 
                   table.nm=" ", 
                   node.data = node.data,
                   edge.data = edge.data
  );
  data.idType <<- idType; 
  overall.graph <<- g;
  if(fileType != "netjson"){
    convertIgraph2JSONFromFile(net.nm, "networkanalyst_0.json", data.idType);
  }

  analSet$ppi.net <- ppi.net;
  analSet$overall.graph <- overall.graph;
  saveSet(paramSet, "paramSet");
  if(json3d){
    
    output <- 2;
  }else{
    output <- 1;
  }
  return(saveSet(analSet, "analSet", 1));
}


# create igraph from the edgelist saved from graph DB
# and decompose into subnets

convertIgraph2JSONFromFile <- function(net.nm, filenm, idType){
  paramSet <- readSet(paramSet, "paramSet");
  anal.type <- paramSet$anal.type;
  if(fileTypeu == "netjson"){
    return(0);
  }
  g <- ppi.comps[[net.nm]];
  
  # annotation
  nms <- V(g)$name;
  lbls <- ppi.net$node.data[which(ppi.net$node.data[,1] %in% nms) ,2]
  
  
  # setup shape (gene circle, other squares)
  shapes <- rep("circle", length(nms));
  
  # get edge data
  edge.mat <- get.edgelist(g);
  
  edge.mat1 = data.frame(edge.mat)
  edge.mat1$color = "target"
  edge.mat1 = as.matrix(edge.mat1)
  
  if(!is.null(E(g)$weight)){
    E(g)$correlation <- E(g)$weight
    E(g)$weight <- abs(E(g)$weight)
    ppi.comps[[net.nm]] <- g
    ppi.comps <<- ppi.comps
  }
  
  if(!is.null(E(g)$correlation)){
    edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], color = edge.mat1[,3], correlation=E(g)$correlation);
  }else{
    edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], color = edge.mat1[,3]);
  }
  
  # get the note data
  
  node.btw <- as.numeric(betweenness(g));
  #node.clo <- as.numeric(closeness(g));
  #node.adh <- as.numeric(adhesion(g));
  node.eig <- eigen_centrality(g);
  node.eig <- as.numeric(node.eig$vector);
  node.tra <- transitivity(g,type=c("local"))
  
  node.dgr <- as.numeric(degree(g));
  node.exp <- as.numeric(get.vertex.attribute(g, name="Expression", index = V(g)));
  
  if(length(node.exp) == 0){
    node.exp <- rep(0,length(node.dgr)); 
  }
  
  # node size to degree values
  if(vcount(g)>500){
    min.size <- 1;
  }else if(vcount(g)>200){
    min.size <- 2;
  }else{
    min.size <- 3;
  }
  
  minval <- min(node.dgr, na.rm=T);
  maxval <- max(node.dgr, na.rm=T);
  result <- maxval-minval;
  
  if(result == 0){
    node.sizes <- rep((log(node.dgr))^2, length(nms));
  }else{
    node.sizes <- as.numeric(rescale2NewRange((log(node.dgr))^2, min.size, 9));
  }
  
  centered <- T;
  notcentered <- F;
  
  # update node color based on betweenness
  require("RColorBrewer");
  topo.val <- log(node.btw+1);
  topo.colsb <- ComputeColorGradient(topo.val, "black", notcentered, FALSE);
  topo.colsw <-  ComputeColorGradient(topo.val, "white", notcentered, FALSE);
  topo.colsc <-  ComputeColorGradient(topo.val, "colorblind", notcentered, TRUE);
  
  # color based on expression
  bad.inx <- is.na(node.exp) | node.exp==0;
  if(!all(bad.inx)){
    exp.val <- node.exp;
    node.colsb.exp <- ComputeColorGradient(exp.val, "black", centered, FALSE); 
    node.colsw.exp <- ComputeColorGradient(exp.val, "white", centered, FALSE);
    node.colsc.exp <- ComputeColorGradient(exp.val, "colorblind", centered, TRUE);
    node.colsb.exp[bad.inx] <- "#d3d3d3"; 
    node.colsw.exp[bad.inx] <- "#c6c6c6"; 
    node.colsc.exp[bad.inx] <- "#99ddff";
  }else{
    node.colsb.exp <- rep("#d3d3d3",length(node.exp)); 
    node.colsw.exp <- rep("#c6c6c6",length(node.exp)); 
    node.colsc.exp <- rep("#99ddff",length(node.exp)); 
  }
  
  # now update for bipartite network
  gene.inx <- nms %in% edge.mat[,"source"];
  mir.inx <- nms %in% edge.mat[,"target"];
  node_attr = list.vertex.attributes(g);
  
  attr=list();
  for(j in 1:length(node_attr)){
    attr[[node_attr[j]]] = vertex_attr(g, node_attr[j])
  }
  attr_names <- names(attr);
  attr_nd <- list();
  arr <- list()
  for(i in 1:length(node.sizes)){
    for(j in 1:length(attr)){
      attr_nd[node_attr[j]] = as.character(unlist(attr[node_attr[j]])[i])
    }
    arr[[i]] = attr_nd;
  }
  network_prop = list();
  for(i in 1:length(node.sizes)){
    network_prop[[i]]  <- list(
      #Closeness = node.clo[i],
      Eigen = node.eig[i],
      Transitivity = node.tra[i]
    )
  }
  pos.xy <- PerformLayOut(net.nm, "Default");
  lblsu <<- nms;
  
  #assign type
  typeVec <- rep("circle",length(node.exp)); 
  if(!is.null(V(g)$type)){
    typeVec <- V(g)$type;
    ut <- unique(typeVec)
    typeVec[which(typeVec == ut[1])] <- "circle"
    if(!is.na(ut[2])){
      typeVec[which(typeVec == ut[2])] <- "square"
    }
    if(!is.na(ut[3])){
      typeVec[which(typeVec == ut[3])] <- "diamond"
    }
  }
  
  if(is.null(V(g)$moltype)){
    mol.types <- rep("Unknown",length(node.exp)); 
  }else{
    mol.types <- V(g)$moltype; 
  }
  
  color.vec <- .gg_color_hue(length(unique(mol.types)))
  
  # now create the json object
  nodes <- vector(mode="list");
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id=nms[i],
      idnb=i,
      x = pos.xy[i,1],
      y = pos.xy[i,2],
      label=lbls[i],
      size=node.sizes[i], 
      type=typeVec[i],
      molType = mol.types[i],
      colorb=topo.colsb[i],
      colorw=topo.colsw[i],
      colorc=topo.colsc[i],
      topocolb=topo.colsb[i],
      topocolw=topo.colsw[i],
      topocolc=topo.colsc[i],
      expcolb=node.colsb.exp[i],
      expcolw=node.colsw.exp[i],
      expcolc=node.colsc.exp[i],
      user=c(arr[[i]], network_prop[[i]]),
      attributes=list(
        expr = node.exp[i],
        degree=node.dgr[i], 
        between=node.btw[i]
      )
    );
  }
  
  # save node table
  nd.tbl <- data.frame(Id=nms, Degree=node.dgr, Betweenness=round(node.btw,2));
  # order 
  ord.inx <- order(nd.tbl[,2], nd.tbl[,3], decreasing = TRUE)
  nd.tbl <- nd.tbl[ord.inx, ];
  fast.write(nd.tbl, file="node_table.csv", row.names=FALSE);
  
  netUploadU <<-1
  if(length(V(g)$name)>100){
    modules <- FindCommunities("walktrap", FALSE);
  }else{
    modules <- "NA"
  }
  
  # covert to json
  #globalProperties <-list();
  #globalProperties[["Diameter"]] <-diameter(g)
  #globalProperties[["Radius"]] <-radius(g)
  #globalProperties[["Average path length"]] <-signif(mean_distance(g), 3);
  #globalProperties[["Clustering coefficient"]] <- signif(transitivity(g, type="global"), 3);

  netData <- list( nodes=nodes, edges=edge.mat, idType=idType, analType=anal.type, org=data.org, naviString = "network", modules=modules , tblNm="NA", nodeTypes= as.vector(unique(mol.types)), nodeColors= color.vec);
  
  if(!is.null(E(g)$correlation)){
    netData[["maxCorrelation"]] <- max(E(g)$correlation)
    netData[["minCorrelation"]] <- min(abs(E(g)$correlation))
  }
  paramSet$jsonNms$network <<- filenm
  if(!filenm %in% partialToBeSaved){
    partialToBeSaved <<- c(partialToBeSaved, c(filenm))
  }
  saveSet(paramSet,"paramSet");
  sink(filenm);
  cat(RJSONIO::toJSON(netData));
  sink();
}

read.sif <- function (sif.file, format = "graphNEL", directed = FALSE, header = TRUE, sep = "\t", ...) {
  
  net <- read.csv(file = sif.file, sep = sep, colClasses = "character", header = header, ...)
  
  # Assume form: node1 linktype node2 side.info..
  if ( ncol(net) > 2 ) { 
    nas <- apply(net, 1, function (x) {any(is.na(x[c(1,3)]))})
    if (any(nas)) {
      net <- net[!nas, ]
      warning("NAs removed from network node list, ", sum(nas), " edges removed.")
    }
    net <- graph.edgelist(as.matrix(net[, -2]), directed = directed);
    
  } else if ( ncol(net) == 2 ) { # assume form: node1 node2
    
    # remove NA nodes 
    nas <- apply(net, 1, function (x) {any(is.na(x))})
    if (any(nas)) {
      net <- net[!nas, ]
      warning("NAs removed from network node list, ", sum(nas), " edges removed.")
    }
    
    net <- graph.edgelist(cbind(net[,1],net[,2]), directed = directed)
  }
  
  if (format == "graphNEL") { net <- igraph.to.graphNEL(net) }
  
  return(net);
}

cbind_dif <- function(x = list()){ ## different length columns as a data.frame
    # Find max length
    max_length <- max(unlist(lapply(x, length)))

    # Set length of each vector as
    res <- lapply(x, function(x){
        length(x) <- max_length
        return(x)
    })

    return(as.data.frame(res))
}