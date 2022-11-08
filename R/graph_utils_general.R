##################################################
## R script for NetworkAnalyst
## Description: General graph manipulation functions 
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

GetColorSchema <- function(my.grps){
  # test if total group number is over 9
  my.grps <- as.factor(my.grps);
  grp.num <- length(levels(my.grps));
  
  if(grp.num > 9){
    pal12 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
               "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
               "#FFFF99", "#B15928");
    dist.cols <- colorRampPalette(pal12)(grp.num);
    lvs <- levels(my.grps);
    colors <- vector(mode="character", length=length(my.grps));
    for(i in 1:length(lvs)){
      colors[my.grps == lvs[i]] <- dist.cols[i];
    }
  }else{
    colors <- as.numeric(my.grps)+1;
  }
  return (colors);
}


DecomposeGraph <- function(gObj,analSet, minNodeNum = 3, jsonBool = F){
  # now decompose to individual connected subnetworks
    if(jsonBool == "netjson"){
        comps <-list(gObj)
    }else{
        comps <-decompose.graph(gObj, min.vertices=minNodeNum);
    }
  if(length(comps) == 0){
    current.msg <<- paste("No subnetwork was identified with at least", minNodeNum, "nodes!");
    return(NULL);
  }
  
  # first compute subnet stats
  net.stats <- ComputeSubnetStats(comps);
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  comps <- comps[ord.inx];
  names(comps) <- rownames(net.stats) <- paste("subnetwork", 1:length(comps), sep="");
  
  # note, we report stats for all nets (at least 3 nodes);
  hit.inx <- net.stats$Node >= minNodeNum;
  comps <- comps[hit.inx];
  #overall <- list();
  #overall[["overall"]] <- g
  #ppi.comps <- append(overall, ppi.comps, after=1);
  
  # now record
  ppi.comps <<- comps;
  net.stats <<- net.stats;
  sub.stats <- unlist(lapply(comps, vcount)); 
  analSet$ppi.comps <- comps;
  analSet$net.stats <- net.stats;
  analSet$substats <- sub.stats;
  return(analSet);
}


ComputeSubnetStats <- function(comps){
  net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
  colnames(net.stats) <- c("Node", "Edge", "Query");
  for(i in 1:length(comps)){
    g <- comps[[i]];
    net.stats[i,] <- c(vcount(g),ecount(g),sum(V(g)$name %in% seed.proteins));
    #print(net.stats[i, ]);
  }
  return(net.stats);
}


# new range [a, b]
rescale2NewRange <- function(qvec, a, b){
  q.min <- min(qvec);
  q.max <- max(qvec);
  if(length(qvec) < 50){
    a <- a*2;
  }
  if(q.max == q.min){
    new.vec <- rep(8, length(qvec));
  }else{
    coef.a <- (b-a)/(q.max-q.min);
    const.b <- b - coef.a*q.max;
    new.vec <- coef.a*qvec + const.b;
  }
  return(new.vec);
}


ComputeColorGradient <- function(nd.vec, background="black", centered, colorblind){
  require("RColorBrewer");
  
  minval <- min(nd.vec, na.rm=TRUE);
  maxval <- max(nd.vec, na.rm=TRUE);
  res <- maxval-minval;
  
  if(res == 0){
    return(rep("#FF0000", length(nd.vec)));
  }
  color <- GetColorGradient(background, centered, colorblind);
  breaks <- generate_breaks(nd.vec, length(color), center = centered);
  return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
}

generate_breaks <- function(x, n, center = F){
  if(center){
    m <- max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res <- seq(-m, m, length.out = n + 1)
  }
  else{
    res <- seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
  }
  return(res)
}

scale_vec_colours <- function(x, col = rainbow(10), breaks = NA){
  breaks <- sort(unique(breaks));
  return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

GetColorGradient <- function(background, center, colorblind=F) {
  if (background == "black") {
    if (center) {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#6baed6", "#bdd7e7", "#eff3ff"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(c(colorRampPalette(c("#31A231", "#5BC85B", "#90EE90", "#C1FFC1"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      }
    } else {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#bbfdff"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(colorRampPalette(rev(heat.colors(9)))(100))
      }
    }
  } else {
    if (center) {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#bbfdff"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(c(colorRampPalette(c("#137B13", "#31A231", "#5BC85B", "#90EE90"))(50), colorRampPalette(c("#FF7783", "#E32636", "#BD0313", "#96000D"))(50)))
      }
    } else {
      return(colorRampPalette(hsv(h = seq(0.72, 1, 0.035), s = 0.72, v = 1))(100))
    }
  }
}

# also save to GraphML
ExportNetwork <- function(fileName){
  current.net <- ppi.comps[[current.net.nm]];
  write.graph(current.net, file=fileName, format="graphml");
}


ExtractModule<- function(nodeids){
  set.seed(8574);
  nodes <- strsplit(nodeids, ";")[[1]];
  
  g <- ppi.comps[[current.net.nm]];
  # try to see if the nodes themselves are already connected
  hit.inx <- V(g)$name %in% nodes; 
  gObj <- induced.subgraph(g, V(g)$name[hit.inx]);
  
  # now find connected components
  comps <-decompose.graph(gObj, min.vertices=1);
  
  if(length(comps) == 1){ # nodes are all connected
    g <- comps[[1]];
  }else{
    # extract modules
    paths.list <-list();
    sd.len <- length(nodes);
    for(pos in 1:sd.len){
      paths.list[[pos]] <- get.shortest.paths(g, nodes[pos], nodes[-(1:pos)])$vpath;
    }
    nds.inxs <- unique(unlist(paths.list));
    nodes2rm <- V(g)$name[-nds.inxs];
    g <- simplify(delete.vertices(g, nodes2rm));
  }
  nodeList <- get.data.frame(g, "vertices");
  if(nrow(nodeList) < 3){
    return ("NA");
  }

  if(ncol(nodeList) == 1){
    nodeList <- data.frame(Id=nodeList[,1], Label=nodeList[,1]);
  }
  
  module.count <- module.count + 1;
  module.nm <- paste("module", module.count, sep="");
  colnames(nodeList) <- c("Id", "Label");
  ndFileNm <- paste(module.nm, "_node_list.csv", sep="");
  fast.write(nodeList, file=ndFileNm, row.names=F);
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  edgFileNm <- paste(module.nm, "_edge_list.csv", sep="");
  fast.write(edgeList, file=edgFileNm, row.names=F);
  
  filenm <- paste(module.nm, ".json", sep="");
  
  # record the module 
  ppi.comps[[module.nm]] <<- g;
  UpdateSubnetStats();
  
  module.count <<- module.count
  
  convertIgraph2JSON(module.nm, filenm);
  return (filenm);
}

PerformLayOut <- function(net.nm, algo, focus=""){
  paramSet <- readSet(paramSet, "paramSet");
  lib.path <- paramSet$lib.path;
  g <- ppi.comps[[net.nm]];
  vc <- vcount(g);
  if(algo == "Default"){
    if(vc > 3000) {
      pos.xy <- layout.lgl(g, maxiter = 100);
    }else if(vc > 2000) {
      pos.xy <- layout.lgl(g, maxiter = 150);
    }else if(vc > 1000) {
      pos.xy <- layout.lgl(g, maxiter = 200);
    }else if(vc < 150){
      pos.xy <- layout.kamada.kawai(g);
    }else{
      pos.xy <- layout.fruchterman.reingold(g);
    }
  }else if(algo == "FrR"){
    pos.xy <- layout.fruchterman.reingold(g);
  }else if(algo == "random"){
    pos.xy <- layout.random(g);
  }else if(algo == "lgl"){
    if(vc > 3000) {
      pos.xy <- layout.lgl(g, maxiter = 100);
    }else if(vc > 2000) {
      pos.xy <- layout.lgl(g, maxiter = 150);
    }else {
      pos.xy <- layout.lgl(g, maxiter = 200);
    }
  }else if(algo == "gopt"){
    # this is a slow one
    if(vc > 3000) {
      maxiter <- 50;
    }else if(vc > 2000) {
      maxiter <- 100;
    }else if(vc > 1000) {
      maxiter <- 200;
    }else{
      maxiter <- 500;
    }
    pos.xy <- layout.graphopt(g, niter=maxiter);
  }else if(algo == "fr"){
    pos.xy <- layout_with_fr(g, dim=3, niter=500)
  }else if(algo == "kk"){
    pos.xy <- layout_with_kk(g, dim=3, maxiter=500)
  }else if(algo == "tree"){
    if(ppi.net$db.type == "signal"){
      g <- ppi.comps[[net.nm]];
      current.net.nm <<- net.nm
      # annotation
      nms <- V(g)$name;
      hit.inx <- match(nms, ppi.net$node.data[,1]);
      lbls <- ppi.net$node.data[hit.inx,2];
      
      ids.arr <- c(as.character(edge.infoU$Source), as.character(edge.infoU$Target));
      type.arr <- c(as.character(edge.infoU$TYPEA), as.character(edge.infoU$TYPEB));
      inx <- !duplicated(ids.arr);
      ids.arr <- ids.arr[inx];
      type.arr <- type.arr[inx];
      names(type.arr) <- ids.arr;
      node.types <- list();
      node.types[[1]] <- type.arr[nms];
      
      path <- paste(lib.path,  "/signor/sig_location.rds", sep="");
      sig.sets <- readRDS(path);
      path2 <- paste(lib.path,  "/signor/sig_types.rds", sep="");
      sig.sets2 <- readRDS(path);
      sig.types <- rep("cytoplasm", length(lbls));
      names(sig.types) <- nms
      hit.inx <- sig.sets[,1] %in% nms
      subset <- sig.sets[hit.inx,]
      sig.types[subset[,1]] = subset[,2]
      
      layout.levels <- sig.types
      
      extra.inx <- layout.levels %in% "extracellular"
      memb.inx <- layout.levels %in% "membrane" | toupper(lbls) %in% sig.sets[["membrane"]]
      cyto.inx <- layout.levels %in% "cytoplasm"
      nucl.inx <- layout.levels %in% "nucleus"
      layout.levels[c(1:length(layout.levels))] <- 0
      layout.levels <- as.numeric(unname(layout.levels))
      
      #assign layers according to cell location
      layout.levels[extra.inx]<-max(layout.levels)+1
      layout.levels[memb.inx]<-max(layout.levels)+1
      layout.levels[cyto.inx]<-max(layout.levels)+1
      if("TRUE" %in% names(table(cyto.inx)) ){
        cyt.num <- max(layout.levels)
      }else{
        cyt.num <- -1
      }
      layout.levels[nucl.inx]<-max(layout.levels) +1
      
      pheno.inx <- node.types[[1]] == "phenotype"
      layout.levels[pheno.inx] <- max(layout.levels)+1    
      
      #assign more vertical space for cytoplasmic proteins
      if(cyt.num != -1){
        layout.levels[which(layout.levels>cyt.num)] <- layout.levels[which(layout.levels>cyt.num)] + 0.5
        layout.levels[which(layout.levels<cyt.num)] <- layout.levels[which(layout.levels<cyt.num)] - 0.5
      }
      
      l <- layout_with_sugiyama(g, layers= layout.levels, vgap=vc/5)
      pos.xy <- -l$layout
    }else{
      l <- layout_with_sugiyama(g, vgap=vc/4)
      pos.xy <- -l$layout
    }
  }else if(algo == "circular_tripartite"){
    library(ggforce)
    l <- layout_with_sugiyama(g, layers = V(g)$layers*(vc/3) +30)
    layout <- l$layout
    
    radial <- radial_trans(
      r.range = rev(range(layout[,2])),
      a.range = range(layout[,1]),
      offset = 0
    )
    coords <- radial$transform(layout[,2], layout[,1])
    layout[,1] <- coords$x
    layout[,2] <- coords$y
    pos.xy <-layout
  }else if(algo == "tripartite"){
    l <- layout_with_sugiyama(g, layers = V(g)$layers*(vc/4))
    pos.xy <- -l$layout[,2:1] 
  }else if(algo == "concentric"){
    library(graphlayouts)
    # the fist element in the list for concentric is the central node.
    if(focus==""){
      inx=1;
    }else{
      inx = which(V(g)$name == focus)
    }
    coords <- layout_with_focus(g,inx)
    pos.xy <- coords$xy
  }else if(algo == "backbone"){
    library(graphlayouts)
    if(length(V(g)$name)<2000){
      coords = layout_with_stress(g)
      pos.xy = coords
    }else{
      coords = layout_with_sparse_stress(g,pivots=100)
      pos.xy = coords
    }
    
  }
  pos.xy;
}

UpdateNetworkLayout <- function(algo, filenm, focus=""){
  current.net <- ppi.comps[[current.net.nm]];
  pos.xy <- PerformLayOut(current.net.nm, algo, focus);
  nms <- V(current.net)$name;
  nodes <- vector(mode="list");
  if(algo %in% c("fr", "kk")){
    for(i in 1:length(nms)){
      nodes[[i]] <- list(
        id=nms[i],
        x=pos.xy[i,1],
        y=pos.xy[i,2],
        z=pos.xy[i,3]
      );
    }
  }else{
    for(i in 1:length(nms)){
      nodes[[i]] <- list(
        id=nms[i], 
        x=pos.xy[i,1], 
        y=pos.xy[i,2]
      );
    }
  }
  # now only save the node pos to json
  netData <- list(nodes=nodes);
  sink(filenm);
  cat(RJSONIO::toJSON(netData));
  sink();
  return(filenm);
}

getGraphStatsFromFile <- function(){
  g <- ppi.comps[[net.nmu]];
  nms <- V(g)$name;
  edge.mat <- get.edgelist(g);
  return(c(length(nms), nrow(edge.mat)));        
}

GetNetUploaded <- function(){
  netUploadU
}

GetNetsNameString <- function(){
  paste(rownames(net.stats), collapse="||");
}

# support walktrap, infomap and lab propagation
FindCommunities <- function(method="walktrap", use.weight=FALSE){
  paramSet <- readSet(paramSet, "paramSet")
  seed.expr <- paramSet$seed.expr;
  # make sure this is the connected
  current.net <- ppi.comps[[current.net.nm]];
  g <- current.net;
  if(!is.connected(g)){
    g <- decompose.graph(current.net, min.vertices=2)[[1]];
  }
  total.size <- length(V(g));
  
  if(use.weight){ # this is only tested for walktrap, should work for other method
    # now need to compute weights for edges
    egs <- get.edges(g, E(g)); #node inx
    nodes <- V(g)$name;
    # conver to node id
    negs <- cbind(nodes[egs[,1]],nodes[egs[,2]]);
    
    # get min FC change
    base.wt <- min(abs(seed.expr))/10;
    
    # check if user only give a gene list without logFC or all same fake value
    if(length(unique(seed.expr)) == 1){
      seed.expr <- rep(1, nrow(negs));
      base.wt <- 0.1; # weight cannot be 0 in walktrap
    }
    
    wts <- matrix(base.wt, ncol=2, nrow = nrow(negs));
    for(i in 1:ncol(negs)){
      nd.ids <- negs[,i];
      hit.inx <- match(names(seed.expr), nd.ids);
      pos.inx <- hit.inx[!is.na(hit.inx)];
      wts[pos.inx,i]<- seed.expr[!is.na(hit.inx)]+0.1;
    }
    nwt <- apply(abs(wts), 1, function(x){mean(x)^2})    
  }
  
  if(method == "walktrap"){
    fc <- walktrap.community(g);
  }else if(method == "infomap"){
    fc <- infomap.community(g);
  }else if(method == "labelprop"){
    fc <- label.propagation.community(g);
  }else{
    print(paste("Unknown method:", method));
    return ("NA||Unknown method!");
  }
  
  if(length(fc) == 0 || modularity(fc) == 0){
    return ("NA||No communities were detected!");
  }
  
  # only get communities
  communities <- communities(fc);
  community.vec <- vector(mode="character", length=length(communities));
  gene.community <- NULL;
  qnum.vec <- NULL;
  pval.vec <- NULL;
  rowcount <- 0;
  nms <- V(g)$name;
  hit.inx <- match(nms, ppi.net$node.data[,1]);
  sybls <- ppi.net$node.data[hit.inx,2];
  names(sybls) <- V(g)$name;
  for(i in 1:length(communities)){
    # update for igraph 1.0.1 
    path.ids <- communities[[i]];
    psize <- length(path.ids);
    if(psize < 5){
      next; # ignore very small community
    }
    if(netUploadU == 1){
      qnums <- psize;
    }else{
      hits <- seed.proteins %in% path.ids;
      qnums <- sum(hits);
    }
    if(qnums == 0){
      next; # ignor community containing no queries
    }
    
    rowcount <- rowcount + 1;
    pids <- paste(path.ids, collapse="->");
    #path.sybls <- V(g)$Label[path.inx];
    path.sybls <- sybls[path.ids];
    com.mat <- cbind(path.ids, path.sybls, rep(i, length(path.ids)));
    gene.community <- rbind(gene.community, com.mat);
    qnum.vec <- c(qnum.vec, qnums);
    
    # calculate p values (comparing in- out- degrees)
    subgraph <- induced.subgraph(g, path.ids);
    in.degrees <- degree(subgraph);
    out.degrees <- degree(g, path.ids) - in.degrees;
    ppval <- suppressWarnings(wilcox.test(in.degrees, out.degrees)$p.value);
    ppval <- signif(ppval, 3);
    pval.vec <- c(pval.vec, ppval);
    
    # calculate community score
    community.vec[rowcount] <- paste(c(psize, qnums, ppval, pids), collapse=";");
  }

  if(length(pval.vec)>1){
    ord.inx <- order(pval.vec, decreasing=F);
    community.vec <- community.vec[ord.inx];
    qnum.vec <- qnum.vec[ord.inx];
    ord.inx <- order(qnum.vec, decreasing=T);
    community.vec <- community.vec[ord.inx];
  }
  
  all.communities <- paste(community.vec, collapse="||");
  if(!is.null(gene.community)){
    colnames(gene.community) <- c("Id", "Label", "Module");
    fast.write(gene.community, file="module_table.csv", row.names=F);
    return(all.communities);
  }else{
    return("NA");
  }
  
}

community.significance.test <- function(graph, vs, ...) {
  subgraph <- induced.subgraph(graph, vs)
  in.degrees <- degree(subgraph)
  out.degrees <- degree(graph, vs) - in.degrees
  wilcox.test(in.degrees, out.degrees, ...)
}


convertIgraph2JSON <- function(net.nm, filenm){
  paramSet <- readSet(paramSet, "paramSet");
  anal.type <- paramSet$anal.type;
  lib.path <- paramSet$lib.path;

  g <- ppi.comps[[net.nm]];
  current.net.nm <<- net.nm
  # annotation
  nms <- V(g)$name;
  hit.inx <- match(nms, ppi.net$node.data[,1]);
  lbls <- ppi.net$node.data[hit.inx,2];
  
  # setup shape (gene circle, other squares)
  shapes <- rep("circle", length(nms));
  
  # get edge data
  edge.mat <- get.edgelist(g);
  if(!is.null(E(g)$mechanism)){
    edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], direction=E(g)$direction, mechanism=E(g)$mechanism);
  }else{
    edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], direction=E(g)$direction);
  }
  
  # now get coords
  #pos.xy <- PerformLayOut_mem(net.nm, "Default");
  if(ppi.net$db.type == "signal"){
    pos.xy <- PerformLayOut(net.nm, "tree");
  }else{
    pos.xy <- PerformLayOut(net.nm, "Default");
  }
  
  # get the note data
  node.btw <- as.numeric(betweenness(g));
  node.dgr <- as.numeric(degree(g));
  node.exp <- as.vector(expr.vec[nms]);
  
  node.exp[is.na(node.exp)] <- 0;
  # node size to degree values
  if(vcount(g)>500){
    min.size <- 2;
  }else if(vcount(g)>200){
    min.size <- 3;
  }else{
    min.size <- 4;
  }
  node.sizes <- as.numeric(rescale2NewRange((log(node.dgr))^2, min.size, 12));
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
    node.colsc.exp[bad.inx] <- "#c6c6c6"; 
    # node.colsw.exp[bad.inx] <- "#b3b3b3";
  }else{
    node.colsb.exp <- rep("#d3d3d3",length(node.exp)); 
    node.colsw.exp <- rep("#c6c6c6",length(node.exp)); 
    node.colsc.exp <- rep("#b3b3b3",length(node.exp)); 
  }
  
  # now update for bipartite network
  if(is.null(V(g)$moltype)){
    mol.types <- rep("NA",length(node.exp)); 
  }else{
    mol.types <- V(g)$moltype; 
  }

  notbipartite <- c("ppi", "tissuecoex", "cellcoex", "tissueppi")
  if(!ppi.net$db.type %in% c("ppi", "tissuecoex", "cellcoex", "tissueppi", "uploaded")){ # the other part miRNA or TF will be in square
    if(ppi.net$db.type == "tfmir"){
      db.path <- paste(lib.path, data.org, "/mirlist.rds", sep="");
      db.map <-  readRDS(db.path);
      
      tf.inx <- nms %in% edge.mat[,"source"];
      mir.inx <- nms %in% db.map[,1];
      
      shapes[tf.inx] <- "diamond";
      shapes[mir.inx] <- "square";
      
      node.sizes[mir.inx] <- node.sizes[mir.inx] + 0.5;
      
      # update mir node color
      
      node.colsw.exp[tf.inx] <- topo.colsw[tf.inx] <- "#00ff00";
      node.colsw.exp[tf.inx] <- topo.colsw[tf.inx] <- "#00ff00";
      node.colsc.exp[tf.inx] <- topo.colsc[tf.inx] <- "#00ff00";
      node.colsw.exp[mir.inx] <- topo.colsw[mir.inx] <- "#306EFF"; # dark blue
      node.colsb.exp[mir.inx] <- topo.colsb[mir.inx] <- "#98F5FF";
      node.colsc.exp[mir.inx] <- topo.colsb[mir.inx] <- "#98F5FF";
      mol.types <- rep("Protein",length(node.exp)); 
      mol.types[tf.inx] <- "TF";
      mol.types[mir.inx] <- "miRNA";
    }else if(ppi.net$db.type == "signal"){
      ids.arr= c(as.character(edge.infoU$Source), as.character(edge.infoU$Target))
      type.arr= c(as.character(edge.infoU$TYPEA), as.character(edge.infoU$TYPEB))
      inx <- !duplicated(ids.arr)
      ids.arr <- ids.arr[inx]
      type.arr <- type.arr[inx]
      names(type.arr) <- ids.arr
      node.types <- list();
      node.types[[1]] <- type.arr[nms]
      
      path <- paste(lib.path,  "/signor/sig_location.rds", sep="");
      sig.sets <- readRDS(path)
      path2 <- paste(lib.path,  "/signor/sig_types.rds", sep="");
      sig.sets2 <- readRDS(path)
      sig.types <- rep("cytoplasm", length(lbls));
      names(sig.types) <- nms
      hit.inx <- sig.sets[,1] %in% nms
      subset <- sig.sets[hit.inx,]
      sig.types[subset[,1]] = subset[,2]
      
      pheno.inx <- node.types[[1]] == "phenotype"
      mem.inx <- toupper(lbls) %in% sig.sets[["membrane"]]
      
      sig.types[mem.inx]<-"membrane";
      sig.types[pheno.inx]<-"phenotype";
      
      
      node.cols <- rep("#ff4500", length(node.dgr));
      ntype <- unique(node.types[[1]])
      color.vec <- .gg_color_hue(length(ntype))
      for(i in 1:length(ntype)){
        node.cols[which(node.types[[1]] ==ntype[i])]=color.vec[i]
      }
      colVec <- unique(node.cols)
      
      mir.inx <- is.na(doEntrez2UniprotMapping(nms)) # if not gene of host org
      shapes[mir.inx] <- "square";
      node.sizes[mir.inx] <- node.sizes[mir.inx] + 0.5;
      
      # update mir node color
      node.colsw.exp <- topo.colsw <- node.cols; # dark blue
      node.colsb.exp <- topo.colsb <- node.cols;
      node.colsc.exp <- topo.colsc <- node.cols;
      
      pos.xy <- PerformLayOut(net.nm, "tree");
      mol.types <- unname(node.types[[1]])
      uni.vec <- doEntrez2UniprotMapping(nms)
    }  else if(ppi.net$db.type == "multi" || ppi.net$db.type == "tf"){
      if(!is.null(V(g)$Types)){
        node.types <-V(g)$Types
        node.cols <- topo.colsb
        ntype <- unique(node.types)
        color.vec <- .gg_color_hue(length(ntype))
        shapes = node.types;
        shape.vec = c("square", "diamond", "equilateral");
        for(i in 1:length(ntype)){
          
          if(ntype[i] %in% c("Protein","Seed")){
            #shapes[shapes == "Protein"] = "circle"
            shapes[shapes == "Seed"] = "circle"
            node.cols[which(node.types =="Seed")]="#BD0313"
            #node.cols[which(node.types =="Protein")]="#BD0313"

          }else{
            shapes[which(shapes ==ntype[i])]=shape.vec[i]
            node.cols[which(node.types ==ntype[i])]=color.vec[i]
          }
        }
        colVec <- unique(node.cols)
        mol.types<-node.types
      }else{
        mol.types <- V(g)$Type; 
        #seed.inx <- nms %in% unique(seed.proteins);
        #mol.types[seed.inx] <- "Seed"
        shapes = nrep("circle", length(nms))
      }
      
      topo.colsw <- node.cols; # dark blue
      topo.colsb <- node.cols;
      topo.colsc <- node.cols;
    }else{
      mir.inx <- nms %in% edge.mat[,"target"];
      shapes[mir.inx] <- "square";
      node.sizes[mir.inx] <- node.sizes[mir.inx] + 0.5;
      
      # update mir node color
      topo.colsw[mir.inx] <- "#306EFF"; # dark blue
      topo.colsb[mir.inx] <- "#98F5FF";
      topo.colsc[mir.inx] <- "#98F5FF";
      dbType <- ppi.net$db.type
      moltype <- "";
      if(dbType == "mir"){
        moltype <- "miRNA";
      }else if(dbType == "tf"){
        moltype <- "Transcription Factor";
      }else if(dbType == "drug"){
        moltype <- "Drug";
      }else if(dbType == "crossppi"){
        mol.types <- rep("Host Protein",length(node.exp)); 
        moltype <- "Foreign Protein";
      }else if(dbType == "disease"){
        moltype <- "Disease";
      }else if(dbType == "chem"){
        moltype <- "Chemical";
      }
      mol.types[mir.inx] = moltype;
    }
  }else{
    seed.inx <- nms %in% unique(seed.proteins);
    mol.types[seed.inx] <- "Seed"
    mol.types[!seed.inx] <- "Protein"
  }

  seed.inx <- nms %in% unique(seed.proteins);
  mol.types[seed.inx] <- "Seed"
  
  newreg.inx <- nms %in% regids;
  mol.types[newreg.inx] <- "regulator"
  reg.inx <- mol.types == "regulator"
  shapes[reg.inx] <- "diamond";
  node.sizes[reg.inx] <- node.sizes[reg.inx] + 1;
  
  # update mir node color
  topo.colsw[reg.inx] <- "#228B22"; # dark blue
  topo.colsb[reg.inx] <- "#00FF00";
  topo.colsc[reg.inx] <- "#00FF00";
  
  freq = table(mol.types)
  
  duplicated.types=mol.types
  for(i in 1:length(unique(mol.types))){
    duplicated.types[duplicated.types == names(freq[i])]=order(freq)[i]
  }
  
  node.cols <- topo.colsb
  ntype <- names(freq)
  color.vec <- .gg_color_hue(length(ntype))
  for(i in 1:length(ntype)){
    if(ntype[i] %in% c("Protein", "Seed", "Gene") ){
      color.vec[i] = "#BD0313"
    }else{
      node.cols[which(mol.types ==ntype[i])]=color.vec[i]
    }
  }
  
  if(length(ntype)==1){
    seed.inx <- nms %in% unique(seed.proteins);
    mol.types[seed.inx] <- "Seed"
    mol.types[!seed.inx] <- "Protein"
    
  }else{
    topo.colsw <- node.cols; # dark blue
    topo.colsb <- node.cols;
    topo.colsc <- node.cols;
  }
  colVec <- color.vec
  
  V(g)$moltype <- mol.types;
  V(g)$layers <- as.numeric(as.factor(mol.types));
  
  ppi.comps[[net.nm]] <<- g;
  
  seed.inx <- nms %in% unique(seed.proteins);
  seed_arr <- rep("notSeed",length(node.dgr));
  seed_arr[seed.inx] <- "seed";
  
  # now create the json object
  nodes <- vector(mode="list");
  for(i in 1:length(node.sizes)){

    nodes[[i]] <- list(
      id=nms[i],
      idnb = i, 
      label=lbls[i],
      x = pos.xy[i,1],
      y = pos.xy[i,2],
      molType = mol.types[i],
      size=node.sizes[i], 
      seedArr = seed_arr[i],
      type=shapes[i],
      colorb=topo.colsb[i],
      colorw=topo.colsw[i],
      colorc=topo.colsc[i],
      expcolb=node.colsb.exp[i],
      expcolw=node.colsw.exp[i],
      expcolc=node.colsc.exp[i],
      highlight = 0,
      attributes=list(
        expr = node.exp[i],
        degree=node.dgr[i], 
        between=node.btw[i])
    );
    if(ppi.net$db.type == "signal"){
      nodes[[i]]["uniprot"] = uni.vec[i]
      nodes[[i]]["molLocation"] = sig.types[i]
    }else if(!is.null(nrow(node.infoU)) && nrow(node.infoU)>0){
        if("domain" %in% colnames(node.infoU)){
            inx <- which(node.infoU[,1] == nms[i])
            bool <- is.integer(inx) && length(inx) == 0
            if(!bool){
            nodes[[i]]["species"] = node.infoU[inx,"species"];
            nodes[[i]]["domain"] = node.infoU[inx,"domain"];
            if( nodes[[i]]["label"] == "hypothetical_protein"){
                nodes[[i]]["label"] = nodes[[i]]["id"]
            }

            }else{
            nodes[[i]]["species"] = "NA";
            nodes[[i]]["domain"] = "NA";
            }
        }
    }
  }

  a = 0;
  
  # save node table
  nd.tbl <- data.frame(Id=nms, Label=lbls, Degree=node.dgr, Betweenness=round(node.btw,2), Expression=node.exp);
  # order 
  ord.inx <- order(nd.tbl[,3], nd.tbl[,4], decreasing = TRUE)
  nd.tbl <- nd.tbl[ord.inx, ];
  fast.write(nd.tbl, file="node_table.csv", row.names=FALSE);
  
  if(length(V(g)$name)>100 && ppi.net$db.type != "uploaded"){
    modules <- FindCommunities("walktrap", FALSE);
  }else{
    modules <- "NA"
  }
  node.res <- ppi.net$node.res
  if("domain" %in% colnames(node.res) || ppi.net$db.type == "signal"){
    node.info <- edge.infoU;
    node.types <- node.types
  }else{
    node.info <- "";
    node.types <- "";
  }
  
  #calculate clustering coefficient

  #globalProperties <-list();
  #globalProperties[["Diameter"]] <-diameter(g)
  #globalProperties[["Radius"]] <-radius(g)
  #globalProperties[["Average path length"]] <-signif(mean_distance(g), 3);
  #globalProperties[["Clustering coefficient"]] <- signif(transitivity(g, type="global"), 3);
  # covert to json
  if(ppi.net$db.type == "signal"){
    netData <- list(nodes=nodes, edges=edge.mat, org=data.org, analType=anal.type, naviString = "network", modules=modules, node.info = node.info, nodeTypes= unique(unname(node.types[[1]])), nodeColors = colVec, tblNm=table.nmu, idType="entrez");
  }else{
    netData <- list(nodes=nodes, edges=edge.mat, org=data.org, analType=anal.type, naviString = "network", modules=modules, tblNm=table.nmu, nodeTypes= ntype, nodeColors = colVec,idType="entrez");
  }
  partialToBeSaved <<- c(partialToBeSaved, c(filenm))
  
  sink(filenm);
  cat(RJSONIO::toJSON(netData));
  sink();
  
  analSet[[filenm]] <- netData;

  #if(!.on.public.web){
  #  library(httr);
  #  r <- POST("localhost:8080/NetworkAnalyst/faces/R_REQUEST_NET", body = list(organism = data.org, idtype = "entrez", network = rjson::toJSON(netData)))
  #}
  return(analSet);
}


# from to should be valid nodeIDs
GetShortestPaths <- function(from, to){
  current.net <- ppi.comps[[current.net.nm]];
  paths <- get.all.shortest.paths(current.net, from, to)$res;
  if(length(paths) == 0){
    return (paste("No connection between the two nodes!"));
  }
  
  path.vec <- vector(mode="character", length=length(paths));
  for(i in 1:length(paths)){
    path.inx <- paths[[i]]; 
    path.ids <- V(current.net)$name[path.inx];
    path.sybls <- path.ids;
    pids <- paste(path.ids, collapse="->");
    psbls <- paste(path.sybls, collapse="->");
    path.vec[i] <- paste(c(pids, psbls), collapse=";")
  }
  
  if(length(path.vec) > 50){
    path.vec <- path.vec[1:50];
  }
  
  all.paths <- paste(path.vec, collapse="||");
  return(all.paths);
}

GetNetsName <- function(){
  rownames(net.stats);
}

GetNetsEdgeNum <- function(){
  as.numeric(net.stats$Edge);
}

GetNetsNodeNum <- function(){
  as.numeric(net.stats$Node);
}

GetNetsQueryNum <- function(){
  as.numeric(net.stats$Query);
}

CheckCor <- function(){
  if(is.null(E(overall.graph)$correlation)){
    return(0);
  }else{
    return(1);
  }
}

GetNetworkTopology <- function(netnm){
  g <- ppi.comps[[netnm]];
  globalProperties <-list();
  globalProperties[["Diameter"]] <-diameter(g);
  globalProperties[["Radius"]] <-radius(g);
  globalProperties[["Average path length"]] <-signif(mean_distance(g), 3);
  globalProperties[["Clustering coefficient"]] <- signif(transitivity(g, type="global"), 3);
  propertiesVector <- c(globalProperties[[1]], globalProperties[[2]], globalProperties[[3]], globalProperties[[4]]);
  return(propertiesVector);
}


PlotDegreeHistogram <- function(imgNm, netNm = "NA", dpi=72, format="png"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
    library(ggplot2)
    if(netNm != "NA"){
        overall.graph <- ppi.comps[[netNm]];
    }
    G.degrees <- degree(overall.graph)

    G.degree.histogram <- as.data.frame(table(G.degrees))
    G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])

    p <- ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
        geom_point() +
        scale_x_continuous("Degree\n(nodes containing that amount of connections)",
                           breaks = c(1, 3, 10, 30, 100, 300),
                           trans = "log10") +
        scale_y_continuous("Frequency\n(number of nodes)",
                           breaks = c(1, 3, 10, 30, 100, 300, 1000),
                           trans = "log10") +
        ggtitle("Degree Distribution (log-log)") +
        theme_bw()  +
        theme(plot.title = element_text(hjust = 0.5))
    print(p)
  dev.off();
}

PlotBetweennessHistogram <- function(imgNm, netNm = "NA",dpi=72, format="png"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
    Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
        library(ggplot2)
        if(netNm != "NA"){
            overall.graph <- ppi.comps[[netNm]];
        }
        G.degrees <- betweenness(overall.graph)

        G.degree.histogram <- as.data.frame(table(G.degrees))
        G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])

        p <- ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
            geom_point() +
            scale_x_continuous("Betweenness\n(nodes with that amount of betweenness)",
                               breaks = c(1, 3, 10, 30, 100, 300,1000,3000,10000,30000),
                               trans = "log10") +
            scale_y_continuous("Frequency\n(number of nodes)",
                               breaks = c(1, 3, 10, 30, 100, 300, 1000),
                               trans = "log10") +
            ggtitle("Betweenness Distribution (log-log)") +
            theme_bw()  +
            theme(plot.title = element_text(hjust = 0.5))
        print(p)
    dev.off();
}


##########################################
############# private utility methods #### 
##########################################

.gg_color_hue <- function(grp.num, filenm=NULL) {
  grp.num <- as.numeric(grp.num)
  pal18 <- c( "#4363d8","#98f5ff" ,  "#3cb44b", "#f032e6", "#ffe119", "#e6194B", "#f58231", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#42d4f4","#000075", "#ff4500");
  if(grp.num <= 18){ # update color and respect default
    colArr <- pal18[1:grp.num];
  }else{
    colArr <- colorRampPalette(pal18)(grp.num);
  }
  if(is.null(filenm)){
    return(colArr);
  }else{
    sink(filenm);
    cat(rjson::toJSON(colArr));
    sink();
    return(filenm);
  }
}