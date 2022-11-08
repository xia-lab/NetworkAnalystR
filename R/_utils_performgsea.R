my.performgsea<- function(dataSetObj=NA, file.nm, fun.type, netNm, mType, selectedFactorInx, mode = "multi"){
    if(!exists("my.reg.enrich")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/networkanalystr/_utils_dogsea.Rc");    
    }
    dataSet <- .get.mSet(dataSetObj);
    current.geneset <- .loadEnrichLib(fun.type)
    require("fgsea");

    if(anal.type == "onedata"){
        datnorm <- dataSet$data.norm;
        inx  <- order(dataSet$fst.cls)
        sampleNms <- colnames(dataSet$data.norm);
        sampleNms <- sampleNms[inx]
        rankedVec<- ComputeRankedVec(dataSet, selectedFactorInx);

    }else{
        if(dataSet$name != selDataNm){
        dataSet <- qs::qread(selDataNm);
        }
        inx  <- rep(T, ncol(dataSet$data))
        datnorm <- dataSet$data
        sampleNms <- colnames(dataSet$data);
        ds <- inmex.ind[selDataNm][[1]]
        rankedVec <- ComputeRankedVec(dataSet, 1);
    }

    if(rankedVec == 0){
    current.msg <<- "Selected ranking method is not suitable. Please select another one!";
    return(0);
    }
  
    if(mode == "simple"){
        fgseaRes <- fgsea(pathways = current.geneset, 
                        stats = rankedVec,
                        minSize=15,
                        maxSize=500,
                        nperm=10000);
    } else {
        fgseaRes <- fgsea(pathways = current.geneset, 
                        stats = rankedVec,
                        minSize=15,
                        maxSize=500);
    }
  
    fgseaRes <- fgseaRes[!duplicated(fgseaRes$pathway),]
  
    rownames(fgseaRes) <- make.names(fgseaRes$pathway, unique=TRUE)
    fgseaRes <- fgseaRes[,c("size","ES", "pval", "pathway", "padj")]

    if(nrow(fgseaRes)<1){
        SetListNms()
        initsbls <- doEntrez2SymbolMapping(list.genes);
        names(initsbls) <- list.genes
        netData <- list(sizes=listSizes, genelist=initsbls);
        netName <- paste0(netNm, ".json");
        sink(netName);
        cat(RJSONIO::toJSON(netData));
        sink();
        return(0);
    }
  
    fgseaRes <- fgseaRes[order(fgseaRes$pval),]
    fgseaRes <<- fgseaRes
    if(nrow(fgseaRes[which(fgseaRes$pval < 0.05),])<20 ){
        if(nrow(fgseaRes)>20){
            fgseaRes <- fgseaRes[c(1:20),]
        }
    }else{
        fgseaRes <- fgseaRes[which(fgseaRes$pval < 0.05),]
    } 

    current.mset <- current.geneset[fgseaRes$pathway]
    current.mset <- current.mset[!duplicated(names(current.mset))]

    ora.vec <- names(rankedVec)
    ora.nms <- doEntrez2SymbolMapping(ora.vec)

    hits.query <- lapply(current.mset, 
                       function(x) {
                         unique(ora.nms[ora.vec%in%unlist(x)]);
                       });
    
    set.num <- unlist(lapply(current.mset, function(x){length(unique(x))}), use.names=TRUE);
    names(hits.query) <- names(current.mset);
    hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=TRUE);
    qs::qsave(hits.query, "hits_query.qs");
    fgseaRes$hits <- hit.num[which(fgseaRes$pathway  %in% names(hit.num))] 
    fgseaRes$total <- set.num[which(fgseaRes$pathway %in% names(set.num))]
  
    fgseaRes <- fgseaRes[which(fgseaRes$hits>1),]
    fgseaRes <- fgseaRes[which(fgseaRes$hits<500),]
    fgseaRes <- fgseaRes[which(fgseaRes$total<2000),]
    if(nrow(fgseaRes)<1){
        SetListNms();
        initsbls <- doEntrez2SymbolMapping(list.genes)
        names(initsbls) <- list.genes
        netData <- list(sizes=listSizes, genelist=initsbls);
        netName <- paste0(netNm, ".json");
        sink(netName);
        cat(RJSONIO::toJSON(netData));
        sink();
        return(0);
    }
  
    fgseaRes=fgseaRes[order(fgseaRes$pval),]
    if(nrow(fgseaRes[which(fgseaRes$pval < 0.05),])<20 ){
        if(nrow(fgseaRes)>20){
            fgseaRes <- fgseaRes[c(1:20),]
        }
    }else{
        fgseaRes <- fgseaRes[which(fgseaRes$padj < 0.05),]
    } 

    fgseaRes <- data.frame(fgseaRes, stringsAsFactors=FALSE)
  
    #get gene symbolsPerformGSEA
    current.msg <<- "Functional enrichment analysis was completed";
  
    # write json
    fun.anot <- hits.query; 
    fun.pval <- fgseaRes[,3]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
    #fun.pval<-signif(fun.pval,5);  
    fun.padj <- fgseaRes[,5]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
    #fun.padj<-signif(fun.padj,5);  
    es.num <- fgseaRes[,2]; if(length(es.num) ==1) { es.num <- matrix(es.num) };
    fun.ids <- as.vector(current.setids[names(fun.anot)]); 
    if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };

    json.res <- list(
        fun.link = current.setlink[1],
        fun.anot = fun.anot,
        fun.ids = fun.ids,
        fun.pval = fun.pval,
        fun.padj = fun.padj,
        pathname = fgseaRes[,"pathway"],
        es.num = es.num,
        hits = fgseaRes[,"hits"],
        total = fgseaRes[,"total"],
        cls = dataSet$meta.info[inx,],
        sample.nms = sampleNms
    );

    if(mType == "network"){
        res.mat<-matrix(0, nrow=length(fun.pval), ncol=4);
        colnames(res.mat)<-c("Total","Hits", "P.Value", "FDR");
        res.mat[,"Total"] <- fgseaRes[,"total"];
        res.mat[,"Hits"] <- fgseaRes[,"hits"];
        res.mat[,"P.Value"] <- fgseaRes[,"pval"];
        res.mat[,"FDR"] <- fgseaRes[,"padj"];
        res.mat <- data.matrix(data.frame(res.mat, stringsAsFactors=FALSE));
        rownames(res.mat) <- fgseaRes[,"pathway"];
        enr.mat <<- res.mat;
        list.genes <<- doEntrez2SymbolMapping(rownames(dataSet$sig.mat));
        SetListNms();
        PrepareEnrichNet(netNm, "meta", "mixed");
        file.nm <- gsub("gsea", "enrichment", file.nm)
        json.res$naviString <- "Enrichment Network"
    }else{
        json.res$org <- data.org
        json.res$analType <- anal.type
        json.res$naviString <- "GSEA";
    }

    json.mat <- RJSONIO::toJSON(json.res);
    json.nm <- paste(file.nm, ".json", sep="");
    partialToBeSaved <<- c(partialToBeSaved, c(json.nm, "current_geneset.qs"))
    sink(json.nm)
    cat(json.mat);
    sink();
  
    # write csv
  
    fgseaRes <<- fgseaRes 
    ftype <- fun.type
    if(fun.type %in% c("bp", "mf", "cc")){
        ftype <- paste0("go_", fun.type);
    }

    csvDf <- data.frame(Name=fgseaRes$pathway, Total=fgseaRes$total, Hits=fgseaRes$hits, EnrichmentScore=fgseaRes$ES, Pval=fgseaRes$pval, Padj=fgseaRes$padj);
    fast.write(csvDf, file=paste0(file.nm, ".csv"));
  
    return(.set.mSet(dataSet));
}  

