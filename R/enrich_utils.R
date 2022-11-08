##################################################
## R script for NetworkAnalyst
## Description: Functions for heatmaps
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

PerformVolcanoEnrichment<-function(dataSet=NA, file.nm, fun.type, IDs, type, inx){
    dataSet <- .get.mSet();
    inx <- as.numeric(inx)
    if(anal.type == "onedata"){
        if(dataSet$type == "array"){
            sigmat <- dataSet$sig.mat
        } else {
            sigmat <- dataSet$sig.mat
        }
    }else{
        sigmat <- inmex.ind[selDataNm][[1]][which(inmex.ind[selDataNm][[1]][,'Pval'] < as.numeric(pvalu)),];
    }
    sigm <<- sigmat

    if(type == "focus"){
        gene.vec <- unlist(strsplit(IDs, "; "));
    }else if(type == "all"){
    gene.vecup <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),]);
    gene.vecdown <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),]);
        gene.vec <- c(gene.vecup, gene.vecdown);
    }else if(type == "up"){
        gene.vec <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),]);
    }else{
        gene.vec <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),]);
    }
    sym.vec <- doEntrez2SymbolMapping(gene.vec);
    names(gene.vec) <- sym.vec;
    res <- .performEnrichAnalysis(dataSet, file.nm, fun.type, gene.vec);
    if(res == 0){
    return(res);
    }else{
    return(.set.mSet(dataSet));
    }
}


PerformGSEA<- function(dataSet=NA, file.nm, fun.type, netNm, mType, selectedFactorInx, mode = "multi"){
    save.image("gsea.RData");
    dataSet <- .get.mSet();
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
  
    #get gene symbols
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

GetQEA.pathNames<-function(){
  current.geneset <- qs::qread("current_geneset.qs")
  hit.inx <- match(rownames(analSet$qea.mat),names(current.geneset));
  return(names(current.geneset)[hit.inx]);
}

PlotGSView <-function(cmpdNm,  format="png", dpi=72, width=NA){
  library("ggplot2");
  library("fgsea");
  current.geneset <- qs::qread("current_geneset.qs")
  imgName <- gsub("\\/", "_",  cmpdNm);
  imgName <- gsub(" ", "_",  imgName);
  imgName <- paste(imgName, "_dpi", dpi, ".", format, sep="");
  #indx<-which(rownames(boxplot_id)==cmpdNm);
  #gene.id <- boxplot_id[indx,1];
  Cairo(file = imgName, dpi=72, width=340, height=300, type="png", bg="transparent");
  g <- plotEnrichment(current.geneset[[cmpdNm]], rankedVec)
  print(g)
  dev.off();
  return(imgName);
}

PlotGShm <-function(cmpdNm, IDs){
  ids <- unlist(strsplit(IDs, "; "));
  cmpdNm <- gsub(" ", "_",  cmpdNm);
  cmpdNm <- gsub("/", "_",  cmpdNm);
  if(anal.type == "onedata"){
    subset <- dataSet$data.norm[which(doEntrez2SymbolMapping(rownames(dataSet$data.norm)) %in% ids),]
    if(length(subset)<1){
      subset <- dataSet$data.norm[which(rownames(dataSet$data.norm) %in% ids),]
    }
  }else{
    if(dataSet$name != selDataNm){
      dataSet <- qs::qread(selDataNm);
    }
    subset <- dataSet$data[which(doEntrez2SymbolMapping(rownames(dataSet$data)) %in% ids),]
    if(length(subset)<1){
      subset <- dataSet$data[which(rownames(dataSet$data) %in% ids),]
    }
  }
  inx <- order(dataSet$fst.cls)
  subset <- subset[,inx]
  
  json.res <- list(
    data=subset,
    ids=doEntrez2SymbolMapping(rownames(subset)),
    entrez = rownames(subset)
  )

  json.mat <- RJSONIO::toJSON(json.res);
  json.nm <- paste(cmpdNm,"_hm", ".json", sep="");
  
  sink(json.nm);
  cat(json.mat);
  sink();
  
  return(json.nm)
}

#For GSEA of AnalOverview page
ComputeRankedVec <- function(data, inx = 1){
    opt <- rankOptGlobal;
    cls <- data$cls
    if(anal.type == "metadata"){
        if(opt %in% c("mwt", "s2n", "wcx", "stu")){
           matr <- as.matrix(data$data)
        }else{
           matr <- as.matrix(qs::qread("meta.resTable.qs"));
        }
    }else{
        if(opt %in% c("mwt", "s2n", "wcx", "stu")){
           matr <- as.matrix(data$data.norm)
        }else{
           matr <- as.matrix(data$resTable);
        }
    }
    if(opt == "mwt"){
        res <- CalculateMWT(matr, cls)
        rankedVec <- res$MWT
        names(rankedVec) <- rownames(matr)
    }else if(opt == "stu"){
        inx1 <- which(data$cls==levels(data$cls)[1]);
        inx2 <- which(data$cls==levels(data$cls)[2]);
        res <- apply(matr, 1, function(x) {
                tmp <- try(t.test(x[inx1], x[inx2], paired = FALSE, var.equal = TRUE));
                if(class(tmp) == "try-error") {
                    return(c(NA, NA));
                }else{
                    return(c(tmp$statistic, tmp$p.value));
                }
        })
      res <- t(res)
      rankedVec <- res[,1]
      posInx <- sign(rankedVec) == 1
      rankedVec[posInx] <- 1000-rankedVec[posInx]
      names(rankedVec) <- rownames(matr)

    }else if(opt == "wcx"){
        inx1 <- which(data$cls==levels(data$cls)[1]);
        inx2 <- which(data$cls==levels(data$cls)[2]);
        res <- apply(matr, 1, function(x) {
                tmp <- try(wilcox.test(x[inx1], x[inx2], paired = FALSE));
                if(class(tmp) == "try-error") {
                    return(c(NA, NA));
                }else{
                    return(c(tmp$statistic, tmp$p.value));
                }
        })
      res <- t(res)
      rankedVec <- res[,2]
      names(rankedVec) <- rownames(matr)

    }else if (opt == "s2n"){
        res <- CalculateS2N(matr, as.numeric(cls)-1)
        rankedVec <- res;
    }else if(opt == "pval"){
        m <- length(colnames(matr))
        if (dataSet$de.method=="limma"){
            if("t" %in% colnames(matr)){
            rankedVec <- as.vector(matr[,"t"]);
            }else{
            rankedVec <- as.vector(matr[,"F"]);
            return(0);
            }
        } else if (dataSet$de.method=="deseq2"){
            rankedVec <- as.vector(matr[,"stat"]);
        } else {
            rankedVec <- as.vector(matr[,"LR"]);
        }
        names(rankedVec) <- rownames(matr);
    }else{
        rankedVec <- as.vector(matr[,inx]);
        names(rankedVec) <- rownames(matr);
    }
        rankedVec <- sort(rankedVec)
        rankedVec <- rankedVec[unique(names(rankedVec))]
        rankedVec <<- rankedVec
        return(rankedVec)
}


CalculateS2N <- function(data, vec = y.vec, ...) {
        
        A <- data + 0.00000001
        
        ind1 <- which(vec==1) # cases
        n1 <- length(ind1)    
        M1 <- rowMeans(A[,ind1])
        A2 <- A*A    
        S1 <- rowMeans(A2[,ind1])   
        S1 <- S1 - M1*M1    
        S1 <- sqrt(abs((n1/(n1-1)) * S1))   
        
        ind2 <- which(vec==0) # controls
        n2 <- length(ind2)
        M2 <- rowMeans(A[,ind2])
        S2 <- rowMeans(A2[,ind2])   
        S2 <- S2 - M2*M2    
        S2 <- sqrt(abs((n2/(n2-1)) * S2))   
        
        # small sigma "fix" as used in GeneCluster
        S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
        S2 <- ifelse(S2 == 0, 0.2, S2) 
        S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
        S1 <- ifelse(S1 == 0, 0.2, S1) 
        M1 <- M1 - M2
        S1 <- S1 + S2
        s2n <- M1/S1
        
        return(s2n)
}

CalculateMWT <- function(xdat,grp,na.rm=TRUE){
    ## basic statistics
    glab = unique(grp)
    n1 = sum(grp==glab[1])
    n2 = sum(grp==glab[2])
    d1 = n1-1
    d2 = n2-1
    m1 = rowMeans(xdat[,grp==glab[1]], na.rm=na.rm)
    m2 = rowMeans(xdat[,grp==glab[2]], na.rm=na.rm)
    
    s2.g1 = rowSums((xdat[,grp==glab[1]]-m1)^2, na.rm=na.rm)/d1
    
    ## We might either have all NA in one group or variance = 0
    ## (e.g. might happen with RMA with small samples)
    ## In this situation we want to remove the gene
    s2.g1[s2.g1 == 0] <- NA
    
    s2.g2 = rowSums((xdat[,grp==glab[2]]-m2)^2, na.rm=na.rm)/d2
    s2.g2[s2.g2 == 0] <- NA
    
    ## If either s2.g1 or s2.g2 are NA this will be NA
    sig2 = (d1*s2.g1 + d2*s2.g2)/(d1+d2)
    fac = 1/n1 + 1/n2
    se2 = (sig2 * fac)
    
    ## F test
    
    lev.test = levene(xdat,grp)
    fFDR = lev.test$FDR
    fStat = lev.test$statistic
    
    
    ## ordinary Welch statistics
    se2.sep = s2.g1/n1 + s2.g2/n2
    df = se2.sep^2/((s2.g1/n1)^2/d1 + (s2.g2/n2)^2/d2)
    
    ## weighted formulas
    df.w  = fFDR*(d1+d2) + (1-fFDR)*df
    se2.w = fFDR*se2 + (1-fFDR)*se2.sep
    ds = est.hyper(z=log(se2.w),D=mean(df.w,na.rm=na.rm),d12=d1+d2)   
    
    ## ....................................... moderated Welch
    se2.com = (ds$d0*ds$s2 + df.w*se2.w)/(ds$d0 + df.w)
    Wm.t = (m1-m2)/sqrt(se2.com) ## Welch t
    df.com = ds$d0 + df.w      ## df
    
    
    Wm.pval = pt(-abs(Wm.t), df= df.com) * 2
    
    ## ................. Compute Global FDR
    
    fdr <- NULL
    
    return(list(MWT= Wm.t, coefficients=cbind((m1-m2)),pvalue = Wm.pval))
    
}

levene <- function (xdat, grp, na.rm = TRUE) 
{
    glab = unique(grp)
    ngr = length(glab)
    n = mn = s2 = NULL
    X0 = NULL
    for (i in 1:ngr) {
        ndx = grp == glab[i]
        mni = rowMeans(xdat[, ndx], na.rm = na.rm)
        x0 = xdat[, ndx] - mni
        X0 = cbind(X0, x0)
    }
    xdat = abs(X0)
    for (i in 1:ngr) {
        ndx = grp == glab[i]
        ni = sum(ndx)
        mni = rowMeans(xdat[, ndx], na.rm = na.rm)
        x0 = xdat[, ndx] - mni
        s2i = rowSums(x0 * x0, na.rm = na.rm)/(ni - 1)
        n = c(n, ni)
        mn = cbind(mn, mni)
        s2 = cbind(s2, s2i)
    }
    N = sum(n)
    mmn = rowSums(xdat, na.rm = na.rm)/N
    mn0 = mn - mmn
    num = rowSums(t(t(mn0 * mn0) * n))/(ngr - 1)
    den = rowSums(t(t(s2 * (n - 1))))/(N - ngr)
    F3 = num/den
    pval = pf(F3, df1 = ngr - 1, df2 = N - ngr)
    pval = ifelse(pval < 0.5, 2 * pval, 2 * (1 - pval))
    lvn.FDR = pval2FDR(pval)
    return(list(statistic = F3, pvalue = pval, FDR = lvn.FDR))
}

pval2FDR <-function (pval, lim = 0.7) 
{
    n1 = length(pval)
    ok.id <- 1:n1
    if (any(is.na(pval))) {
        ok.id <- which(!is.na(pval))
        pval <- na.omit(pval)
    }
    n = length(pval)
    Fp = rank(pval)/length(pval)
    p0 = sum(pval > lim)/((1 - lim) * n)
    p0 = min(p0, 1)
    FDRp = p0 * pmin(pval/Fp, 1)
    ord = order(pval)
    FDR.o = FDRp[ord]
    b = rev(cummin(rev(FDR.o)))
    FDR = rep(0, n)
    FDR[ord] = b
    out.FDR <- rep(NA, n1)
    out.FDR[ok.id] <- FDR
    attr(out.FDR, "p0") <- p0
    return(out.FDR)
}

est.hyper <- function (z, D, d12) 
{
    f2 <- function(d0, D) {
        var(z, na.rm = TRUE) - trigamma(D/2) - trigamma(d0/2)
    }
    lim = f2(100, D)
    if (lim < 0) 
        d0.est <- 100
    if (lim > 0) 
        d0.est <- uniroot(f2, c(1, 100), D = D, extendInt = "yes")$root
    s2.est <- exp(mean(z, na.rm = TRUE) - digamma(D/2) + digamma(d0.est/2) - 
        log(d0.est/D))
    return(list(d0 = d0.est, s2 = s2.est))
}

# note: hit.query, resTable must synchronize
# ora.vec should contains entrez ids, named by their gene symbols
PerformVolcanoOneEnrichment <- function(file.nm, fun.type, IDs, inx){
    # prepare lib
    inx <- as.numeric(inx);
    if(anal.type == "onedata"){
        if(dataSet$type == "array"){
            sigmat <- dataSet$sig.mat
        } else {
            sigmat <- dataSet$sig.mat
        }
    }else{
        sigmat <- inmex.ind[selDataNm][[1]][which(inmex.ind[selDataNm][[1]][,'Pval'] < as.numeric(pvalu)),];
    }
    sigm <<- sigmat

    one.path.vec <- unlist(strsplit(IDs, "; "));

    gene.vecup <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),]);
    gene.vecdown <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),]);
    ora.vec <- c(gene.vecup, gene.vecdown);


    sym.vec <- doEntrez2SymbolMapping(ora.vec);
    names(ora.vec) <- sym.vec;
    
    current.geneset <- list()
    current.geneset[["Set"]] <- one.path.vec
    current.geneset[["Set2"]] <- one.path.vec

    # prepare query
    ora.nms <- names(ora.vec);

    # prepare for the result table
    set.size<-length(current.geneset);
    res.mat<-matrix(0, nrow=set.size, ncol=5);
    rownames(res.mat)<-names(current.geneset);
    colnames(res.mat)<-c("Total", "Expected", "Hits", "P.Value", "FDR");

    # need to cut to the universe covered by the pathways, not all genes
    current.universe <- unique(unlist(current.geneset)); 
    hits.inx <- ora.vec %in% current.universe;
    ora.vec <- ora.vec[hits.inx];
    ora.nms <- ora.nms[hits.inx];

    q.size<-length(ora.vec);

    # get the matched query for each pathway
    hits.query <- lapply(current.geneset, 
        function(x) {
            ora.nms[ora.vec%in%unlist(x)];
        }
    );

    qs::qsave(hits.query, "hits_query.qs");

    names(hits.query) <- names(current.geneset);
    hit.num<-unlist(lapply(hits.query, function(x){length(unique(x))}), use.names=FALSE);

    # total unique gene number
    #uniq.count <- length(current.universe);
    uniq.count <- nrow(dataSet$data.norm);
    
    # unique gene count in each pathway
    set.size <- unlist(lapply(current.geneset, length));

    res.mat[,1]<-set.size;
    res.mat[,2]<-q.size*(set.size/uniq.count);
    res.mat[,3]<-hit.num;

    # use lower.tail = F for P(X>x)
    raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
    res.mat[,4]<- raw.pvals;
    res.mat[,5] <- raw.pvals;

    # now, clean up result, synchronize with hit.query
    res.mat <- res.mat[hit.num>0,,drop = F];
    hits.query <- hits.query[hit.num>0];
    if(nrow(res.mat)> 0){
        # order by p value
        ord.inx<-order(res.mat[,4]);
        #res.mat <- signif(res.mat[ord.inx,],3);
        hits.query <- hits.query[ord.inx];

        imp.inx <- res.mat[,4] <= 0.05;
        if(sum(imp.inx) < 10){ # too little left, give the top ones
            topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
            res.mat <- res.mat[1:topn,];
            hits.query <- hits.query[1:topn];
        }else{
            res.mat <- res.mat[imp.inx,];
            hits.query <- hits.query[imp.inx];
            if(sum(imp.inx) > 120){
                # now, clean up result, synchronize with hit.query
                res.mat <- res.mat[1:120,];
                hits.query <- hits.query[1:120];
            }
        }
    }else{
        return(0);
    }

    #get gene symbols
    resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
    res.mat[,"Hits"] <- res.mat[,"Hits"]
    enr.mat <<- res.mat
    current.msg <<- "Functional enrichment analysis was completed";

    # write json
    fun.anot <- hits.query; 
    total <- resTable[,2]; if(length(total) ==1) { total <- matrix(total) };
    fun.pval <- resTable[,5]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
    fun.padj <- resTable[,6]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
    hit.num <- resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    fun.ids <- as.vector(current.setids[names(fun.anot)]); 
    if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
        json.res <- list(
                    fun.link = current.setlink[1],
                    fun.anot = fun.anot,
                    fun.ids = fun.ids,
                    fun.pval = fun.pval,
                    fun.padj = fun.padj,
                    hit.num = hit.num,
                    total= total
        );
     json.mat <- rjson::toJSON(json.res);
     json.nm <- paste(file.nm, ".json", sep="");
     
     sink(json.nm)
     cat(json.mat);
     sink();

    # write csv
    fun.hits <<- hits.query;
    fun.pval <<- resTable[,5];
    hit.num <<- resTable[,4];
    csv.nm <- paste(file.nm, ".csv", sep="");    
    fast.write(resTable, file=csv.nm, row.names=F);

    return(.set.mSet(dataSet));
}