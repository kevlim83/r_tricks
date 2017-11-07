#hack reactomepa for arabidopsis
library(igraph)
library(reshape2)
library(randomcoloR)

clusterPathways <- function(g,thres){

    g2<-delete.edges(g,E(g)[E(g)$weight<thres])
    cf<-cluster_fast_greedy(g2)
    n<-length(cf)
    palette <- distinctColorPalette(n)
    V(g2)$color<-palette(membership(g2))    
    plot(g2,vertex.size=3,vertex.label=NA)
    list(g2,cf)

}

constructPathwayNetwork <- function(x){

    geneSets <- x@geneSets
    y <- as.data.frame(x)
    n <- nrow(y)
    id <- y[,1]
    geneSets <- geneSets[id]

    w <- matrix(NA, nrow=n, ncol=n)
    colnames(w) <- rownames(w) <- y$Description

    for (i in 1:n) {
        for (j in i:n) {
            w[i,j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
        }
    }

    wd <- melt(w)
    wd <- wd[wd[,1] != wd[,2],]
    wd <- wd[!is.na(wd[,3]),]
    g <- graph.data.frame(wd[,-3], directed=F)
    E(g)$weight<-wd[,3]
    g

}
                         
overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}


myenrichPathway <- function(gene,
                          organism="arabidopsis",
                          pvalueCutoff = 0.05,
                          pAdjustMethod="BH",
                          qvalueCutoff = 0.2,
                          universe,
                          minGSSize=10,
                          maxGSSize=500,
                          readable=FALSE) {

    Reactome_DATA <- get_myReactome_DATA(organism)
    
    res <- enricher_internal(gene,
                           pvalueCutoff=pvalueCutoff,
                           pAdjustMethod=pAdjustMethod,
                           qvalueCutoff=qvalueCutoff,
                           universe = universe,
                           minGSSize = minGSSize,
                           maxGSSize = maxGSSize,
                           USER_DATA = Reactome_DATA)

    if (is.null(res))
        return(res)

    res@keytype <- "ENTREZID"
    res@organism <- organism
    OrgDb <- getDb(organism)
    if (readable) {
        res <- setReadable(res, OrgDb)
    }
    res@ontology <- "Reactome"
    return(res)
}


mygsePathway <- function(geneList,
                       organism      = "arabidopsis",
                       exponent      = 1,
                       nPerm         = 1000,
                       minGSSize     = 10,
                       maxGSSize     = 500,
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH",
                       verbose       = TRUE,
                       seed          = FALSE,
                       by = 'fgsea') {

    Reactome_DATA <- get_myReactome_DATA(organism)
    
    res <- GSEA_internal(geneList      = geneList,
                         exponent      = exponent,
                         nPerm         = nPerm,
                         minGSSize     = minGSSize,
                         maxGSSize     = maxGSSize,
                         pvalueCutoff  = pvalueCutoff,
                         pAdjustMethod = pAdjustMethod,
                         verbose       = verbose,
                         USER_DATA     = Reactome_DATA,
                         seed          = seed,
                         by = by)

    if (is.null(res))
        return(res)
    
    res@organism <- organism
    res@setType <- "Reactome"
    res@keytype <- "ENTREZID"
    
    return(res)
}

get_myReactome_DATA <- function(organism = "arabidopsis") {
    ReactomePA_Env <- get_Reactome_Env()
    
    if (exists("organism", envir=ReactomePA_Env, inherits = FALSE)) {
        org <- get("organism", envir=ReactomePA_Env)
        if (org == organism &&
            exists("PATHID2EXTID", envir = ReactomePA_Env) &&
            exists("EXTID2PATHID", envir = ReactomePA_Env) &&
            exists("PATHID2NAME",  envir = ReactomePA_Env)) {
            
            ## use_cached
            return(ReactomePA_Env)
        }
    }
    
    #this should return ath locus ids
    ALLEG <- getmyALLEG(organism)
    
    #this should use the gmt file
    EXTID2PATHID <- athid2path
    PATHID2EXTID <- path2athid
    #EXTID2PATHID <- as.list(reactomeEXTID2PATHID)
    #EXTID2PATHID <- EXTID2PATHID[names(EXTID2PATHID) %in% ALLEG]
    
    #PATHID2EXTID <- as.list(reactomePATHID2EXTID) ## also contains reactions
    
    #PATHID2NAME <- as.list(reactomePATHID2NAME)
    #PI <- names(PATHID2NAME)
    PATHID2NAME <- pathid2name
    PI <- names(pathid2name)
    ## > PATHID2NAME[['68877']]
    ## [1] "Homo sapiens: Mitotic Prometaphase" "Homo sapiens: Mitotic Prometaphase"
    PATHID2NAME <- lapply(PATHID2NAME, function(x) x[1])
    names(PATHID2NAME) <- PI
    
    PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% names(PATHID2NAME)]
    PATHID2EXTID <- PATHID2EXTID[names(PATHID2EXTID) %in% unique(unlist(EXTID2PATHID))]
    PATHID2EXTID <- lapply(PATHID2EXTID, function(x) intersect(x, ALLEG))
    
    PATHID2NAME <- PATHID2NAME[names(PATHID2NAME) %in% names(PATHID2EXTID)]

    PATHID2NAME <- unlist(PATHID2NAME)
    PATHID2NAME <- gsub("^\\w+\\s\\w+:\\s+", "", PATHID2NAME) # remove leading spaces
    
    assign("PATHID2EXTID", PATHID2EXTID, envir=ReactomePA_Env)
    assign("EXTID2PATHID", EXTID2PATHID, envir=ReactomePA_Env)
    assign("PATHID2NAME", PATHID2NAME, envir=ReactomePA_Env)
    return(ReactomePA_Env)
}

getmyALLEG <- function(organism) {
    annoDb <- getDb(organism)
    require(annoDb, character.only = TRUE)
    annoDb <- eval(parse(text=annoDb))
    eg=keys(annoDb)
    return(eg)
}

