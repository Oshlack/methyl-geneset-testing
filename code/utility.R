loadAnnotation <- function(arrayType = c("450k", "EPIC")) {
    # Function to load 450K or EPIC array annotation data and save
    # locally as RData file for faster repeat execution.
    # Jovana Maksimovic
    # Modified: 14 August 2020

    arrayType <- match.arg(arrayType)
    outDir <- here::here("data/annotations")
    if (!dir.exists(outDir)) dir.create(outDir)
    annFile <- here::here(glue::glue("data/annotations/ann{arrayType}.RData"))

    if(file.exists(annFile)){
        # load annotation from file for faster execution
        load(annFile)
    } else {
        if(arrayType == "EPIC") {
            ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
        } else {
            ann <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
        }
        # save annotation to file to speed up subsequent execution
        save(ann, file=annFile)
    }

    ann
}

loadFlatAnnotation <- function(ann) {
    # Function to load 450K or EPIC flattened array annotation that associates
    # CpGs with Entrez IDs and save as RData file locally for faster repeat
    # execution.
    # Jovana Maksimovic
    # Modified: 14 August 2020

    arrayType <- ifelse(nrow(ann) > 500000, "EPIC", "450k")
    outDir <- here::here("data/annotations")
    if (!dir.exists(outDir)) dir.create(outDir)
    flatFile <- here::here(glue::glue("data/annotations/flatAnn{arrayType}.RData"))

    if(file.exists(flatFile)){
        # load flat annotation from file for faster execution
        load(flatFile)
    } else {
        flatAnn <- missMethyl:::.getFlatAnnotation(anno=ann)
        # save flat annotation to file to speed up subsequent execution
        save(flatAnn, file=flatFile)
    }

    flatAnn
}

getMode <- function(x) {
    # Function determine the mode.
    # Jovana Maksimovic

    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

readData <- function(saveAs, #dataSrc = c("BC","RA"),
                     saveRG = FALSE){
    # Function to load, normalise and filter the flow-sorted blood cell data
    # and save relevant objects as RData file for downstream use.
    # Jovana Maksimovic
    # Modified: 14 August 2020

    #dataSrc <- match.arg(dataSrc)

    #if(dataSrc == "BC"){
    dat <- getBloodCellData()
    rawRg <- dat$rawRg
    targets <- dat$targets

    # } else if (dataSrc == "RA"){
    #   dat <- getRmArthritisData()
    #   rawRg <- dat$rawRg
    #   targets <- dat$targets
    # }

    # calculate the detection p-values
    detP <- minfi::detectionP(rawRg)

    # normalise
    normGr <- minfi::preprocessQuantile(rawRg)

    # filtering
    fltGr <- filterQual(normGr, detP)
    fltGr <- filterProbes(fltGr)

    save(detP, normGr, fltGr, targets, file = saveAs)
    if(saveRG) save(rawRg, file = glue("{dirname(saveAs)}/RG.{basename(saveAs)}"))
}

getBloodCellData <- function(){
    # Function download the flow-sorted blood cell data from ExperimentHub at
    # Bioconductor and return RGset and sample information table as a list.
    # Jovana Maksimovic

    hub <- ExperimentHub::ExperimentHub()
    AnnotationHub::query(hub, "FlowSorted.Blood.EPIC")
    FlowSorted.Blood.EPIC <- hub[["EH1136"]]

    # get all the sorted cell samples
    rawRg <- FlowSorted.Blood.EPIC[,FlowSorted.Blood.EPIC$CellType != "MIX"]

    targets <- colData(rawRg)

    # give the samples descriptive names
    sampleNames(rawRg) <- targets$Sample_Name

    return(list(rawRg = rawRg, targets = targets))
}

filterQual <- function(normGr, detP, pval=0.01){
    # Function to remove probes for which the detection p-value > 0.01 in at
    # least ONE sample.
    # Jovana Maksimovic

    # ensure probes are in the same order in both objects
    detP <- detP[match(rownames(normGr),rownames(detP)),]

    # remove any probes that have failed in one or more samples
    keep <- rowSums(detP < pval) == ncol(normGr)
    fltGr <- normGr[keep,]
}

filterProbes <- function(datGr, dist=2, mafcut=0, and=TRUE, rmcrosshyb = TRUE,
                         rmXY=TRUE){
    # Function to remove SNP, cross-hybridising and sex chromosome probes.
    # Jovana Maksimovic

    # filter out SNP probes and multi-mapping probes
    keep <- rownames(DMRcate::rmSNPandCH(getM(datGr), dist = dist,
                                         mafcut = mafcut, and = and,
                                         rmcrosshyb = rmcrosshyb, rmXY=rmXY))
    fltGr <- datGr[rownames(datGr) %in% keep,]
    fltGr
}

cpgsEgGoFreqs <- function(flatAnn){
    # Function count the number of CpGs per gene per GO category.
    # Jovana Maksimovic

    cpgsEg <- table(flatAnn$entrezid)
    egGo <- reshape2::melt(missMethyl:::.getGO()$idList)
    colnames(egGo) <- c("ENTREZID","GO")

    m1 <- match(egGo$ENTREZID, names(cpgsEg))
    cpgsEgGo <- data.frame(egGo[!is.na(m1),], cpgsEg[m1[!is.na(m1)]])
    colnames(cpgsEgGo)[3:4] <- c("ENTREZID.","Freq")

    gene.go <- paste(cpgsEgGo$ENTREZID, cpgsEgGo$GO, sep=".")
    d <- duplicated(gene.go)
    cpgsEgGo <- cpgsEgGo[!d,]
    cpgsEgGo
}

gsaseq <- function(sig.de, universe, collection, plot.bias=FALSE,
                   gene.length=NULL, sort = TRUE)
    # Generalised version of goana with user-specified gene sets
    # Gene sets collections must be Entrez Gene ID
    # Can take into account gene length bias
    # Belinda Phipson
    # 4 May 2020
{

    if(!is.vector(sig.de))
        stop("Input DE list is not a character vector")

    if(is.null(universe)){
        stop("Please supply the universe: all genes analysed in the experiment")
    }

    if(!is.null(gene.length)){
        if(length(gene.length)!=length(universe)){
            stop("Gene length and universe must be the same length")
        }
    }

    # Check collection is a list with character vectors
    if(!is.list(collection))
        collection <- list(collection=collection)
    collection <- lapply(collection, as.character)
    # Make sure gene set collections don't have any NAs
    collection <- lapply(collection, function(x) x[!is.na(x)])
    # Remove genes that are not in the universe from collections
    collection <- lapply(collection, function(x) x[x %in% universe])
    # Remove collections with no genes left after universe filter
    inUniv <- sapply(collection, function(x) length(x) > 0)
    collection <- collection[inUniv]

    test.de <- as.numeric(universe %in% sig.de)

    # Estimate prior probabilities
    if(!is.null(gene.length)){
        pwf <- .estimatePWF(D=test.de,bias=gene.length)
        if(plot.bias)
            .plotBiasSeq(D=test.de,bias=as.vector(gene.length))
    }

    results <- matrix(NA,ncol=4,nrow=length(collection))
    colnames(results) <- c("N","DE","P.DE","FDR")
    rownames(results) <- names(collection)
    results[,"N"] <- unlist(lapply(collection,length))
    results[,"DE"] <- unlist(lapply(collection, function(x) sum(sig.de %in% x)))
    Nuniverse <- length(universe)
    m <- length(sig.de)

    # Hypergeometric test with prior probabilities
    if(!is.null(gene.length)){
        for(i in 1:length(collection)){
            InSet <- universe %in% collection[[i]]
            pw.red <- sum(pwf[InSet])/results[i,"N"]
            pw.white <- sum(pwf[!InSet])/(Nuniverse-results[i,"N"])
            odds <- pw.red/pw.white
            results[i,"P.DE"] <- BiasedUrn::pWNCHypergeo(results[i,"DE"],
                                                         results[i,"N"],
                                                         Nuniverse-results[i,"N"],
                                                         m,
                                                         odds,
                                                         lower.tail=FALSE) +
                BiasedUrn::dWNCHypergeo(results[i,"DE"],
                                        results[i,"N"],
                                        Nuniverse-results[i,"N"],
                                        m, odds)
        }
    }
    # Hypergeometric test without prior probabilities
    else{
        for(i in 1:length(collection)){
            results[i,"P.DE"] <- stats::phyper(q=results[i,"DE"]-0.5, m=m,
                                               n=Nuniverse-m,
                                               k=results[i,"N"],
                                               lower.tail=FALSE)
        }
    }
    results[,"FDR"] <- stats::p.adjust(results[,"P.DE"],method="BH")
    if(sort){
        o <- order(results[,"P.DE"])
        results[o,]
    }
    else
        data.frame(results)
}

.estimatePWF <- function(D,bias)
    # An alternative to goseq function nullp, which is transformation invariant
    # Belinda Phipson and Gordon Smyth
    # 6 March 2015
{
    prior.prob <- bias
    o <- order(bias)
    prior.prob[o] <- limma::tricubeMovingAverage(D[o],span=0.5)
    prior.prob
}

.plotBiasSeq <- function(D,bias)
    # Plotting function to show gene level CpG density bias
    # Belinda Phipson
    # 5 March 2015
{
    o <- order(bias)
    splitf <- rep(1:100,each=200)[1:length(bias)]
    avgbias <- tapply(bias[o],factor(splitf),mean)
    sumDM <- tapply(D[o],factor(splitf),sum)
    propDM <- sumDM/table(splitf)
    graphics::par(mar=c(5,5,2,2))
    graphics::plot(avgbias,as.vector(propDM),
                   xlab="Gene length (binned)",
                   ylab="Proportion Differential Expression",cex.lab=1.5,
                   cex.axis=1.2)
    graphics::lines(stats::lowess(avgbias,propDM),col=4,lwd=2)
}

getBiasDat <- function (sig.cpg, all.cpg = NULL, collection,
                        array.type = c("450K", "EPIC"), plot.bias = FALSE,
                        prior.prob = TRUE, anno = NULL, equiv.cpg = TRUE,
                        fract.counts = TRUE)
{
    # Function to get the probe-bias data for plotting using ggplot2.
    # Jovana Maksimovic

    if (!is.vector(sig.cpg))
        stop("Input CpG list is not a character vector")
    array.type <- match.arg(toupper(array.type), c("450K", "EPIC"))
    if (!is.null(anno)) {
        out <- getMappedEntrezIDs(sig.cpg = sig.cpg, all.cpg = all.cpg,
                                  array.type = array.type, anno)
    }
    else {
        out <- getMappedEntrezIDs(sig.cpg = sig.cpg, all.cpg = all.cpg,
                                  array.type = array.type)
    }
    sorted.eg.sig <- out$sig.eg
    eg.universe <- out$universe
    freq_genes <- out$freq
    test.de <- out$de
    frac <- out$fract.counts
    equiv <- out$equiv

    # data for bias plot
    D <- test.de
    bias <- as.vector(equiv)
    o <- order(bias)
    splitf <- rep(1:100, each = 200)[1:length(bias)]
    avgbias <- tapply(bias[o], factor(splitf), mean)
    sumDM <- tapply(D[o], factor(splitf), sum)
    propDM <- sumDM/table(splitf)

    # save data to make bias ggplot
    data.frame(avgbias = avgbias, propDM = as.vector(propDM),
               stringsAsFactors = FALSE)
}

shift_legend <- function(p, pos = "center", plot = FALSE) {
    # Function to move ggplot2 legend into empty panel.
    # Jovana Maksimovic

    pnls <- cowplot::plot_to_gtable(p) %>%
        gtable::gtable_filter("panel") %>%
        with(setNames(grobs, layout$name)) %>%
        purrr::keep(~identical(.x,zeroGrob()))

    if(length(pnls) == 0){
        p
    } else {
        lemon::reposition_legend(p, pos, panel=names(pnls), plot = plot)

    }

}

# colour palette used for manuscript figures
methodPal <- c("#a0e85b",
               "#154e56",
               "#61cab8",
               "#218841",
               "#bfd6fa",
               "#333a9e",
               "#5794d7",
               "#996ddb",
               "#a113b2",
               "#e9b4f5",
               "#a1085c",
               "#f75ef0")

# method names used for manuscript figures
dict <- c("mgsa.glm" = "mGLM",
          "mgsa.ora" = "mRRA (ORA)",
          "mgsa.gsea" = "mRRA (GSEA)",
          "champ.wt" = "ebGSEA (WT)",
          "champ.kpmt" = "ebGSEA (KPMT)",
          "hgt" = "HGT",
          "hgt.cpg" = "HGT-mod",
          "hgt.cpg.fc.ec" = "GOmeth",
          "goregion-gometh" = "GOregion",
          "goana" = "HGT",
          "mmethyl.gm1000" = "GOmeth (1000)",
          "mmethyl.gm5000" = "GOmeth (5000)",
          "mmethyl.gometh" = "GOmeth",
          "mmethyl.hgt" = "HGT",
          "gometh" = "GOmeth",
          "champ.ebgsea" = "ebGSEA",
          "mgsa.glm.3" = "mGLM",
          "mgsa.glm.6" = "mGLM",
          "mgsa.glm.3" = "mGLM",
          "gometh-probe-top" = "GOmeth (5000)",
          "gometh-probe-fdr" = "GOmeth (FDR < 0.05)",
          "glm" = "mGLM",
          "ora" = "mRRA (ORA)",
          "gsea" = "mRRA (GSEA)")

# colour labelling dictionary to use with ggplot2
methodCols <- methodPal[c(1,2,3,5,6,7,8,12,9,7,11,12,12,7,
                          12,5,1,1,1,12,10,1,2,3)]
names(methodCols) <- unname(dict)

find_modes<- function(x,adjust=1)
    # modified from http://stackoverflow.com/questions/27418461/calculate-the-modes-in-a-multimodal-distribution-in-r
    # Written by Dr. Belinda Phipson
{
    dens <- density(x,adjust=adjust)
    dy <- dens$y
    modes <- NULL
    for ( i in 2:(length(dy)-1) ){
        if ( (dy[i] > dy[i-1]) & (dy[i] > dy[i+1]) ) {
            modes <- c(modes,i)
        }
    }
    if ( length(modes) == 0 ) {
        return(NA)
    }
    return(dens$x[modes])
}

find_num_modes <- function(x,adjust=1)
    # for a given set of data, returns the number of modes in the data
    # Written by Dr. Belinda Phipson
{
    dy <- density(x,adjust=adjust)$y
    modes <- NULL
    for ( i in 2:(length(dy)-1) ){
        if ( (dy[i] > dy[i-1]) & (dy[i] > dy[i+1]) ) {
            modes <- c(modes,i)
        }
    }
    return(length(modes))
}
