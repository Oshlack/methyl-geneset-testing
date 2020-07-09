loadAnnotation <- function(arrayType = c("450k", "EPIC")) {

  arrayType <- match.arg(arrayType)
  annFile <- here::here(glue::glue("data/ann{arrayType}.RData"))

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

  arrayType <- ifelse(nrow(ann) > 500000, "EPIC", "450k")

  flatFile <- here::here(glue::glue("data/flatAnn{arrayType}.RData"))

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
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

readData <- function(saveAs, dataSrc = c("BC","RA"), saveRG = FALSE){

  dataSrc <- match.arg(dataSrc)

  if(dataSrc == "BC"){
    dat <- getBloodCellData()
    rawRg <- dat$rawRg
    targets <- dat$targets

  } else if (dataSrc == "RA"){
    dat <- getRmArthritisData()
    rawRg <- dat$rawRg
    targets <- dat$targets
  }

  # calculate the detection p-values
  detP <- detectionP(rawRg)

  # normalise
  normGr <- minfi::preprocessQuantile(rawRg)

  # filtering
  fltGr <- filterQual(normGr, detP)
  fltGr <- filterProbes(fltGr)

  save(detP, normGr, fltGr, targets, file = saveAs)
  if(saveRG) save(rawRg, file = glue("{dirname(saveAs)}/RG.{basename(saveAs)}"))
}

cpgsEgGoFreqs <- function(flatAnn){

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

getBloodCellData <- function(){

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

getRheumArthritisData <- function(){
  require(GEOquery)
  require(here)
  require(minfi)

  cat("Downloading data...\n")
  dataFiles <- getGEOSuppFiles(GEO = "GSE42861", baseDir = here("data"),
                               filter_regex = "RAW")
  outPath <- here("data/GSE42861")

  cat("Extracting TAR archive...\n")
  untar(tarfile = rownames(dataFiles), exdir = outPath)

  cat("Downloading series matrix...\n")
  seriesMatrix <- getGEO(GEO = "GSE42861", destdir = outPath,
                         GSEMatrix = TRUE, parseCharacteristics = FALSE)
  pd <- pData(phenoData(seriesMatrix[[1]]))
  tmp <- as.character(pd$supplementary_file)
  tmp <- strsplit2(tmp,"suppl/")[,2]
  pd$id <- strsplit2(tmp,"_Grn")[,1]
  targets <- pd[pd$`subject:ch1` == "Normal",]
  targets$Basename <- paste0(outPath,"/",targets$id)

  cat("Unzipping GZ files...\n")
  idatgz <- c(paste0(targets$Basename,"_Red.idat.gz"),
              paste0(targets$Basename,"_Grn.idat.gz"))
  lapply(idatgz, function(f) gunzip(f))

  cat("Reading IDAT files...\n")
  rawRg <- read.metharray.exp(targets = targets)

  return(list(rawRg = rawRg, targets = targets))
}

filterQual <- function(normGr, detP, pval=0.01){
  # ensure probes are in the same order in both objects
  detP <- detP[match(rownames(normGr),rownames(detP)),]

  # remove any probes that have failed in one or more samples
  keep <- rowSums(detP < pval) == ncol(normGr)
  fltGr <- normGr[keep,]
}

filterProbes <- function(datGr, dist=2, mafcut=0, and=TRUE, rmcrosshyb = TRUE,
                         rmXY=TRUE){
  require(DMRcate)

  # filter out SNP probes and multi-mapping probes
  keep <- rownames(rmSNPandCH(getM(datGr), dist = dist, mafcut = mafcut, and = and,
                              rmcrosshyb = rmcrosshyb, rmXY=rmXY))
  fltGr <- datGr[rownames(datGr) %in% keep,]
  fltGr
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

# .plotBias <- function (D, bias) {
#     o <- order(bias)
#     splitf <- rep(1:100, each = 200)[1:length(bias)]
#     avgbias <- tapply(bias[o], factor(splitf), mean)
#     sumDM <- tapply(D[o], factor(splitf), sum)
#     propDM <- sumDM/table(splitf)
#     par(mar = c(5, 5, 2, 2))
#     plot(avgbias, as.vector(propDM), xlab = "Number of CpGs per gene (binned)",
#          ylab = "Proportion Differential Methylation", cex.lab = 1.5,
#          cex.axis = 1.2)
#     lines(lowess(avgbias, propDM), col = 4, lwd = 2)
# }

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
          "gometh-probe-fdr" = "GOmeth (FDR < 0.05)")

methodCols <- methodPal[c(1,2,3,5,6,7,8,12,9,7,11,12,12,7,
                          12,5,1,1,1,12,10)]
names(methodCols) <- unname(dict)
