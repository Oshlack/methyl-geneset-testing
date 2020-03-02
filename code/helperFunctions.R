readData <- function(saveAs, dataSrc = c("BC","RA"), saveRG = FALSE){

  require(here)
  require(minfi)
  #require(ExperimentHub)

  dataSrc <- match.arg(dataSrc)

  #hub <- ExperimentHub()
  #query(hub, "FlowSorted.Blood.EPIC")
  #FlowSorted.Blood.EPIC <- hub[["EH1136"]]

  # get all the sorted cell samples
  #rawRg <- FlowSorted.Blood.EPIC[,FlowSorted.Blood.EPIC$CellType != "MIX"]

  #targets <- colData(rawRg)
  # give the samples descriptive names
  #sampleNames(rawRg) <- targets$Sample_Name

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
  normGr <- preprocessQuantile(rawRg)

  # filtering
  fltGr <- filterQual(normGr, detP)
  fltGr <- filterProbes(fltGr)

  save(detP, normGr, fltGr, targets, file = saveAs)
  if(saveRG) save(rawRg, file = glue("{dirname(saveAs)}/RG.{basename(saveAs)}"))
}

getBloodCellData <- function(){
  require(minfi)
  require(ExperimentHub)

  hub <- ExperimentHub()
  query(hub, "FlowSorted.Blood.EPIC")
  FlowSorted.Blood.EPIC <- hub[["EH1136"]]

  # get all the sorted cell samples
  rawRg <- FlowSorted.Blood.EPIC[,FlowSorted.Blood.EPIC$CellType != "MIX"]

  targets <- colData(rawRg)

  # give the samples descriptive names
  sampleNames(rawRg) <- targets$Sample_Name

  return(list(rawRg = rawRg), targets = targets)
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

loadAnnotation <- function(arrayType = c("450k", "EPIC")) {
  require(here)
  require(glue)

  arrayType <- match.arg(arrayType)
  annFile <- here(glue("data/ann{arrayType}.RData"))

  if(file.exists(annFile)){
    # load annotation from file for faster execution
    load(annFile)
  } else {
    if(arrayType == "EPIC") {
      require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    } else {
      require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      ann <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    }
    # save annotation to file to speed up subsequent execution
    save(ann, file=annFile)
  }

  ann
}

loadFlatAnnotation <- function(ann) {
  require(here)
  require(glue)
  require(missMethyl)

  arrayType <- ifelse(nrow(ann) > 500000, "EPIC", "450k")

  flatFile <- here(glue("data/flatAnn{arrayType}.RData"))

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

cpgsByGO <- function(flatAnn){

  require(AnnotationDbi)
  require(org.Hs.eg.db)

  cpgsEg <- table(flatAnn$entrezid)
  egGo <- AnnotationDbi::select(org.Hs.eg.db, keys = names(cpgsEg),
               columns = c("ENTREZID","GO"))
  egGo <- unique(egGo[!is.na(egGo$GO),c("ENTREZID","GO")])

  m1 <- match(egGo$ENTREZID, names(cpgsEg))
  cpgsGo <- data.frame(egGo[!is.na(m1),], cpgsEg[m1[!is.na(m1)]])

  cpgsGo$dup <- paste(cpgsGo$ENTREZID, cpgsGo$GO, sep=".")
  d <- duplicated(cpgsGo$dup)
  cpgsGo <- cpgsGo[!d,]
  cpgsGo
}

makeLong <- function(nullSim){

  simLong <- numeric(0)

  for(cpgNum in names(nullSim)){

    cpgSim <- nullSim[[cpgNum]]

    for(i in 1:length(cpgSim)){
      sim <- cpgSim[[i]]
      tmp <- sapply(sim, function(x) x$P.DE)
      rownames(tmp) <- rownames(sim[[1]])
      long <- melt(tmp)
      colnames(long) <- c("id","adj","pval")
      long$sim <- rep(i, nrow(long))
      simLong <- rbind(simLong, long)
    }
  }
  each <- nrow(simLong)/length(nullSim)
  simLong$cpg <- rep(as.numeric(names(nullSim)), each=each)
  simLong
}

postProcess <- function(sampleNums, prefix = "nullDat", subset = NULL){

  require(here)
  require(glue)

  pnames <- c("P.DE","P.DE","pvalue","pvalue","pvalue","P(WT)")
  longDat <- numeric(0)
  samps <- numeric(0)

  for(s in sampleNums){
    load(here(glue("output/{prefix}.{s}.RData")))
    allSims <- get(paste0(prefix, ".", s))
    simCount <- 0

    for(sim in allSims){
      simCount <- simCount + 1
      if(is.null(subset)) subset <- 1:length(sim)

      pvals <- lapply(subset, function(j){
        if(!grepl("(",pnames[j], fixed = TRUE)){
          sim[[j]][,pnames[j]]
        } else {
          unname(sim[[j]][[1]][,pnames[j]])
        }
      })

      names(pvals) <- names(sim)[subset]
      methods <- lapply(1:length(pvals), function(i){
        rep(names(pvals)[i], length(pvals[[i]]))
      })

      tmp <- data.frame(pval=unlist(pvals), method=unlist(methods))
      tmp$sim <- rep(simCount, nrow(tmp))
      samps <- c(samps, rep(s, nrow(tmp)))
      longDat <- rbind(longDat, tmp)
    }
  }

  longDat$samps <- samps
  save(longDat, file = here(glue("output/{prefix}.long.RData")))
}

getGOByKey <- function(keytype = c("ENTREZID","SYMBOL"),
                       minsize=NULL, maxsize=NULL){
  require(GO.db)
  require(org.Hs.eg.db)

  keytype <- match.arg(keytype)
  keys <- keys(org.Hs.eg.db, keytype = keytype)
  GeneID.PathID <- suppressMessages(select(org.Hs.eg.db, keys=keys,
                                           columns=c(keytype,"GO","ONTOLOGY"),
                                           keytype=keytype))
  d <- !duplicated(GeneID.PathID[, c(keytype, "GO")])
  GeneID.PathID <- GeneID.PathID[d, c(1,2,4)]
  GOID.TERM <- suppressMessages(select(GO.db, keys=unique(GeneID.PathID$GO),
                                       columns=c("GOID","ONTOLOGY","TERM"),
                                       keytype="GOID"))
  go <- tapply(GeneID.PathID[,keytype], GeneID.PathID$GO, list)

  len <- sapply(go,length)
  minsize <- ifelse(is.null(minsize), 1, minsize)
  mazsize <- ifelse(is.null(maxsize), max(len), maxsize)
  go <- go[which(len >= minsize & len <= maxsize)]

  go
}

getCollection <- function(collection=c("GO","KEGG"), keytype=c("ENTREZID","SYMBOL")){

  keytype <- match.args(keytype)

  if(keytype == "ENTREZID"){
    goEntrez <- getGOAnyKey(keytype = "ENTREZID")

  } else {
    goSymbol <- getGOAnyKey(keytype = "SYMBOL")

  }
}

addDetails <- function(result, collection=c("GO","KEGG"), pname=c("P.DE","pvalue")){

  collection <- match.arg(toupper(collection),c("GO","KEGG"))
  pname <- match.arg(pname)

  if(collection == "GO"){
    go <- missMethyl:::.getGO()
    result <- merge(go$idTable,result,by.x="GOID",by.y="row.names")
    rownames(result) <- result$GOID

  } else if(collection == "KEGG"){
    kegg <- missMethyl:::.getKEGG()
    result <- merge(kegg$idTable,result,by.x="PathwayID",by.y="row.names")
    rownames(result) <- result$PathwayID
  }

  result <- result[order(result[,pname]),]
  result[,-1]
}

# runBgAnalysis <- function(scriptFile, method, arrayType, set, inputFile){
#   out <- glue("{method}.{arrayType}.{set}")
#   system2(command = glue("{Sys.getenv('R_HOME')}/bin/Rscript"),
#           arg=c(scriptFile, method, arrayType, set, inputFile), wait = FALSE,
#           stdout = here(glue("output/{out}.out")),
#           stderr = here(glue("output/{out}.err")))
# }

runBgAnalysis <- function(scriptFile, outPrefix, args){

  if(is.null(args)) {
    arg <- scriptFile
  } else {
    arg <- c(scriptFile, args)
  }

  system2(command = glue("{Sys.getenv('R_HOME')}/bin/Rscript"),
          arg = arg, wait = FALSE,
          stdout = here(glue("output/{outPrefix}.out")),
          stderr = here(glue("output/{outPrefix}.err")))
}

plotTopSets <- function(pval, term, n, method, col = NULL, cex.names = 0.8, main = NULL){

  pvals <- ifelse(pval == 0, min(pval[pval != 0]), pval)

  # par(mar = c(5, 20, 1, 3) + 0.1, oma = c(0, 0, 0, 0))
  mp <- graphics::barplot(rev(-log10(pvals)), names.arg = glue("{rev(term)} (N = {rev(n)})"),
                horiz = TRUE, las = 2, cex.names = cex.names, col = col,
                xlim = c(0, max(-log10(pvals)) + 0.1*max(-log10(pvals))),
                xlab = expression(-log[10](FDR)),
                main = main)
  mtext(method, side = 4, outer = FALSE, line = 1)
  abline(v = -log10(0.05), lty = 2)
  sapply(1:length(pval), function(i){
    if(pval[i] == 0) text(-log10(pvals[i]), rev(mp)[i], labels = "*", pos = 4)
  })
}

