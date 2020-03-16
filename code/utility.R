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
  colnames(cpgsEgGo)[3:4] <- c("ENTREZID.","NoCpgs")

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
