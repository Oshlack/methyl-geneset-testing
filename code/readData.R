
loadData <- function(saveAs){

  require(here)
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

  # calculate the detection p-values
  detP <- detectionP(rawRg)

  # normalise
  normGr <- preprocessQuantile(rawRg)

  # filtering
  fltGr <- filterQual(normGr, detP)
  fltGr <- filterProbes(fltGr)

  save(detP, normGr, fltGr, targets, file = saveAs)
}

filterQual <- function(normGr, detP, pval=0.01){
  # ensure probes are in the same order in both objects
  detP <- detP[match(rownames(normGr),rownames(detP)),]

  # remove any probes that have failed in one or more samples
  keep <- rowSums(detP < pval) == ncol(normGr)
  fltGr <- normGr[keep,]
}

filterProbes <- function(datGr, dist=2, mafcut=0.05, and=TRUE, rmcrosshyb = TRUE,
                       rmXY=FALSE){
  require(DMRcate)

  # filter out SNP probes and multi-mapping probes
  keep <- rownames(rmSNPandCH(getM(datGr), dist = dist, mafcut = mafcut, and = and,
                             rmcrosshyb = rmcrosshyb, rmXY=rmXY))
  fltGr <- datGr[rownames(datGr) %in% keep,]
  fltGr
}
