#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

args

library(ChAMP)
library(here)
library(glue)
source(here("code/helperFunctions.R"))

if (length(args) != 7) {
  stop("Usage: Rscript runDmrMethods.R {maxGap} {minProbes} {cores} {B} {adjP} {output file} {input file}\n",
       call.=FALSE)
}

maxGap <- as.numeric(args[1]) # max distance between probes in same cluster
minProbes <- as.numeric(args[2]) # valid cluster must have > minProbes
cores <- as.numeric(args[3]) # num. cores to use for bootstraps
B <- as.numeric(args[4]) # num. bootstraps
adjP <- as.numeric(args[5]) # dmr adjusted p-value cutoff
dmrFile <- args[6] # name of results file
load(args[7]) # saved input objects

dmrList <- vector("list", ncol(tfit))

for(i in 1:ncol(tfit)){

  cellType <- names(tfit$contrasts[,i])[tfit$contrasts[,i] != 0]
  cat(paste(cellType, collapse = " vs. "))
  cat("\nDMRcate\n")

  cpgAnn <- cpg.annotate("array", mVals, what = "M", arraytype = "EPIC",
                         analysis.type = "differential", design = tfit$design,
                         contrasts = TRUE, cont.matrix = tfit$contrasts,
                         coef = colnames(tfit$contrasts)[i])
  dmrOut <- dmrcate(cpgAnn, lambda = 1000, mc.cores = cores, betacutoff = 0.1,
                    min.cpgs = 3)
  dmrc <- as.data.frame(extractRanges(dmrOut))

  cat("\nbumphunter\n")
  bumph <- champ.DMR(bVals[,targets$CellType %in% cellType],
                     pheno = targets$CellType[targets$CellType %in% cellType],
                     arraytype = "EPIC", method = "Bumphunter", #minProbes = minProbes,
                     cores = cores, maxGap = maxGap, B = B,
                     adjPvalDmr = adjP)$BumphunterDMR

  cat("\nProbeLasso\n")
  probel <- champ.DMR(bVals[,targets$CellType %in% cellType],
                     pheno = targets$CellType[targets$CellType %in% cellType],
                     arraytype = "EPIC", method = "ProbeLasso", PDFplot = FALSE,
                     adjPvalDmr = adjP)$ProbeLassoDMR

  dmrList[[i]] <- list(dmrc = dmrc, bumph = bumph, probel = probel)
}

save(dmrList, file = dmrFile)
