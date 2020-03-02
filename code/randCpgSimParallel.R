#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(here)
library(glue)

gsetTestRandCpgs <- function(cpgNum, anno, simNum){

  randTest <- vector("list", simNum)

  for(i in 1:simNum){
    cat(glue("{i} "),file = here(glue("output/sim{cpgNum}.{simNum}.log")), append = TRUE)

    randSet <- sample(x = anno$Name, size = cpgNum) # randomly sample cpgNum sites
    # GO testing without ANY adjustment
    none <- gometh(sig.cpg = randSet, array.type = "EPIC", anno = anno,
                    prior.prob = FALSE, frac.count = FALSE, equiv.cpg = FALSE)
    # GO testing with adjustment for no. cpgs
    cpg <- gometh(sig.cpg = randSet, array.type = "EPIC", anno = anno,
                  prior.prob = TRUE, frac.count = FALSE, equiv.cpg = FALSE)
    cpg.fc <- gometh(sig.cpg = randSet, array.type = "EPIC", anno = anno,
                  prior.prob = TRUE, frac.count = TRUE, equiv.cpg = FALSE)
    cpg.ec <- gometh(sig.cpg = randSet, array.type = "EPIC", anno = anno,
                        prior.prob = TRUE, frac.count = FALSE, equiv.cpg = TRUE)
    cpg.fc.ec <- gometh(sig.cpg = randSet, array.type = "EPIC", anno = anno,
                    prior.prob = TRUE, frac.count = TRUE, equiv.cpg = TRUE)
    randTest[[i]] <- list(none=none,
                          cpg=cpg,
                          cpg.fc=cpg.fc,
                          cpg.ec=cpg.ec,
                          cpg.fc.ec=cpg.fc.ec)
  }

  randTest
}

if (length(args) != 2) {
  stop("Usage: Rscript randCpgSim.R {no. CpGs} {no. sims}\n", call.=FALSE)

} else if (length(args) == 2) {

  library(missMethyl)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(BiocParallel)
  annoEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

  # run sims on several different length CpG sets
  cpgNumSet <- c(100,500,1000,5000,10000)
  
  # set seed for reproducibility
  set.seed(42)
  #randCpgSim <- gsetTestRandCpgs(cpgNum = args[1], anno = annoEPIC, simNum = args[2])
  bpparam <- MulticoreParam(length(cpgNumSet), RNGseed = 42)
  nullSim <- bplapply(cpgNumSet, gsetTestRandCpgs, anno=ann, simNum=100, BPPARAM = bpparam)
  
  # save results to file
  save(nullSim, file = here(glue("output/nullSim{args[1]}.{args[2]}.RData")))

}
