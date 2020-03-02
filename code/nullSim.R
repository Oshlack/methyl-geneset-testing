#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(here)
library(glue)

gsetTestRandCpgs <- function(cpgNum, anno, simNum){

  require(missMethyl)

  results <- vector("list", simNum)

  if(nrow(anno) > 500000){
    array.type <- "EPIC"
  } else {
    array.type <- "450k"
  }

  for(i in 1:simNum){
    cat(glue("{i} "),file = here(glue("output/nullSim.{array.type}.{cpgNum}.{simNum}.log")),
        append = TRUE)

    randSet <- sample(x = anno$Name, size = cpgNum) # randomly sample cpgNum sites
    # GO testing without ANY adjustment
    none <- gometh(sig.cpg = randSet, array.type = array.type, anno = anno,
                    prior.prob = FALSE, frac.count = FALSE, equiv.cpg = FALSE)
    # GO testing with adjustment for no. cpgs
    cpg <- gometh(sig.cpg = randSet, array.type = array.type, anno = anno,
                  prior.prob = TRUE, frac.count = FALSE, equiv.cpg = FALSE)
    cpg.fc <- gometh(sig.cpg = randSet, array.type = array.type, anno = anno,
                  prior.prob = TRUE, frac.count = TRUE, equiv.cpg = FALSE)
    cpg.ec <- gometh(sig.cpg = randSet, array.type = array.type, anno = anno,
                        prior.prob = TRUE, frac.count = FALSE, equiv.cpg = TRUE)
    cpg.fc.ec <- gometh(sig.cpg = randSet, array.type = array.type, anno = anno,
                    prior.prob = TRUE, frac.count = TRUE, equiv.cpg = TRUE)

    results[[i]] <- list(none=none,
                          cpg=cpg,
                          cpg.fc=cpg.fc,
                          cpg.ec=cpg.ec,
                          cpg.fc.ec=cpg.fc.ec)
  }

  results
}

if (length(args) != 2) {
  stop("Usage: Rscript nullSim.R {450k or EPIC} {no. sims}\n", call.=FALSE)

} else if (length(args) == 2) {

  if(args[1] %in% c("450k","EPIC")){

    #library(missMethyl)
    #library(minfi)
    library(BiocParallel)

    if(args[1] == "EPIC"){
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

    } else {
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      ann <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

    }

    # run sims on several different length CpG sets
    cpgNums <- c(50,100,500,1000,5000,10000)

    # set seed for reproducibility
    set.seed(42)

    bpparam <- MulticoreParam(length(cpgNums), RNGseed = 42)
    nullSim <- bplapply(cpgNums, gsetTestRandCpgs, anno=ann, simNum=args[2],
                        BPPARAM = bpparam)
    names(nullSim) <- cpgNums

    # save results to file
    save(nullSim, file = here(glue("output/nullSim{args[1]}.{args[2]}.RData")))

  } else {
    stop("Usage: Rscript nullSim.R {450k or EPIC} {no. sims}\n", call.=FALSE)

  }
}
