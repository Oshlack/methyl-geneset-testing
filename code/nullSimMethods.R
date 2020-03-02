#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(here)
library(glue)

gsetTestNull <- function(sampleNum, anno, simNum, data){

  require(missMethyl)

  results <- vector("list", simNum)

  array.type <- "450k"

  for(i in 1:simNum){
    cat(glue("{i} "),file = here(glue("output/nsMethods.{sampleNum}.{simNum}.log")),
        append = TRUE)

    set1 <- sample(x = 1:ncol(data), size = sampleNum)
    set2 <- sample(x = 1:ncol(data)[-set1], size = sampleNum)
    subSet <- c(set1,set2)

    design <- model.matrix(~factor(c(rep(1,sampleNum), rep(0, sampleNum))))
    fit <- lmFit(data[,subSet], design)
    fit <- eBayes(fit, robust = TRUE)

    gm1 <- gometh(sig.cpg = rownames(topTable(fit, n = 500, coef = 2)),
                  array.type = array.type, anno = anno,
                  prior.prob = TRUE, frac.count = TRUE, equiv.cpg = TRUE)
    gm2 <- gometh(sig.cpg = sig.cpg = rownames(topTable(fit, n = 1000, coef = 2)),
                  array.type = array.type, anno = anno,
                  prior.prob = TRUE, frac.count = TRUE, equiv.cpg = TRUE)
    glm <- methylglm(cpg.pval = fit$p.value[,2], FullAnnot = ann, minsize = 5,
                     maxsize = 5000, GS.list = NULL, GS.idtype = "SYMBOL",
                     parallel = TRUE, BPPARAM = MulticoreParam(10))
    ora <- methylRRA(cpg.pval = fit$p.value[,2], method = "ORA", FullAnnot = ann,
                     minsize = 5, maxsize = 5000, GS.list = NULL,
                     GS.idtype = "SYMBOL")
    gsea <- methylRRA(cpg.pval = fit$p.value[,2], method = "GSEA", FullAnnot = ann,
                      minsize = 5, maxsize = 5000, GS.list = NULL,
                      GS.idtype = "SYMBOL")
    ebgs <- champ.ebGSEA(beta = data[,subSet], pheno=design[,2], minN=5, adjPval=1,
                         arraytype="450k")

    results[[i]] <- list(gm1=gm1,
                         gm2=gm2,
                         glm=glm,
                         ora=ora,
                         gsea=gsea,
                         ebgs=ebgs)
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

    # run sims on different numbers of samples
    sampleNums <- c(5, 10, 25, 50, 100)

    # set seed for reproducibility
    set.seed(42)

    bpparam <- MulticoreParam(length(sampleNums), RNGseed = 42)
    nullSim <- bplapply(cpgNums, gsetTestNull, anno=ann, simNum=args[2],
                        BPPARAM = bpparam)
    names(nullSim) <- sampleNums

    # save results to file
    save(nullSim, file = here(glue("output/nullSim{args[1]}.{args[2]}.RData")))

  } else {
    stop("Usage: Rscript nullSim.R {450k or EPIC} {no. sims}\n", call.=FALSE)

  }
}
