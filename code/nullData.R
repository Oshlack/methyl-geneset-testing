#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(here)
library(glue)
library(missMethyl)
library(ChAMP)
library(methylGSA)
library(minfi)
source(here("code/helperFunctions.R"))

# set seed for reproducibility
set.seed(42)

# load(here("data/inputNullDat.RData"))
load(args[1]) # load processed data
sampleNum <- args[2] # get sample number
simNum <- args[3] # number of simulations
#cores <- args[4] # number of cores to use

cat("Cleaning up BROAD lists for analysis...\n")
# set up BROAD gene lists for testing
data("PathwayList")
k <- keys(org.Hs.eg.db, keytype = "SYMBOL")
keep <- sapply(PathwayList,function(x) any(x %in% k))
PathwayList <- PathwayList[keep]
cat("Converting BROAD list Symbols to Entrez IDs...\n")
# convert SYMBOLS to ENTREZID for gometh
collection <- suppressMessages(lapply(PathwayList, function(x){
  tmp <- select(org.Hs.eg.db, x, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
  tmp[!is.na(tmp)]
}))

cat("Starting analysis...\n")
ann <- loadAnnotation(arrayType = "450k")

tmp <- vector("list", simNum)

for(i in 1:simNum){
    cat(paste0("Sim. no. ", i, " of ", simNum, "...\n"))
    set1 <- sample(x = 1:ncol(res), size = sampleNum)
    set2 <- sample(x = (1:ncol(res))[-set1], size = sampleNum)
    subSet <- c(set1,set2)

    design <- model.matrix(~factor(c(rep(1, sampleNum), rep(0, sampleNum))))
    fit <- lmFit(res[,subSet], design)
    fit <- eBayes(fit, robust = TRUE)

    gm1 <- topGSA(gsameth(sig.cpg = rownames(topTable(fit, n = 500, coef = 2)),
                          all.cpg = rownames(res), collection = collection,
                           array.type = "450k", anno=ann, prior.prob = TRUE,
                           equiv.cpg = FALSE), number = Inf)
    gm2 <- topGSA(gsameth(sig.cpg = rownames(topTable(fit, n = 1000, coef = 2)),
                          all.cpg = rownames(res), collection = collection,
                          array.type = "450k", anno=ann, prior.prob = TRUE,
                          equiv.cpg = TRUE), number = Inf)

    glm <- methylglm(cpg.pval = fit$p.value[,2], FullAnnot = ann, minsize = 3,
                      maxsize = 3000, GS.list = PathwayList, GS.idtype = "SYMBOL",
                      parallel = TRUE, BPPARAM = MulticoreParam(5))
    ora <- methylRRA(cpg.pval = fit$p.value[,2], method = "ORA", FullAnnot = ann,
                      minsize = 3, maxsize = 3000, GS.list = PathwayList,
                      GS.idtype = "SYMBOL")
    gsea <- methylRRA(cpg.pval = fit$p.value[,2], method = "GSEA", FullAnnot = ann,
                       minsize = 3, maxsize = 3000, GS.list = PathwayList,
                       GS.idtype = "SYMBOL")

    ebgs <- champ.ebGSEA(beta = res[,subSet], pheno=design[,2], minN=5, adjPval=1,
                          arraytype="450k")

    tmp[[i]] <- list(gm1=gm1,
                      gm2=gm2,
                      glm=glm,
                      ora=ora,
                      gsea=gsea,
                      ebgs=ebgs)
    #tmp[[i]] <- list(gm1=gm1,
    #                 gm2=gm2)
}

# save results to file
cat("Saving results...\n")
out <- paste("nullDat", sampleNum, sep=".")
#out <- paste("nullDat.gm", sampleNum, sep=".")
outFile <- here(glue("output/{out}.RData"))

assign(out, tmp)
save(list=c(eval(out)), file = outFile)
