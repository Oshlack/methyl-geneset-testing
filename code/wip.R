#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
package <- args[1]
set <- args[2]
input <- args[3]

library(org.Hs.eg.db)
library(ChAMP)
source(here::here("code/utility.R"))

obj <- readRDS(input)
obj$tfit -> tfit
obj$maxsize -> maxsize
obj$minsize -> minsize
obj$mVals -> mVals
obj$targets -> targets

pvalList <- lapply(colnames(tfit$contrasts), function(coef){
  top <- limma:::topTreat(tfit, coef = coef, num = Inf)
  pval <- top$P.Value
  names(pval) <- rownames(top)
  pval
})
names(pvalList) <- colnames(tfit$contrasts)

ann <- loadAnnotation(arrayType="EPIC")

if(set == "KEGG"){
  collection <- missMethyl:::.getKEGG()$idList
  collection <- suppressMessages(lapply(collection, function(x){
    AnnotationDbi::select(org.Hs.eg.db, x, columns = "SYMBOL",
                          keytype = "ENTREZID")$SYMBOL
  }))

} else if (set == "GO"){
  collection <- missMethyl:::.getGO()$idList
  collection <- suppressMessages(lapply(collection, function(x){
    AnnotationDbi::select(org.Hs.eg.db, x, columns = "SYMBOL",
                          keytype = "ENTREZID")$SYMBOL
  }))

} else if (set == "BROAD"){
  data("PathwayList")
  k <- keys(org.Hs.eg.db, keytype = "SYMBOL")
  keep <- sapply(PathwayList,function(x) any(x %in% k))
  collection <- PathwayList[keep]

}

res <- lapply(pvalList, function(pval){

  tmp <- vector("list", 3)
  names(tmp) <- c("mgsa.glm", "mgsa.ora", "mgsa.gsea")

  tmp[[1]] <- tibble::as_tibble(methylGSA::methylglm(cpg.pval = pval,
              FullAnnot = ann, minsize = minsize, maxsize = maxsize,
              GS.list = collection, GS.idtype = "SYMBOL", parallel = TRUE,
              BPPARAM = BiocParallel::MulticoreParam(min(BiocParallel::bpworkers(),
                                                         10)))[, c("ID","pvalue")])
  tmp[[2]] <- tibble::as_tibble(methylGSA::methylRRA(cpg.pval = pval,
              method = "ORA", FullAnnot = ann, minsize = minsize,
              maxsize = maxsize, GS.list = collection,
              GS.idtype = "SYMBOL")[, c("ID","pvalue")])
  tmp[[3]] <- tibble::as_tibble(methylGSA::methylRRA(cpg.pval = pval,
              method = "GSEA", FullAnnot = ann, minsize = minsize,
              maxsize = maxsize, GS.list = collection,
              GS.idtype = "SYMBOL")[, c("ID","pvalue")])
  dplyr::bind_rows(tmp, .id = "method")

})

res <- dplyr::bind_rows(res, .id = "contrast")

out <- glue::glue("{outDir}/{package}.{set}.rds")
saveRDS(res, out)
