#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
package <- args[1]
set <- args[2]
input <- args[3]
outDir <- args[4]

#library(org.Hs.eg.db)
source(here::here("code/utility.R"))
library(EnsDb.Hsapiens.v75)

#  set <- "GO"
#  package <- "MethylGSA"
#  input <- here::here("data/cache-intermediates/blood.contrasts.rds")
# # #input <- here::here("data/cache-intermediates/cancer.cellline.contrasts.rds")
#  outDir <- here::here(glue::glue("output/compare-methods/BLOOD"))
# # #outDir <- here::here(glue::glue("output/compare-methods/CANCER"))

obj <- readRDS(input)
obj$tfit -> tfit
obj$maxsize -> maxsize
obj$minsize -> minsize
obj$mVals -> mVals
obj$targets -> targets

edb <- EnsDb.Hsapiens.v75
ensGenes <- ensembldb::genes(edb, columns = c("gene_id", "symbol", "entrezid"),
                             return.type = "DataFrame")
ensGenes$entrezid <- sapply(ensGenes$entrezid, function(x) x[1])

pvalList <- lapply(colnames(tfit$contrasts), function(coef){
  top <- limma:::topTreat(tfit, coef = coef, num = Inf)
  pval <- top$P.Value
  names(pval) <- rownames(top)
  pval
})
names(pvalList) <- colnames(tfit$contrasts)

if(nrow(mVals) > 500000) arrayType <- "EPIC" else arrayType <- "450K"
ann <- loadAnnotation(arrayType=arrayType)

if(set == "KEGG"){
  collection <- missMethyl:::.getKEGG()$idList
  # collection <- suppressMessages(lapply(collection, function(x){
  #   AnnotationDbi::select(org.Hs.eg.db, x, columns = "SYMBOL",
  #                         keytype = "ENTREZID")$SYMBOL
  # }))
  collection <- suppressMessages(lapply(collection, function(x){
    tmp <- ensGenes$symbol[ensGenes$entrezid %in% x]
    tmp[!is.na(tmp)]
  }))

} else if (set == "GO"){
  collection <- missMethyl:::.getGO()$idList
  # collection <- suppressMessages(lapply(collection, function(x){
  #   AnnotationDbi::select(org.Hs.eg.db, x, columns = "SYMBOL",
  #                         keytype = "ENTREZID")$SYMBOL
  # }))
  collection <- suppressMessages(lapply(collection, function(x){
    tmp <- ensGenes$symbol[ensGenes$entrezid %in% x]
    tmp[!is.na(tmp)]
  }))

} else if (set == "BROAD"){
  library(ChAMPdata)
  data("PathwayList")
  #k <- keys(org.Hs.eg.db, keytype = "SYMBOL")
  keep <- sapply(PathwayList, function(x) any(x %in% ensGenes$symbol))
  collection <- PathwayList[keep]

}

res <- lapply(pvalList, function(pval){

  tmp <- vector("list", 3)
  names(tmp) <- c("mgsa.glm", "mgsa.ora", "mgsa.gsea")

  tmp[[1]] <- methylGSA::methylglm(cpg.pval = pval, #parallel = TRUE, BPPARAM = MulticoreParam(20),
              FullAnnot = ann, minsize = minsize, maxsize = maxsize,
              GS.list = collection, GS.idtype = "SYMBOL")[, c("ID","pvalue")]
  tmp[[2]] <- methylGSA::methylRRA(cpg.pval = pval,
              method = "ORA", FullAnnot = ann, minsize = minsize,
              maxsize = maxsize, GS.list = collection,
              GS.idtype = "SYMBOL")[, c("ID","pvalue")]
  tmp[[3]] <- methylGSA::methylRRA(cpg.pval = pval,
              method = "GSEA", FullAnnot = ann, minsize = minsize,
              maxsize = maxsize, GS.list = collection,
              GS.idtype = "SYMBOL")[, c("ID","pvalue")]

  dplyr::bind_rows(tmp, .id = "method")

})
res <- dplyr::bind_rows(res, .id = "contrast")
res <- tibble::add_column(res, sub = "n", .after = 1)
res$set <- set

out <- glue::glue("{outDir}/{package}.{set}.rds")
saveRDS(res, out)
