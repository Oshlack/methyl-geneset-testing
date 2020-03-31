#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
package <- args[1]
set <- args[2]
input <- args[3]
outDir <- args[4]

library(org.Hs.eg.db)
library(ChAMP)
source(here::here("code/utility.R"))

obj <- readRDS(input)
obj$tfit -> tfit
obj$maxsize -> maxsize
obj$minsize -> minsize
obj$mVals -> mVals
obj$targets -> targets

contList <- lapply(colnames(tfit$contrasts), function(coef){

    tmp <- NULL
    tmp$c1 <- rownames(limma::topTreat(tfit, coef = coef, num = 5000))
    tmp$c2 <- rownames(limma::topTreat(tfit, coef = coef, num = 10000))
    tmp$p1 <- rownames(limma::topTreat(tfit, coef = coef, p.value = 0.01, num = Inf))
    tmp$p2 <- rownames(limma::topTreat(tfit, coef = coef, p.value = 0.05, num = Inf))
    tmp
})
names(contList) <- colnames(tfit$contrasts)

ann <- loadAnnotation(arrayType="EPIC")

if(set == "KEGG"){
  collection <- missMethyl:::.getKEGG()$idList

} else if (set == "GO"){
  collection <- missMethyl:::.getGO()$idList

} else if (set == "BROAD"){
  data("PathwayList")
  k <- keys(org.Hs.eg.db, keytype = "SYMBOL")
  keep <- sapply(PathwayList,function(x) any(x %in% k))
  collection <- suppressMessages(lapply(PathwayList[keep], function(x){
      tmp <- select(org.Hs.eg.db, x, columns = "ENTREZID",
                    keytype = "SYMBOL")$ENTREZID
      tmp[!is.na(tmp)]
  }))
}

universe <- rownames(mVals)

res <- lapply(contList, function(cont){
    tmp <- lapply(cont, function(sig){
        tmp <- vector("list", 2)
        names(tmp) <- c("mmethyl.hgt", "mmethyl.gometh")

        tmp[[1]] <- tibble::rownames_to_column(missMethyl::topGSA(missMethyl::gsameth(sig.cpg = sig,
              all.cpg = universe, collection = collection,
              array.type = "EPIC", anno=ann, prior.prob = FALSE,
              equiv.cpg = FALSE), number = Inf), var = "ID")[, c("ID","P.DE")]
        tmp[[2]] <- tibble::rownames_to_column(missMethyl::topGSA(missMethyl::gsameth(sig.cpg = sig,
              all.cpg = universe, collection = collection, fract.counts = TRUE,
              array.type = "EPIC", anno=ann, prior.prob = TRUE,
              equiv.cpg = TRUE), number = Inf), var = "ID")[, c("ID","P.DE")]

        dplyr::bind_rows(tmp, .id = "method")

    })
    dplyr::bind_rows(tmp, .id = "sub")

})
res <- dplyr::bind_rows(res, .id = "contrast")
colnames(res)[5] <- "pvalue"
res$set <- set

out <- glue::glue("{outDir}/{package}.{set}.rds")
saveRDS(res, out)
