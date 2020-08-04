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

res <- lapply(1:ncol(tfit), function(i){

    tmp <- vector("list", 2)
    names(tmp) <- c("champ.wt", "champ.kpmt")

    cellType <- names(tfit$contrasts[,i])[tfit$contrasts[,i] != 0]
    ebgs <- data.frame(champ.ebGSEA(beta = mVals[,targets$CellType %in% cellType],
                       pheno = targets$CellType[targets$CellType %in% cellType],
                       minN = 5, adjPval=1, arraytype = "EPIC")[[1]])

    tmp[[1]] <- tibble::rownames_to_column(ebgs, var = "ID")[, c("ID", "P.WT.")]
    colnames(tmp[[1]])[2] <- "pvalue"
    tmp[[2]] <- tibble::rownames_to_column(ebgs, var = "ID")[, c("ID", "P.KPMT.")]
    colnames(tmp[[2]])[2] <- "pvalue"
    dplyr::bind_rows(tmp, .id = "method")

})

names(res) <- colnames(tfit$contrasts)
res <- dplyr::bind_rows(res, .id = "contrast")
res <- tibble::add_column(res, sub = "n", .after = 1)
res$set <- "BROAD"

out <- glue::glue("{outDir}/{package}.{set}.rds")
saveRDS(res, out)
