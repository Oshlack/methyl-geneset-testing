#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
package <- args[1]
set <- args[2]
input <- args[3]
outDir <- args[4]

#library(org.Hs.eg.db)
library(ebGSEA)
source(here::here("code/utility.R"))

obj <- readRDS(input)
obj$tfit -> tfit
obj$maxsize -> maxsize
obj$minsize -> minsize
#obj$mVals -> mVals
minfi::ilogit2(obj$mVals) -> bVals
obj$targets -> targets

if(nrow(bVals) > 500000) arrayType <- "EPIC" else arrayType <- "450k"

if(set == "KEGG"){
    collection <- missMethyl:::.getKEGG()$idList

} else if (set == "GO"){
    collection <- missMethyl:::.getGO()$idList

} else if (set == "BROAD"){
    library(ChAMPdata)
    library(EnsDb.Hsapiens.v75)

    edb <- EnsDb.Hsapiens.v75
    ensGenes <- ensembldb::genes(edb, columns = c("gene_id", "symbol", "entrezid"),
                                 return.type = "DataFrame")
    ensGenes$entrezid <- sapply(ensGenes$entrezid, function(x) x[1])

    data("PathwayList")
    #k <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    keep <- sapply(PathwayList, function(x) any(x %in% ensGenes$symbol))
    # collection <- suppressMessages(lapply(PathwayList[keep], function(x){
    #     tmp <- select(org.Hs.eg.db, x, columns = "ENTREZID",
    #                   keytype = "SYMBOL")$ENTREZID
    #     tmp[!is.na(tmp)]
    # }))
    collection <- suppressMessages(lapply(PathwayList[keep], function(x){
        tmp <- ensGenes$entrezid[ensGenes$symbol %in% x]
        tmp[!is.na(tmp)]
    }))
}

res <- lapply(1:ncol(tfit), function(i){

    tmp <- vector("list", 2)
    #names(tmp) <- c("champ.wt", "champ.kpmt")
    names(tmp) <- c("ebgsea.wt", "ebgsea.kpmt")

    groups <- names(tfit$contrasts[,i])[tfit$contrasts[,i] != 0]
    # ebgs <- data.frame(champ.ebGSEA(beta = mVals[,targets$CellType %in% cellType],
    #                    pheno = targets$CellType[targets$CellType %in% cellType],
    #                    minN = 5, adjPval=1, arraytype = "EPIC")[[1]])

    samps <- as.logical(unname(rowSums(tfit$design[,groups])))
    pheno <- unname(tfit$design[samps, groups][,1])

    gtRanks <- doGT(pheno.v = pheno,
                    data.m = bVals[, samps],
                    array = arrayType,
                    ncores = 1)
    ebgs <- data.frame(doGSEAwt(rankEID.m = gtRanks, ptw.ls = collection,
                               ncores = 1, minN = minsize, adjPVth = 1)$`Rank(P)`)

    tmp[[1]] <- tibble::rownames_to_column(ebgs, var = "ID")[, c("ID", "P.WT.")]
    colnames(tmp[[1]])[2] <- "pvalue"
    tmp[[2]] <- tibble::rownames_to_column(ebgs, var = "ID")[, c("ID", "P.KPMT.")]
    colnames(tmp[[2]])[2] <- "pvalue"
    dplyr::bind_rows(tmp, .id = "method")

})

names(res) <- colnames(tfit$contrasts)
res <- dplyr::bind_rows(res, .id = "contrast")
res <- tibble::add_column(res, sub = "n", .after = 1)
res$set <- set

out <- glue::glue("{outDir}/{package}.{set}.rds")
saveRDS(res, out)
