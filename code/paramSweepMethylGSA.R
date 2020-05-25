#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
method <- args[1]
input <- args[2]
min <- args[3]
max <- args[4]
outDir <- args[5]

library(EnsDb.Hsapiens.v75)
library(methylGSA)
library(tidyverse)
library(tictoc)
library(BiocParallel)
source(here::here("code/utility.R"))

cat("Convert GO entrez IDs to SYMBOLs\n")
edb <- EnsDb.Hsapiens.v75
trans <- genes(edb, columns = c("symbol", "gene_id", "entrezid"),
               return.type = "data.frame")
ens <- reshape2::melt(trans$entrezid, value.name = "entrezid")
colnames(ens)[2] <- "gene_id"
goentrez <- missMethyl:::.getGO()$idList
go <- reshape2::melt(goentrez, value.name = "entrezid")
colnames(go)[2] <- "GO"
go %>% mutate(entrezid = as.character(entrezid)) -> go
ens %>% mutate(entrezid = as.character(entrezid)) -> ens

go %>% inner_join(ens) %>%
    inner_join(data.frame(trans[, c("symbol", "gene_id")])) -> go
go %>% dplyr::select(GO, symbol) %>%
    distinct() -> go
symbol <- split(go$symbol, f = factor(go$GO))

cat("Loading annotation\n")
ann <- loadAnnotation(arrayType = "EPIC")
cat(glue::glue("Loading input: {input}"), "\n")
load(input)
result <- NULL

cat("Start loop\n")
for(i in 1:ncol(tfit$contrasts)){
    cat("Method starting\n")
    tic(method)

    if(method == "glm"){
        out <- try(methylglm(cpg.pval = tfit$p.value[,i],
                             FullAnnot = ann, minsize = min,
                             maxsize = max, GS.list = symbol,
                             GS.idtype = "SYMBOL", parallel = TRUE,
                             BPPARAM = MulticoreParam(workers = 10)),
                   silent = TRUE)

    } else if (method == "ora"){
        out <- try(methylRRA(cpg.pval = tfit$p.value[,i],
                             method = "ORA", FullAnnot = ann, minsize = min,
                             maxsize = max, GS.list = symbol,
                             GS.idtype = "SYMBOL"), silent = TRUE)

    } else if (method == "gsea") {
        out <- try(methylRRA(cpg.pval = tfit$p.value[,i],
                             method = "GSEA", FullAnnot = ann, minsize = min,
                             maxsize = max, GS.list = symbol,
                             GS.idtype = "SYMBOL"), silent = TRUE)

    }

    cat("Method done\n")

    toc(log = TRUE, quiet = TRUE)
    runtime <- limma::strsplit2(tic.log(format = TRUE)[[1]], " ")[2]
    tic.clearlog()

    if(class(out) == "try-error"){
        tmp <- data.frame(ID = NA, Size = NA, pvalue = NA, padj = NA,
                          method = method, time = as.numeric(runtime),
                          contrast = colnames(tfit$contrasts)[i],
                          minsize = min, maxsize = max,
                          stringsAsFactors = FALSE)

    } else {
        tmp <- remove_rownames(out[, c("ID","Size","pvalue","padj")])
        tmp$method <- method
        tmp$time <- as.numeric(runtime)
        tmp$contrast <- colnames(tfit$contrasts)[i]
        tmp$minsize <- min
        tmp$maxsize <- max
        tmp <- data.frame(tmp, stringsAsFactors = FALSE)

    }
    result <- bind_rows(result, tmp)
    cat("Results added\n")

}

cat("Saving results\n")
outFile <- glue::glue("{outDir}/{method}.{min}.{max}.rds")
saveRDS(result, file = outFile)
cat("Results saved\n")
