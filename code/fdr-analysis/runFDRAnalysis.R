#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

sampleNo <- args[1]
simNo <- args[2]
input <- args[3]
broadFile <- args[4]
outDir <- args[5]

library(here)
library(glue)
library(missMethyl)
library(ChAMP)
library(methylGSA)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
source(here::here("code/utility.R"))

ann <- loadAnnotation(arrayType = "450k")
dat <- readRDS(input) # read in the KIRC data

# set a seed for each simulation using sim no. indicated by slurm
set.seed(simNo)

set1 <- sample(x = 1:ncol(dat), size = sampleNo)
set2 <- sample(x = (1:ncol(dat))[-set1], size = sampleNo)
subSet <- c(set1, set2)

design <- model.matrix(~factor(c(rep(1, sampleNo), rep(0, sampleNo))))
fit <- limma::lmFit(minfi::logit2(dat[, subSet]), design)
fit <- limma::eBayes(fit, robust = TRUE)

broad <- readRDS(broadFile)

# run gometh
tmp <- topGSA(gsameth(sig.cpg = rownames(limma::topTable(fit, n = 1000, coef = 2)),
                      all.cpg = rownames(dat), collection = broad$entrez,
                      array.type = "450k", anno=ann, prior.prob = TRUE,
                      equiv.cpg = TRUE, fract.counts = TRUE), number = Inf)
gm1k <- tibble::rownames_to_column(tmp, var = "ID")[, c("ID", "P.DE")]
colnames(gm1k)[2] <- "pvalue"
tmp <- topGSA(gsameth(sig.cpg = rownames(limma::topTable(fit, n = 5000, coef = 2)),
                      all.cpg = rownames(dat), collection = broad$entrez,
                      array.type = "450k", anno = ann, prior.prob = TRUE,
                      equiv.cpg = TRUE, fract.counts = TRUE), number = Inf)
gm5k <- tibble::rownames_to_column(tmp, var = "ID")[, c("ID", "P.DE")]
colnames(gm5k)[2] <- "pvalue"

# run methylGSA
tmp <- methylglm(cpg.pval = fit$p.value[, 2], FullAnnot = ann, minsize = 5,
                 maxsize = 5000, GS.list = broad$symbol, GS.idtype = "SYMBOL")
glm <- tmp[,c("ID", "pvalue")]
tmp <- methylRRA(cpg.pval = fit$p.value[, 2], method = "ORA", FullAnnot = ann,
                 minsize = 5, maxsize = 5000, GS.list = broad$symbol,
                 GS.idtype = "SYMBOL")
ora <- tmp[,c("ID", "pvalue")]
tmp <- methylRRA(cpg.pval = fit$p.value[, 2], method = "GSEA", FullAnnot = ann,
                  minsize = 5, maxsize = 5000, GS.list = broad$symbol,
                  GS.idtype = "SYMBOL")
gsea <- tmp[,c("ID", "pvalue")]

# run ebGSEA
tmp <- data.frame(champ.ebGSEA(beta = dat[,subSet], pheno=design[,2], minN=5, adjPval=1,
                     arraytype="450k")[[1]])
wt <- tibble::rownames_to_column(tmp, var = "ID")[,c("ID","P.WT.")]
kpmt <- tibble::rownames_to_column(tmp, var = "ID")[,c("ID","P.KPMT.")]
colnames(wt)[2] <- "pvalue"
colnames(kpmt)[2] <- "pvalue"

res <- list(mmethyl.gm1000 = gm1k,
            mmethyl.gm5000 = gm5k,
            mgsa.glm = glm,
            mgsa.ora = ora,
            mgsa.gsea = gsea,
            champ.wt = wt,
            champ.kpmt = kpmt)
res <- dplyr::bind_rows(res, .id = "method")
res$simNo <- simNo
res$sampleNo <- sampleNo

# save results to file
out <- glue::glue("{outDir}/FDR.{sampleNo}.{simNo}.rds")
saveRDS(res, out)

