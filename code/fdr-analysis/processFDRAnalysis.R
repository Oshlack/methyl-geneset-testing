#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
sampleType <- args[1]

dir <- here::here("output/FDR-analysis")

inFiles <- list.files(dir, pattern = glue::glue("FDR.{sampleType}.*.*.rds"),
                      full.names = TRUE)
outFiles <- list.files(dir, pattern = glue::glue("FDR.{sampleType}.*.out"),
                       full.names = TRUE)
errFiles <- list.files(dir, pattern = glue::glue("FDR.{sampleType}.*.err"),
                       full.names = TRUE)

bpparam <- BiocParallel::MulticoreParam(min(20, BiocParallel::multicoreWorkers()))
tmp <- BiocParallel::bplapply(inFiles, readRDS, BPPARAM = bpparam)
obj <- dplyr::bind_rows(tmp)

saveRDS(obj, glue::glue("{dir}/FDR.{sampleType}.analysis.rds"))

binDir <- glue::glue("{dir}/.bin")
if (!dir.exists(binDir)) dir.create(binDir)

filesstrings::file.move(inFiles, binDir)
filesstrings::file.move(outFiles, binDir)
filesstrings::file.move(errFiles, binDir)
