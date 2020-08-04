#!/usr/bin/env Rscript

sampleNos <- c(5, 10, 20, 40, 80)

dir <- here::here("output/FDR-analysis")
rdsFiles <- numeric(0)
outFiles <- numeric(0)

obj <- vector("list", length(sampleNos))
i = 1

for (sampleNo in sampleNos) {
    inFiles <- sapply(1:100, function(simNo){
        glue::glue("{dir}/FDR.{sampleNo}.{simNo}.rds")
    })

    bpparam <- BiocParallel::MulticoreParam(min(20, BiocParallel::multicoreWorkers()))
    tmp <- BiocParallel::bplapply(inFiles, readRDS, BPPARAM = bpparam)
    obj[[i]] <- dplyr::bind_rows(tmp)

    rdsFiles <- c(rdsFiles, inFiles)
    outFiles <- c(outFiles, glue::glue("{dir}/FDR.{sampleNo}.err"),
                  glue::glue("{dir}/FDR.{sampleNo}.out"))

    i = i + 1
}
obj <- dplyr::bind_rows(obj)

saveRDS(obj, glue::glue("{dir}/FDR-analysis.rds"))

binDir <- glue::glue("{dir}/.bin")
if (!dir.exists(binDir)) dir.create(binDir)

filesstrings::file.move(rdsFiles, binDir)
filesstrings::file.move(outFiles, binDir)
