#!/usr/bin/env Rscript

noCpgs <- c(50, 100, 500, 1000, 5000, 10000)
arrayType <- c("EPIC","450K")

dir <- here::here("output/random-cpg-sims/")
rdsFiles <- numeric(0)
outFiles <- numeric(0)
noSim <- numeric(0)

for (array in arrayType) {
  obj <- vector("list", length(noCpgs))
  i = 1

  for (cpgs in noCpgs) {
    inFiles <- sapply(1:100, function(sim){
      glue::glue("{dir}/{array}.{cpgs}.{sim}.rds")
    })

    bpparam <- BiocParallel::MulticoreParam(min(20, BiocParallel::multicoreWorkers()))
    tmp <- BiocParallel::bplapply(inFiles, readRDS, BPPARAM = bpparam)
    obj[[i]] <- dplyr::bind_rows(tmp)

    rdsFiles <- c(rdsFiles, inFiles)
    outFiles <- c(outFiles, glue::glue("{dir}/{array}.{cpgs}.err"),
                  glue::glue("{dir}/{array}.{cpgs}.out"))
    i = i + 1

  }
  obj <- dplyr::bind_rows(obj)
  obj <- dplyr::mutate(obj, noCpgs = factor(noCpgs, levels = c(50, 100, 500,
                                                               1000, 5000,
                                                               10000)))
  saveRDS(obj, glue::glue("{dir}/{array}.rds"))

}

binDir <- glue::glue("{dir}/.bin")
if (!dir.exists(binDir)) dir.create(binDir)

filesstrings::file.move(rdsFiles, binDir)
filesstrings::file.move(outFiles, binDir)
