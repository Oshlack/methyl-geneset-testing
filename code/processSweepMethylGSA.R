#!/usr/bin/env Rscript

methods <- c("glm", "gsea", "ora")
params <- data.frame(minsize = c(1:5, rep(5, 5)),
                     maxsize = c(rep(5000, 5), seq(7000, 18000, by = 2750)))

dir <- here::here("output/methylgsa-params")
rdsFiles <- numeric(0)
outFiles <- numeric(0)

obj <- vector("list", length(methods) * nrow(params))
i = 1

for (method in methods) {
        inFiles <- sapply(1:nrow(params), function(i){
            minsz <- params$minsize[i]
            maxsz <- params$maxsize[i]
            glue::glue("{dir}/{method}.{minsz}.{maxsz}.rds")
        })

        bpparam <- BiocParallel::MulticoreParam(min(length(inFiles),
                                                    BiocParallel::multicoreWorkers()))
        tmp <- BiocParallel::bplapply(inFiles, readRDS, BPPARAM = bpparam)
        obj[[i]] <- dplyr::bind_rows(tmp)

        rdsFiles <- c(rdsFiles, inFiles)
        outFiles <- c(outFiles, gsub("rds", "err", inFiles),
                      gsub("rds", "out", inFiles))

        i = i + 1
}
obj <- dplyr::bind_rows(obj)

saveRDS(obj, glue::glue("{dir}/methylGSA-param-sweep.rds"))

binDir <- glue::glue("{dir}/.bin")
if (!dir.exists(binDir)) dir.create(binDir)

filesstrings::file.move(rdsFiles, binDir)
filesstrings::file.move(outFiles, binDir)
