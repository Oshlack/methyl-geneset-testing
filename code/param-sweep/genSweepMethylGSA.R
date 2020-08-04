#!/usr/bin/env Rscript

dir <- here::here("code")
jobDir <- glue::glue("{dir}/.job")
outDir <- here::here("output/methylgsa-params")
input <- here::here("data/cache-intermediates/input.RData")

if (!dir.exists(outDir)) dir.create(outDir)
if (!dir.exists(jobDir)) dir.create(jobDir)

methods <- c("glm", "ora", "gsea")
params <- data.frame(minsize = c(1:5, rep(5, 5)),
                     maxsize = c(rep(5000, 5), seq(7000, 18000, by = 2750)))

for (method in methods) {
    for(i in 1:nrow(params)){
        min <- params$minsize[i]
        max <- params$maxsize[i]

        # Start writing to this file
        jobFile <- glue::glue("{jobDir}/params.{method}.{min}.{max}.job")
        sink(file = jobFile)

        # the basic job submission script is a bash script
        cat("#!/bin/bash\n")
        cat(glue::glue("#SBATCH --job-name={method}.{min}.{max}.job"), "\n")
        cat(glue::glue("#SBATCH --output={outDir}/{method}.{min}.{max}.out"), "\n")
        cat(glue::glue("#SBATCH --error={outDir}/{method}.{min}.{max}.err"), "\n")
        cat("#SBATCH --time=12:00:00\n")
        if(method == "gsea"){
            cat("#SBATCH --mem=32767\n")

        } else {
            cat("#SBATCH --mem=16384\n")

        }
        cat("#\n")
        cat("module load R/3.6.0\n")
        cat("#\n")
        cat(glue::glue("Rscript {dir}/paramSweepMethylGSA.R {method} {input} {min} {max} {outDir}"),
            "\n")

        # Close the sink!
        sink()

        # Submit to run on cluster
        system(glue::glue("sbatch {jobFile}"))

    }

}
