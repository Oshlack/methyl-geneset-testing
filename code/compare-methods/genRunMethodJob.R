#!/usr/bin/env Rscript

dir <- here::here("code")
jobDir <- glue::glue("{dir}/.job")
outDir <- here::here("output/compare-methods")
input <- here::here("data/cache-intermediates/blood.contrasts.rds")

if (!dir.exists(outDir)) dir.create(outDir)
if (!dir.exists(jobDir)) dir.create(jobDir)

packages <- c("MissMethyl", "MethylGSA", "ebGSEA")
sets <- c("GO", "KEGG", "BROAD")

for (package in packages) {
    for (set in sets) {
        # if(package == "ChAMP" & set != "BROAD"){
        #     next
        #
        # }

        # Start writing to this file
        jobFile <- glue::glue("{jobDir}/{package}.{set}.job")
        sink(file = jobFile)

        # the basic job submission script is a bash script
        cat("#!/bin/bash\n")
        cat(glue::glue("#SBATCH --job-name={package}.{set}.job"), "\n")
        cat(glue::glue("#SBATCH --output={outDir}/{package}.{set}.out"), "\n")
        cat(glue::glue("#SBATCH --error={outDir}/{package}.{set}.err"), "\n")
        cat("#SBATCH --time=12:00:00\n")
        cat("#SBATCH --mem=32GB\n")
        cat("#\n")
        #cat("module load R/3.6.0\n")
        cat("module load R/4.0.2\n")
        cat("#\n")
        cat(glue::glue("Rscript {dir}/compare-methods/run{package}.R {package} {set} {input} {outDir}"),
            "\n")

        # Close the sink!
        sink()

        # Submit to run on cluster
        system(glue::glue("sbatch {jobFile}"))

    }
}
