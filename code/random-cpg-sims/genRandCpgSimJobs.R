#!/usr/bin/env Rscript

noCpgs <- c(50, 100, 500, 1000, 5000, 10000)
arrayType <- c("EPIC","450K")

dir <- here::here("code")
jobDir <- glue::glue("{dir}/.job")
outDir <- here::here("output/random-cpg-sims")

if (!dir.exists(outDir)) dir.create(outDir)
if (!dir.exists(jobDir)) dir.create(jobDir)

for (cpgs in noCpgs) {
  for (array in arrayType) {
    # Start writing to this file
    jobFile <- glue::glue("{jobDir}/{array}.{cpgs}.job")
    sink(file = jobFile)

    # the basic job submission script is a bash script
    cat("#!/bin/bash\n")
    cat(glue::glue("#SBATCH --job-name={array}.{cpgs}.job"), "\n")
    cat(glue::glue("#SBATCH --output={outDir}/{array}.{cpgs}.out"), "\n")
    cat(glue::glue("#SBATCH --error={outDir}/{array}.{cpgs}.err"), "\n")
    cat("#SBATCH --time=00:15:00\n")
    cat("#SBATCH --mem=6144\n")
    cat("#SBATCH --array=1-100\n")
    cat("#\n")
    cat("module load R/3.6.0\n")
    cat("#\n")
    cat(glue::glue("Rscript {dir}/randomCpgSim.R {array} {cpgs} $SLURM_ARRAY_TASK_ID {outDir}"),
        "\n")

    # Close the sink!
    sink()

    # Submit to run on cluster
    system(glue::glue("sbatch {jobFile}"))
  }
}
