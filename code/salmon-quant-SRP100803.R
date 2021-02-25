#!/usr/bin/env Rscript

dir <- here::here("data/datasets/SRP100803")
jobDir <- glue::glue("{dir}/.job")
if (!dir.exists(jobDir)) dir.create(jobDir)

files <- list.files(here::here("data/datasets/SRP100803"), recursive = TRUE,
                    pattern = "sra", full.names = TRUE)
index <- "/oshlack_lab/shared/genomes/hg19/gencode_v37/indices/salmon/salmon_index_v1.3.0"

for (file in files){

    name <- gsub(".sra", "", basename(file))
    # Start writing to this file
    jobFile <- glue::glue("{jobDir}/{name}.job")
    sink(file = jobFile)

    # the basic job submission script is a bash script
    cat("#!/bin/bash\n")
    cat(glue::glue("#SBATCH --job-name={name}.job"), "\n")
    cat(glue::glue("#SBATCH --output={dirname(file)}/{name}.out"), "\n")
    cat(glue::glue("#SBATCH --error={dirname(file)}/{name}.err"), "\n")
    cat("#SBATCH --time=8:00:00\n")
    cat("#SBATCH --mem=24GB\n")
    cat("#SBATCH --ntasks=1\n")
    cat("#SBATCH --cpus-per-task=8\n")
    cat("#\n")
    cat("module load salmon/1.3.0\n")
    cat("module load sratoolkit/2.9.0\n")
    cat("#\n")
    cat(glue::glue("fastq-dump -I --split-files -O {dirname(file)} {dirname(file)}/{name}.sra\n\n"))
    cat("#\n")
    cat(glue::glue("salmon quant -i {index} -l A -1 {dirname(file)}/{name}_1.fastq -2 {dirname(file)}/{name}_2.fastq -p 8 --validateMappings -o {dir}/quants/{name}_quant\n\n"))
    cat("#\n")

    # Close the sink!
    sink()

    # Submit to run on cluster
    system(glue::glue("sbatch {jobFile}"))

}
