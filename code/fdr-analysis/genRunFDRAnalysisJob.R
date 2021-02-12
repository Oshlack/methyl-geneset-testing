#!/usr/bin/env Rscript

dir <- here::here("code")
jobDir <- glue::glue("{dir}/.job")
outDir <- here::here("output/FDR-analysis")
inputs <- list(KIRC = here::here("data/datasets/TCGA.KIRC.rds"),
               BRCA = here::here("data/datasets/TCGA.BRCA.rds"))

if (!dir.exists(outDir)) dir.create(outDir)
if (!dir.exists(jobDir)) dir.create(jobDir)

library(ChAMP)
library(org.Hs.eg.db)
library(AnnotationDbi)

broadFile <- here::here("output/FDR-analysis/BROAD-sets.rds")

if (!file.exists(broadFile)){
    data("PathwayList")
    k <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    keep <- sapply(PathwayList,function(x) any(x %in% k))
    PathwayList <- PathwayList[keep]

    # convert SYMBOLS to ENTREZID for gometh
    collection <- suppressMessages(lapply(PathwayList, function(x){
        tmp <- select(org.Hs.eg.db, x, columns = "ENTREZID",
                      keytype = "SYMBOL")$ENTREZID
        tmp[!is.na(tmp)]
    }))

    broad <- list(symbol = PathwayList, entrez = collection)
    saveRDS(broad, file = broadFile)

}

for (sampleType in names(inputs)){

    sampleNos <- c(5, 10, 20, 40, 80)
    if (sampleType == "BRCA") sampleNos <- sampleNos[-length(sampleNos)]

    for (sampleNo in sampleNos) {

        # Start writing to this file
        jobFile <- glue::glue("{jobDir}/FDR.{sampleNo}.job")
        sink(file = jobFile)

        # the basic job submission script is a bash script
        cat("#!/bin/bash\n")
        cat(glue::glue("#SBATCH --job-name=FDR.{sampleNo}.job"), "\n")
        cat(glue::glue("#SBATCH --output={outDir}/FDR.{sampleNo}.out"), "\n")
        cat(glue::glue("#SBATCH --error={outDir}/FDR.{sampleNo}.err"), "\n")
        cat("#SBATCH --time=4:00:00\n")
        cat("#SBATCH --mem=16384\n")
        cat("#SBATCH --array=1-100\n")
        cat("#\n")
        cat("module load R/3.6.0\n")
        cat("#\n")
        cat(glue::glue("Rscript {dir}/runFDRAnalysis.R {sampleNo} $SLURM_ARRAY_TASK_ID {input} {broadFile} {outDir} {sampleType}"),
            "\n")

        # Close the sink!
        sink()

        # Submit to run on cluster
        system(glue::glue("sbatch {jobFile}"))

    }
}
