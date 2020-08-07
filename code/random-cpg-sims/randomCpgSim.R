#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (args[1] == "EPIC"){
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

} else if (args[1] == "450K") {
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

}

arrayType <- args[1]
noCpgs <- args[2]
simNo <- args[3]
outDir <- args[4]

# set a seed for each simulation using sim no. indicated by slurm
# this ensures different CpGs are sampled for each of the 100 simulations
set.seed(simNo)
randSet <- sample(x = ann$Name, size = noCpgs) # randomly sample cpgNum sites

res <- vector("list", 5)
names(res) <- c("hgt", "hgt.cpg", "hgt.cpg.fc", "hgt.cpg.ec", "hgt.cpg.fc.ec")

res[[1]] <- tibble::rownames_to_column(missMethyl::gometh(sig.cpg = randSet,
              array.type = arrayType, anno = ann, prior.prob = FALSE,
              fract.counts = FALSE, equiv.cpg = FALSE),
              var = "GO")[, c("GO","P.DE")]

res[[2]] <- tibble::rownames_to_column(missMethyl::gometh(sig.cpg = randSet,
              array.type = arrayType, anno = ann, prior.prob = TRUE,
              fract.count = FALSE, equiv.cpg = FALSE),
              var = "GO")[, c("GO","P.DE")]

res[[3]] <- tibble::rownames_to_column(missMethyl::gometh(sig.cpg = randSet,
              array.type = arrayType, anno = ann, prior.prob = TRUE,
              fract.counts = TRUE, equiv.cpg = FALSE),
              var = "GO")[, c("GO","P.DE")]

res[[4]] <- tibble::rownames_to_column(missMethyl::gometh(sig.cpg = randSet,
              array.type = arrayType, anno = ann, prior.prob = TRUE,
              fract.counts = FALSE, equiv.cpg = TRUE),
              var = "GO")[, c("GO","P.DE")]

res[[5]] <- tibble::rownames_to_column(missMethyl::gometh(sig.cpg = randSet,
                     array.type = arrayType, anno = ann, prior.prob = TRUE,
                     fract.counts = TRUE, equiv.cpg = TRUE),
                     var = "GO")[, c("GO","P.DE")]

res <- dplyr::bind_rows(res, .id = "method")
res <- dplyr::bind_cols(res, simNo = rep(simNo, nrow(res)),
                        noCpgs = rep(noCpgs, nrow(res)))

out <- glue::glue("{outDir}/{arrayType}.{noCpgs}.{simNo}.rds")
saveRDS(res, out)
