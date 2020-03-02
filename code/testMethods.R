#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

library(methylGSA)
library(missMethyl)
library(ChAMP)
library(here)
library(glue)
library(BiocParallel)
library(org.Hs.eg.db)
library(AnnotationDbi)
source(here("code/helperFunctions.R"))

if (length(args) != 4) {
  stop("Usage: Rscript testMethods.R {method} {array} {gene set} {input file}\n",
       call.=FALSE)

}

method <- args[1] # gometh, mRRA, ebGSEA
arrayType <- args[2] # EPIC, 450k
set <- args[3] # GO, KEGG, BROAD
load(args[4]) # saved input objects

if(set == "GO"){
  if(method == "gometh"){
    collection <- getGOByKey(keytype = "ENTREZID", minsize, maxsize)
  } else if (method == "mRRA"){
    collection <- getGOByKey(keytype = "SYMBOL", minsize, maxsize)
  } else {
    stop("GO categories are only compatible with gometh and mRRA methods.")
  }

} else if(set == "KEGG"){

  if(method %in% c("gometh","mRRA")){
    collection <- missMethyl:::.getKEGG()$idList

    if(method == "mRRA"){
      collection <- suppressMessages(lapply(collection, function(x){
        select(org.Hs.eg.db, x, columns = "SYMBOL", keytype = "ENTREZID")$SYMBOL
      }))
    }
  } else {
    stop("KEGG pathways are only compatible with gometh and mRRA methods.")
  }
} else if (set == "BROAD"){

  data("PathwayList")
  k <- keys(org.Hs.eg.db, keytype = "SYMBOL")
  keep <- sapply(PathwayList,function(x) any(x %in% k))
  PathwayList <- PathwayList[keep]

  if(method == "mRRA"){
    collection <- PathwayList
  } else if (method == "gometh"){
    collection <- suppressMessages(lapply(PathwayList, function(x){
      tmp <- select(org.Hs.eg.db, x, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
      tmp[!is.na(tmp)]
    }))
  }
}

if(method == "gometh"){

  sigList <- lapply(colnames(tfit$contrasts), function(coef){
    rownames(topTreat(tfit, coef = coef, num = 5000))
  })

  universe <- rownames(mVals)
  ann <- loadAnnotation(arrayType=arrayType)

  tmp <- lapply(sigList, function(sig){

    none <- topGSA(gsameth(sig.cpg = sig, all.cpg = universe, collection = collection,
                        array.type = arrayType, anno=ann, prior.prob = FALSE,
                        equiv.cpg = FALSE), number = Inf)
    adj <- topGSA(gsameth(sig.cpg = sig, all.cpg = universe, collection = collection,
                       array.type = arrayType, anno=ann, prior.prob = TRUE,
                       equiv.cpg = TRUE), number = Inf)

    if(set != "BROAD"){
      none <- addDetails(none, collection = set, pname = "P.DE")
      adj <- addDetails(adj, collection = set, pname = "P.DE")
    }

    list(none = none, adj = adj)
  })

} else if (method == "mRRA"){

  pvalList <- lapply(colnames(tfit$contrasts), function(coef){
    top <- limma:::topTreat(tfit, coef = coef, num = Inf)
    pval <- top$P.Value
    names(pval) <- rownames(top)
    pval
  })

  ann <- loadAnnotation(arrayType=arrayType)

  tmp <- lapply(pvalList, function(pval){

    glm <- methylglm(cpg.pval = pval, FullAnnot = ann, minsize = minsize,
                    maxsize = maxsize, GS.list = collection, GS.idtype = "SYMBOL",
                    parallel = TRUE, BPPARAM = MulticoreParam(10))
    ora <- methylRRA(cpg.pval = pval, method = "ORA", FullAnnot = ann,
                    minsize = minsize, maxsize = maxsize, GS.list = collection,
                    GS.idtype = "SYMBOL")
    gsea <- methylRRA(cpg.pval = pval, method = "GSEA", FullAnnot = ann,
                     minsize = minsize, maxsize = maxsize, GS.list = collection,
                     GS.idtype = "SYMBOL")

    if(set != "BROAD"){
      glm <- addDetails(glm, collection = set, pname = "pvalue")
      ora <- addDetails(ora, collection = set, pname = "pvalue")
      gsea <- addDetails(gsea, collection = set, pname = "pvalue")
    }

    list(glm = glm, ora = ora, gsea = gsea)
  })

} else if (method == "ebGSEA"){

  tmp <- lapply(1:ncol(tfit), function(i){
    cellType <- names(tfit$contrasts[,i])[tfit$contrasts[,i] != 0]

    list(ebgs=champ.ebGSEA(beta=mVals[,targets$CellType %in% cellType],
                           pheno=targets$CellType[targets$CellType %in% cellType],
                           minN=5, adjPval=1, arraytype=arrayType))
  })
}

out <- paste(method,arrayType,set,sep=".")
outFile <- here(glue("data/{out}.RData"))

assign(out, tmp)
save(list=c(eval(out)), file = outFile)

































