library(here)

# extract columns from simulation results
getSimResults <- function(file,column){

  load(file)
  a <- data.frame(Method="Adj", sapply(randCpgSim,function(sim){
                    sim$Adj[,column]}))
  u <- data.frame(Method="Unadj", sapply(randCpgSim,function(sim){
                    sim$Unadj[,column]}))
  rbind(a,u)
}

# get names of simulation results files
files <- list.files(here("output/test"), pattern = "randCpgSim",full.names = TRUE)
# order results files by no. of Cpgs sampled
n <- sapply(files, function(file){
  as.numeric(gsub("randCpgSim","",
                  unlist(strsplit(unlist(strsplit(file,"/"))[9],".",fixed=TRUE))[1]))
}, USE.NAMES = FALSE)
o <- order(n)
# get raw p-values from all simulation results
simPvals <- lapply(files[o],getSimResults,"P.DE")
names(simPvals) <- n[o]
save(simPvals, file=here("output/allRandCpgSimPvals.RData"))

