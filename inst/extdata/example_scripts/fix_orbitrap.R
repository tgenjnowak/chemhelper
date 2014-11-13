library(chemhelper)
library(parallel)

files=list.files("***********FULL INPUT PATH THAT CONTAINS FILES*************",pattern='.mzXML',full.names=T,recursive=T)
outdir = "***********FULL OUTPUT PATH*************"


# make files and outdir into a lsit
input=cbind(files,outdir)
input <- split(input, 1:NROW(input))


#now make a cluster and convert in parallel
cl <- makeCluster(detectCores())
clusterExport(cl,c("input","rem_satellite_peaks")          )
parLapply(cl,input,function(x) rem_satellite_peaks(x[1],x[2]))
stopCluster(cl)