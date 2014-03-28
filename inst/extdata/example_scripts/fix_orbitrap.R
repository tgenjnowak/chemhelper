library(chemhelper)
library(parallel)

files=list.files('D:/FEM/Other_projects/Marynka',pattern='.mzXML',full.names=T,recursive=T)
outdir = dirname(files)
outdir = gsub('Marynka/','Marynka/fixedfiles/',outdir)

rst

# make files and outdir into a lsit
input=cbind(files,outdir)
input <- split(input, 1:NROW(input))


#now make a cluster and convert in parallel
cl <- makeCluster(getOption("cl.cores", 2))
clusterExport(cl,c("input","rem_satellite_peaks")          )
parLapply(cl,input,function(x) rem_satellite_peaks(x[1],x[2]))
stopCluster(cl)