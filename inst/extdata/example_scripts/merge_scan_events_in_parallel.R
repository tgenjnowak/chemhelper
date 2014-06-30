require(chemhelper)
require(parallel)
require(stringr)

infiles = list.files(path="***********FULL INPUT PATH*************",pattern=".mzML",recursive = T,full.names = T)
outfiles = str_replace(infiles,"org_data","scans_merged")
outfiles = str_replace(outfiles,".mzML",".mzData")




# combine infiles and outfiles in a list so it can be passed to apply functions
arg_list=list()
for (i in 1:length(infiles)){
  arg_list[[i]]=c(infiles[i],outfiles[i])
}



# Do parallel conversion
cl <- makeCluster(detectCores())
clusterExport(cl, "merge_scan_events")
parLapplyLB(cl,arg_list,function(x) merge_scan_events(x[1],x[2],oddrange = c(0,500),evenrange = c(500,2000)))  # remember to change to your mzranges
stopCluster(cl)
