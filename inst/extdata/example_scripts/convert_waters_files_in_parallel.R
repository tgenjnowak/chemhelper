require(stringr)
require(chemhelper)
require(parallel)


#list of folders with files to convert
indir_list=c("Standards_2013_Apr (synthesized peptides).pro/Data",
        "Standards_2013_Feb (cyclic dipeptides etc).pro/Data",
        "Standards_2013_june.pro/Data",
        "Standards_2013_Mar (confirmations for 95B).pro/Data",
        "Standards_2013_sep (95A).pro/Data",
        "Standards_2013-oct.pro/Data",
        "Standards_AgustJuly_2013.pro/Data",
        "Standards_August_2013.pro/Data",
        "Standards_Dec -12.pro/Data",
        "Standards_Feb-13/data",
        "Standards_jan2013.pro/Data",
        "Standards_JanFeb-12.pro/Data",
        "Standards_May2013/Standards4/Data",
        "Standards_Sep2013",
        "Standards_Sept-12/Data"
        )

indir=paste("***********FULL INPUT PATH TO ABOVE FOLDERS*************",indir_list,"/",sep="")

# Now we make a list of matching output floders
outdir=paste("***********FULL OUTPUT PATH *************",indir_list,"/",sep="")
outdir=str_replace(outdir,"/Data","")
outdir=str_replace(outdir,"/data","")


# combine indir and outdir in a list so it can be passed to apply functions
arg_list=list()
for (i in 1:length(indir)){
arg_list[[i]]=c(indir[i],outdir[i])
}



# Do parallel conversion
cl <- makeCluster(getOption("cl.cores", 2))
clusterExport(cl, "convert.waters")
parLapplyLB(cl,arg_list,function(x) convert.waters(x[1],x[2]))
stopCluster(cl)
