merge_scan_events=function(infile,outfile,oddrange,evenrange){


xraw = xcmsRaw(infile,profstep = 0)

# Get all scans seperately
scans=list()
scanrange = matrix(nrow = length(xraw@scanindex),ncol=2)
colnames(scanrange)=c("mzmin","mzmax")
for (i in 1:length(xraw@scanindex)){
  scans[[i]] = getScan(xraw,i)
  scanrange[i,"mzmin"]=min(scans[[i]][,"mz"])
  scanrange[i,"mzmax"]=max(scans[[i]][,"mz"])
}


odd_scans = seq(from=1, to=nrow(scanrange)-1,by=2)
even_scans = seq(from=2, to=nrow(scanrange),by=2)


# Check if larger mass scan are always after lower mass scan
scan_order_check = matrix(nrow = 0,ncol=2)
colnames(scanrange)=c("mzmin","mzmax")

for( i in odd_scans){
  temp =   scanrange[i,"mzmin"]<scanrange[i+1,"mzmin"]
  temp[2] =   scanrange[i,"mzmax"]<scanrange[i+1,"mzmax"]
  names(temp)=c("mzmin","mzmax")
  scan_order_check = rbind(scan_order_check,temp)
}

# If not then stop
if(!(all(scan_order_check)))
  stop("Scans are not in order. No implementation to handle this yet.")



# Enforce scan ranges
if(!missing(oddrange)){
  for( i in odd_scans){
    select = scans[[i]][,"mz"]>oddrange[1]    &    scans[[i]][,"mz"]<oddrange[2]
    scans[[i]]   =    scans[[i]][select,]
  }
}


if(!missing(evenrange)){
  for( i in even_scans){
    select = scans[[i]][,"mz"]>evenrange[1]    &    scans[[i]][,"mz"]<evenrange[2]
    scans[[i]]   =    scans[[i]][select,]
  }
}



#Merge the scans
scans_new=list()
for( i in 1:length(odd_scans)){
  scans_new[[i]]=rbind(  scans[[odd_scans[i]]],  scans[[odd_scans[i]+1]]             )
}


#make sure mz are ordered
for( i in 1:length(scans_new)){
  order = order(scans_new[[i]][,"mz"])
  scans_new[[i]]=scans_new[[i]][order,,drop=FALSE]
}


# Create new scanindex
scanindex_new=numeric()
for (i in 1:length(scans_new)){
  if (i==1){
    scanindex_new[i]=0
  }else{
    scanindex_new[i] = scanindex_new[i-1] + nrow(scans_new[[(i-1)]])
  }
}






# Put data into new object
ob<-new("xcmsRaw")
ob@env = new.env(parent=.GlobalEnv)

scans_new=do.call(rbind,scans_new)
ob@env$mz<-scans_new[,'mz']
ob@env$intensity<-scans_new[,'intensity']

ob@scantime <-xraw@scantime[odd_scans]
ob@scanindex <-as.integer(scanindex_new)
ob@polarity   <- xraw@polarity[odd_scans]
ob@acquisitionNum <- 1:length(ob@scanindex)
ob@profmethod     <- xraw@profmethod
ob@mzrange        <- xraw@mzrange
ob@filepath       <- xraw@filepath
profStep(ob) =   profStep(xraw)
ob=xcms:::remakeTIC(ob)


dir.create(dirname(outfile),showWarnings=F,recursive = T)
write.mzdata(ob,outfile)


}





