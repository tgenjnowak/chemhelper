analyze.xcms.group=function(xset_grouped,mz,rt,rt_tol_sample=60,mz_tol_sample=0.01,rt_tol_group=30,mz_tol_group=0.05){
  
  ## Ungrouped data. In which samples is it found?
  select = abs(xset_grouped@peaks[,"rt"]-rt)<rt_tol_sample & abs(xset_grouped@peaks[,"mz"]-mz)<mz_tol_sample
  if(sum(select)==0){stop("No peaks found.")}
  
  
  # xset_grouped@filepaths[       unique(xset_grouped@peaks[select,"sample"])       ]
  # sum(select)
  # 
  # min(xset_grouped@peaks[select,"rt"])/60
  # max(xset_grouped@peaks[select,"rt"])/60
  # min(xset_grouped@peaks[select,"mz"])
  # max(xset_grouped@peaks[select,"mz"])
  
  
  # grouped data
  select_g = abs(xset_grouped@groups[,"rtmed"]-rt)<rt_tol_group & abs(xset_grouped@groups[,"mzmed"]-mz)<mz_tol_group
  select_g = which(select_g)
  group_report=cbind(rtmed=xset_grouped@groups[select_g,"rtmed"]/60,rtmin=xset_grouped@groups[select_g,"rtmin"]/60,rtmax=xset_grouped@groups[select_g,"rtmax"]/60,mzmed=xset_grouped@groups[select_g,"mzmed"],mzmin=xset_grouped@groups[select_g,"mzmin"],mzmax=xset_grouped@groups[select_g,"mzmax"])
  
  
  #plot
  plot(xset_grouped@peaks[select,"rt"]/60,xset_grouped@peaks[select,"mz"],pch=0,cex=0.5,             ylim=c(min(xset_grouped@groups[select_g,"mzmin"]),max(xset_grouped@groups[select_g,"mzmax"])),                     ,xlim=c(min(xset_grouped@groups[select_g,"rtmin"])/60,max(xset_grouped@groups[select_g,"rtmax"])/60)    ,xlab="rt (min)",ylab="m/z", main=paste("Grouping of m/z=",mz,", rt=",rt/60,"min",sep=""))
    
  for (i in 1:nrow(group_report)){
    points(    xset_grouped@peaks[  xset_grouped@groupidx[[select_g[i]]]   ,"rt"] /60 ,  xset_grouped@peaks[  xset_grouped@groupidx[[select_g[i]]]   ,"mz"]    ,col=i,pch=20)
    rect(group_report[i,"rtmin"], group_report[i,"mzmin"], group_report[i,"rtmax"], group_report[i,"mzmax"],border=NA,col=alpha(i, 0.5))
  }
  
  points(rt/60,mz,col="red",pch=20,cex=3)
  
  return(group_report)
  
}