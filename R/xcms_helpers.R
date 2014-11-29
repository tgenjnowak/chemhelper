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







xsAnnotate_stats <- function(object){

# get the peaklist
peaklist=getPeaklist(object)


# setup output
annotation_statistics=matrix(ncol=2,nrow=10)

annotation_statistics[1,1]='# Feature'
annotation_statistics[2,1]='# Feature groups'
annotation_statistics[3,1]='# Features groups with several features'
annotation_statistics[4,1]='# Features in feature groups with several features'

annotation_statistics[5,1]='# Features annotated as isotopes'
annotation_statistics[6,1]='# Feature with adducts/fragments annotation'
annotation_statistics[7,1]='# Feature groups with fragment/adduct annotation'
annotation_statistics[8,1]='# Annotation groups'
annotation_statistics[9,1]='Mean number of annotation groups per feature group with annotations'
annotation_statistics[10,1]='Median number of annotation groups per feature group with annotations'



# Feature
annotation_statistics[1,2]=nrow(peaklist['adduct'])


# Feature groups
annotation_statistics[2,2]=length(unique(peaklist[,'pcgroup']))


# Features groups with several features
annotation_statistics[3,2]=
  paste(
    sum(table(peaklist[,'pcgroup'])>1),
    ' (',
    round(  sum(table(peaklist[,'pcgroup'])>1)           /       length(unique(peaklist[,'pcgroup'])) *100   ,0)    ,
    ' %)', sep='')


# Features in feature groups with several features
annotation_statistics[4,2]=
  paste(
    sum(table(peaklist[,'pcgroup'])    [table(peaklist[,'pcgroup'])>1]),
    ' (',
    round(  sum(table(peaklist[,'pcgroup'])    [table(peaklist[,'pcgroup'])>1])         /          nrow(peaklist['adduct'])*100      ,)    ,
    ' %)', sep='')



# Features annotated as isotopes
annotation_statistics[5,2]=
  paste(
    nrow(object@isoID),
    ' (',
    round(  nrow(object@isoID)   /  sum(table(peaklist[,'pcgroup'])    [table(peaklist[,'pcgroup'])>1])*100      ,)    ,
    ' %)', sep='')



# Feature with adducts/fragments annotation
annotation_statistics[6,2]=
  paste(
    sum(peaklist['adduct']!=''),
    ' (',
    round(  sum(peaklist['adduct']!='') / (nrow(peaklist['adduct']) - nrow(object@isoID))*100      ,)    ,
    ' %)', sep='')



# Feature groups with fragment/adduct annotation
annotation_statistics[7,2]=
  paste(
    length(table(object@annoGrp[,'psgrp'])),
    ' (',
    round(  length(table(object@annoGrp[,'psgrp']))  / sum(table(peaklist[,'pcgroup'])>1)*100      ,)    ,
    ' %)', sep='')



# Annotation groups
annotation_statistics[8,2]=nrow((object@annoGrp))

#Mean number of annotation groups per feature group with annotations
annotation_statistics[9,2]=round(mean(table(object@annoGrp[,'psgrp']))  ,1)

#Median number of annotation groups per feature group with annotations
annotation_statistics[10,2]=median(table(object@annoGrp[,'psgrp']))


return(annotation_statistics)

}