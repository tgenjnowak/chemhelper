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
  group_report <- cbind(`rtmed (min)` = xset_grouped@groups[select_g,"rtmed"]/60,
                        `rtmin (min)` = xset_grouped@groups[select_g,"rtmin"]/60,
                        `rtmax (min)` = xset_grouped@groups[select_g,"rtmax"]/60,
                        mzmed         = xset_grouped@groups[select_g,"mzmed"],
                        mzmin         = xset_grouped@groups[select_g,"mzmin"],
                        mzmax         = xset_grouped@groups[select_g,"mzmax"])
  
  
  #plot
  plot(xset_grouped@peaks[select,"rt"]/60,xset_grouped@peaks[select,"mz"],pch=0,cex=0.5,             ylim=c(min(xset_grouped@groups[select_g,"mzmin"]),max(xset_grouped@groups[select_g,"mzmax"])),                     ,xlim=c(min(xset_grouped@groups[select_g,"rtmin"])/60,max(xset_grouped@groups[select_g,"rtmax"])/60)    ,xlab="rt (min)",ylab="m/z", main=paste("Grouping of m/z=",mz,", rt=",rt/60,"min",sep=""))
    
  rt_sd <- c()
  for (i in 1:nrow(group_report)){
    rt_points <- xset_grouped@peaks[  xset_grouped@groupidx[[select_g[i]]]   ,"rt"] /60
    rt_sd[i]  <- sd(rt_points*60)
    
    points(   rt_points  ,  xset_grouped@peaks[  xset_grouped@groupidx[[select_g[i]]]   ,"mz"]    ,col=i,pch=20)
    rect(group_report[i,"rtmin (min)"], group_report[i,"mzmin"], group_report[i,"rtmax (min)"], group_report[i,"mzmax"],border=NA,col=alpha(i, 0.5))
  }
  
  points(rt/60,mz,col="red",pch=20,cex=3)
  
  
  group_report <- cbind.data.frame(group_report,`rt_sd (s)` = rt_sd)
  
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






mzmine2xcmsSet <- function(peaklist_mzmine, filepath = NULL, get.scantimes = FALSE){
 
  # to make package check happy about non-standard evaluation.
  type <- 
  value <- 
  `row ID` <- 
  `row m/z` <- 
  `row retention time` <- 
  `All identity` <- 
  egauss <- 
  mu <-
  sigma <-
  h <-
  f <-
  dppm <-
  scpos <-
  scmin <-
  scmax <-lmin <-
  lmax <-
  `All identity elements` <-
  ID <-
  `Molecular formula` <-
  Name <-
  `Identification method` <-
  sample_name <-
  `Peak m/z` <-
  distinct <-
  `.` <-
  `Peak status` <-
  `Peak RT` <-
  `Peak RT end` <-
  `Peak RT start` <-
  mz <-
  rtmax <-
  rtmin <-
  `Peak area` <-
  `Peak height` <-
  mzmin <-
  mzmax <-
  into <-
  intb <-
  maxo <-
  sn <-
  rt <-
  NULL 
  
  

  # Check that filepath was set
  if(is.null(filepath)){stop("filepath has to be filled and give the path to the data files files")}
  
  
  
  # Make an empty xcmsSet
  mzmine_xset <- new("xcmsSet")
  
  
  # completely made up fields
  mzmine_xset@profinfo  <- list(method="bin",step=0.005)
  mzmine_xset@scanrange <- numeric(0)
  
  if(.hasSlot(mzmine_xset, "mslevel")){ mzmine_xset@mslevel   <- numeric(0) } # some versions seem not to create this slot
  
  
  # Convert table to long format
  peaklist_mzmine_long <- peaklist_mzmine %>% gather(type,value,-matches("row ID|row m/z|row retention time|All identity elements|ID|Molecular formula|Name|Identification method")) %>% 
                          separate(type,c("sample_name","type"),sep = " ",extra = "merge") %>% 
                          mutate(type = str_replace_all(type, "filtered ", ""))
  
  
  # Filter, do some fixes, and convert back to wide format
  peaklist_mzmine_peaks <-  peaklist_mzmine_long %>% 
                            filter(type %in% c("Peak m/z","Peak RT","Peak RT start","Peak RT end","Peak height","Peak area","Peak FWHM","Peak tailing factor","Peak asymmetry factor","Peak status")) %>% 
                            select(`row ID`,`sample_name`,`type`,`value`) %>% 
                            mutate(type = str_replace_all(type, "sample_name", "sample")) %>% 
                            spread(type, value) %>% 
                            mutate(sample = as.numeric(factor(sample_name,levels=unique(sample_name)))) %>%
                            arrange(sample) %>% 
                            filter(`Peak m/z`!=0)
  
  
  # Get which peaks have been filled
  mzmine_xset@filled <- which(peaklist_mzmine_peaks$`Peak status` != "DETECTED")
  
  
  # take out filename so we have ordering
  mzmine_xset@filepaths  <- peaklist_mzmine_peaks %>% 
                            select(sample,sample_name) %>% 
                            distinct %>% 
                            arrange(sample) %>% 
                            select(sample_name) %>% 
                            as.matrix %>% as.character %>%
                            paste(collapse="|") %>% 
                            grep(list.files(filepath,full.names = TRUE), value=TRUE)
  
  
  # Make dummy sample grouping
  mzmine_xset@phenoData <- data.frame(class = factor(rep("dummy_class",length(mzmine_xset@filepaths))),row.names = sub("\\.[[:alnum:]]+$", "", basename(as.character(mzmine_xset@filepaths))))
  
  
  # get rt of each scan for each file to fill into @rt
  # Need to check if this can be disabled to save time if CAMERA doesn't need it
  # we are putting the same in "raw" and "corrected"
  if(get.scantimes){
    
    mzmine_xset@rt <- list(raw = vector("list",length(mzmine_xset@filepaths)),corrected = vector("list",length(mzmine_xset@filepaths)))
    
    for(i in 1:length(mzmine_xset@filepaths)){
      xraw <- xcmsRaw(mzmine_xset@filepaths[i],profstep = 0)
      mzmine_xset@rt$raw[[i]] <- mzmine_xset@rt$corrected[[i]] <- xraw@scantime
      names(mzmine_xset@rt$raw[[i]]) <- names(mzmine_xset@rt$corrected[[i]]) <- 1:length(xraw@scantime)
    }
    
  }
  
  
  # Take out group indexes while we still have the grouping info
  mzmine_xset@groupidx <- peaklist_mzmine_peaks %>% select(`row ID`) %>% 
                          mutate(row_nr = 1:n()) %>% 
                          as.matrix %>% 
                          split(x = .[,"row_nr"],f = .[,"row ID"])
  
  names(mzmine_xset@groupidx) <- NULL
  
  
  # get number of detected peaks per peak group
  npeaks <-   peaklist_mzmine_peaks %>% 
              select(`row ID`,`Peak status`) %>% 
              group_by(`row ID`) %>% 
              summarise(npeaks = sum(`Peak status`=="DETECTED")) %>% 
              select(npeaks) %>% as.matrix %>% as.integer
  
  
  # Create the final list of peaks
  as.numeric_quietly <- function(x){temp_fun <- quietly(as.numeric); temp_fun(x)$result }
  
  peaklist_mzmine_peaks %<>% 
                              select(-sample_name,-`row ID`,-`Peak status`) %>%
                              rename(mz = `Peak m/z`) %>%
                              rename(rt = `Peak RT`, rtmax = `Peak RT end`, rtmin = `Peak RT start`) %>% 
                              mutate_each(funs(as.numeric_quietly)) %>%
                              mutate(mzmin = mz-0.005,mzmax = mz+0.005) %>%
                              mutate(rt = rt*60, rtmax = rtmax*60, rtmin = rtmin*60) %>% 
                              mutate(into = `Peak area`, intb = `Peak area`, maxo = `Peak height`) %>% 
                              mutate(sn = NA, egauss = NA, mu = NA, sigma = NA, h = NA, f = NA, dppm = NA, scale = NA, scpos = NA, scmin = NA, scmax = NA, lmin = NA, lmax = NA) %>% 
                              select(mz,mzmin,mzmax,rt,rtmin,rtmax,into,intb,maxo,sn,egauss,mu,sigma,h,f,dppm,scale,scpos,scmin,scmax,lmin,lmax,sample) %>% 
                              mutate_each(funs(as.numeric)) %>% as.matrix
  
  mzmine_xset@peaks     <- peaklist_mzmine_peaks
  
  
  # make peak group stats
  mzmine_xset@groups <- lapply(mzmine_xset@groupidx, function(x){ data.frame( mzmed  = median(mzmine_xset@peaks[x,"mz"]),
                                                                              mzmin  = min(   mzmine_xset@peaks[x,"mzmin"]),
                                                                              mzmax  = max(   mzmine_xset@peaks[x,"mzmax"]),
                                                                              rtmed  = median(mzmine_xset@peaks[x,"rt"]),
                                                                              rtmin  = min(   mzmine_xset@peaks[x,"rtmin"]),
                                                                              rtmax  = max(   mzmine_xset@peaks[x,"rtmax"])
                                                                             )
                                                                }
                              ) %>% do.call(rbind,.) %>%
                                    mutate(npeaks = npeaks) %>% 
                                    mutate(dummy_class = npeaks) %>%  # I don't know what is actually supposed to go here
                                    as.matrix
  
  
  return(mzmine_xset)
}

