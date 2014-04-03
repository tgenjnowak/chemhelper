
# Load data ###########################
#--------------------------------------#
load.camera.rules=function(mode){
  
  if (tolower(mode)=='pos') {ionmode=1}
  if (tolower(mode)=='neg') {ionmode=2}
  if (tolower(mode)=='ei')  {ionmode=3}
  
  
  rules <- read.xlsx2(system.file("extdata",paste('CAMERA_rules_',c('pos','neg','EI')[ionmode],'.xlsx',sep=""), package="chemhelper"),1)  
  
  
  
  
  # Check to find adduct/fragment combinations that would give the same difference. This can course bad annotation: for example +NH4 getting annotated instead of -NH3
  # Better to remove the most unlikely or lower the score for the least likely assignment
  bad_pairs=matrix(ncol=2,nrow=0)
  massdiff=as.numeric(rules[,"massdiff"])
  
  for (i in 1:nrow(rules)){
    for (i2 in 1:nrow(rules)){ 
      if (i<i2) {next}
      diff=abs(abs((abs(massdiff[i])-abs(massdiff[i2]))/2)-1.007276) < 0.001
      pos_neg = (massdiff[i] < 0 & massdiff[i2] > 0) | (massdiff[i] > 0 & massdiff[i2] < 0)
      
      if (diff & pos_neg){
        bad_pairs = rbind(bad_pairs,cbind(i,i2)) 
      }
      
    }
  }
  
  
  if (nrow(bad_pairs)!=0){
  bad_pairs_names=cbind(as.character(rules[bad_pairs[,1],"name"]),as.character(rules[bad_pairs[,2],"name"]))
  warning(
    paste("The following adducts/fragments seem to collide.","\n",
      paste(apply(bad_pairs_names,1, function(x) paste(x,collapse="\t")),collapse="\n"),"\n\n",
      "Consider removing one of them. Example:","\n",
      'rules=rules[            !grepl("[M+NH4]+",rules[,"name"],fixed=T)         ,]'
    ,sep="")
    )
  }
  
  
  
  
  
  
  
  return(rules)
}







#--------------------------------------#