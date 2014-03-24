
# Load data ###########################
#--------------------------------------#
load.camera.rules=function(mode){
  
  if (tolower(mode)=='pos') {ionmode=1}
  if (tolower(mode)=='neg') {ionmode=2}
  if (tolower(mode)=='ei')  {ionmode=3}
  
  
  rules <- read.xlsx(system.file("extdata",paste('CAMERA_rules_',c('pos','neg','EI')[ionmode],'.xlsx',sep=""), package="chemhelper"),1)  
  
}
#--------------------------------------#