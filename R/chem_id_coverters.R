# smile2inchi ##################################
smile2inchi =function(smile,verbose=F){  
  # calls obabel binary in system. openbabel need to be installed
  
  string=paste(as.character(smile),collapse='" -:"')
  output=system(      paste('obabel -ismi -:"',string,'" -oinchi',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  return(output)
}





# inchi2sdf ##################################
inchi2sdf =function(inchi,verbose=F){
  #require(ChemmineR)
  # calls obabel binary in system. openbabel need to be installed
  
  string=paste(as.character(inchi),collapse='" -:"')
  output=system(      paste('obabel -iinchi -:"',string,'" -osdf --gen2D',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  output=read.SDFset(output)
  
  return(output)
}






# inchi.keep.cont ##################################
inchi.keep.cont=function(inchi,verbose=F){
  # calls obabel binary in system. openbabel need to be installed
  
  string=paste(as.character(inchi),collapse='" -:"')
  output=system(      paste('obabel -iinchi -:"',string,'" -r -xr -oinchi',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  return(output) 
}





# load.camera.rules ##################################
load.camera.rules=function(mode){
    
  if (tolower(mode)=='pos') {ionmode=1}
  if (tolower(mode)=='neg') {ionmode=2}
  if (tolower(mode)=='ei')  {ionmode=3}
  
  
  rules <- read.xlsx(system.file("extdata",paste('CAMERA_rules_',c('pos','neg','EI')[ionmode],'.xlsx',sep=""), package="chemhelper"),1)  
  
}


