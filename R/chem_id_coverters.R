smile2inchi =function(smile,verbose=F){  
  # calls obabel binary in system. openbabel need to be installed
  
  string=paste(as.character(smile),collapse='" -:"')
  output=system(      paste('obabel -ismi -:"',string,'" -oinchi',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  return(output)
}




inchi2sdf =function(inchi,verbose=F){
  #require(ChemmineR)
  # calls obabel binary in system. openbabel need to be installed
  
  string=paste(as.character(inchi),collapse='" -:"')
  output=system(      paste('obabel -iinchi -:"',string,'" -osdf --gen2D',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  output=read.SDFset(output)

  return(output)
}


inchi.keep.cont=function(inchi,verbose=F){
  # calls obabel binary in system. openbabel need to be installed
  
  string=paste(as.character(inchi),collapse='" -:"')
  output=system(      paste('obabel -iinchi -:"',string,'" -r -xr -oinchi',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  return(output) 
}
