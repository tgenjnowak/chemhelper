smile2inchi =function(smile){
  # calls obabel binary in system. openbabel need to be installed
  
  string=paste(as.character(smile),collapse='" -:"')
  output=system(      paste('obabel -ismi -:"',string,'" -oinchi --unique /nostereo',sep='')      ,intern=T,ignore.stderr = T)
  
  return(output)
}