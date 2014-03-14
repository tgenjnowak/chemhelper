smile2inchi =function(smile){
  
  string=paste(as.character(smile),collapse='" -:"')
  output=system(      paste('obabel -ismi -:"',string,'" -oinchi --unique /nostereo',sep='')      ,intern=T,ignore.stderr = T)
  
  return(output)
}