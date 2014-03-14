smile2inchi =function(smile){
  output = system(paste('obabel -ismi -:"',smile,'" -oinchi --unique /nostereo',sep=''),intern=T,ignore.stderr = T)
  return(output)
}