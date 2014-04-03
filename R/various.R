clipboard <- function(x, sep="\t", row.names=FALSE, col.names=TRUE){
  con <- pipe("xclip -selection clipboard -i", open="w")
  write.table(x, con, sep=sep, row.names=row.names, col.names=col.names)
  close(con)
}


getExt=function (path){
  parts <- strsplit(path, '\\.')[[1]]
  last <- parts[length(parts)]
  last
}




chemhelper.update=function(){ 
  detach("package:chemhelper", unload=TRUE)
  install_github(repo = "chemhelper", username = "stanstrup")
  library(chemhelper)
}