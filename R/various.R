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
  install_github("stanstrup/chemhelper")
  library(chemhelper)
}




google_ss <- function(gid = NA, key = NA) 
{
  if (is.na(gid)) {stop("\nWorksheetnumber (gid) is missing\n")}
  if (is.na(key)) {stop("\nDocumentkey (key) is missing\n")}
  
  url <- getURL(paste("https://docs.google.com/spreadsheet/pub?key=", key,
                      "&single=true&gid=", gid, "&output=csv", sep = ""),
                cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))
  read.csv(textConnection(url), header = T, sep = ",")
}