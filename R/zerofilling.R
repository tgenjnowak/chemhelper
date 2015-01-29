
half.min    <- function(x){
  if(all(is.nan(x))) return(NaN)
  if(all(is.na(x))) return(NA)
  
  out <- min(x, na.rm = TRUE)/2
  return(out)
}

rand.to.min <- function(x){
  if(all(is.nan(x))) return(NaN)
  if(all(is.na(x))) return(NA)
  
  out <- runif(1, min=0, max=min(x, na.rm = TRUE))
  return(out)
}



rep.na <- function(mat,zero.as.na=TRUE,rep.FUN=half.min) {
  
  mat_new <- mat
  
  if(zero.as.na) {
    mat_new[mat_new==0] <- NA
  }
  
  
  mat_new <- apply(mat_new,2,function(x) {
    out <- x  
    to_rep <- is.na(out)
    out[to_rep] <- replicate(   sum(to_rep)  ,  rep.FUN(out)   )
    return(as.numeric(out))
  })
  
  rownames(mat_new) <- rownames(mat)
  colnames(mat_new) <- colnames(mat)
  
  return(mat_new)
}
