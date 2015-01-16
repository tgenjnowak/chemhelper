
half.min    <- function(x) min(x, na.rm = TRUE)/2
rand.to.min <- function(x) runif(1, min=0, max=min(x, na.rm = TRUE))



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
  
  
  return(mat_new)
}
