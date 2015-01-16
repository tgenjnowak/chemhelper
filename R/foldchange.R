
# Function to make string for combination of factors
pasteC <- function(x,y) paste(x,y,sep=" - ")


# Function to calculate combination of the values of each factor
ratio.w.invert <- function(x,y){
  v <- x/y
  to_invert <- v<1
  
  if(sum(to_invert)>0){
    v[to_invert] <- -1/v[to_invert]
  }
  
  return(v)
}





fold.change <- function(MAT,f,aggr_FUN=mean,combi_FUN="/"){
  
  temp <- aggregate(. ~ class, data = cbind.data.frame(class=f,MAT), aggr_FUN)
  class_computed <- as.character(temp[,1])
  temp <- as.matrix(temp[,-1])
  combs <- t(combn(class_computed,2))
  rownames(temp) <- class_computed
  
  combi_FUN_vec <- Vectorize(combi_FUN)
  
  out <- apply(temp,2,function(x){
    v <- combi_FUN_vec(    x[combs[,1]]      ,        x[combs[,2]]    )
    names(v) <- pasteC(combs[,1],combs[,2])
    return(v)
  })
  
  
}
