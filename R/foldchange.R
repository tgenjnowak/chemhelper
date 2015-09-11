
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






fold.change <- function(MAT,f,aggr_FUN=mean,combi_FUN = function(x,y) "/"(x,y)    ,method=5){
  
  if(method==1){  out <-  fold.change1(MAT=MAT,f=f,aggr_FUN=aggr_FUN,combi_FUN=combi_FUN)       }
  if(method==2){  out <-  fold.change2(MAT=MAT,f=f,aggr_FUN=aggr_FUN,combi_FUN=combi_FUN)       }
  if(method==3){  out <-  fold.change3(MAT=MAT,f=f,aggr_FUN=aggr_FUN,combi_FUN=combi_FUN)       }
  if(method==4){  out <-  fold.change4(MAT=MAT,f=f,aggr_FUN=aggr_FUN,combi_FUN=combi_FUN)       }
  if(method==5){  out <-  fold.change5(MAT=MAT,f=f,aggr_FUN=aggr_FUN,combi_FUN=combi_FUN)       }
  if(method==6){  out <-  fold.change6(MAT=MAT,f=f,aggr_FUN=aggr_FUN,combi_FUN=combi_FUN)       }
  
  return(out)
}










fold.change1 <- function(MAT,f,aggr_FUN=mean,combi_FUN=function(x,y) "/"(x,y)   ){
  
  f_un <- unique(f)
  temp <- matrix(NA,nrow = length(f_un),ncol=ncol(MAT))
  rownames(temp) <- f_un
  
  for(i in 1:length(f_un)){
    temp[i,] <- apply(MAT[f_un[i] == f,,drop=FALSE],2,aggr_FUN)
  }
  
  
  combs <- t(combn(as.character(f_un),2))
  combi_FUN_vec <- Vectorize(combi_FUN)
  
  out <- matrix(NA,nrow = nrow(combs),ncol=ncol(temp))
  rownames(out) <- pasteC(combs[,1],combs[,2])
  colnames(out) <- 1:ncol(temp)
  
  for( i in 1:nrow(combs)){
    out[i,] <- combi_FUN_vec(    temp[combs[i,1],]      ,         temp[combs[i,2],]   )
  }
  
  
  colnames(out) <- colnames(MAT)
  return(out)
}





fold.change2 <- function(MAT,f,aggr_FUN=mean,combi_FUN=function(x,y) "/"(x,y)    ){
  
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
  
  
  colnames(out) <- colnames(MAT)
  return(out)
}




fold.change3 <- function(MAT,f,aggr_FUN=mean,combi_FUN=function(x,y) "/"(x,y)     ){
  
  temp <- recast(data.frame(class=f,MAT),class ~ variable,id.var="class",aggr_FUN)
  
  
  class_computed <- as.character(temp[,1])
  temp <- as.matrix(temp[,-1])
  combs <- t(combn(class_computed,2))
  rownames(temp) <- class_computed
  
  combi_FUN_vec <- Vectorize(combi_FUN)
  
  out <- matrix(NA,nrow = nrow(combs),ncol=ncol(temp))
  rownames(out) <- pasteC(combs[,1],combs[,2])
  colnames(out) <- 1:ncol(temp)
  
  for( i in 1:nrow(combs)){
    out[i,] <- combi_FUN_vec(    temp[combs[i,1],]      ,         temp[combs[i,2],]   )
  }
  
  
  colnames(out) <- colnames(MAT)
  return(out)
}









fold.change4 <- function(MAT,f,aggr_FUN=mean,combi_FUN=function(x,y) "/"(x,y)     ){
  
  # mean by purrr package
  temp <- data.frame(class = f, MAT) %>%
    slice_rows("class") %>%
    by_slice(map, aggr_FUN)
  
  
  class_computed <- as.character(as.matrix(temp[,1]))
  temp <- as.matrix(temp[,-1])
  combs <- t(combn(class_computed,2))
  rownames(temp) <- class_computed
  
  combi_FUN_vec <- Vectorize(combi_FUN)
  
  out <- matrix(NA,nrow = nrow(combs),ncol=ncol(temp))
  rownames(out) <- pasteC(combs[,1],combs[,2])
  colnames(out) <- 1:ncol(temp)
  
  for( i in 1:nrow(combs)){
    out[i,] <- combi_FUN_vec(    temp[combs[i,1],]      ,         temp[combs[i,2],]   )
  }
  
  colnames(out) <- colnames(MAT)
  return(out)
}






fold.change5 <- function(MAT,f,aggr_FUN=mean,combi_FUN=function(x,y) "/"(x,y)   ){
  
  # mean using purrr
  x <- data.frame(class = f, MAT) %>%  slice_rows("class") %>% by_slice(map, aggr_FUN)
  rownames <- as.character(as.data.frame(x[,1])[,1])
  x <- as.matrix(x[,-1])
  rownames(x) <- rownames
  
  # calculate changes between all rows
  i <- combn(unique(f), 2)
  ret <- combi_FUN(x[i[1,],] , x[i[2,],])
  rownames(ret) <- pasteC(i[1,], i[2,])
  
  # Put original colnames
  colnames(ret) <- colnames(MAT)
  
  return(ret)
}







fold.change6 <- function(MAT,f,aggr_FUN=mean,combi_FUN=function(x,y) "/"(x,y)   ){
  
  f_un <- unique(f)
  combs <- t(combn(as.character(f_un),2))
  
  out3 <- data.frame(class = f, MAT) %>%  slice_rows("class") %>% by_slice(map, aggr_FUN) %>% 
                                          do(combi_FUN( slice(.,match(combs[,1], class))[,-1]  ,slice(.,match(combs[,2], class))[,-1]      )) %>% 
                                          as.matrix
  
  # Put original dimnames
  rownames(out3) <- pasteC(combs[,1],combs[,2])
  colnames(out3) <- colnames(MAT)
  
  return(out3)
}

