group_correct <- function (dataset,f,zero.is.na=F){
  
  # copy dataset
  dataset_fixed <- dataset
  
   # Replace na with zero
  if (zero.is.na) { 
    is.zero = dataset_fixed==0
    dataset_fixed[is.zero]=NA
  }
  
  
  split_order_reverse <- order(unlist(split(1:length(f),f)))
  
  
  new_matrix <- apply(dataset_fixed,2,function(x) {
    
    overall_mean=mean(x,na.rm=T)
    split_data = split(x,f)
    new_matrix=lapply(split_data,function(y){     (y/mean(y,na.rm=T))*overall_mean     })
    unlist(new_matrix)
    
  })
  
  new_matrix <- new_matrix[split_order_reverse,]
  
  colnames(new_matrix) <- colnames(dataset)
  rownames(new_matrix) <- rownames(dataset)

  if (zero.is.na) { 
    new_matrix[is.zero]=0
  }
  
  
  if(is.data.frame(dataset)) {new_matrix <- as.data.frame(new_matrix)}
  
  return(new_matrix)
  
}