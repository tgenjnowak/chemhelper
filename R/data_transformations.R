group_correct=function (dataset,f,zero.is.na=F){
  
  dataset_fixed=dataset
  
  if (zero.is.na) { 
    is.zero = dataset_fixed==0
    dataset_fixed[is.zero]=NA
  }
  
  
  new_matrix = apply(dataset_fixed,2,function(x) {
    
    overall_mean=mean(x,na.rm=T)
    split_data = split(x,f)
    new_matrix=lapply(split_data,function(y){     (y/mean(y,na.rm=T))*overall_mean     })
    unlist(new_matrix)
    
  })
  
  colnames(new_matrix) = colnames(dataset)
  rownames(new_matrix) = rownames(dataset)

  if (zero.is.na) { 
    new_matrix[is.zero]=0
  }
  
  
  if(is.data.frame(dataset)) {new_matrix=as.data.frame(new_matrix)}
  
  return(new_matrix)
  
}