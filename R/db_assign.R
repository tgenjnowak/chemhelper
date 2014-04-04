db.comp.assign=function(mz,rt,comp_name_db,mz_db,rt_db,mzabs=0.01,ppm=15,ret_tol=30){



#Make vector to put IDed compounds into
ID_list = vector("list", length(mz))


result_summary = vector("numeric", length(comp_name_db))
mul_assign = vector("numeric", length(comp_name_db))


# Find the marker and putting in the name in the dataset
for (i in 1:length(comp_name_db)) {
  # Finding where the compounds are
  hit = (mz > (mz_db[i] - mzabs)) & 
    (mz < (mz_db[i] + mzabs)) & 
    (abs((mz-mz_db[i])/mz_db[i]) < mzabs) &
    (rt > (rt_db[i] - ret_tol)) & 
    (rt < (rt_db[i] + ret_tol))
  
  
  # Did we find it and is it unique?
  result = sum(hit)
  result_summary[i] = result
  
  
  # Make string. Add "non-unique" if we cannot get a unique hit for the
  # compound.
  if (result > 1){ # More than one hit. We don't know which is right
    ID = paste(comp_name_db[i],'(non-unique hit)')
  }
  
  if (result == 1){ # YES! we know where the marker is
    ID = comp_name_db[i]
  }
  
  
  if (result != 0){
    
    
    # Where should we place the marker?
    hit_loc = which(hit)
    
    for (i_2 in 1:length(hit_loc)){
      current_hit_loc = hit_loc[i_2]
      
      if (    is.null(ID_list[[current_hit_loc]])     ) {
        ID_list[[current_hit_loc]] = ID
        mul_assign[i] = 0
      }else{
        ID_list[[current_hit_loc]] = paste(ID_list[[current_hit_loc]],' OR ',ID)
        mul_assign[i] = 1
      }
    }
    
    
  }
  
  
  
}
rm(i, ID, comp_name_db, current_hit_loc,hit,hit_loc,i_2,mz_db,result)





## Give a report 
assign_summary = matrix(data=0,ncol=1, nrow=4)
assign_summary[1,1] = sum(result_summary == 1)
assign_summary[2,1] = sum(result_summary > 1)
assign_summary[3,1] = sum(result_summary == 0)
assign_summary[4,1] = sum(mul_assign)
rm(mul_assign, result_summary, rt_db)

rownames(assign_summary) = c('Unique hits'   ,   'Non-unique hits'  ,   'Compounds not found'    ,   'Markers had multiple compounds assigned' )
print(assign_summary)
rm(assign_summary)





## Give a list of the IDed compounds and fix ID list
choose = sapply(ID_list,function(x) !is.null(x))
hit_list = cbind(mz[choose],rt[choose],ID_list[choose])

ID_list = as.character(ID_list)
ID_list[!choose]=""
rm(choose)

return(ID_list)





}