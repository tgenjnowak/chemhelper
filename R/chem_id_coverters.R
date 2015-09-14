name2struc =function(input_names, input_pubchem=as.numeric(matrix(data=NA,nrow=length(input_names)))     ){
  
  if(!(packageVersion("rpubchem")>="1.5.0.3")){
  stop('At least version 1.5.0.3 of rpubchem is required for this function to work.
       Latest bleading edge version can be install with:
       install_github(repo="cdkr",username = "rajarshi",subdir = "rpubchem")
       ')
  }
  
  
  nas_logi = is.na(input_pubchem)
  nas_idx = which(nas_logi)
  output=matrix(nrow=0,ncol=5)
  colnames(output)=c('org_row','input_name','output_name','pubchem_CID','inchi')
  
  if (length(nas_idx)==0){
    warning('No compound names left that doesn\'t already have pubchem ids',immediate. =T)
  }
  
  
  nas_idx = nas_idx[    !duplicated(input_names[nas_idx])     ]
  
  
  
  for (i in nas_idx){
    
    org_row            = i
    input_name         = as.character(input_names[i])
    pubchem_CID        = as.character(        CTSgetR(input_name,from='Chemical Name',to='PubChem CID',limit.values=F)[,'PubChem CID']          )
    pubchem_CID        = as.numeric(unlist(str_split(pubchem_CID,',')))
    output_name        = as.character(        CTSgetR(pubchem_CID,to='Chemical Name',from='PubChem CID',limit.values=F)[,'Chemical Name']      )
    
    smiles             = try(     get.cid(pubchem_CID)  , silent = TRUE)
    
    if (   inherits(smiles, "try-error") )   {
      cat("\n Lookup failed. Re-trying in 30sec... \n")
      Sys.sleep(30)
      smiles        =     get.cid(pubchem_CID)  
    }
    
    
    smiles=smiles[,'CanonicalSmile']
    inchi              = smile2inchi(     smiles          )
    
    
    
    if (length(inchi)==0){inchi=''}
    
    output = rbind(output,           cbind(org_row,input_name,output_name,pubchem_CID,inchi)            ) 
  }
  
  
  return(output)
}







pubchem2inchi <- function(cid,skip=NULL,silent=T){  
  
  output=character(length=length(cid))
  
  cid_org =   suppressWarnings(     as.numeric(cid)     )
  cid = cid_org
  
  if (   !missing('skip') & !(length(skip)==0)   ) {cid=cid[-skip]}
  
  cid_unique = unique(cid)
  cid_unique=cid_unique[!is.na(cid_unique)]
  
  
  for (i in cid_unique){
    
    if (!(silent==T)){
      cat(paste("Looking up pubchem CID ",i," (",which(i==cid_unique)," of ",length(cid_unique),")","\n",sep=""))
    }
    
    
    
    smiles             = try(     get.cid(i)  , silent = TRUE)
    
    if (   inherits(smiles, "try-error") )   {
      cat("\n Lookup failed. Re-trying in 30sec... \n")
      Sys.sleep(30)
      smiles             = try(     get.cid(i)  , silent = TRUE)
    }
    
    if (   inherits(smiles, "try-error") )   {
      cat("\n Lookup failed again. Re-trying in 60sec... \n")
      Sys.sleep(60)
      smiles             = try(     get.cid(i)  , silent = TRUE)     
    }
    
    if (   inherits(smiles, "try-error") )   {
      cat("\n Lookup failed again. Re-trying in 90sec... \n")
      Sys.sleep(90)
      smiles           =   get.cid(i)      
    }
    
    smiles=smiles[,'CanonicalSmile']
    inchi              = smile2inchi(     smiles          )
      
        
    
    idx = i==cid_org
    if (     !missing('skip') & !(length(skip)==0)    ) {idx[skip]=F}  
    output[idx]=inchi
  }
  
  
  
  return(output)
}








