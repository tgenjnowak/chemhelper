# Converters ###########################
#--------------------------------------#

smile2inchi =function(smile,verbose=F){  
  # calls obabel binary in system. openbabel need to be installed
  
  
  string=paste(as.character(smile),collapse='" -:"')
  output=system(      paste('obabel -ismi -:"',string,'" -oinchi',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  return(output)
}


inchi2smile =function(inchi,verbose=F){  
  # calls obabel binary in system. openbabel need to be installed
  
  string=paste(as.character(inchi),collapse='" -:"')
  output=system(      paste('obabel -iinchi -:"',string,'" -osmi',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  return(output)
}




inchi2sdf =function(inchi,verbose=F){
  #require(ChemmineR)
  # calls obabel binary in system. openbabel need to be installed
  # Cannot handle an infinite number of inchiùs in one chunk. So we split it.
    
  
  inchi = split(inchi, ceiling(seq_along(inchi)/200))
  
  output = sapply(inchi,function(x) {
    string=paste(as.character(x),collapse='" -:"')
    system(      paste('obabel -iinchi -:"',string,'" -osdf --gen2D',sep='')      ,intern=T,ignore.stderr = !verbose)
  })
  
  output = as.character(unlist(output))
  output=read.SDFset(output)
  
  return(output)
}






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













#--------------------------------------#






# Modify molecules ###########################
#--------------------------------------#

inchi.keep.cont=function(inchi,verbose=F){
  # calls obabel binary in system. openbabel need to be installed
  
  output=character(length=length(inchi))
  loc_ok= grepl(pattern="InChI=",x=inchi)
  inchi_good=inchi[loc_ok]
  
  
  string=paste(as.character(inchi_good),collapse='" -:"')
  inchi_cleaned=system(      paste('obabel -iinchi -:"',string,'" -r -xr -oinchi',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  output[which(loc_ok)]=inchi_cleaned
  output[which(!loc_ok)]=inchi[which(!loc_ok)]
  
  
  return(output) 
}







inchi.rem.stereo =function(inchi,verbose=F){  
  # calls obabel binary in system. openbabel need to be installed
  # Cannot handle an infinite number of inchiùs in one chunk. So we split it.
  
  inchi = split(inchi, ceiling(seq_along(inchi)/200))
  
  
  output = sapply(inchi,function(x) {
  string=paste(as.character(x),collapse='" -:"')
  system(      paste('obabel -iinchi -:"',string,'" -oinchi -xT/nostereo',sep='')      ,intern=T,ignore.stderr = !verbose)
  })
  
  output = as.character(unlist(output))
  return(output)
}



inchi.rem.charges =function(inchi,verbose=F){  
  # calls obabel binary in system. openbabel need to be installed
  # Cannot handle an infinite number of inchiùs in one chunk. So we split it.
  
  inchi = split(inchi, ceiling(seq_along(inchi)/200))
  
  output = sapply(inchi,function(x) {
    string=paste(as.character(x),collapse='" -:"')
    system(      paste('obabel -iinchi -:"',string,'" -oinchi -xT/nochg',sep='')      ,intern=T,ignore.stderr = !verbose)
  })
  
  output = as.character(unlist(output))
  return(output)
  
}




#--------------------------------------#






