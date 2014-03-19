# Converters ###########################
#--------------------------------------#

smile2inchi =function(smile,verbose=F){  
  # calls obabel binary in system. openbabel need to be installed
  
  string=paste(as.character(smile),collapse='" -:"')
  output=system(      paste('obabel -ismi -:"',string,'" -oinchi',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  return(output)
}






inchi2sdf =function(inchi,verbose=F){
  #require(ChemmineR)
  # calls obabel binary in system. openbabel need to be installed
  
  string=paste(as.character(inchi),collapse='" -:"')
  output=system(      paste('obabel -iinchi -:"',string,'" -osdf --gen2D',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  output=read.SDFset(output)
  
  return(output)
}






name2struc =function(input_names,    input_pubchem=as.numeric(matrix(data=NA,nrow=length(input_names)))     ){
  
  
  nas_logi = is.na(input_pubchem)
  nas_idx = which(nas_logi)
  output=matrix(nrow=0,ncol=5)
  colnames(output)=c('org_row','input_name','output_name','pubchem_CID','inchi')
  
  if (length(nas_idx)==0){
    warning('No compound names left that doesn\'t already have pubchem ids',immediate. =T)
  }
  
  
  
  for (i in nas_idx){
      
    org_row            = i
    input_name         = as.character(input_names[i])
    pubchem_CID        = as.character(        CTSgetR(input_name,from='Chemical Name',to='PubChem CID',limit.values=F)[,'PubChem.CID']          )
    pubchem_CID        = as.numeric(unlist(str_split(pubchem_CID,',')))
    output_name        = as.character(        CTSgetR(pubchem_CID,to='Chemical Name',from='PubChem CID',limit.values=F)[,'Chemical.Name']      )
    inchi              = smile2inchi(         get.cid(pubchem_CID)[,'CanonicalSmiles']      )
    
    if (length(inchi)==0){inchi=''}
    
    output = rbind(output,           cbind(org_row,input_name,output_name,pubchem_CID,inchi)            )
    
    
  }
  
  return(output)
}







pubchem2inchi =function(cid,skip){  
    
  output=character(length=length(cid))
    
  cid_org =   suppressWarnings(     as.numeric(cid)     )
  cid = cid_org
  
  if (!missing('skip')) {cid[-skip]}
  
  cid_unique = unique(cid)
  cid_unique=cid_unique[!is.na(cid_unique)]
  
    
  for (i in cid_unique){
  inchi =      smile2inchi(         get.cid(i)[,'CanonicalSmiles']      )
  
  idx = i==cid_org
  if (!missing('skip')) {idx[skip]=F}  
  output[idx]=inchi
  }
    
  
  
  return(output)
}













#--------------------------------------#






# Modify molecules ###########################
#--------------------------------------#

inchi.keep.cont=function(inchi,verbose=F){
  # calls obabel binary in system. openbabel need to be installed
  
  string=paste(as.character(inchi),collapse='" -:"')
  output=system(      paste('obabel -iinchi -:"',string,'" -r -xr -oinchi',sep='')      ,intern=T,ignore.stderr = !verbose)
  
  return(output) 
}
#--------------------------------------#





# Load data ###########################
#--------------------------------------#
load.camera.rules=function(mode){
    
  if (tolower(mode)=='pos') {ionmode=1}
  if (tolower(mode)=='neg') {ionmode=2}
  if (tolower(mode)=='ei')  {ionmode=3}
  
  
  rules <- read.xlsx(system.file("extdata",paste('CAMERA_rules_',c('pos','neg','EI')[ionmode],'.xlsx',sep=""), package="chemhelper"),1)  
  
}
#--------------------------------------#

