mass_decompose = function(mass,mode="neutral",elements=c('C','H','N','O','P','S'),ppm=30, filter=TRUE) 
{
  
  
  
  # Get the neutral mass
  if(tolower(mode)=="pos"){mode=1}
  if(tolower(mode)=="neg"){mode=2}
  if(tolower(mode)=="neutral"){mode=3}
  
  mass=mass+1.0073*c(-1,1,0)[mode]
  
  
  # Shortcut to allow more masses
  if(tolower(elements)[1]=="expand"){
    elements=c('C','H','N','O','P','S','Na','K','Cl','Br','Fe','Al','Cu','Zn','As')
  }
  
  
  # Do the decomposing
  elements=initializeElements(elements)
  formulas=decomposeMass(mass,ppm=30,elements=elements)
  
  # Make nice table of results
  formulas=cbind(getFormula(formulas),getValid(formulas),formulas$DBE,signif(getMass(formulas)+1.0073*c(1,-1,0)[mode],digits=7), signif((getMass(formulas)-(mass))/(mass)*1e6,digits=3), signif(getScore(formulas),digits=3))
  
  colnames(formulas)=c("Formula","Valid","DBE","Calc. m/z","ppm","Score")
  
  
  if(filter){
    formulas=formulas[   floor(as.numeric(formulas[,"DBE"])) == as.numeric(formulas[,"DBE"]) & as.numeric(formulas[,"DBE"])>-1 & formulas[,"Valid"]=="Valid"       ,,drop=F]
  }
  
  
  
  return(formulas)
  
  
  
  
  
}