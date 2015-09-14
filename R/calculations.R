


formula2adduct = function(formula,add,sub){
  
  elements=c(  initializePSE(),initializeCharges()   )
  formulas=vector("character",length(formula))
  
  for (i in 1:length(formula)){
  formulas[i] = addMolecules(formula[i],add,elements=elements)$formula
  formulas[i] = subMolecules(formulas[i],sub,elements=elements)$formula
  }
  
  return(formulas)

}




mass_decompose = function(mass,intensities=NULL,ion="neutral",elements=c('C','H','N','O','P','S'),ppm=30, filter.DBE=TRUE,filter.nitrogen=TRUE,simplify=TRUE, minElements="C0", maxElements="C999999") 
{
  
  elements=initializeElements(elements)
  formulas_list=vector("list", length(ion))
  
  for(i in 1:length(ion)){
    
    
    offset <- 
      switch(tolower(ion[i]),
             "neutral" = 0,
             "pos" = 1.007276455,
             "neg" = -1.007276455,
             "[m+h]+" = 1.007276455,
             "[m-h]-" = -1.007276455,
             
             "[m+na]+" = 22.989218,
             "[m+k]+" = 38.963158,
             
             "[m-2h+na]-" = 20.974666,
             "[m-2h+k]-" = 36.948606,
             "[m+cl]-" = 34.969402,
             "[m-h+hcoona]-" = 66.980164,
             "[m-h+hcooh]-" = 44.998194
      )
    
    
    
    add <- 
      switch(tolower(ion[i]),
             "neutral" = "H0",
             "pos" = "H+",
             "neg" = "-",
             "[m+h]+" = "H+",
             "[m-h]-" = "-",
             
             "[m+na]+" = "Na+",
             "[m+k]+" = "K+",
             
             "[m-2h+na]-" = "Na-",
             "[m-2h+k]-" = "K-",
             "[m+cl]-" = "Cl-",
             "[m-h+hcoona]-" = "HCOONa-",
             "[m-h+hcooh]-" = "HCOOH-"
      )
    
    
    
    sub <- 
      switch(tolower(ion[i]),
             "neutral" = "H0",
             "pos" = "H0",
             "neg" = "H",
             "[m+h]+" = "H0",
             "[m-h]-" = "H",
             
             "[m+na]+" = "H0",
             "[m+k]+" = "H0",
             
             "[m-2h+na]-" = "H2",
             "[m-2h+k]-" = "H2",
             "[m+cl]-" = "H0",
             "[m-h+hcoona]-" = "H",
             "[m-h+hcooh]-" = "H"
      )
    
    
    
    
    
    
    if(is.null(offset)){stop("The selected ion is not supported")}
    
    # Get the neutral mass
    mass_offsetted=mass-offset
    
    
    # Shortcut to allow more masses
    if(tolower(elements)[1]=="expand"){
      elements=initializeElements(c('C','H','N','O','P','S','Na','K','Cl','Br','Fe','Al','Cu','Zn','As'))
    }
    
    
    # Do the decomposing
    if(is.null(intensities)){
    formulas=decomposeMass(mass_offsetted,ppm=ppm,elements=elements,mzabs = 0, minElements=minElements, maxElements=maxElements)
    }else{
    formulas=decomposeIsotopes(mass_offsetted,intensities,ppm=ppm,elements=elements,mzabs = 0, minElements=minElements, maxElements=maxElements)      
    }
    
    # is nothing found then continue
    if (is.null(formulas)){next}
    
    
    # Make nice table of results
    formulas=cbind(getFormula(formulas),getValid(formulas),formulas$DBE,signif(getMass(formulas)+offset,digits=7), signif((getMass(formulas)-(mass_offsetted[1]))/(mass_offsetted[1])*1e6,digits=3), signif(getScore(formulas),digits=3))
    
    colnames(formulas)=c("Formula","Nitrogen rule","DBE","Calc. m/z","ppm","Rdisop score")
    
    
    if(filter.DBE){
      formulas=formulas[   floor(as.numeric(formulas[,"DBE"])) == as.numeric(formulas[,"DBE"]) & as.numeric(formulas[,"DBE"])>-1        ,,drop=F]
    }

    if(filter.nitrogen){
      formulas=formulas[  formulas[,"Nitrogen rule"]=="Valid"       ,,drop=F]
    }
    
    
    # Add adduct formulas
    if(nrow(formulas)>0){
      formulas = cbind(formulas            ,    formula2adduct(formulas[,"Formula"],add,sub)     )
      colnames(formulas)=c("Formula","Nitrogen rule","DBE","Calc. m/z","ppm","Rdisop score","Adduct formula")
    }
    formulas_list[[i]]=formulas
  }
  
  
  if(simplify){
    if(length(formulas_list)==1){formulas_list=formulas_list[[1]]}
  }
  
  return(formulas_list)
  
}





