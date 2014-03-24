getExt=function (path){
  parts <- strsplit(path, '\\.')[[1]]
  last <- parts[length(parts)]
  last
}







rem_satellite_peaks=function(files,outdir){
  

  for (i in 1:length(files)){
   
    file=files[i]
    
    # openMS needs mzML files
    if (!(getExt(files[i])=='mzML')){
      tempfile = tempfile(pattern = basename(sub("\\.[^.]*$", "", files[i]) ), tmpdir = outdir, fileext = ".mzXML")
      
      system(paste('msconvert --mzML "',files[i],'" -o "',outdir,'"',sep=""),intern=T)

      
    }
    
    
    
    
  }  
}





#C:\Program Files\OpenMS-1.11\bin>SpectraFilterWindowMower.exe -in "D:\temp\QC_Bo
#x6_AFTER_1.mzML" -out "D:\temp\test.mzML"  -algorithm:windowsize 0.3 -algorithm:
  #peakcount 1