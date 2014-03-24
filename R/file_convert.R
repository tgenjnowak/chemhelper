getExt=function (path){
  parts <- strsplit(path, '\\.')[[1]]
  last <- parts[length(parts)]
  last
}







rem_satellite_peaks=function(files,outdir){
  

  for (i in 1:length(files)){
   
    file=files[i]
    
    # openMS needs mzML files so we use msconvert to convert first of we don't have mzML
    if (!(getExt(files[i])=='mzML')){
      tempfile = tempfile(pattern = basename(sub("\\.[^.]*$", "", files[i]) ), tmpdir = outdir, fileext = ".mzXML")
      system(paste('msconvert --mzML "',files[i],'" -o "',outdir,'" --outfile ',basename(tempfile),'"',sep=""),intern=T)
      file=tempfile
    }
    
    
        # do the fixing    
    system(paste('SpectraFilterWindowMower -in "',file,'" -out "',outdir,'/',basename(sub("\\.[^.]*$", "", files[i])),'.mzML','" -algorithm:windowsize 0.3 -algorithm:peakcount 1',sep=""),intern=T)
    
    
    # delete temp mzML
    if (!(getExt(files[i])=='mzML')){
      unlink(tempfile)
    }
  }  
}


