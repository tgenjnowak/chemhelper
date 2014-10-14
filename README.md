chemhelper
==========

R functions helpful in working with chemical data. A number of functions help analyze metabolomics data processed with xcms/CAMERA.
* Converts files from Waters .raw format to mzData (convert.waters(indir,outdir)). MassLynx need to be installed and masswolf need to be in path.
* Extended list of adducts/fragments for CAMERA (load.camera.rules)
* Annotated mz/rt from dataset by database of known compounds' mz/rt (db.comp.assign(mz,rt,comp_name_db,mz_db,rt_db,mzabs=0.01,ppm=15,ret_tol=30))
* For XCMS data: Why is my feature not listed? Try analyze.xcms.group too see if it picked in the files and if it was added to a group.
* pubchem2inchi
* inchi2sdf
* smile2inchi


### Installation
First install all dependencies as needed.
```R
install.packages("devtools")
install.packages("jsonlite")

library(devtools)
install_github("dgrapov/CTSgetR")
install_github("rajarshi/cdkr",subdir = "rpubchem") # Version from CRAN is currently out of date.

install.packages("xlsx")
install.packages("stringr")
install.packages("RCurl")
install.packages("scales")

source("http://bioconductor.org/biocLite.R")
biocLite()
library(BiocInstaller)

biocLite("ChemmineR")
biocLite("xcms")
biocLite("CAMERA")
biocLite("Rdisop")

install_github("sneumann/MetShot")
```


Then chemhelper.
```R
install_github("stanstrup/chemhelper")
library(chemhelper)
```


### Updating
```R
chemhelper.update()
```
