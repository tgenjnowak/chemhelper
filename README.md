chemhelper
==========

R functions helpful in working with chemical data. A number of functions help analyze metabolomics data processed with xcms/CAMERA.

#### Data conversion and fixes
* `convert.waters`: Converts files from Waters .raw format to mzData. MassLynx need to be installed and masswolf need to be in path. (this works around the problem of properly converting Waters data described in the supplementary of dx.doi.org/10.1007/s00216-013-6954-6).
* `merge_scan_events`: Merge even numbered scans with preceding odd numbered scan. This is used if for example the m/z range was split between different traces ("functions" in the lingo of some vendors).
* `rem_satellite_peaks`: Remove satellite/shoulder peaks (orbitrap artifacts).

#### XCMS/CAMERA helpers
* `load.camera.rules`: Extended list of adducts/fragments for CAMERA.
* `analyze.xcms.group`: For XCMS data. Why is my feature not listed? This function plots all peaks (in all samples) in a given rt / m/z area. The peak group boundaries are draw as rectangles. This can be used to debug peak groups that are erroneously split or merged.
* `db.comp.assign`: Have a database of compounds' rt and m/z you want to match to your dataset? This function does exactly that.


#### Convert compound identifiers
* `name2struc`,`pubchem2inchi`,`inchi2sdf`,`smile2inchi`,`inchi2smile`
* `inchi.keep.cont`: Keep only largest continues part of a molecule InChI (AKA remove salts).
* `inchi.rem.charges`: Remove charges from a molecule InChI.
* `inchi.rem.stereo`: Remove stereochemistry from a molecule InChI.


#### Data massaging
* `group_correct`: Batch-correct dataset feature-wise.

#### Data analysis
* `mass_decompose`: Mass decomposition. Calculate the elementary compositions from an exact mass. This function is s wrapper for Rdisop's decomposeMass function, but gives more output and filtering options. The output is a table with diagnostic information for each possible formula.


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
