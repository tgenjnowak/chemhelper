chemhelper
==========

R functions helpful in working with chemical data. A number of functions to help analyze metabolomics data processed with xcms/CAMERA.
*Many functions have been moved to the packages:* [*convert.waters.raw*](https://github.com/stanstrup/convert.waters.raw), [*obabel2R*](https://github.com/stanstrup/obabel2R), and [*massageR*](https://github.com/stanstrup/massageR).







#### Data conversion and fixes
* `merge_scan_events`: Merge even numbered scans with preceding odd numbered scan. This is used if for example the m/z range was split between different traces ("functions" in the lingo of some vendors).
* `rem_satellite_peaks`: Remove satellite/shoulder peaks (orbitrap artifacts).
* `orbifilter`, `xcmsRaw.orbifilter`: Remove satellite/shoulder peaks (orbitrap artifacts). Newer and better than the approach in `rem_satellite_peaks`.


#### XCMS/CAMERA helpers
* `load.camera.rules`: Extended list of adducts/fragments for CAMERA.
* `analyze.xcms.group`: For XCMS data. Why is my feature not listed? This function plots all peaks (in all samples) in a given rt / m/z area. The peak group boundaries are draw as rectangles. This can be used to debug peak groups that are erroneously split or merged.
* `db.comp.assign`: Have a database of compounds' rt and m/z you want to match to your dataset? This function does exactly that.
* `xsAnnotate_stats`: Gives some statistics from an xsAnnotate object annotated by CAMERA.


#### Convert compound identifiers
* `name2struc`,`pubchem2inchi`


#### Chemical calculations
* `mass_decompose`: Mass decomposition. Calculate the elementary compositions from an exact mass. This function is s wrapper for Rdisop's decomposeMass function, but gives more output and filtering options. The output is a table with diagnostic information for each possible formula.




### Installation
First install all dependencies as needed.
```R
install.packages("devtools")
install.packages("jsonlite")

library(devtools)
install_github("dgrapov/CTSgetR")
install_github("rajarshi/cdkr",subdir = "rpubchem") # Version from CRAN is currently outdated.

install.packages("xlsx")
install.packages("stringr")
install.packages("RCurl")
install.packages("scales")

source("http://bioconductor.org/biocLite.R")
biocLite()
library(BiocInstaller)

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
