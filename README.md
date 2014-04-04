chemhelper
==========

R functions helpful in working with chemical data.
* Converts files from Waters .raw format to mzData (convert.waters(indir,outdir)). MassLynx need to be installed and masswolf need to be in path.
* Extended list of adducts/fragments for CAMERA (load.camera.rules)
* Annotated mz/rt from dataset by database of known compounds' mz/rt (db.comp.assign(mz,rt,comp_name_db,mz_db,rt_db,mzabs=0.01,ppm=15,ret_tol=30))
* pubchem2inchi
* inchi2sdf
* smile2inchi


### Installation
```R
install.packages("devtools")
library(devtools)
install_github(repo = "chemhelper", username = "stanstrup")
library(chemhelper)
```



### Updating
```R
chemhelper.update()
```