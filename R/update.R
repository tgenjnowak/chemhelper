

chemhelper.update=function(){ 
require(devtools)
detach("package:chemhelper", unload=TRUE)
install_github(repo = "chemhelper", username = "stanstrup")
library(chemhelper)
}