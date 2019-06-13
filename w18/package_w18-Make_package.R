############################################################
## This script guides you through the process of building ##
##           a custom R package from scratch              ##
############################################################
## Author: Erick Sousa (29/07/2018)
## More info: https://github.com/klutometis/roxygen#roxygen2

## Packages needed
library(devtools)
library(roxygen2)

## Create package directory
parentDir = "/media/erick/OS/Bioinformatics/R" #yourDir
setwd(parentDir)
if (!dir.exists("./w18")){
  create("w18")
  } else { print ("El directorio ya existe.")
    }

## Put your code in the R folder that has been generated
## and use roxygen2 package to write the documentation
## in the same R script that contains the function(s)

## Process documentation
setwd("./w18")
document()

## Installation
# You need to run this from the parentDir that contains the w18 folder
setwd(parentDir)
install("w18")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(w18)

# Some testings that the package will effectively work
filterMatches()
?getWindowCond()
makeSuperwindows()


#########################################################################
## In case you need to do a clean re-install of the package run code below:
remove.packages('w18')
.rs.restartR()
setwd(parentDir)
install("w18")
pacman::p_load(w18)

