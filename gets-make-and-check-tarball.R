####################################################
## This file makes and checks the tarball of the 
## gets package.
##
## First created 2 July 2019 by Genaro Sucarrat
##
## Requirements:
##
## - reading and writing permissions to the working
##   directory that is used
## - that the work directory contains the gets-devel
##   folder, and the gets-base-source-devel.R and
##   gets-isat-source-devel.R files
##
####################################################

##initiate:
##=========
 
##clean workspace:
rm(list = ls())

##Set working directory:
setwd(choose.dir())
#Examples in Windows:
#setwd("C:/Program files/R/R-devel/bin/")
#setwd("C:/Program files/R/R-3.5.3/bin/")

##Load sources
source("gets-base-source-devel.R")
source("gets-isat-source-devel.R")

##Note: On Windows, make sure you have writing and modifying
##permission in the working directory: First (right-click
##mouse in folder, choose properties -> security ...etc.).
##Subsequent steps will not work withoug writing and modyfing
##permissions.

##gets, gets-devel, gets-skeleton:
##================================

##make skeleton (i.e. folder structure w/files):
package.skeleton(name="gets")

##copy files, make gets folder, etc.:
file.rename("gets", "gets-skeleton")
dir.create("gets")
fileNames <- dir("./gets-devel/")
for(i in 1:length(fileNames)){
  file.copy( paste0("./gets-devel/", fileNames[i]), "gets", recursive=TRUE)
}

fileNames <- dir("./gets-skeleton/R")
fileNames <- setdiff(fileNames, "gets-internal.R") #remove gets-internal.R
for(i in 1:length(fileNames)){
  file.copy( paste0("./gets-skeleton/R/", fileNames[i]), "./gets/R/")
}

## Note: Check that version is correct in gets-internal.R,
## DESCRIPTION, start-up message and gets-package.RD

##build and check tarball:
##========================

##build tarball:
system("R CMD build gets --resave-data")
## The --resave-data option is recommended for better compression.
## However, it is not obligatory. In principle, the latest
## development version of R should be used for the build.

##check tarball (needs internet):
fileNames <- dir()
tarballWhere <- grep(".tar.gz", fileNames)
tarballName <- fileNames[tarballWhere]
system( paste0("R CMD check ", tarballName, " --as-cran") )
## Note: The --as-cran option is obligatory according to cran policy.
## Check also the user manual for line exceedances in the .Rd files.
## This is not detected by the tarball check. An alternative way to
## check tarball, see: http://win-builder.r-project.org/

##If desirable, install package:
##==============================

system( paste0("R CMD INSTALL ", tarballName) )
#system("R CMD INSTALL gets_0.19.tar.gz")
#system("R CMD INSTALL --build gets")

##Upload to CRAN:
##===============

##When all the testing is done, upload the tarball, that is, the
##?.tar.gz? file, see "Writing R extensions" for the current
## submission guidelines. Upload package via:
## https://xmpalantir.wu.ac.at/cransubmit/
