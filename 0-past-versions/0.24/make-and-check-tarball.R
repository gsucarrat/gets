####################################################
## This file makes and checks the tarball of the 
## gets package, and installs it
##
## First created 2 July 2019 by Genaro Sucarrat
##
## Requirements:
##
## - Reading and writing permissions to the working
##   directory that is used. To ensure this in
##   Windows, first right-click mouse in folder,
##   and then choose properties -> security ...etc..
##
## - That the work directory contains: the gets-devel
##   folder, and the gets-base-source.R and gets-
##   isat-source.R files
##
## Contents:
##
## 1 SET WORK DIRECTORY
## 2 CLEAN WORK-DIRECTORY AND WORKSPACE
## 3 BUILD THE GETS FOLDER
## 4 BUILD AND CHECK THE TARBALL
## 5 INSTALL PACKAGE
## 6 UPLOAD TO CRAN
##
####################################################


####################################################
## 1 SET WORK DIRECTORY
####################################################
  
##Set working directory:
#setwd(choose.dir())
#Examples in Windows:
#setwd("C:/Program files/R/R-devel/bin/")
#setwd("C:/Program files/R/R-3.5.3/bin/")
#setwd("C:/Program files/R/R-3.6.3/bin/")
#setwd("C:/Users/sucarrat/Documents/R/R-devel/bin")
setwd("C:/Users/sucarrat/Documents/R/gs/gets/devel/")


####################################################
## 2 CLEAN WORK-DIRECTORY AND WORKSPACE
####################################################

##delete files and folders from previous builds?
doDelete <- TRUE #TRUE or FALSE?

##clean work-directory:
if(doDelete){

  ##files and folders of the work-directory:
  fileNames <- dir()
  
  ##delete tarball, if it already exists:             
  toBeDeleted <- fileNames[ grep(".tar.gz", fileNames) ] 
  if(length(toBeDeleted)>0){ file.remove(toBeDeleted) }
  
  ##delete folders "gets", "gets-skeleton" and "gets.Rcheck",
  ##if they already exist:
  toBeDeleted <- intersect(c("gets","gets-skeleton","gets.Rcheck"),
    fileNames)
  if(length(toBeDeleted)>0){
    for(i in toBeDeleted){
      unlink(i, recursive=TRUE) #delete folder+its content
    }
  }

} #end if(doDelete)

##clean workspace:
rm(list = ls())


####################################################
## 3 BUILD THE GETS FOLDER
####################################################

##Load sources:
source("gets-base-source.R")
source("gets-isat-source.R")

##make skeleton (i.e. folder structure w/files):
package.skeleton(name="gets")

##rename "gets" folder to "gets-skeleton":
file.rename("gets", "gets-skeleton")

##create new gets folder:
dir.create("gets")

##copy files from gets-devel to gets:
fileNames <- dir("./gets-devel/")
for(i in 1:length(fileNames)){
  file.copy( paste0("./gets-devel/", fileNames[i]), "gets", recursive=TRUE)
}

##copy files from skeleton to gets:
fileNames <- dir("./gets-skeleton/R")
fileNames <- setdiff(fileNames, "gets-internal.R") #remove gets-internal.R
for(i in 1:length(fileNames)){
  file.copy( paste0("./gets-skeleton/R/", fileNames[i]), "./gets/R/")
}

## Remember to check, manually, that version and date are correct in:
## - DESCRIPTION
## - NEWS
## - /R/gets-internal.R (start-up message)
## - /man/gets-package.Rd


####################################################
## 4 BUILD AND CHECK THE TARBALL
####################################################

##build tarball:
system("R CMD build gets --resave-data")
## - The --resave-data option is recommended for better compression.
##   However, it is not obligatory.
## - In principle, the latest development version of R should be
##   used for the build, but this sometimes leads to spurious errors.

##check tarball (needs internet):
fileNames <- dir()
tarballWhere <- grep(".tar.gz", fileNames)
tarballName <- fileNames[tarballWhere]
system( paste0("R CMD check ", tarballName, " --as-cran") )
#system( paste0("R CMD check ", tarballName) )
## Note: The --as-cran option is obligatory according to cran policy.
## Check also the user manual for line exceedances in the .Rd files.
## This is not detected by the tarball check. An alternative way to
## check tarball, see: http://win-builder.r-project.org/


####################################################
## 5 INSTALL PACKAGE
####################################################

##remove old version first:
remove.packages("gets")

##install new version:
system( paste0("R CMD INSTALL ", tarballName) )
#system("R CMD INSTALL gets_0.24.tar.gz")
#system("R CMD INSTALL --build gets")


####################################################
## 6 UPLOAD TO CRAN
####################################################

## This must be done by the maintainer of the package.
##
## When all the testing is done, upload the tarball, that is, the
## ?.tar.gz? file, see "Writing R extensions" for the current
## submission guidelines. Upload package via:
## https://xmpalantir.wu.ac.at/cransubmit/
