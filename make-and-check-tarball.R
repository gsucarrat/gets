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
## - The work directory contains: the gets folder
##
## Contents:
##
## 1 SET WORK DIRECTORY
## 2 CLEAN WORK-DIRECTORY AND WORKSPACE
## 3 BUILD AND CHECK THE TARBALL
## 4 INSTALL PACKAGE
## 5 UPLOAD TO CRAN
##
####################################################


####################################################
## 1 SET WORK DIRECTORY
####################################################
  
##Set working directory:
#setwd(choose.dir())
#Examples in Windows:
#setwd("C:/Program files/R/R-devel/bin/")
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
  if( length(toBeDeleted)>0 ){ file.remove(toBeDeleted) }
  
  ##delete "gets.Rcheck" folder, if it already exists:
  toBeDeleted <- intersect("gets.Rcheck", fileNames)
  if( length(toBeDeleted)>0 ){
    for(i in toBeDeleted){
      unlink(i, recursive=TRUE) #delete folder+its content
    }
  }

} #end if(doDelete)


####################################################
## 3 BUILD AND CHECK THE TARBALL
####################################################

## Remember to check, manually, that version and date are correct in:
## - DESCRIPTION
## - NEWS
## - /R/gets-internal.R (start-up message)
## - /man/gets-package.Rd

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
#system( paste0("R CMD check ", tarballName, " --as-cran") )
system( paste0("R CMD check ", tarballName) )
## Note: The --as-cran option is obligatory according to cran policy.
## Check also the user manual for line exceedances in the .Rd files.
## This is not detected by the tarball check. An alternative way to
## check tarball, see: http://win-builder.r-project.org/


####################################################
## 4 INSTALL PACKAGE
####################################################

##remove old version first:
remove.packages("gets")

##install new version:
system( paste0("R CMD INSTALL ", tarballName) )
#system("R CMD INSTALL gets_0.25.tar.gz")
#system("R CMD INSTALL --build gets")


####################################################
## 5 UPLOAD TO CRAN
####################################################

## This must be done by the maintainer of the package.
##
## When all the testing is done, upload the tarball, that is, the
## ?.tar.gz? file, see "Writing R extensions" for the current
## submission guidelines. Upload package via:
## https://xmpalantir.wu.ac.at/cransubmit/
