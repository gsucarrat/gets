####################################################
## This file makes and checks the tarball of the 
## gets package, and installs it
##
## First created 2 July 2019 by Genaro Sucarrat
##
## Requirements:
##
## - writing permission to the working directory that
##   is used.
##
## Contents:
##
## 1 SET DIRECTORIES
## 2 CLEAN WORK-DIRECTORY AND WORKSPACE
## 3 BUILD AND CHECK THE TARBALL
## 4 INSTALL PACKAGE
## 5 UPLOAD TO CRAN
##
####################################################


####################################################
## 1 SET DIRECTORIES
####################################################
  
##set working directory:
##======================

##this directory will be used to store and test the
##tarball, i.e. the auxiliary files and folders
##that are created during the tests will appear in
##this folder.

##set working directory:
setwd("C:/Users/sucarrat/Documents/R/gs/gets/github/")
#setwd(choose.dir()) #interactively

##where will the 'gets' folder be located?:
##=========================================

whereFolder <- paste0(getwd(), "/gets")


####################################################
## 2 CLEAN WORK-DIRECTORY AND WORKSPACE
####################################################

##delete files and folders from previous builds?:
## - 'gets' folder
## - '*.Rcheck' folder(s)
## - '*.tar.gz' file(s)

doDelete <- TRUE #TRUE or FALSE?

##clean work-directory:
if( doDelete ){

  ##files and folders of the work-directory:
  fileNames <- dir()

  ##the deleted files and folders:
  deletedItems <- NULL

  ##delete "gets" folder, if it already exists:
  toBeDeleted <- intersect("gets", fileNames)
  if( length(toBeDeleted)>0 ){
    for(i in toBeDeleted){
      unlink(i, recursive=TRUE) #delete folder+its content
    }
  }
  deletedItems <- c(deletedItems, toBeDeleted)
  
  ##delete tarball, if it already exists:             
  toBeDeleted <- fileNames[ grep(".tar.gz", fileNames) ] 
  if( length(toBeDeleted)>0 ){ file.remove(toBeDeleted) }
  deletedItems <- c(deletedItems, toBeDeleted)
  
  ##delete "gets.Rcheck" folder, if it already exists:
  toBeDeleted <- intersect("gets.Rcheck", fileNames)
  if( length(toBeDeleted)>0 ){
    for(i in toBeDeleted){
      unlink(i, recursive=TRUE) #delete folder+its content
    }
  }
  deletedItems <- c(deletedItems, toBeDeleted)

  ##print deleted items:
  if( length(deletedItems) > 0){
    cat("\n")
    cat("The following was deleted: \n") 
    cat(deletedItems, sep="\n")
    cat("\n")
  }else{
    cat("\n")
    cat("No items were deleted\n")
    cat("\n")
  }
  
} #end if(doDelete)


####################################################
## 3 OBTAIN AND MODIFY 'gets' FOLDER
####################################################

##copy the 'gets' folder from Github:
fromFolder <- paste0( getwd(), "/contents/gets" )
file.copy(fromFolder, to=getwd(), recursive = TRUE, overwrite = TRUE)

##delete 'tests' folder:
unlink("./gets/tests", recursive = TRUE)

##delete some files:
toBeDeleted <- paste0("./gets/", c(".gitignore",".Rbuildignore","gets.Rproj"))
file.remove(toBeDeleted)

##modify DESCRIPTION file:
fileName <- paste0(getwd(), "/gets/DESCRIPTION")
fileTxt <- readChar(fileName, file.info(fileName)$size)
fileTxt <- gsub("lgarch, xtable, Matrix, microbenchmark, testthat",
  "lgarch, xtable, Matrix, microbenchmark", fileTxt)
fileTxt <- gsub("\nConfig/testthat/edition: 3\r\n|\r", "", fileTxt)
writeChar(fileTxt, fileName)


####################################################
## 4 BUILD AND CHECK THE TARBALL
####################################################

## Remember to check, manually, that version and date are correct in:
## - DESCRIPTION
## - NEWS
## - /R/gets-internal.R (start-up message)
## - /man/gets-package.Rd

##compress the vignette:
##======================

##check whether qpdf works:
##system( "qpdf.exe" ) #should work!
##system( paste0( getwd(), "/qpdf-10.1.0/bin/qpdf.exe") )

##non-automatic approach no. 1:
#library(tools)
#compactPDF( paste0(getwd(), "/introduction.pdf"), gs_quality = "ebook")
#compactPDF( paste0(getwd(), "/introduction.pdf") )
#compactPDF( paste0(getwd(), "/introduction.pdf"),
#  qpdf=paste0( getwd(), "/qpdf-10.1.0/bin/qpdf.exe" ))

##non-automatic approach no. 2:
## - compress introduction.pdf
## - put it, together with introduction.Rnw, in the gets/inst/doc folder
## - make sure the --no-build-vignettes is used when invoking R CMD build

##build tarball:
##==============

##note: the following command assumes the gets package,
##i.e. the folder 'gets' with the source, is contained
##in whereFolder.

#system( paste0("R CMD build ", whereFolder, " --resave-data") )
#system( paste0("R CMD build ", whereFolder, " --no-build-vignettes") )
system( paste0("R CMD build ", whereFolder, " --no-build-vignettes --resave-data") )
#system( paste0("R CMD build ", whereFolder, " --compact-vignettes") )

## - The --resave-data option is recommended by CRAN for
##   better compression, but it is not obligatory.
## - In principle, the latest development version of R should be
##   used for the build, but this sometimes leads to spurious errors.
## - I have not been able to make the "--compact-vignettes" option
##   work yet, since it seems I have not been able to properly enable
##   qpdf for R

##check tarball (needs internet):
##===============================

fileNames <- dir()
tarballWhere <- grep(".tar.gz", fileNames)
tarballName <- fileNames[ tarballWhere ]
system( paste0("R CMD check ", tarballName, " --as-cran") )
#system( paste0("R CMD check ", tarballName) )
## Note: The --as-cran option is obligatory according to cran policy.
## Check also the PDF user manual for line exceedances in the .Rd files.
## This is not detected by the tarball check. An alternative way to
## check tarball, see: http://win-builder.r-project.org/


####################################################
## 4 INSTALL PACKAGE
####################################################

##remove old version first:
remove.packages("gets")

##install new version:
system( paste0("R CMD INSTALL ", tarballName) )
#system("R CMD INSTALL gets_0.27.tar.gz")
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
