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
setwd("C:/Users/sucarrat/Documents/R/gs/gets/devel/")
# setwd("~/GitHub/gets/")
#setwd(choose.dir()) #interactively

##where is the 'gets' package located?:
##=====================================

whereFolder <- "C:/Users/sucarrat/Documents/R/gs/gets/devel/gets"

# whereFolder <- "~/GitHub/gets/gets"

####################################################
## 2 CLEAN WORK-DIRECTORY AND WORKSPACE
####################################################

##delete files and folders from previous builds?:
## - '*.Rcheck' folder(s)
## - '*.tar.gz' file(s)

doDelete <- TRUE #TRUE or FALSE?

##clean work-directory:
if( doDelete ){

  ##files and folders of the work-directory:
  fileNames <- dir()

  ##the deleted files and folders:
  deletedItems <- NULL
  
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
## 2.5 MOVE CODE TO TEMPORARY FOLDER - M-orca Feb 2021
####################################################
# Find the folder above our source code (i.e. the folder above the package structure)
folder_above <- dirname(whereFolder)

# Create a temporary folder there
if(dir.exists(paste0(folder_above,"/temp"))){
  unlink(dir.exists(paste0(folder_above,"/temp")))
} else{
  dir.create(paste0(folder_above,"/temp"))
}

# Copy the package to the temporary folder
file.copy(from = whereFolder, to = paste0(folder_above,"/temp"), recursive = TRUE,overwrite = TRUE)

### Remove Unwanted Files ###
unlink(paste0(folder_above,"/temp/gets/tests"),recursive = TRUE)
file.remove(c(paste0(folder_above,"/temp/gets/",c(".gitignore",".Rbuildignore","gets.Rproj"))))
  

# Remove Unwanted Text from the DESCRIPTION
fileName <- paste0(folder_above,"/temp/gets/DESCRIPTION")
text <- readChar(fileName, file.info(fileName)$size)
# Remove testthat package from DESCRIPTION
text <- gsub("lgarch, xtable, Matrix, microbenchmark, testthat","lgarch, xtable, Matrix, microbenchmark", text)
# Remove Edition and line breaks from the DESCRIPTIONS
text <- gsub("\nConfig/testthat/edition: 3\r\n|\r","", text)

# Write the changes to the DESCRIPTION
fileConn<-file(fileName)
writeLines(text, fileConn)
close(fileConn)


####################################################
## 3 BUILD AND CHECK THE TARBALL
####################################################

## Remember to check, manually, that version and date are correct in:
## - DESCRIPTION
## - NEWS
## - /R/gets-internal.R (start-up message)
## - /man/gets-package.Rd

##build tarball:
##==============

##note: the following command assumes the gets package,
##i.e. the folder 'gets' with the source, is contained
##in whereFolder.

## M-orca, Feb 2021:
system( paste0("R CMD build ", paste0(folder_above,"/temp/gets"), " --resave-data") )

#system( paste0("R CMD build ", whereFolder, " --resave-data") )
#system( paste0("R CMD build ", whereFolder, " --compact-vignettes") )

## - The --resave-data option is recommended by CRAN for
##   better compression, but it is not obligatory.
## - In principle, the latest development version of R should be
##   used for the build, but this sometimes leads to spurious errors.

## M-orca, Feb 2021:
### Remove temporary directory ### 
unlink(paste0(folder_above,"/temp"),recursive = TRUE, force = TRUE)



##compress the vignette:
##======================
#
#system(
# "/qpdf-10.1.0/bin/qpdf.exe introduction.pdf introduction_NEW.pdf"
#)
#
#compactPDF("C:/Users/sucarrat/Documents/R/gs/gets/devel")

##check tarball (needs internet):
##===============================

fileNames <- dir()
tarballWhere <- grep(".tar.gz", fileNames)
tarballName <- fileNames[ tarballWhere ]
#system( paste0("R CMD check ", tarballName, " --as-cran") )
system( paste0("R CMD check ", tarballName) )
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
