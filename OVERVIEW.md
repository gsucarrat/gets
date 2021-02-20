# Overview of the files and folders on Github
The folder *gets* contains the 'raw-material' of the current development version of the package. The file *gets_devel.tar.gz* contains a tested tarball of the package. See further below for instructions on how to install this version.

## Source code
The source code of the package is contained in the following files in the *gets/R/* folder:

    gets-base-source.R
    gets-isat-source.R
    gets-dlogitx-source.R
    
The first file contains the base functions of the package, and the functions related to the *arx()* function. The second file contains the functions related to the *isat()* function. The third file contains the functions related to the *dlogitx()* function. The associated help-manual files are contained in the folder named *gets/man/*.

## Create and check the tarball
The file in the *0-make-and-check-tarball-by-G* folder is used by the maintainer to build and CRAN-check the tarball. The files in *0-test-files-by-G* are used by the maintainer to test the code before a new version is released. A subset of these tests is automatically invoked by the test-files in the *gets/tests/* folder.

## Other files
The remaining files on Github are auxiliary files, or prototypes of new ideas and functionality.
