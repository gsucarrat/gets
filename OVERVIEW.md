# Overview of the files and folders here on Github
The folder *gets* contains the 'raw-material' of the current development version of the package. Whenever the contents changes, the package changes.

## Source code
The source code of the package is contained in two files in the *gets/R/* folder:

    gets-base-source.R
    gets-isat-source.R

The latter contains the functions related to the *isat()* function, whereas the former contains the rest of the functions. The associated help-manual files are contained in the folder named *man*.

## Create and check the tarball
To build and CRAN-check the tarball, the file *make-and-check-tarball.R* can be used. Its result is a tarball named *gets_XXXX.tar.gz*, which is the file that is uploaded to CRAN once a new version is ready. The folder *0-past-versions* contains (recent) past versions, and *0-test-files* contains most of the files used to test and check the code before a new version is released.

## Other files
The remaining files on Github are auxiliary files, or prototypes of new ideas and functionality.
