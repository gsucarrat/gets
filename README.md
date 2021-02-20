# gets

![R-CMD-check](https://github.com/gsucarrat/gets/workflows/R-CMD-check/badge.svg)

The R Package *gets* provides General-to-Specific (GETS) modelling of the mean and variance of a regression, and Indicator Saturation (ISAT) methods for detecting and testing for structural breaks in the mean. Facilities for *user-specified* GETS and ISAT methods are also provided.

The CRAN webpage of the package is: [https://CRAN.R-project.org/package=gets]( https://CRAN.R-project.org/package=gets). Here, on Github, we provide the development version of the package. For an overview of the files and folders, see below.

# Installation
The following R command installs the stable version from CRAN:

    install.packages('gets', dependencies = TRUE)

To install the current development version available here at Github, first download the tarball (i.e. the file named gets_devel.tar.gz). Next, run:

    system("R CMD INSTALL --build gets")

Alternatively, you can try installing the development version directly from GitHub with the following code:

    install.packages(
      "https://github.com/gsucarrat/gets/raw/master/gets_devel.tar.gz",
      repos = NULL, type = "source"
    )
    
# Resources
* An introduction (PDF): [https://doi.org/10.18637/jss.v086.i03](https://doi.org/10.18637/jss.v086.i03)
* User-specified GETS and ISAT (PDF): [https://journal.r-project.org/archive/2021/RJ-2021-024/](https://journal.r-project.org/archive/2021/RJ-2021-024/)
* Webpage: [http://www.sucarrat.net/R/gets](http://www.sucarrat.net/R/gets)
* CRAN webpage: [https://CRAN.R-project.org/package=gets]( https://CRAN.R-project.org/package=gets)
* [https://felixpretis.climateeconometrics.org/software/](https://felixpretis.climateeconometrics.org/software/)

# License
This package is free and open source software, licensed under GPL (>= 2)

# Overview of files and folders
The folder *gets* contains the 'raw-material' of the current development version of the package. The file *gets_devel.tar.gz* contains a tested tarball of the package. See further below for instructions on how to install this version.

### Source code
The source code of the package is contained in the following files in the *gets/R/* folder:

    gets-base-source.R
    gets-isat-source.R
    gets-dlogitx-source.R
    
The first file contains the base functions of the package, and the functions related to the *arx()* function. The second file contains the functions related to the *isat()* function. The third file contains the functions related to the *dlogitx()* function. The associated help-manual files are contained in the folder named *gets/man/*.

### Create and check the tarball
The file in *0-make-and-check-tarball-by-G* is used by the maintainer to build and CRAN-check the tarball. 

### Test-files
The files in *0-test-files-by-G* are used by the maintainer to test the code before a new version is released. A subset of these tests is automatically invoked by the test-files in the *gets/tests/* folder.

### Other files
The remaining files on Github are auxiliary files, or prototypes of new ideas and functionality.
