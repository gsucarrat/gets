# gets
The R Package *gets* provides General-to-Specific (GETS) modelling of the mean and variance of a regression, and Indicator Saturation (ISAT) methods for detecting and testing for structural breaks in the mean. Facilities for *user-specified* GETS and ISAT methods are also provided.

The CRAN webpage of the package is: [https://CRAN.R-project.org/package=gets]( https://CRAN.R-project.org/package=gets). Here, on Github, we provide the development version of the package. For an overview of the files and folders, see OVERVIEW.md.

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
* User-specified GETS and ISAT (PDF): [https://mpra.ub.uni-muenchen.de/96653/](https://mpra.ub.uni-muenchen.de/96653/)
* Webpage: [http://www.sucarrat.net/R/gets](http://www.sucarrat.net/R/gets)
* CRAN webpage: [https://CRAN.R-project.org/package=gets]( https://CRAN.R-project.org/package=gets)
* [https://felixpretis.climateeconometrics.org/software/](https://felixpretis.climateeconometrics.org/software/)

# License
This package is free and open source software, licensed under GPL (>= 2)
