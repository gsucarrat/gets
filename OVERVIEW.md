# The R package 'gets': Overview of files and folders on Github
The folder *gets* contains the current development version of the package. Whenever its content changes, the package changes.

The source code of the package is contained in two files in the *gets/R/* folder:

    gets-base-source.R
    gets-isat-source.R

The latter contains the functions related to the *isat()* function, whereas the former contains the rest of the functions. The associated help-manual files are contained in the folder named *man*.

To create or build the tarball, the file *make-and-check-tarball.R* can be used. Its result is the tarball named *gets_XXXX.tar.gz*, which is the file that is uploaded to CRAN once a new version is ready. The folder *0-past-versions* contains (recent) past versions, and *0-test-files* contains most of the files used to test and check the code before a new version is released.

The remaining files on Github are auxiliary files, or prototypes of new ideas and functionality.
