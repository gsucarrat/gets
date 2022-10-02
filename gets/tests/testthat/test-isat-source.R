##################################################
## Test file for the 'isat' source of the 'gets'
## package. First created 23 September 2014, Oslo.
##
##
## 1 INITIATE
## 2 TEST iim(), sim() AND tim()
## 3 TEST MAIN isat() ARGUMENTS
## 4 TEST USER-DEFINED DIAGNOSTICS
## 5 TEST USER-DEFINED ESTIMATION
## 6 TEST USER-DEFINED GOF FUNCTION
## 7 TEST PARALLEL COMPUTING
## 8 TEST ROBUST COEFFICIENT COVARIANCES
## 9 SIMULATIONS
##
##################################################

##################################################
##1 INITIATE
##################################################



test_environment <- FALSE


##################################################
## 2 TEST iim(), sim() AND tim()
##################################################

# x <- 6
# iim(x); sim(x); tim(x)
# which.ones <- c(2,5)
# iim(x, which.ones=which.ones)
# sim(x, which.ones=which.ones)
# tim(x, which.ones=which.ones)
#
# set.seed(123)
# x <- ts(rnorm(5), start=2010)
# iim(x); sim(x); tim(x)
#
# set.seed(123)
# x <- ts(rnorm(6), frequency=4, end=c(2015,4))
# iim(x); sim(x); tim(x)
# which.ones <- c(2014.25,2014.75,2015,2015.50)
# iim(x, which.ones=which.ones)
# sim(x, which.ones=which.ones)
# tim(x, which.ones=which.ones)
#
# ##the following commands should return the
# ##error-message: 'which.ones' not in index
# which.ones <- c(2001.25,2002.50,2003.75,2004)
# iim(x, which.ones=which.ones)
# sim(x, which.ones=which.ones)
# tim(x, which.ones=which.ones)
#
# ##used to yield error:
# set.seed(123); y <- rnorm(30); x <- as.vector(sim(y, which.ones=15))
# isat(y, mxreg=x, LAPACK=FALSE)

##some ideas for the future?:
#y <- ymon <- yqtr <- rnorm(24)
###zooreg:
#y <- zooreg(y, frequency=4, start=c(2000,1))
#index(y)
#
###quarterly:
#yqtr <- yearqtr(yqtr, start=c(2000,1))
#
#index(y)
#tmp <- as.yearqtr(index(y))
#format(tmp, format="%YQ%q")
#iim(y)


##################################################
## 3 TEST MAIN isat() ARGUMENTS
##################################################

##control the plotting:
##=====================

getOption("plot")
options(plot=FALSE)
options(print.searchinfo=FALSE)
#options(plot=FALSE)
#options(plot=NULL)

test_that("TEST MAIN isat() ARGUMENTS - very basic", {

  ##control the plotting:
  ##=====================

  getOption("plot")
  options(plot = FALSE)
  #options(plot=FALSE)
  #options(plot=NULL)


  ##generate some data:
  ##===================

  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n) #or: y <- rt(dgp.n, df=4.1)
  mX <- matrix(rnorm(dgp.n * 3), dgp.n, 3)
  y[1:10] <- y[1:10] + 4 #step-shift
  #plot(as.zoo(y), col="blue") #plot


  ##some basic tests:
  ##=================

  ##for visual inspection:
  expect_silent(isat(y, sis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, iis=TRUE, sis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, iis=TRUE, sis=FALSE, print.searchinfo = FALSE))
  expect_silent(isat(y, sis=FALSE, tis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, sis=TRUE, tis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, iis=TRUE, sis=FALSE, tis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, iis=TRUE, sis=TRUE, tis=TRUE, print.searchinfo = FALSE))

  expect_error(isat(y, sis=FALSE))

  expect_silent(isat(y, ar=0, sis=TRUE, print.searchinfo = FALSE)) # checking that ar=0 works
  expect_silent(isat(y, ar=1:2, sis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, iis=TRUE, sis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, iis=TRUE, sis=FALSE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, iis=FALSE, sis=FALSE, tis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, iis=TRUE, sis=FALSE, tis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, iis=FALSE, sis=TRUE, tis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, iis=TRUE, sis=TRUE, tis=TRUE, print.searchinfo = FALSE))

  expect_silent(isat(y, ar=1:2, mxreg=mX, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, mxreg=mX, iis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, mxreg=mX, iis=FALSE, sis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, mxreg=mX, iis=TRUE, sis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, mxreg=mX, iis=FALSE, tis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, mxreg=mX, iis=TRUE, tis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, mxreg=mX, iis=FALSE, sis=TRUE, tis=TRUE, print.searchinfo = FALSE))
  expect_silent(isat(y, ar=1:2, mxreg=mX, iis=TRUE, sis=TRUE, tis=TRUE, print.searchinfo = FALSE))

  # Check the messages
  expect_message(isat(y, ar=1:2, mxreg=mX, iis=TRUE, sis=TRUE, tis=TRUE, print.searchinfo = TRUE))


  ##yielded error in version 0.9 to 0.23:
  expect_silent(isat(y, ar=1:2, mxreg=as.data.frame(mX), print.searchinfo = FALSE))


})

test_that("TEST MAIN isat() ARGUMENTS - slightly more advanced",{
  ##issue reported by F-bear regarding version 0.14 (email 30/3-2018).
  ##in version 0.14 the impulse dummies were not detected:
  set.seed(123)
  y <- rnorm(100, 0, 1)
  y[30] <- y[30]+10
  y[40] <- y[40]+10
  y[60] <- y[60]-10
  # plot(as.zoo(y))

  expect_silent(a <- isat(y, iis=TRUE, sis=FALSE, t.pval=0.05, plot=FALSE, print.searchinfo = FALSE))

  # check that any iis are identified
  expect_true(grepl("iis",a$ISfinalmodels))

  # check that the correct iis are identified
  expect_true(all(c("iis30", "iis40", "iis60") %in% a$ISfinalmodels[[1]]))

  z <- rnorm(100, 0, 1)
  z[1:30] <- z[1:30] + 10
  b <- isat(z, iis=TRUE,  sis=FALSE, t.pval=0.05, plot=FALSE, print.searchinfo = FALSE)

  # check that any iis are identified
  expect_true(grepl("iis",b$ISfinalmodels))

  # check that no sis are identified
  expect_true(!grepl("sis",b$ISfinalmodels))




  ##issue reported by F-bear regarding version 0.20 (email 26/9-2019).
  ## "...there seems to be a serious bug in the latest version of the package.
  ## ISnames seems to be null, even if there are impulses retained. This
  ## seems to break a lot of other functions building on the ISnames element.
  ## For example, this does not work and returns null: [fixed by G in 0.21]"
  set.seed(123)
  y <- rnorm(100, 0, 1)
  my <- isat(y, iis=TRUE, sis=FALSE, t.pval=0.05,plot=FALSE, print.searchinfo = FALSE)
  my

  expect_true(!is.null(my$ISnames))


})

test_that("TEST MAIN isat() ARGUMENTS - Testing the extraction functions", {


  ##test the extraction functions:
  ##==============================

  ##same data as earlier:
  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n) #or: y <- rt(dgp.n, df=4.1)
  mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
  y[1:10] <- y[1:10] + 4 #step-shift
  # plot(as.zoo(y), col="blue") #plot

  isatmod <- isat(y, ar=1:2, mxreg=mX, iis=TRUE, sis=TRUE, tis=TRUE, print.searchinfo = FALSE)

  expect_output(print(isatmod))
  expect_type(summary(isatmod), "character")
  expect_equal(length(summary(isatmod)),90)


  coef(isatmod)
  plot(cbind(fitted(isatmod),
             fitted(isatmod, spec="m"),
             fitted(isatmod, spec="v")))
  logLik(isatmod)
  plot(cbind(residuals(isatmod),
             residuals(isatmod, std=FALSE),
             residuals(isatmod, std=TRUE)))
  paths(isatmod)
  expect_error(paths(mod01)) #should return the error-message: object 'mod01' not found
  plot(isatmod)
  expect_error(predict(isatmod),"'newmxreg' is NULL") #should return the error-message: 'newmxreg' is NULL
  predict(isatmod, newmxreg=matrix(0,12,5))
  predict(isatmod, n.ahead=1, newmxreg=matrix(0,1,5)) #used to yield error
  predict(isatmod, newmxreg=matrix(0,12,5), newindex=13:24)
  predict(isatmod, newmxreg=matrix(0,12,5), return=FALSE)
  predict(isatmod, newmxreg=matrix(0,12,5), plot=FALSE)
  predict(isatmod, newmxreg=matrix(0,12,5), return=FALSE, plot=FALSE)
  terminals(isatmod)
  expect_error(terminals(mod01), "object 'mod01' not found") #should return the error-message: object 'mod01' not found
  vcov(isatmod)

})



test_that("TEST MAIN isat() ARGUMENTS - test further arguments",{

  options(print.searchinfo = FALSE)

  set.seed(123)
  y <- rnorm(30)

  ##test further arguments:
  ##=======================

  expect_silent(isat(y, print.searchinfo = FALSE, t.pval = 0.1)) #default: 0.001
  expect_silent(isat(y, print.searchinfo = FALSE, do.pet = TRUE)) #default: FALSE
  expect_silent(isat(y, print.searchinfo = FALSE, wald.pval = 0.1)) #default: 0.001
  expect_silent(isat(y, print.searchinfo = FALSE, ar.LjungB = list(lag = NULL, pval = 0.01))) #default: NULL

  #no search, because the gum does not pass arch-diagnostics:
  expect_silent(isat(y, print.searchinfo = FALSE, arch.LjungB = list(lag = NULL, pval = 0.01))) #default: NULL
  expect_silent(isat(y, print.searchinfo = FALSE, normality.JarqueB = 0.025)) #default: NULL
  expect_silent(isat(y, print.searchinfo = FALSE, info.method = "aic")) #default: sc

  expect_silent(isat(y, print.searchinfo = FALSE, arch.LjungB = list(lag = NULL, pval = 0.9))) #default: NULL

  #should return warning:
  expect_warning(isat(y, print.searchinfo = FALSE, include.gum = TRUE)) #default: NULL

  expect_silent(isat(y, print.searchinfo = FALSE, include.1cut = TRUE)) #default: FALSE
  expect_silent(isat(y, print.searchinfo = FALSE, include.empty = TRUE)) #default: FALSE
  expect_silent(isat(y, print.searchinfo = FALSE, max.paths = 3)) #default: NULL (i.e. "multi-path")
  expect_silent(isat(y, print.searchinfo = FALSE, max.paths = 5)) #default: NULL (i.e. "multi-path")
  expect_silent(isat(y, print.searchinfo = FALSE, parallel.options = NULL)) #default: NULL
  expect_silent(isat(y, print.searchinfo = FALSE, parallel.options = 2)) #default: NULL
  if(parallel::detectCores()<5){expect_error(isat(y, print.searchinfo = FALSE, parallel.options = 5))} #default: NULL
  #isat(y, print.searchinfo = FALSE, parallel.options = 4) #default: NULL
  expect_silent(isat(y, print.searchinfo = FALSE, parallel.options = 2, max.paths = 2))
  expect_silent(isat(y, print.searchinfo = FALSE, turbo = TRUE))
  expect_silent(isat(y, print.searchinfo = FALSE, parallel.options = 2, max.paths = 2, turbo = TRUE)) #default: NULL

})



options(plot = TRUE)


# You'd then also provide a helper that skips tests where you can't
# be sure of producing exactly the same output
expect_snapshot_plot <- function(name, code) {
  # Other packages might affect results
  skip_if_not_installed("ggplot2", "2.0.0")
  # Or maybe the output is different on some operation systems
  skip_on_ci()
  # You'll need to carefully think about and experiment with these skips

  name <- paste0(name, ".png")

  # Announce the file before touching `code`. This way, if `code`
  # unexpectedly fails or skips, testthat will not auto-delete the
  # corresponding snapshot file.
  announce_snapshot_file(name = name)

  # To use expect_snapshot_file() you'll typically need to start by writing
  # a helper function that creates a file from your code, returning a path
  save_png <- function(code, width = 400, height = 400) {
    path <- tempfile(fileext = ".png")
    png(path, width = width, height = height)
    on.exit(dev.off())
    code

    path
  }

  path <- save_png(code)
  expect_snapshot_file(path, name)
}

test_that("TEST MAIN isat() ARGUMENTS - additional test of predict.isat",{

  ##same data as earlier:
  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n) #or: y <- rt(dgp.n, df=4.1)
  mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
  y[1:10] <- y[1:10] + 4 #step-shift
  # plot(as.zoo(y), col="blue") #plot

  isatmod <- isat(y, ar = 1:2, mxreg = mX, iis = TRUE, sis = TRUE, tis = TRUE, print.searchinfo = FALSE)


  ##additional test of predict.isat():
  ##==================================

  expect_snapshot_plot("Test 1",predict(isatmod, newmxreg = matrix(0,12,5),  ci.levels=seq(0.20,0.95,by=0.05),n.sim=20000))
  expect_snapshot_plot("Test 2",predict(isatmod, newmxreg = matrix(0,12,5),  ci.levels=seq(0.20,0.95,by=0.05),n.sim=20000, plot.options=list(shades=seq(20,95,by=5))))
  expect_snapshot_plot("Test 3",predict(isatmod, newmxreg = matrix(0,12,5),  ci.levels=seq(0.20,0.95,by=0.05),n.sim=20000, plot.options=list(shades=seq(95,20,by=-5))))
  expect_snapshot_plot("Test 4",predict(isatmod, newmxreg = matrix(0,12,5),  ci.levels=seq(0.20,0.95,by=0.05),n.sim=20000, plot.options=list(shades=seq(100,25,by=-5))))
  expect_snapshot_plot("Test 5",predict(isatmod, newmxreg = matrix(0,12,5),  plot.options = list(keep=1)))
  expect_snapshot_plot("Test 6",predict(isatmod, newmxreg = matrix(0,12,5),  plot.options = list(line.at.origin=FALSE)))
  expect_snapshot_plot("Test 7",predict(isatmod, newmxreg = matrix(0,12,5),  plot.options = list(start.at.origin=TRUE)))
  expect_snapshot_plot("Test 8",predict(isatmod, newmxreg = matrix(0,12,5),  plot.options = list(start.at.origin=TRUE, fitted=FALSE)))
  expect_snapshot_plot("Test 9",predict(isatmod, newmxreg = matrix(0,12,5),  plot.options = list(dot.at.origin=FALSE)))
  expect_snapshot_plot("Test 10",predict(isatmod, newmxreg = matrix(0,12,5), plot.options = list(hlines=c(-2,0,2,4,6,8))))
  expect_snapshot_plot("Test 11",predict(isatmod, newmxreg = matrix(0,12,5), plot.options = list(col=c("darkred","green"))))
  expect_snapshot_plot("Test 12",predict(isatmod, newmxreg = matrix(0,12,5), plot.options = list(lty=c(3,2))))
  expect_snapshot_plot("Test 13",predict(isatmod, newmxreg = matrix(0,12,5), plot.options = list(lwd=3)))
  expect_snapshot_plot("Test 14",predict(isatmod, newmxreg = matrix(0,12,5), plot.options = list(ylim=c(-8,16))))
  expect_snapshot_plot("Test 15",predict(isatmod, newmxreg = matrix(0,12,5), plot.options = list(ylab="G-values")))
  expect_snapshot_plot("Test 16",predict(isatmod, newmxreg = matrix(0,12,5), plot.options = list(main="Plot slightly lower when 'main' is specified")))
  expect_snapshot_plot("Test 17",predict(isatmod, newmxreg = matrix(0,12,5), plot.options = list(legend.text=c("Prognose","Faktisk"))))
  expect_snapshot_plot("Test 18",predict(isatmod, newmxreg = matrix(0,12,5), plot.options = list(fitted=FALSE)))
  expect_snapshot_plot("Test 19",predict(isatmod, newmxreg=matrix(0,12,5),   plot.options = list(newmactual=rep(0,6))))
  expect_snapshot_plot("Test 20",predict(isatmod, newmxreg=matrix(0,12,5),   plot.options = list(shades=c(95,50))))
  expect_snapshot_plot("Test 21",predict(isatmod, newmxreg = matrix(0, 12, 5), plot.options = list(shades = c(50, 95)))) #invert shades

})

test_that("TEST MAIN isat() ARGUMENTS - Test that mconst is named correctly",{

  ##In the following model (isatmod), the constant was not correctly
  ##named 'mconst' at one point. Instead, it was named 'mxreg',
  ##which created problems for predict.isat: predict(isatmod). Issue
  ##solved 31/7/2019.

  set.seed(123)
  y <- rnorm(30)
  isatmod <- isat(y, print.searchinfo = FALSE)
  expect_identical(names(coef(isatmod))[1], "mconst")
  expect_snapshot_plot("Test 22",predict(isatmod, plot=TRUE))

  ##same y, but slightly different model:
  isatmod <- isat(y, ar=1, print.searchinfo = FALSE)
  expect_identical(names(coef(isatmod))[1], "mconst")
  expect_snapshot_plot("Test 23",predict(isatmod, plot=TRUE))


})






test_that("TEST MAIN isat() ARGUMENTS - test uis argument",{

  ##test uis argument:
  ##==================
  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n)
  ##uis as matrix:
  uis <- iim(dgp.n)
  expect_silent(isat(y, sis=FALSE, uis=uis, print.searchinfo = FALSE))
  expect_silent(isat(y, sis=FALSE, uis=uis, max.paths=1, print.searchinfo = FALSE))
  uis <- sim(dgp.n)

  ##as of August 2020, these did not work;
  ##but as of 2 October 2020 they do!:
  expect_silent(isat(y, print.searchinfo = FALSE, sis=FALSE, uis=uis))
  expect_silent(isat(y, print.searchinfo = FALSE, sis=FALSE, uis=uis[,seq(1,dgp.n,2)]))
  expect_silent(isat(y, print.searchinfo = FALSE, sis=FALSE, uis=uis, max.paths=1))

  uis <- tim(dgp.n)
  expect_silent(isat(y, print.searchinfo = FALSE, sis=FALSE, uis=uis))
  expect_silent(isat(y, print.searchinfo = FALSE, sis=FALSE, uis=uis[,seq(1,dgp.n,2)]))
  expect_silent(isat(y, print.searchinfo = FALSE, sis=FALSE, uis=uis, max.paths=1))

  ##used to crash (uis is a matrix):
  set.seed(123)
  y <- rnorm(30)
  z <- rnorm(30)
  mX <- matrix(rnorm(1*30), 30, 1)
  expect_silent(isat(y, mxreg = z, iis = FALSE, sis = TRUE, uis = mX, print.searchinfo = FALSE))

  ##used to crash (uis is a matrix):
  ##as of August 2020, these did not work;
  ##but as of 2 October 2020 they do!:
  dgpN <- 50
  set.seed(123)
  y <- rnorm(dgpN)
  x <- rnorm(dgpN)
  x_mis <- sim(dgpN)*x
  colnames(x_mis) <- paste("mis", seq(2:(NCOL(x_mis)+1)), sep="")

  # The following will still produce messages (despite print.searchinfo = FALSE): "Warning: uis specified but no mxfull variable given. Using mconst instead."
  # this officially is not a warning but just a message
  expect_message(isat(y, print.searchinfo = FALSE, ar=1, mxreg=x, sis=FALSE, uis = x_mis, t.pval=0.05))
  expect_message(isat(y, print.searchinfo = FALSE, ar=1, mxreg=x, sis=FALSE, uis = x_mis, t.pval=0.05, max.paths=1))

  expect_snapshot_plot("Test 24",isat(y, ar = 1, mxreg = x, sis = FALSE, uis = x_mis, t.pval = 0.05, plot = TRUE, print.searchinfo = FALSE))
  expect_snapshot_plot("Test 25",isat(y, ar = 1, mxreg = x, sis = FALSE, uis = x_mis, t.pval = 0.05, max.paths = 1, plot =
                                        TRUE, print.searchinfo = FALSE))

  ##used to yield error because NCOL(mX) > length(y):
  set.seed(123)
  y <- rnorm(20)
  mX <- matrix(rnorm(20*40), 20, 40)
  expect_silent(isat(y, sis=FALSE, uis = mX, print.searchinfo = FALSE))

  ##uis as list:
  uis <- list(sis = sim(y) , tis = tim(y))
  ##as of August 2020, these do not work (why?):
  ##but as of 14 October 2020 they do!:
  expect_silent(isat(y, print.searchinfo = FALSE, sis = FALSE, uis = uis))
  expect_silent(isat(y, print.searchinfo = FALSE, sis = FALSE, uis = uis, max.paths = 1))
  expect_silent(isat(y, print.searchinfo = FALSE, sis = FALSE, uis = uis, max.paths = 2))

  ##uis as data.frame:
  ##in an email 5/10-2020, F-bear reported an error produced by
  ##the following code:
  ## UPDATE 08/2022: potentially fixed by M-orca
  # before, when passing a data.frame for uis, each column was split up into a separate list
  # now a data.frame is treated like a matrix
  set.seed(123)
  mx <- data.frame(matrix(runif(500),ncol=5))
  colnames(mx) <- paste("x", seq(1:5), sep="")
  mx <- as.data.frame(mx)
  #mx <- as.zoo(mx) #conversion to zoo solves the issue partially
  yx <- 2*mx[,1] + rnorm(100, 0, 1)

  # The following will still produce messages (despite print.searchinfo = FALSE): "Warning: uis specified but no mxfull variable given. Using mconst instead."
  # this officially is not a warning but just a message
  expect_message(isat(yx, iis=TRUE, sis=FALSE, uis=mx, t.pval=0.01, print.searchinfo = FALSE))
  expect_message(isat(yx, iis=TRUE, sis=FALSE, uis=mx, t.pval=0.01, print.searchinfo = TRUE))

})

test_that("TEST MAIN isat() ARGUMENTS - test blocks argument",{

  ##test blocks argument:
  ##=====================

  ##same data as earlier:
  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n)

  myblocks <- list()
  myblocks[[1]] <- list(10:20, 30:40)
  expect_silent(isat(y, print.searchinfo = FALSE, iis=TRUE, sis=FALSE, tis=FALSE, uis=FALSE, blocks=myblocks))
  expect_silent(isat(y, print.searchinfo = FALSE, iis=TRUE, sis=FALSE, tis=FALSE, uis=FALSE, blocks=myblocks, max.paths=1))
  expect_silent(isat(y, print.searchinfo = FALSE, iis=FALSE, sis=TRUE, tis=FALSE, uis=FALSE, blocks=myblocks))
  expect_silent(isat(y, print.searchinfo = FALSE, iis=FALSE, sis=TRUE, tis=FALSE, uis=FALSE, blocks=myblocks, max.paths=1))

  uis <- list(sim(dgp.n), tim(dgp.n))
  myblocks[[2]] <- list(7:19, 27:34, 40:45)

  ##as of August 2020, these did not work;
  ##but as of 2 October 2020 they do!:
  expect_silent(isat(y, print.searchinfo = FALSE, iis=FALSE, sis=FALSE, tis=FALSE, uis=uis, blocks=myblocks))
  expect_silent(isat(y, print.searchinfo = FALSE, iis=FALSE, sis=FALSE, tis=FALSE, uis=uis, blocks=myblocks, max.paths=1))

})

test_that("TEST MAIN isat() ARGUMENTS - further tests of predict.isat",{

  ##further tests of predict.isat:
  ##==============================

  ##issue reported by Steven Sabol (email 27/01-2017). The
  ##code should work as of 0.11:
  set.seed(123)
  mxreg <- zooreg(matrix(rnorm(700), ncol = 7), start = 2002 ,frequency = 12)
  colnames(mxreg) <- c("c1","c2","c3","c4","c5","c6","c7")
  y = zooreg(rnorm(88),start = 2002 ,frequency = 12)

  expect_silent(isat_mod <- isat(y, mxreg = mxreg, mc = TRUE, ar = 4, sis = TRUE, t.pval = 0.01, vcov.type = "white", print.searchinfo = FALSE))
  newmxreg  <- tail(na.trim(mxreg),12)
  new_index <- index(tail(na.trim(mxreg),12))

  # TODO!!!!
  ##as of 17 July 2019, does not work:
  #prediction_isat <- predict.isat(isat_mod, newmxreg = newmxreg,
  #                                n.ahead = 12, newindex = new_index, return = TRUE)
})

test_that("TEST MAIN isat() ARGUMENTS - tests of biascorr, isattest, isatvar",{

  ##tests of biascorr, isattest, isatvar, ...etc.:
  ##==============================================

  ##issue reported by Gareth Thomas (EViews, email 22/5-2017):
  ##"The following code produces a blank graphics output.  It seems as
  ##though isattest calls the R graphics even if plot=FALSE is set".
  ##Should be solved as of 0.13.
  set.seed(123)
  d <- matrix(0,100,1)
  d[35:55] <- 1
  e <- rnorm(100, 0, 1)
  y <- d*2 +e

  ##Static Test against hnull=0 using bias-correction:
  ys <- isat(y, sis=TRUE, iis=FALSE, tis=FALSE, t.pval=0.01, plot=FALSE, print.searchinfo = FALSE)
  expect_silent(isattest(ys, hnull = 0, lr = FALSE, ci.pval = 0.99, plot = FALSE, biascorr = TRUE))
  expect_identical(ncol(isattest(ys, hnull = 0, lr = FALSE, ci.pval = 0.99, plot = FALSE, biascorr = TRUE)),as.integer(4))
  expect_identical(names(isattest(ys, hnull = 0, lr = FALSE, ci.pval = 0.99, plot = FALSE, biascorr = TRUE)),c("ci.low","ci.high", "bias.high","bias.low"))


})


# ##################################################
# ## 4 TEST USER-DEFINED DIAGNOSTICS
# ##################################################

test_that("TEST USER-DEFINED DIAGNOSTICS",{

  ##generate some data:
  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n)
  y[1:10] <- y[1:10] + 4
  mX <- matrix(rnorm(dgp.n*3),dgp.n,3)

  ##user-defined Shapiro-Wilks test for normality in the residuals:
  SWtest <- function(x, ...){
    tmp <- shapiro.test(x$residuals)
    result <- c(tmp$statistic, NA, tmp$p.value)
    return(result)
  }

  # TODO this is not able to be tested using automatic testing
  # for manual testing, just uncomment the below - these should execute normally
  #expect_silent(isat(y, user.diagnostics = list(name = "SWtest", pval = 1e-10), print.searchinfo = FALSE))
  #expect_identical(nrow(isat(y, user.diagnostics = list(name = "SWtest", pval = 1e-10), print.searchinfo = FALSE)$diagnostics), as.integer(3))
  #expect_identical(row.names(isat(y, user.diagnostics = list(name = "SWtest", pval = 1e-10), print.searchinfo = FALSE)$diagnostics)[3], "SWtest")

  ##test the envir entry:
  rm("SWtest") #make sure SWtest is not defined in the global environment
  myenv <- new.env()
  assign("SWtest",
         function(x, ...){
           tmp <- shapiro.test(x$residuals)
           result <- c(tmp$statistic, NA, tmp$p.value)
           result <- rbind( as.numeric(c(tmp$statistic, NA, tmp$p.value)) )
           rownames(result) <- "SWtest"
           return(result)
         },
         envir = myenv) #close assign
  expect_error(isat(y, user.diagnostics = list(name = "SWtest", pval = 0.025), print.searchinfo = FALSE)) #should not work

  expect_silent(isat(y, print.searchinfo = FALSE, user.diagnostics = list(name = "SWtest", pval = 1e-05, envir = myenv)))
  expect_silent(isat(y, print.searchinfo = FALSE, user.diagnostics = list(name = "SWtest", pval = 0.025,  envir = myenv)))

})


##################################################
## 5 TEST USER-DEFINED ESTIMATION
##################################################

test_that("TEST USER-DEFINED ESTIMATION",{
  
  skip_if(!test_environment)
  
  ##generate some data:
  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n)
  y[1:10] <- y[1:10] + 4
  mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
  
  ##define estimator:
  myEstimator <- function(y, x){ ols(y,x) }
  
  ##isat w/user-defined estimator:
  expect_silent(isat(y, user.estimator=list(name="myEstimator"), print.searchinfo = FALSE))
  
  ##isat w/user-defined estimator and test for normality:
  SWtest <- function(x, ...){
    tmp <- shapiro.test(x$residuals)
    result <- c(tmp$statistic, NA, tmp$p.value)
    return(result)
  }
  expect_silent(isat(y, user.estimator=list(name="myEstimator"),
                     user.diagnostics=list(name="SWtest", pval=1e-10), print.searchinfo = FALSE))
  
  expect_identical(row.names(isat(y, user.estimator=list(name="myEstimator"),
                                  user.diagnostics=list(name="SWtest", pval=1e-10), print.searchinfo = FALSE)$diagnostics)[3], "SWtest")
  
  expect_true(is.numeric(isat(y, user.estimator=list(name="myEstimator"),
                              user.diagnostics=list(name="SWtest", pval=1e-10), print.searchinfo = FALSE)$diagnostics[3,3]))
  
})


test_that("TEST USER-DEFINED ESTIMATION",{
  
  skip_if(!test_environment)
  
  ##generate some data:
  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n)
  y[1:10] <- y[1:10] + 4
  mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
  
  ##faster ols?:
  ##There are packages and routines that make OLS faster in
  ##certain situations, e.g. the Matrix package. The code below
  ##creates a new function, ols2, which is essentially a copy
  ##of ols(y, x, method=3), but based on routines from the Matrix
  ##package.
  library(Matrix)
  ols2 <- function(y, x){
    out <- list()
    out$n <- length(y)
    if (is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if (out$k > 0) {
      x <- as(x, "dgeMatrix")
      out$xpy <- crossprod(x, y)
      out$xtx <- crossprod(x)
      out$coefficients <- as.numeric(solve(out$xtx,out$xpy))
      out$xtxinv <- solve(out$xtx)
      out$fit <- out$fit <- as.vector(x %*% out$coefficients)
    }else{
      out$fit <- rep(0, out$n)
    }
    out$residuals <- y - out$fit
    out$residuals2 <- out$residuals^2
    out$rss <- sum(out$residuals2)
    out$sigma2 <- out$rss/out$df
    if(out$k > 0){ out$vcov <- as.matrix(out$sigma2 * out$xtxinv) }
    out$logl <-
      -out$n * log(2 * out$sigma2 * pi)/2 - out$rss/(2 * out$sigma2)
    return(out)
  }
  
  ##isat w/ols2:
  expect_silent(isat(y, user.estimator=list(name="ols2"), print.searchinfo = FALSE))
  # checking that coefficients are identical across estimators
  expect_identical(round(isat(y, user.estimator=list(name="ols2"), print.searchinfo = FALSE)$coefficients,7),
                   round(isat(y, print.searchinfo = FALSE)$coefficients,7))
  
  expect_identical(round(isat(y, ar = 1, user.estimator=list(name="ols2"), print.searchinfo = FALSE)$coefficients,7),
                   round(isat(y, ar = 1, print.searchinfo = FALSE)$coefficients,7))
  
})

# M-orca 02/10/22: commented out the speed comparison section
# all of the below should work, but this is isn't needed each time that we do automatic testing
# ##compare speed 1:
# system.time(isat(y))
# system.time(isat(y, user.estimator=list(name="ols2")))
# ##Conclusion: here, ols is faster than ols2
# 
# ##comparisons 2: w/microbenchmark, see
# ##https://nelsonareal.net/blog/2017/06/speeding_up_ols.html
# library(microbenchmark)
# microbenchmark( ols(y,mX), ols2(y,mX), times=10)
# microbenchmark( isat(y),
#                 isat(y, user.estimator=list(name="ols2")),
#                 times=10)
# ##Conclusion ("small T"): ols is faster than ols2
# 
# ##compare speed 2:
# set.seed(123); dgp.n <- 1000; y <- rnorm(dgp.n)
# system.time(isat(y))
# system.time(isat(y, user.estimator=list(name="ols2")))
# ##Conclusion: sample size matters, additional experiments
# ##suggests the speed increase is increasing in sample size


##################################################
## 6 TEST USER-DEFINED GOF FUNCTION
##################################################

test_that("TEST USER-DEFINED GOF FUNCTION", {
  
  skip_if(!test_environment)
  
  ##generate some data:
  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n)
  y[1:10] <- y[1:10] + 4
  mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
  
  
  ##define estimator:
  myEstimator <- function(y, x){ ols(y,x) }
  
  ##user-defined gof-function:
  myGof <- function(object){ infocrit(object) }
  
  ##do isat:
  expect_silent(isat(y, gof.function=list(name="myGof"), print.searchinfo = FALSE))
  expect_silent(isat(y, user.estimator=list(name="myEstimator"), gof.function=list(name="myGof"), print.searchinfo = FALSE))
  expect_silent(isat(y, user.diagnostics=list(name="SWtest", pval=1e-10), user.estimator=list(name="myEstimator"),gof.function=list(name="myGof"), print.searchinfo = FALSE))
  
  ##adjusted R-squared:
  myGof <- function(object){
    yvar <- object$fit + object$residuals
    TSS <- sum( (yvar - mean(yvar))^2 )
    RSS <- sum(object$residuals^2)
    Rsquared <- 1 - RSS/TSS
    result <- 1 - (1-Rsquared)*(object$n-1)/(object$n-object$k)
    return(result)
  }
  
  ##do isat while maximising R-squared:
  expect_silent(isat(y, gof.function=list(name="myGof"), gof.method="max", print.searchinfo = FALSE))
  
  ##minimise the number of parameters/regressors:
  myGof <- function(x, ...){ return( x$k ) }
  expect_silent(isat(y, gof.function=list(name="myGof"), gof.method="min", print.searchinfo = FALSE))
  
  
})



##################################################
## 7 TEST PARALLEL COMPUTING
##################################################

test_that("TEST PARALLEL COMPUTING",{
  
  ##generate some data:
  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n)
  y[1:10] <- y[1:10] + 4
  mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
  
  expect_silent(isat(y, print.searchinfo = FALSE))
  expect_silent(isat(y, parallel.options=2, print.searchinfo = FALSE))
  # check that parallel and non-parallel equal same result
  expect_identical(isat(y, print.searchinfo = FALSE)$coefficients,
                   isat(y, parallel.options=2, print.searchinfo = FALSE)$coefficients)
  
  
  expect_silent(isat(y, mxreg = mX, print.searchinfo = FALSE))
  expect_silent(isat(y, mxreg = mX, parallel.options=2, print.searchinfo = FALSE))
  # check that parallel and non-parallel equal same result
  expect_identical(isat(y, mxreg = mX, print.searchinfo = FALSE)$coefficients,
                   isat(y, mxreg = mX, parallel.options=2, print.searchinfo = FALSE)$coefficients)
  
  skip_if(!test_environment)
  
  ##user-defined Shapiro-Wilks test for normality in the residuals:
  SWtest <- function(x, ...){
    tmp <- shapiro.test(x$residuals)
    result <- c(tmp$statistic, NA, tmp$p.value)
    return(result)
  }
  
  expect_silent(isat(y, user.diagnostics=list(name="SWtest", pval=1e-10), print.searchinfo = FALSE))
  expect_silent(isat(y, user.diagnostics=list(name="SWtest", pval=1e-10),parallel.options=2, print.searchinfo = FALSE))
  # check that parallel and non-parallel equal same result
  expect_identical(isat(y, user.diagnostics=list(name="SWtest", pval=1e-10), print.searchinfo = FALSE)$coefficients,
                   isat(y, user.diagnostics=list(name="SWtest", pval=1e-10),parallel.options=2, print.searchinfo = FALSE)$coefficients)
  
  ##user-defined estimator:
  myEstimator <- function(y, x){ ols(y,x) }
  expect_silent(isat(y, user.estimator = list(name = "myEstimator"), print.searchinfo = FALSE))
  expect_silent(isat(y, user.estimator = list(name = "myEstimator"), parallel.options = 2, print.searchinfo = FALSE))
  expect_identical(isat(y, user.estimator = list(name = "myEstimator"), print.searchinfo = FALSE)$coefficients,
                   isat(y, user.estimator = list(name = "myEstimator"), parallel.options = 2, print.searchinfo = FALSE)$coefficients)
  ##user-defined gof:
  myGof <- function(x, ...){ return( x$k ) }
  expect_silent(isat(y, gof.function=list(name="myGof"), gof.method="min", print.searchinfo = FALSE))
  expect_silent(isat(y, gof.function=list(name="myGof"), gof.method="min", parallel.options=2, print.searchinfo = FALSE))
  expect_identical(isat(y, gof.function=list(name="myGof"), gof.method="min", print.searchinfo = FALSE)$coefficients,
                   isat(y, gof.function=list(name="myGof"), gof.method="min", parallel.options=2, print.searchinfo = FALSE)$coefficients)
  
  ##all three user-defined:
  expect_silent(isat(y, user.diagnostics = list(name = "SWtest", pval = 1e-10),
                     user.estimator = list(name = "myEstimator"),
                     gof.function = list(name = "myGof"), gof.method = "min",
                     parallel.options = 2, print.searchinfo = FALSE))
  
  
})



##################################################
## 8 TEST ROBUST COEFFICIENT COVARIANCES
##################################################

# test_that("TEST ROBUST COEFFICIENT COVARIANCES",{
# 
# M-orca 02/10/2022 disabled as continues to throw errors - we need to fix this overall 
# TODO fix this
# skip("Testing of robust coefficient covariances currently skipped")
# 
## "white" and "newey-west" coefficient covariances generally lead to
## either outright errors, or strange results. Currently, therefore,
## until more numerically stable versions are derived, users are
## discouraged to use these robust coefficient covariances. The tests
## here, therefore, only check whether the arguments work or not. The
## last set of checks suggest the source of the problem is near-zero
## residuals.
# 
# set.seed(123) 
# dgp.n <- 100
# y <- rnorm(dgp.n)
# plotarg <- FALSE
# 
# ##test "white" vcov (all yield errors?):
# expect_warning(isat(y, iis=TRUE, sis=FALSE, vcov.type="w", print.searchinfo = FALSE))
# expect_warning(expect_error(isat(y, iis=FALSE, sis=TRUE, vcov.type="w", print.searchinfo = FALSE)))
# expect_warning(expect_error(isat(y, iis=TRUE, sis=TRUE, vcov.type="w", print.searchinfo = FALSE)))
# expect_warning(expect_error(isat(y, iis=FALSE, tis=TRUE, vcov.type="w", print.searchinfo = FALSE)))
# expect_warning(expect_error(isat(y, iis=TRUE, tis=TRUE, vcov.type="w", print.searchinfo = FALSE)))
# expect_warning(expect_error(isat(y, iis=FALSE, sis=TRUE, vcov.type="w", tis=TRUE, plot=plotarg, print.searchinfo = FALSE)))
# expect_warning(expect_error(isat(y, iis=TRUE, sis=TRUE, vcov.type="w", tis=TRUE, plot=plotarg, print.searchinfo = FALSE)))
# 
# ##test "newey-west" vcov:
# set.seed(123) 
# dgp.n <- 100
# y <- rnorm(dgp.n)
# #y <- rt(dgp.n, df=4.1)
# mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
# plotarg <- FALSE
# 
# ##test "newey-west" vcov:
# isat(y, iis=TRUE, sis=FALSE, vcov.type="n", plot=plotarg)
# #these ones yields errors!:
# isat(y, iis=FALSE, sis=TRUE, vcov.type="n")
# isat(y, iis=TRUE, sis=TRUE, vcov.type="n")
# isat(y, iis=FALSE, tis=TRUE, vcov.type="n")
# isat(y, iis=TRUE, tis=TRUE, vcov.type="n")
# isat(y, iis=FALSE, sis=TRUE, vcov.type="n", tis=TRUE, plot=plotarg)
# isat(y, iis=TRUE, sis=TRUE, vcov.type="w", tis=TRUE, plot=plotarg)
# 
# ##code that sheds light on what the possible source of
# ##the problem is (near-zero residuals)
# set.seed(123)
# y <- rnorm(100, 0, 1)
# isat(y, vcov.type=c("newey-west"))
# isat(y, sis=FALSE, uis=sim(y, which.ones=1:10))
# isat(y, sis=FALSE, uis=sim(y, which.ones=1:10), vcov.type="newey-west")
# isat(y, sis=FALSE, uis=sim(y, which.ones=seq(1,20,2)), vcov.type="newey-west")
# 
# })




##################################################
## 9 SIMULATIONS
##################################################

##For the future...




# ##################################################
# ## 10 Real world applications/bugs found by users
# ##################################################



test_that("TEST MAIN isat() ARGUMENTS - perfect linearity works", {
  data(hpdata)
  y <- zooreg(hpdata$GCQ, 1959, frequency=4)
  dlogy <- diff(log(y))
  x <- zooreg(hpdata$GYDQ, 1959, frequency=4)
  xmatrix <- as.zoo(ts(data.frame(GYDQ = diff(log(x))), start = index(dlogy)[1], end = index(dlogy)[144], frequency = 4))
  
  a <- isat(dlogy, mxreg=xmatrix[,"GYDQ",drop = FALSE], iis=TRUE, sis = FALSE, t.pval =0.05, plot = FALSE, print.searchinfo = FALSE)
  
  xmatrix$GYDQ_lincom <- xmatrix$GYDQ * 4
  b <- isat(dlogy, mxreg=xmatrix, iis=TRUE, sis = FALSE, t.pval =0.05, plot = FALSE, print.searchinfo = FALSE)
  
  ##should be identical, bu were not in (until?) version 0.36:
  expect_identical(coef(a), coef(b))
  
  # second identified issue by M-orca 25/09/2022
  set.seed(123)
  y <- rnorm(100)
  
  mxreg_other <- data.frame(z = rep(2,100), # perfectly collinear with the intercept
                            x = rnorm(100))
  
  a <- isat(y, mc = TRUE, mxreg = mxreg_other, plot = FALSE, print.searchinfo = FALSE)
  b <- isat(y, mc = TRUE, mxreg = mxreg_other["x"],plot = FALSE, print.searchinfo = FALSE)
  
  expect_identical(a$coefficients, b$coefficients)
})
