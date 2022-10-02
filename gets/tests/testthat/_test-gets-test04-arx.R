##################################################
## Test file for gets package. First created
## 23 September 2014, Oslo.
##
## 1 INITIATE
## 2 TEST MAIN arx() ARGUMENTS
## 3 MORE TESTS OF predict.arx
## 3 TEST USER-DEFINED DIAGNOSTICS
## 4 TEST USER-DEFINED ESTIMATION
## 5 SIMULATIONS (FOR THE FUTURE)
##
##################################################

##################################################
##1 INITIATE
##################################################

# not needed. Simply click "Run tests" in RStudio in the top right or devtools::test_file("test-filename.R")
# to test this manually (line by line), first load the gets package using:
# devtools::load_all()


##################################################
## 2 TEST MAIN arx() ARGUMENTS
##################################################

getOption("plot")
options(plot=TRUE)
#options(plot=FALSE)
#options(plot=NULL)

##generate some data:
##===================

set.seed(123)
y.n <- 60
y <- arima.sim(list(ar=0.7),y.n)
y <- ts(y, frequency=4, end=c(2015,4))
mX <- matrix(rnorm(10*y.n), y.n, 10)
colnames(mX) <- paste("xvar", 1:NCOL(mX), sep="")
mX <- ts(mX, frequency=4, end=c(2015,4))
mX[1:5,2] <- NA
stepdum1 <- sim(y, which.ones=index(y)[floor(y.n/2)])
stepdum2 <- sim(y, which.ones=index(y)[floor(y.n/1.5)])
mX <- cbind(as.zoo(mX), as.zoo(stepdum1), as.zoo(stepdum2))
y[1] <- NA
y[y.n] <- NA

##test each argument separately and together:
arx(y) #should return "Warning message: In plot.arx(out) : No estimated...etc."
arx(y, normality.JarqueB=TRUE)
arx(y, mc=TRUE)
arx(y, ar=c(1,3))
arx(y, ewma=list(length=c(2,4)))
arx(y, mxreg=mX)
arx(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)), mxreg=mX)
arx(y, vc=TRUE)
arx(y, arch=c(2,4))
arx(y, asym=c(1,3))
arx(y, log.ewma=list(length=c(3,5)))
arx(y, vxreg=cbind(log(mX[,1:2]^2), mX[,3:4]))
arx(y, vc=TRUE, arch=c(2,4), asym=c(1,3),log.ewma=list(length=c(3,5)),vxreg=cbind(log(mX[,1:2]^2), mX[,3:4]))
arx(y, mc=TRUE, ar=c(1,3), vcov.type="o")
arx(y, mc=TRUE, ar=c(1,3), vcov.type="w")
arx(y, mc=TRUE, ar=c(1,3), vcov.type="n")
arx(y, mc=TRUE, ar=c(1,3), qstat.options=c(5,5))
arx(y, mc=TRUE, ar=c(1,3), tol=1e-15)
# arx(y, mc=TRUE, ar=c(1,3), tol=1, LAPACK=FALSE) #should crash
arx(y, mc=TRUE, ar=c(1,3), tol=1, LAPACK=TRUE)


# Unit Testing

test_that("arx test each argument",{
  expect_message(arx(y)) #should return "Warning message: In plot.arx(out) : No estimated...etc."
  expect_message(arx(y, normality.JarqueB=TRUE))
  
  expect_silent(arx(y, mc=TRUE))
  expect_silent(arx(y, ar=c(1,3)))
  expect_silent(arx(y, ewma=list(length=c(2,4))))
  expect_silent(arx(y, mxreg=mX))
  expect_silent(arx(y, mc=TRUE, ar=c(1,3), ewma=list(length=c(2,4)), mxreg=mX))
  expect_silent(arx(y, vc=TRUE))
  expect_silent(arx(y, arch=c(2,4)))
  expect_silent(arx(y, asym=c(1,3)))
  expect_silent(arx(y, log.ewma=list(length=c(3,5))))
  expect_silent(arx(y, vxreg=cbind(log(mX[,1:2]^2), mX[,3:4])))
  expect_silent(arx(y, vc=TRUE, arch=c(2,4), asym=c(1,3),log.ewma=list(length=c(3,5)),vxreg=cbind(log(mX[,1:2]^2), mX[,3:4])))
  expect_silent(arx(y, mc=TRUE, ar=c(1,3), vcov.type="o"))
  expect_silent(arx(y, mc=TRUE, ar=c(1,3), vcov.type="w"))
  expect_silent(arx(y, mc=TRUE, ar=c(1,3), vcov.type="n"))
  expect_silent(arx(y, mc=TRUE, ar=c(1,3), qstat.options=c(5,5)))
  expect_silent(arx(y, mc=TRUE, ar=c(1,3), tol=1e-15))
  
  expect_error(arx(y, mc=TRUE, ar=c(1,3), tol=1, LAPACK=FALSE)) #should crash
  
  expect_silent(arx(y, mc=TRUE, ar=c(1,3), tol=1, LAPACK=TRUE))
})

##only mean specification:
mod01 <- arx(y, ar=1:4, mxreg=mX)
print(mod01)
print(mod01, signif.stars=TRUE) # Nov 2020 Moritz: don't think we need this anymore, our default is signif.stars=TRUE anyway
coef(mod01)
coef(mod01, spec="m")
coef(mod01, spec="v") #should be NULL
coef(mod01, spec="b") #should be the same as "m"

test_that("Test arx Mean specification output",{
  expect_output(print(mod01))
  expect_output(print(mod01, signif.stars=TRUE)) # Nov 2020 M-orca: don't think we need this anymore, our default is signif.stars=TRUE anyway
  expect_vector(coef(mod01))
  expect_vector(coef(mod01, spec="m"))
  expect_null(coef(mod01, spec="v")) #should be NULL
  expect_vector(coef(mod01, spec="b")) #should be the same as "m"
  
})



plot(ES(mod01, level=c(0.99,0.95,0.9))) #expected shortfall
plot(cbind(fitted(mod01),
           fitted(mod01, spec="m"),
           fitted(mod01, spec="v")))

save_png <- function(code, width = 1000, height = 600) {
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = width, height = height)
  on.exit(dev.off())
  code
  
  path
}


test_that("Test the arx plots (snapshot - skipped on CI)",{
  skip_on_ci()
  
  expect_snapshot_file(cran = FALSE, path = save_png(plot(ES(mod01, level=c(0.99,0.95,0.9)))), name = "arx_expectedshotfall_mean_1.png")
  expect_snapshot_file(cran = FALSE, path = save_png(plot(cbind(fitted(mod01),
                                                                fitted(mod01, spec="m"),
                                                                fitted(mod01, spec="v")))), name = "arx_fittedplots.png")
})

fitted(mod01, spec="b") #should be NULL
logLik(mod01)
plot(mod01)
#predict(mod01) #should return the error-message: 'newmxreg' is NULL

test_that("arx fitted and loglik and simple plot",{
  expect_null(fitted(mod01, spec="b")) #should be NULL
  expect_true(isClass("logLik",logLik(mod01)))
  skip_on_ci()
  expect_snapshot_file(cran = FALSE, path = save_png(plot(mod01)), name = "arx_standardplot.png")
})


# Predict

# Visual Inspection
# predict(mod01) #should return the error-message: 'newmxreg' is NULL)
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)))
predict(mod01, n.ahead=30, newmxreg=matrix(0,30,NCOL(mX)))
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)), spec="mean")
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),plot=FALSE, return=FALSE) #no plot, return nothing
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)), plot=FALSE)
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),return=FALSE, plot=TRUE) #plot only
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),plot.options=list(fitted=TRUE))
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),plot.options=list(newmactual=rep(0,5)))
predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),spec="variance") #should return "Set 'vc=TRUE'..etc."
predict(mod01, n.ahead=5, spec="variance") #should return "Set 'vc=TRUE'..etc."

# Unit testing

test_that("Testing arx predictions mean",{
  expect_error(predict(mod01)) #should return the error-message: 'newmxreg' is NULL)
  expect_silent(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX))))
  expect_silent(predict(mod01, n.ahead=30, newmxreg=matrix(0,30,NCOL(mX))))
  expect_silent(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)), spec="mean"))
  expect_silent(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),plot=FALSE, return=FALSE)) #no plot, return nothing
  expect_silent(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)), plot=FALSE))
  expect_silent(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),return=FALSE, plot=TRUE)) #plot only
  expect_silent(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),plot.options=list(fitted=TRUE)))
  expect_silent(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),plot.options=list(newmactual=rep(0,5))))
  expect_message(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),spec="variance"),regexp = "Set 'vc = TRUE' to enable a plot of the variance predictions") #should return "Set 'vc=TRUE'..etc."
  expect_message(predict(mod01, n.ahead=5, spec="variance"),regexp = "Set 'vc = TRUE' to enable a plot of the variance predictions") #should return "Set 'vc=TRUE'..etc."
  
  # Snapshot testing
  skip_on_ci()
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)))), name = "arx_predict_mean_1.png")
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod01, n.ahead=30, newmxreg=matrix(0,30,NCOL(mX)))), name = "arx_predict_mean_2.png")
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)), spec="mean")), name = "arx_predict_mean_3.png")
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),return=FALSE, plot=TRUE)), name = "arx_predict_mean_4.png") #plot only
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),plot.options=list(fitted=TRUE))), name = "arx_predict_mean_5.png")
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod01, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),plot.options=list(newmactual=rep(0,5)))), name = "arx_predict_mean_6.png")
  
})


# Recursive testing

# Visual Inspection
recursive(mod01)
recursive(mod01, return=FALSE) #only plot
recursive(mod01, plot=FALSE) #only return (values)
recursive(mod01, plot=FALSE, return=FALSE) #return nothing
recursive(mod01, std.errors=FALSE)
# recursive(arx(y)) #should return the error-message: No mean-equation
# recursive(arx(y), spec="variance") #should return the error-message No variance-equation
recursive(arx(y, mc=TRUE, plot=FALSE))
recursive(arx(y, mc=TRUE, plot=FALSE), return=FALSE)
recursive(arx(y, ar=1, plot=FALSE))
recursive(arx(y, mxreg=mX, plot=FALSE))
plot(cbind(residuals(mod01),
           residuals(mod01, std=FALSE),
           residuals(mod01, std=TRUE)))
sigma(mod01)
rsquared(mod01)
summary(mod01)
plot(VaR(mod01, level=c(0.99,0.95,0.9))) #value-at-risk
vcov(mod01)
vcov(mod01, spec="m")
vcov(mod01, spec="v") #should return NULL

# Unit testing

test_that("Recursive testing",{
  # Output tests
  
  expect_silent(recursive(mod01))
  expect_silent(recursive(mod01, plot=FALSE)) #only return (values)
  expect_silent(recursive(mod01, plot=FALSE, return=FALSE)) #return nothing
  expect_silent(recursive(mod01, std.errors=FALSE))
  
  expect_error(recursive(arx(y))) #should return the error-message: No mean-equation
  expect_error(recursive(arx(y), spec="variance")) #should return the error-message No variance-equation
  
  expect_silent(recursive(arx(y, mc=TRUE, plot=FALSE)))
  expect_silent(recursive(arx(y, mc=TRUE, plot=FALSE), return=FALSE))
  expect_silent(recursive(arx(y, ar=1, plot=FALSE)))
  expect_silent(recursive(arx(y, mxreg=mX, plot=FALSE)))
  
  expect_vector(sigma(mod01))
  expect_vector(rsquared(mod01))
  expect_silent(summary(mod01))
  # expect_output(plot(VaR(mod01, level=c(0.99,0.95,0.9)))) #value-at-risk
  expect_vector(vcov(mod01))
  expect_vector(vcov(mod01, spec="m"))
  expect_null(vcov(mod01, spec="v")) #should return NULL
  
  
  # Snapshot tests
  skip_on_ci()
  
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(mod01)), name = "arx_recursive_mean_1.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(mod01, return=FALSE)), name = "arx_recursive_mean_2.png") #only plot
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(mod01, std.errors=FALSE)), name = "arx_recursive_mean_3.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(arx(y, mc=TRUE, plot=FALSE))), name = "arx_recursive_mean_4.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(arx(y, mc=TRUE, plot=FALSE), return=FALSE)), name = "arx_recursive_mean_5.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(arx(y, ar=1, plot=FALSE))), name = "arx_recursive_mean_6.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(arx(y, mxreg=mX, plot=FALSE))), name = "arx_recursive_mean_7.png")
  expect_snapshot_file(cran = FALSE, path = save_png(plot(cbind(residuals(mod01),
                                                                residuals(mod01, std=FALSE),
                                                                residuals(mod01, std=TRUE)))), name = "arx_recursive_mean_8.png")
  expect_snapshot_file(cran = FALSE, path = save_png(plot(VaR(mod01, level=c(0.99,0.95,0.9)))), name = "arx_recursive_mean_9.png") #value-at-risk
})


##both mean and variance specifications:
# Visual Inspection
mX <- mX[,1:2]
mod02 <- arx(y, mc=TRUE, ar=1:3, mxreg=mX,arch=1:4,asym=1:2, vxreg=log(mX^2))
print(mod02)
print(mod02, signif.stars=TRUE)
coef(mod02)
coef(mod02, spec="m")
coef(mod02, spec="v")
coef(mod02, spec="b")
plot(ES(mod02, level=c(0.99,0.95,0.9))) #expected shortfall
plot(cbind(fitted(mod02),
           fitted(mod02, spec="m"),
           fitted(mod02, spec="v"),
           fitted(mod02, spec="b")))
logLik(mod02)
plot(mod02)

# Unit Testing

test_that("arx Mean and Variance Specification",{
  
  expect_output(print(mod02))
  expect_output(print(mod02, signif.stars=TRUE))
  expect_vector(coef(mod02))
  expect_vector(coef(mod02, spec="m"))
  expect_vector(coef(mod02, spec="v"))
  expect_vector(coef(mod02, spec="b"))
  expect_true(isClass("logLik",logLik(mod02)))
  
  
  # snapshot testing - skipped on ci
  skip_on_ci()
  expect_snapshot_file(cran = FALSE, path = save_png(plot(ES(mod02, level=c(0.99,0.95,0.9)))),
                       name = "arx_expectedshotfall_meanvariance_1.png") #expected shortfall
  expect_snapshot_file(cran = FALSE, path = save_png(plot(cbind(fitted(mod02),
                                                                fitted(mod02, spec="m"),
                                                                fitted(mod02, spec="v"),
                                                                fitted(mod02, spec="b")))),
                       name = "arx_fitted_meanvariance_1.png")
  
  expect_snapshot_file(cran = FALSE, path = save_png(plot(mod02)),name = "arx_meanvariance_1.png")
})


# visual inspection

# predict(mod02) #should return the error-message 'newvxreg' is NULL
# predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX))) #should return the error-message 'newmxreg' is NULL
predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),newmxreg=matrix(0,5,NCOL(mX)))
predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),newmxreg=matrix(0,5,NCOL(mX)), plot=FALSE)
predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)))
predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),plot.options=list(fitted=TRUE))
predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),plot.options=list(newvactual=rep(1,5)))
predict(mod02, n.ahead=30, newvxreg=matrix(0,30,NCOL(mX)),newmxreg=matrix(0,30,NCOL(mX)))
predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),newmxreg=matrix(0,5,NCOL(mX)), spec="mean")
# predict(mod02, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),spec="variance") #should return the error-message 'newvxreg' is NULL


test_that("Testing arx predictions mean and variance",{
  expect_error(predict(mod02)) #should return the error-message 'newvxreg' is NULL
  expect_error(predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)))) #should return the error-message 'newmxreg' is NULL
  
  expect_vector(predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),newmxreg=matrix(0,5,NCOL(mX))))
  expect_vector(predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),newmxreg=matrix(0,5,NCOL(mX)), plot=FALSE))
  expect_vector(predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX))))
  expect_vector(predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),plot.options=list(fitted=TRUE)))
  expect_vector(predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),plot.options=list(newvactual=rep(1,5))))
  expect_vector(predict(mod02, n.ahead=30, newvxreg=matrix(0,30,NCOL(mX)),newmxreg=matrix(0,30,NCOL(mX))))
  expect_vector(predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),newmxreg=matrix(0,5,NCOL(mX)), spec="mean"))
  expect_error(predict(mod02, n.ahead=5, newmxreg=matrix(0,5,NCOL(mX)),spec="variance")) #should return the error-message 'newvxreg' is NULL
  
  
  # Snapshot - skip on ci
  skip_on_ci()
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),newmxreg=matrix(0,5,NCOL(mX)))), name = "arx_predict_meanvariance_1.png")
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)))), name = "arx_predict_meanvariance_2.png")
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),plot.options=list(fitted=TRUE))), name = "arx_predict_meanvariance_3.png")
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod02, spec="variance", n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),plot.options=list(newvactual=rep(1,5)))), name = "arx_predict_meanvariance_4.png")
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod02, n.ahead=30, newvxreg=matrix(0,30,NCOL(mX)),newmxreg=matrix(0,30,NCOL(mX)))), name = "arx_predict_meanvariance_5.png")
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod02, n.ahead=5, newvxreg=matrix(0,5,NCOL(mX)),newmxreg=matrix(0,5,NCOL(mX)), spec="mean")), name = "arx_predict_meanvariance_6.png")
})


# Recursive

# Visual Inspection
recursive(mod02)
recursive(mod02, spec="variance")
recursive(mod02, spec="variance", return=FALSE)
recursive(mod02, spec="variance", plot=FALSE, return=FALSE) #return nothing
recursive(mod02, spec="variance", std.errors=FALSE)
recursive(arx(y, vc=TRUE), spec="variance")
recursive(arx(y, arch=1), spec="variance")
recursive(arx(y, vxreg=mX), spec="variance")
recursive(arx(y, vc=TRUE), spec="variance", return=FALSE)
plot(cbind(residuals(mod02),
           residuals(mod02, std=FALSE),
           residuals(mod02, std=TRUE)))
summary(mod02)
plot(VaR(mod02, level=c(0.99,0.95,0.9))) #value-at-risk
vcov(mod02)
vcov(mod02, spec="m")
vcov(mod02, spec="v")


# Unit Testing

test_that("Recursive Tests arx Mean and Variance specification",{
  
  expect_silent(recursive(mod02))
  expect_silent(recursive(mod02, spec="variance"))
  expect_silent(recursive(mod02, spec="variance", return=FALSE))
  expect_silent(recursive(mod02, spec="variance", plot=FALSE, return=FALSE)) #return nothing
  expect_silent(recursive(mod02, spec="variance", std.errors=FALSE))
  expect_silent(recursive(arx(y, vc=TRUE), spec="variance"))
  expect_silent(recursive(arx(y, arch=1), spec="variance"))
  expect_silent(recursive(arx(y, vxreg=mX), spec="variance"))
  expect_silent(recursive(arx(y, vc=TRUE), spec="variance", return=FALSE))
  expect_true(isClass("summaryDefault",summary(mod02)))
  expect_vector(vcov(mod02))
  expect_vector(vcov(mod02, spec="m"))
  expect_vector(vcov(mod02, spec="v"))
  
  
  # snapshot testing
  skip_on_ci()
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(mod02)), name = "arx_recursive_meanvariance_1.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(mod02, spec="variance")), name = "arx_recursive_meanvariance_2.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(mod02, spec="variance", return=FALSE)), name = "arx_recursive_meanvariance_3.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(mod02, spec="variance", std.errors=FALSE)), name = "arx_recursive_meanvariance_4.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(arx(y, vc=TRUE), spec="variance")), name = "arx_recursive_meanvariance_5.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(arx(y, arch=1), spec="variance")), name = "arx_recursive_meanvariance_6.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(arx(y, vxreg=mX), spec="variance")), name = "arx_recursive_meanvariance_7.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(arx(y, vc=TRUE), spec="variance", return=FALSE)), name = "arx_recursive_meanvariance_8.png")
  expect_snapshot_file(cran = FALSE, path = save_png(plot(cbind(residuals(mod02),
                                                                residuals(mod02, std=FALSE),
                                                                residuals(mod02, std=TRUE)))), name = "arx_recursive_meanvariance_9.png")
  expect_snapshot_file(cran = FALSE, path = save_png(plot(VaR(mod02, level=c(0.99,0.95,0.9)))), name = "arx_recursive_meanvariance_10.png") #value-at-risk
  
})


##only variance specification:

# Visual Inspection
mod03 <- arx(y, arch=1:4, asym=1:2, log.ewma=3, vxreg=log(mX^2))
print(mod03)
print(mod03, signif.stars=TRUE)
coef(mod03)
coef(mod03, spec="m") #should be NULL
coef(mod03, spec="v")
coef(mod03, spec="b")
plot(ES(mod01, level=c(0.99,0.95,0.9)), plot.type="single", col=c("blue","red","green4")) #expected shortfall
plot(cbind(fitted(mod03),
           fitted(mod03, spec="m"),
           fitted(mod03, spec="v")))
fitted(mod03, spec="b") #should be NULL
logLik(mod03)
plot(mod03)
plot(cbind(residuals(mod03),
           residuals(mod03, std=FALSE),
           residuals(mod03, std=TRUE)))
predict(mod03, n.ahead=24, newvxreg=matrix(0,24,2))
predict(mod03, n.ahead=24, newvxreg=matrix(0,24,2), plot.options=list(fitted=TRUE, newvactual=rep(1,24)))
# recursive(mod03) #should return the error-message No mean-equation
recursive(mod03, spec="variance")
recursive(mod03, spec="variance", return=FALSE)
recursive(mod03, spec="variance", plot=FALSE, return=FALSE) #return nothing
recursive(mod03, spec="variance", std.errors=FALSE)
summary(mod03)
plot(VaR(mod03, level=c(0.99,0.95,0.9)), plot.type="single", col=c("blue","red","green4")) #value-at-risk
vcov(mod03)

# Unit testing

test_that("arx only variance specification",{
  expect_silent(mod03 <- arx(y, arch=1:4, asym=1:2, log.ewma=3, vxreg=log(mX^2)))
  expect_output(print(mod03))
  expect_output(print(mod03, signif.stars=TRUE))
  expect_vector(coef(mod03))
  expect_null(coef(mod03, spec="m")) #should be NULL
  expect_vector(coef(mod03, spec="v"))
  expect_vector(coef(mod03, spec="b"))
  expect_null(fitted(mod03, spec="b")) #should be NULL
  expect_vector(logLik(mod03))
  expect_vector(predict(mod03, n.ahead=24, newvxreg=matrix(0,24,2)))
  expect_vector(predict(mod03, n.ahead=24, newvxreg=matrix(0,24,2), plot.options=list(fitted=TRUE, newvactual=rep(1,24))))
  expect_error(recursive(mod03)) #should return the error-message No mean-equation
  expect_silent(recursive(mod03, spec="variance"))
  expect_silent(recursive(mod03, spec="variance", plot=FALSE, return=FALSE))#return nothing
  expect_silent(recursive(mod03, spec="variance", std.errors=FALSE))
  expect_silent(summary(mod03))
  expect_silent(vcov(mod03))
  
  # Snapshot testing
  skip_on_ci()
  
  expect_snapshot_file(cran = FALSE, path = save_png(mod03 <- arx(y, arch=1:4, asym=1:2, log.ewma=3, vxreg=log(mX^2))), name = "arx_variance_1.png")
  expect_snapshot_file(cran = FALSE, path = save_png(plot(ES(mod01, level=c(0.99,0.95,0.9)), plot.type="single", col=c("blue","red","green4"))), name = "arx_variance_2.png") #expected shortfall
  expect_snapshot_file(cran = FALSE, path = save_png(plot(cbind(fitted(mod03),
                                                                fitted(mod03, spec="m"),
                                                                fitted(mod03, spec="v")))), name = "arx_variance_3.png")
  expect_snapshot_file(cran = FALSE, path = save_png(plot(mod03)), name = "arx_variance_4.png")
  expect_snapshot_file(cran = FALSE, path = save_png(plot(cbind(residuals(mod03),
                                                                residuals(mod03, std=FALSE),
                                                                residuals(mod03, std=TRUE)))), name = "arx_variance_5.png")
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod03, n.ahead=24, newvxreg=matrix(0,24,2))), name = "arx_variance_6.png")
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mod03, n.ahead=24, newvxreg=matrix(0,24,2), plot.options=list(fitted=TRUE, newvactual=rep(1,24)))), name = "arx_variance_7.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(mod03, spec="variance", return=FALSE)), name = "arx_variance_8.png")
  expect_snapshot_file(cran = FALSE, path = save_png(recursive(mod03, spec="variance", std.errors=FALSE)), name = "arx_variance_9.png")
  expect_snapshot_file(cran = FALSE, path = save_png(plot(VaR(mod03, level=c(0.99,0.95,0.9)), plot.type="single", col=c("blue","red","green4"))), name = "arx_variance_10.png") #value-at-risk
})

##ISSUE!!! DISCOVERED 11/8-2020 BY G.:
## M-orca: not included in unit testing for now (Nov 2020)
vcov(mod03, spec="m") #should return NULL
vcov(mod03, spec="v")


##################################################
## 3 TEST USER-DEFINED DIAGNOSTICS
##################################################

##user-defined Shapiro-Wilks test for normality in the residuals:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  result <- c(tmp$statistic, NA, tmp$p.value)
  return(result)
}

# mod06 <- arx(y, ar=1:4, mxreg=mX, user.diagnostics=list(name="SWtest"))
#print(mod06)
#print(mod06, signif.stars=TRUE)
#mod06 <- arx(y, ar=1:4, mxreg=mX,user.diagnostics=list(name="SWtest", pval=0.025))
# the pval argument is ignored (as it should), I think
#print(mod06)


test_that("User defined diagnostics",{
  mod06 <- arx(y, ar=1:4, mxreg=mX, user.diagnostics=list(name="SWtest", envir = environment(SWtest)))
  expect_output(print(mod06))
  expect_output(print(mod06, signif.stars=TRUE))
  mod06 <- arx(y, ar=1:4, mxreg=mX,user.diagnostics=list(name="SWtest", pval=0.025, envir = environment(SWtest)))
  expect_output(print(mod06))
  
})

##test the envir entry:
## M-Orca: this does work
rm("SWtest") #make sure SWtest is not defined in the global environment
myenv <- new.env()
assign("SWtest",
       function(x, ...){
         tmp <- shapiro.test(x$residuals)
         result <- c(tmp$statistic, NA, tmp$p.value)
         result <- rbind( as.numeric(c(tmp$statistic, NA, tmp$p.value)) )
         rownames(result) <- "myenv-test"
         return(result)
       },
       envir=myenv) #end function
# mod06 <- arx(y, ar=1:4, mxreg=mX, user.diagnostics=list(name="SWtest")) #should not work
mod06 <- arx(y, ar=1:4, mxreg=mX, user.diagnostics=list(name="SWtest", envir=myenv)) #should work
print(mod06)

test_that("Testing the user-defined environment",{
  # expect_error(arx(y, ar=1:4, mxreg=mX, user.diagnostics=list(name="SWtest"))) #should not work
  mod06 <- arx(y, ar=1:4, mxreg=mX, user.diagnostics=list(name="SWtest", envir=myenv)) #should work
  expect_output(print(mod06))
  
})

##################################################
## 4 TEST USER-DEFINED ESTIMATION
##################################################

## Rules arx: The returned result should be a list with at least three items
## named "coefficients", "df" and "vcov". The item named "df" is used to
## compute the p-values associated with the t-statistics: coef/std.err.
## The returned result should be a list with at least three items:
## "coefficients", "df" and "vcov".

##user-defined estimator (minimal):
Gfun <- function(y, x, method=3){
  tmp <- ols(y, x, method=method)
  tmp <- list(coefficients=tmp$coefficients, df=tmp$df, vcov=tmp$vcov)
  return(tmp)
}

# Visual Inspection

mod07 <- arx(y, ar=1:4, mxreg=mX,user.estimator=list(name="Gfun", envir=environment(Gfun)), plot=FALSE)
summary(mod07)
print(mod07)
# mod07 <- arx(y, ar=1:4, mxreg=mX, user.estimator=list(name="Gfun"), plot=TRUE) #should work but produce warning
summary(mod07)
print(mod07)
coef(mod07)
coef(mod07, spec="m")
coef(mod07, spec="v") #should be null
coef(mod07, spec="b") #mean coefs only
residuals(mod07) #should be null
residuals(mod07, std=FALSE) #should be null
residuals(mod07, std=TRUE) #should be null
fitted(mod07) #should be null
fitted(mod07, spec="m")
fitted(mod07, spec="v")
fitted(mod07, spec="b") #should be NULL
# logLik(mod07) #should produce warning
# plot(mod07) #should produce warning
# recursive(mod07) #should return the error-message "Not available..."
vcov(mod07)
vcov(mod07, spec="m")
vcov(mod07, spec="v")


# Unit testing

test_that("TEST USER-DEFINED ESTIMATION",{
  expect_silent(mod07 <- arx(y, ar=1:4, mxreg=mX,user.estimator=list(name="Gfun",envir=environment(Gfun)), plot=FALSE))
  expect_silent(summary(mod07))
  expect_output(print(mod07))
  expect_message(mod07 <- arx(y, ar=1:4, mxreg=mX, user.estimator=list(name="Gfun",envir=environment(Gfun)), plot=TRUE)) #should work but produce warning
  expect_silent(summary(mod07))
  expect_output(print(mod07))
  expect_vector(coef(mod07))
  expect_vector(coef(mod07, spec="m"))
  expect_null(coef(mod07, spec="v")) #should be null
  expect_vector(coef(mod07, spec="b")) #mean coefs only
  expect_null(residuals(mod07)) #should be null
  expect_null(residuals(mod07, std=FALSE)) #should be null
  expect_null(residuals(mod07, std=TRUE)) #should be null
  expect_null(fitted(mod07)) #should be null
  expect_null(fitted(mod07, spec="m"))
  expect_null(fitted(mod07, spec="v"))
  expect_null(fitted(mod07, spec="b")) #should be NULL
  expect_warning(logLik(mod07)) #should produce warning
  expect_message(plot(mod07)) #should produce warning
  expect_error(recursive(mod07)) #should return the error-message "Not available..."
  expect_silent(vcov(mod07))
  expect_null(vcov(mod07, spec="m"))
  expect_null(vcov(mod07, spec="v"))
})




##user-defined estimator (usual):
Gfun <- function(y, x, ...){
  tmp <- ols(y, x, method=3)
  tmp$mean.fit <- tmp$fit
  return(tmp)
}
mod08 <- arx(y, ar=1:4, mxreg=mX,user.estimator=list(name="Gfun", envir=environment(Gfun)), plot=FALSE)
summary(mod08)
print(mod08)
# mod08 <- arx(y, ar=1:4, mxreg=mX,user.estimator=list(name="Gfun"), plot=TRUE) #should produce warning
summary(mod08)
print(mod08)
coef(mod08)
coef(mod08, spec="m")
coef(mod08, spec="v") #should be null
coef(mod08, spec="b")
residuals(mod08)
residuals(mod08, std=FALSE)
residuals(mod08, std=TRUE) ##should be NULL
fitted(mod08)
fitted(mod08, spec="m")
fitted(mod08, spec="v") #should be NULL
fitted(mod08, spec="b") #should be NULL
logLik(mod08)
plot(mod08) #should produce warning

test_that("Test user-defined estimator G-fun usual",{
  expect_silent(mod08 <- arx(y, ar=1:4, mxreg=mX,user.estimator=list(name="Gfun",envir= environment(Gfun)), plot=FALSE))
  expect_silent(summary(mod08))
  expect_output(print(mod08))
  expect_message(mod08 <- arx(y, ar=1:4, mxreg=mX,user.estimator=list(name="Gfun",envir= environment(Gfun)), plot=TRUE)) #should produce warning
  expect_silent(summary(mod08))
  expect_output(print(mod08))
  expect_vector(coef(mod08))
  expect_vector(coef(mod08, spec="m"))
  expect_null(coef(mod08, spec="v")) #should be null
  expect_vector(coef(mod08, spec="b"))
  expect_vector(residuals(mod08))
  expect_vector(residuals(mod08, std=FALSE))
  expect_null(residuals(mod08, std=TRUE)) ##should be NULL
  expect_vector(fitted(mod08))
  expect_vector(fitted(mod08, spec="m"))
  expect_null(fitted(mod08, spec="v")) #should be NULL
  expect_null(fitted(mod08, spec="b")) #should be NULL
  expect_silent(logLik(mod08))
  expect_message(plot(mod08)) #should produce warning
})


## IMPORTANT ## 
## G-Man this does not work, fixing it (i.e. changing Gfun) requires some work!:
## M-orca: excluded for now from Unit-testing (Nov 2020)

# predict(mod08, n.ahead=24, newmxreg=matrix(0,24,2))
# recursive(mod08) #should return the error-message "Not available..."
# vcov(mod08)
# vcov(mod08, spec="m")
# vcov(mod08, spec="v") #should return NULL

##user-defined estimator (fast ols):
library(Matrix)
ols2 <- function(y, x){
  out <- list()
  out$n <- length(y)
  if (is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
  out$df <- out$n - out$k
  if (out$k > 0) {
    x <- as(x, "dgeMatrix")
    out$xpy <- Matrix::crossprod(x, y) # change by M-orca Nov 2020 Matrix:: necessary for ci
    out$xtx <- Matrix::crossprod(x) # change by M-orca Nov 2020 Matrix:: necessary for ci
    out$coefficients <- as.numeric(solve(out$xtx,as.matrix(out$xpy))) # change by M-orca Nov 2020 as.matrix necessary for ci
    out$xtxinv <- solve(out$xtx)
    out$fit <- out$fit <- as.vector(x %*% out$coefficients)
  }else{ out$fit <- rep(0, out$n)	}
  out$residuals <- y - out$fit
  out$residuals2 <- out$residuals^2
  out$rss <- sum(out$residuals2)
  out$sigma2 <- out$rss/out$df
  if (out$k > 0) { out$vcov <- as.matrix(out$sigma2 * out$xtxinv) }
  out$logl <-
    -out$n * log(2 * out$sigma2 * pi)/2 - out$rss/(2 * out$sigma2)
  return(out)
}

# Visual inspection

mod09 <- arx(y, ar=1:4, mxreg=mX,user.estimator=list(name="ols2", envir=environment(ols2)), plot=FALSE)
summary(mod09)
print(mod09)
# mod09 <- arx(y, ar=1:4, mxreg=mX,user.estimator=list(name="ols2"), plot=TRUE) #should produce warning
summary(mod09)
print(mod09)
coef(mod09)
coef(mod09, spec="m")
coef(mod09, spec="v") #should be null
coef(mod09, spec="b")
residuals(mod09)
residuals(mod09, std=FALSE)
residuals(mod09, std=TRUE) #should be null
fitted(mod09) #should be null
fitted(mod09, spec="m") #should be null
fitted(mod09, spec="v") #should be null
fitted(mod09, spec="b") #should be NULL
logLik(mod09)
plot(mod09) #should produce warning
# recursive(mod09) #should return the error-message "Not available..."
vcov(mod09)
vcov(mod09, spec="m")
vcov(mod09, spec="v")

# Unit testing

test_that("ols2 user-defined testing (fast ols)",{
  expect_silent(mod09 <- arx(y, ar=1:4, mxreg=mX,user.estimator=list(name="ols2", envir=environment(ols2)), plot=FALSE))
  expect_silent(summary(mod09))
  expect_output(print(mod09))
  expect_message(mod09 <- arx(y, ar=1:4, mxreg=mX,user.estimator=list(name="ols2", envir=environment(ols2)), plot=TRUE)) #should produce warning
  expect_silent(summary(mod09))
  expect_output(print(mod09))
  expect_vector(coef(mod09))
  expect_vector(coef(mod09, spec="m"))
  expect_null(coef(mod09, spec="v")) #should be null
  expect_vector(coef(mod09, spec="b"))
  expect_vector(residuals(mod09))
  expect_vector(residuals(mod09, std=FALSE))
  expect_null(residuals(mod09, std=TRUE)) #should be null
  expect_null(fitted(mod09)) #should be null
  expect_null(fitted(mod09, spec="m")) #should be null
  expect_null(fitted(mod09, spec="v")) #should be null
  expect_null(fitted(mod09, spec="b")) #should be NULL
  expect_vector(logLik(mod09))
  expect_message(plot(mod09)) #should produce warning
  expect_error(recursive(mod09)) #should return the error-message "Not available..."
  expect_vector(vcov(mod09))
  expect_null(vcov(mod09, spec="m"))
  expect_null(vcov(mod09, spec="v"))
})


##################################################
## 5 SIMULATIONS (FOR THE FUTURE)
##################################################
