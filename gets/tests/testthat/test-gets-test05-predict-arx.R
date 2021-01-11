##################################################
## 2 TESTS OF MEAN PREDICTIONS
##################################################

##generate some data:
##===================

##small dgp:
set.seed(123)
vY <- rnorm(20)

test_that("Mean Predictions: ar(0) model without constant",{
  
  mymodel <- arx(vY, mc=FALSE)
  
  ##predictions of the mean:
  functionVals <- predict(mymodel, spec="mean", n.ahead=3)
  
  ##correct predictions:
  yhat1 <- yhat2 <- yhat3 <- 0
  correctVals <- c(yhat1,yhat2,yhat3)
  
  
  ##do they correspond?:
  expect_identical(as.vector(functionVals),correctVals)
  
})
##=============================

test_that("Mean Predictions: ar(0) model w/constant",{
  mymodel <- arx(vY, mc=TRUE)
  
  ##predictions of the mean:
  functionVals <- predict(mymodel, spec="mean", n.ahead=3)
  
  ##correct predictions:
  yhat1 <- yhat2 <- yhat3 <- coef(mymodel)[1]
  correctVals <- c(yhat1,yhat2,yhat3)
  
  ##do they correspond?:
  expect_identical(as.vector(functionVals),as.vector(correctVals))
})


##============
test_that("Mean Predictions: ar(1) model",{
  mymodel <- arx(vY, mc=TRUE, ar=1)
  
  ##predictions of the mean:
  functionVals <- predict(mymodel, spec="mean", n.ahead=12)
  
  ##correct predictions:
  yhat <- rep(NA,13)
  yhat[1] <- vY[length(vY)] #actual value at forecast origin
  for(i in 2:13){
    yhat[i] <- coef(mymodel)[1] + coef(mymodel)[2]*yhat[i-1]
  }
  correctVals <- yhat[-1]
  
  ##do they correspond?:
  expect_identical(as.vector(functionVals),correctVals)
  
})


##==============
test_that("Mean Predictions: ar(1)-x model",{
  mX <- rnorm(length(vY))
  mymodel <- arx(vY, mc=TRUE, ar=1, mxreg=mX)
  
  ##predictions of the mean:
  functionVals <- predict(mymodel, spec="mean", n.ahead=3,
                          newmxreg=matrix(1,3,1))
  
  ##correct predictions:
  yhat0 <- vY[length(vY)] #actual value at forecast origin
  yhat1 <- coef(mymodel)[1] + coef(mymodel)[2]*yhat0 + coef(mymodel)[3]*1
  yhat2 <- coef(mymodel)[1] + coef(mymodel)[2]*yhat1 + coef(mymodel)[3]*1
  yhat3 <- coef(mymodel)[1] + coef(mymodel)[2]*yhat2 + coef(mymodel)[3]*1
  correctVals <- c(yhat1,yhat2,yhat3)
  
  ##do they correspond?:
  expect_identical( as.vector(functionVals) , as.vector(correctVals) )
  
})




##===============
test_that("Mean Predictions: EqWMA(2) model",{
  mymodel <- arx(vY, mc=FALSE, ewma=list(length=2))
  regressorsMean(vY, mc = FALSE, ewma=list(length=2))
  mean( c(vY[1],vY[2]) ) #obs no. 3
  mean( c(vY[2],vY[3]) ) #obs no. 4
  mean( c(vY[18],vY[19]) ) #obs no. 20
  mean( c(tail(vY,n=2),tail(vY,n=3)) )
  
  ##predictions of the mean:
  functionVals <- predict(mymodel, spec="mean", n.ahead=3)
  
  ##correct predictions:
  yhat1 <- coef(mymodel)[1]*mean( c(vY[20],vY[19]) )
  yhat2 <- coef(mymodel)[1]*mean(c(yhat1,vY[20]))
  yhat3 <- coef(mymodel)[1]*mean(c(yhat2,yhat1))
  correctVals <- c(yhat1,yhat2,yhat3)
  
  ##do they correspond?:
  expect_identical( as.vector(functionVals) , as.vector(correctVals) )
})



##################################################
## 3 TESTS OF VARIANCE PREDICTIONS
##################################################

##generate some data:
##===================

##small dgp:
set.seed(123)
vY <- rnorm(20)

##ar(0) model without constant:
##=============================

test_that("Variance Predictions: ar(0) model without constant",{
  mymodel <- arx(vY, mc=FALSE)
  
  ##predictions of the variance:
  functionVals <- predict(mymodel, spec="variance", n.ahead=3)
  
  ##correct predictions:
  sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
  correctVals <- c(sd2hat1,sd2hat2,sd2hat3)
  
  ##do they correspond?:
  expect_identical( as.vector(functionVals) , correctVals )
  
})
# 
# 
##============
test_that("ar(0) model",{
  mymodel <- arx(vY, mc=TRUE)
  
  ##predictions of the variance:
  functionVals <- predict(mymodel, spec="variance", n.ahead=3)
  
  ##correct predictions:
  sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
  correctVals <- c(sd2hat1,sd2hat2,sd2hat3)
  
  ##do they correspond?:
  expect_identical( as.vector(functionVals) , correctVals )
  
})





##============

test_that("ar(1) model",{
  mymodel <- arx(vY, mc=TRUE, ar=1)
  
  ##predictions of the variance:
  functionVals <- predict(mymodel, spec="variance", n.ahead=3)
  
  ##correct predictions:
  sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
  correctVals <- c(sd2hat1,sd2hat2,sd2hat3)
  
  ##do they correspond?:
  expect_identical( as.vector(functionVals) , correctVals )
  #expect_equal( functionVals, correctVals )
  
})


##============================

test_that("arch(0) model with constant",{
  mymodel <- arx(vY, mc=FALSE, vc=TRUE)
  
  ##predictions of the variance:
  functionVals <- predict(mymodel, spec="variance", n.ahead=3)
  
  ##correct predictions:
  sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
  correctVals <- c(sd2hat1,sd2hat2,sd2hat3)
  
  ##do they correspond?:
  #expect_identical( as.vector(functionVals) , correctVals )
  expect_equal(as.vector(functionVals),correctVals,tolerance = 1e-10)
  
})

##arch(1) models:
##===============

save_png <- function(code, width = 1000, height = 600) {
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = width, height = height)
  on.exit(dev.off())
  code
  
  path
}

test_that("arch(1) models graphs correct", {

  mymodel1 <- arx(vY, mc=FALSE, arch=1)
  mymodel2 <- arx(vY, mc=TRUE, arch=1)
  mymodel3 <- arx(vY, mc=TRUE, ar=1, arch=1)
  
  skip_on_ci()
  
  # saves the plots as png files, they coincide with what is displayed from plot
  # when changes, will get error
  expect_snapshot_file(cran = FALSE, path = save_png(plot(mymodel1)), name = "arch1_1.png")
  expect_snapshot_file(cran = FALSE, path = save_png(plot(mymodel2)), name = "arch1_2.png")
  expect_snapshot_file(cran = FALSE, path = save_png(plot(mymodel3)), name = "arch1_3.png")
  
  # build an error
  # save model3 as error.png
  # expect_snapshot_file(cran = FALSE, path = save_png(plot(mymodel3)), name = "error.png")
  # then manipulate model3 so won't coincide with that file
  mymodel3$std.residuals[[2]] <- 1 # from 1.3517 to 1
  expect_snapshot_file(cran = FALSE, path = save_png(plot(mymodel3)), name = "error.png")
  
  # build another error
  # expect_snapshot_file(cran = FALSE, path = save_png(plot(mymodel3)), name = "error2.png")
  # change color of one of the lines (fitted)
  expect_snapshot_file(cran = FALSE, path = save_png(plot(mymodel3,col = c("green", "blue"))),
                       name = "error2.png")
  
})

test_that("arch(1) models predictions correct", {
  
  skip_on_ci()
  
  # snapshot_output saves the output, so don't have to compare to a manually
  # created vector
  mymodel <- arx(vY, mc=FALSE, arch=1)
  set.seed(123)
  expect_snapshot_output(predict(mymodel, spec="variance"))
  
  mymodel <- arx(vY, mc=TRUE, arch=1)
  set.seed(123)
  expect_snapshot_output(predict(mymodel, spec="variance"))
  
  mymodel <- arx(vY, mc=TRUE, ar=1, arch=1)
  set.seed(123)
  expect_snapshot_output(predict(mymodel, spec="variance"))
  
})

##ar(1)-x model:
##==============
mX <- rnorm(length(vY))

test_that("ar(1)-x model", {
  
  mymodel <- arx(vY, mc=TRUE, ar=1, mxreg=mX)
  
  ##predictions of the variance:
  set.seed(123)
  functionVals <- predict(mymodel, spec="variance", n.ahead=3)
  
  ##correct predictions:
  sd2hat1 <- sd2hat2 <- sd2hat3 <- sigma(mymodel)^2
  correctVals <- c(sd2hat1,sd2hat2,sd2hat3)
  
  ##do they correspond?:
  expect_identical(as.vector(functionVals), correctVals)
  
})

##arch(1)-x model:
##================

test_that("arch(1)-x model", {
  skip_on_ci()
  
  mymodel <- arx(vY, mc=FALSE, arch=1, vxreg=mX)
  expect_snapshot_file(cran = FALSE, path = save_png(plot(mymodel)), name = "arch1x_1.png")
  
  set.seed(123)
  expect_snapshot_output(
    predict(mymodel, spec="variance", n.ahead=3, newvxreg=matrix(1,3,1)))
  
})


##################################################
## 4 TESTS OF plot.options ARGUMENTS
##################################################

##ar(1) model:
##============

test_that("ar(1) plot options work", {
  
  skip_on_ci()
  
  # problem is that predict() creates a plot as a by-product / in addition to
  # the predictions
  # cannot save result of predict() and then plot that since it is only a vector
  # solution: set plot = TRUE for the predict() commands
  
  mymodel <- arx(vY, mc=TRUE, ar=1)
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE)), name = "pred_ar1_1.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(keep=1))), name = "pred_ar1_2.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(line.at.origin=TRUE))), name = "pred_ar1_3.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(start.at.origin=FALSE))), name = "pred_ar1_4.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(start.at.origin=FALSE, fitted=TRUE))), name = "pred_ar1_5.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(dot.at.origin=FALSE))), name = "pred_ar1_6.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(hlines=c(-2,-1,0,1,2)))), name = "pred_ar1_7.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(col=c("darkred","green")))), name = "pred_ar1_8.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(lty=c(3,2)))), name = "pred_ar1_9.png")
  ##only the forecast is lwd=3, should both be?:
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(lwd=3))), name = "pred_ar1_10.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(lwd=c(1,3)))), name = "pred_ar1_11.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(ylim=c(-8,8)))), name = "pred_ar1_12.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(ylab="G-values"))), name = "pred_ar1_13.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(main="Plot is slightly lower when 'main' is specified"))), name = "pred_ar1_14.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(legend.text=c("Prognose","Faktisk")))), name = "pred_ar1_15.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(fitted=TRUE))), name = "pred_ar1_16.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(newmactual=rep(0,6)))), name = "pred_ar1_17.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(shades=c(95,50)))), name = "pred_ar1_18.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(shades=c(50,95)))), name = "pred_ar1_19.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(shades=c(95,50)))), name = "pred_ar1_20.png")

  ##does not work, but should it?: jbat: no because have two lines, so need lty = c(3,3)
  #predict(mymodel, plot.options=list(lty=3))

})

##arch(1) model:
##==============

test_that("arch(1) plot options work", {
  
  skip_on_ci()
  
  mymodel <- arx(vY, mc=FALSE, vc=TRUE, arch=1)
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE)), name = "pred_arch1_1.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(keep=1))), name = "pred_arch1_2.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(line.at.origin=TRUE))), name = "pred_arch1_3.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(start.at.origin=FALSE))), name = "pred_arch1_4.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(start.at.origin=FALSE, fitted=TRUE))), name = "pred_arch1_5.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(dot.at.origin=FALSE))), name = "pred_arch1_6.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(hlines=0:4))), name = "pred_arch1_7.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(col=c("darkred","green")))), name = "pred_arch1_8.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(lty=c(3,2)))), name = "pred_arch1_9.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(lwd=3))), name = "pred_arch1_10.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(ylim=c(-6,8)))), name = "pred_arch1_11.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(ylab="G-values"))), name = "pred_arch1_12.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(main="Plot is slightly lower when 'main' is specified"))), name = "pred_arch1_13.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(legend.text=c("Prognose","Residualene kvadrert")))), name = "pred_arch1_14.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(fitted=TRUE))), name = "pred_arch1_15.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(newvactual=rep(1,6)))), name = "pred_arch1_16.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(shades=c(95,50)))), name = "pred_arch1_17.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, plot = TRUE, plot.options=list(shades=c(50,95)))), name = "pred_arch1_18.png")

})

##################################################
## 5 FURTHER TESTS
##################################################

getOption("plot")
options(plot=TRUE)
#options(plot=FALSE)
#options(plot=NULL)

##used to yield error:
set.seed(123)
y <- 2+rnorm(20)
mymodel <- arx(y, mc=TRUE)
predict(mymodel, n.ahead=1)

##used to yield error in the plot:
predict(mymodel, ci.levels=c(0.99,0.80, 0.50))

##used to yield error:
set.seed(123)
y <- 25+rnorm(100)
mymodel <- arx(y, mc=TRUE, ar=1:2) ##y has mean 25
pred <- predict(mymodel, n.ahead=20, plot.options=list(keep=50))
pred
##provide some actual values (needn't be same length as n.ahead)
y.actual <- 25+rnorm(12)
preda <- predict(mymodel, n.ahead=20,
                 plot.options=list(newmactual=y.actual))
preda

##used to yield error in the plotting:
set.seed(123)
dgp.n <- 50
y <- rnorm(dgp.n) #or: y <- rt(dgp.n, df=4.1)
mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
y[1:10] <- y[1:10] + 4 #step-shift
arxmod <- arx(y, mc=TRUE, ar=1:2,
              mxreg=cbind(mX, tim(y, which.ones=c(7,16))))
predict(arxmod, n.ahead=1, newmxreg=matrix(0,1,5),
        plot.options=list(start.at.origin=FALSE,
                          line.at.origin=TRUE, fitted=TRUE))

##used to produce graphical error in the plot; since version 0.25
##the following message is returned to the user: "The values of
##'newindex' are not entirely out-of-sample, so no plot produced"

test_that("plot error is reported", {
  
  set.seed(123)
  y <- rnorm(20)
  arxmod <- arx(y, mc=TRUE, ar=1)
  ##generates predictions with user-specified index, but no plot:
  # not sure why expect_message() does not work
  set.seed(123)
  expect_message(predict(arxmod, newindex=19:30))

})


##################################################
## 6 SIMULATIONS VS. ANALYTICAL FORMULA
##################################################

##TO DO: COMPARE THE SIMULATED QUANTILES AND THE
##ANALYTICAL ONES IN THE SPECIAL CASE OF 1-STEP
##AHEAD, WHEN THE MODEL IS AN AR(1) AND THE ERROR
##IS N(0,1).

##true ar(1) parameter:
phi1 <- 0.95

##simulate:
set.seed(123)
y <- arima.sim(list(ar=phi1), 1000)

##large sample (T=1000):
##======================

##estimate ar(1):
mymodel <- arx(y, mc=FALSE, ar=1)

##predictions of the mean:
predict(mymodel, plot=TRUE)
predict(mymodel, n.ahead=24, plot=TRUE)

test_that("Check that bootstrap and innov=rnorm simulations are similar in large sample",{
  skip_on_ci()
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, n.ahead=24, n.sim=10000, plot=TRUE)), name = "pred_bootstr_largen.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, n.ahead=24, n.sim=10000, plot=TRUE,innov=rnorm(10000*24))), name = "pred_innovrnorm_largen.png")
})

##conclusion ("large sample"): there is not much difference
##in the fans produced by the bootstrap and innov=rnorm.

##small sample (T=20):
##====================

##estimate ar(1):
mymodel <- arx(y[1:20], mc=FALSE, ar=1)

##predictions of the mean:
predict(mymodel, n.ahead=24, plot=TRUE)

test_that("Check that bootstrap and innov=rnorm simulations are similar in small sample",{
  skip_on_ci()
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, n.ahead=24, n.sim=10000, plot=TRUE)), name = "pred_bootstr_smalln.png")
  set.seed(123)
  expect_snapshot_file(cran = FALSE, path = save_png(predict(mymodel, n.ahead=24, n.sim=10000, plot=TRUE,innov=rnorm(10000*24))), name = "pred_innovrnorm_smalln.png")
})


##conclusion (small sample): there can be a large difference
##in the fans produced by the bootstrap and innov=rnorm.
