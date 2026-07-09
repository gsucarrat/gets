test_that("predict.isat",{
  
  ##simulate data:
  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n) #or: y <- rt(dgp.n, df=4.1)
  mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
  y[1:10] <- y[1:10] + 4 #step-shift
  # plot(as.zoo(y), col="blue") #plot
  
  isatmod <- isat(y, ar=1:2, mxreg=mX, iis=TRUE, sis=TRUE, tis=TRUE, print.searchinfo = FALSE)
  
  mX.pred <- cbind(matrix(rnorm(10*3),10,3))
  mx.pred.tis <- as.matrix(tail(tim(x = dgp.n + 10, which.ones = c(7,16)), 10))
  mX.pred.final <- cbind(mX.pred, mx.pred.tis)
  
  expect_silent(predict(isatmod, newmxreg = mX.pred.final, n.ahead = nrow(mX.pred), quiet = TRUE))
  expect_message(predict(isatmod, newmxreg = mX.pred.final, n.ahead = nrow(mX.pred), quiet = FALSE), regexp = "You have provided new data for the following IIS, SIS, or TIS indicators")
  # only provide one tim
  expect_message(predict(isatmod, newmxreg = mX.pred.final[,c(1,2,3,4)], n.ahead = nrow(mX.pred), quiet = FALSE), regexp = "You have provided new data for the following IIS, SIS, or TIS indicators")
  expect_message(pred_fully_spec <- predict(isatmod, newmxreg = mX.pred.final, n.ahead = nrow(mX.pred)), regexp = "You have provided new data for the following IIS, SIS, or TIS indicators")
  expect_silent(pred_onlymx_spec <- predict(isatmod, newmxreg = mX.pred, n.ahead = nrow(mX.pred)))
  # check that those are all equal
  expect_equal(pred_fully_spec, pred_onlymx_spec)
  
  
  # names are different but numbers should be the same
  uis_tis <- tim(dgp.n-2)
  colnames(uis_tis) <- gsub("tis","uis",colnames(uis_tis))
  isatmod.uis <- isat(y, ar=1:2, mxreg=mX, iis=TRUE, sis=TRUE, tis=FALSE, uis = list(uis_tis), print.searchinfo = FALSE)
  
  expect_error(predict(isatmod.uis, newmxreg = mX.pred, n.ahead = 10), "Your object isat contains UIS indicators but not all UIS indicators are in newmxreg")
  mX.pred.final.noname <- mX.pred.final
  colnames(mX.pred.final.noname) <- NULL
  expect_silent(predict(isatmod.uis, newmxreg = mX.pred.final.noname, n.ahead = 10))
  
  
  # now check all
  
  set.seed(123)
  xregs <- matrix(rnorm(4*70), 70, 4)
  y <- arima.sim(list(ar=0.4), 70) + 100 + (xregs[,2]*0.3) + 100
  y[33] <- y[33] + 30
  y[1:20] <- y[1:20] - 5
  y[50:70] <- y[50:70] + 1:21
  
  suppressMessages(isat_result <- isat(y, mc = TRUE, mxreg = xregs, iis = TRUE, sis = TRUE, tis = TRUE, print.searchinfo = FALSE, plot = FALSE))
  
  expect_silent(full_isat_result <- predict(isat_result, newmxreg = matrix(rnorm(4*12), 12, 4)))
  expect_message(predict(isat_result, newmxreg = cbind(matrix(rnorm(4*12), 12, 4),matrix(0,nrow = 12,dimnames = list(NULL, "iis33")))), regexp = "You have provided new data for the following IIS, SIS, or TIS indicators")
  expect_silent(predict(isat_result, newmxreg = cbind(matrix(rnorm(4*12), 12, 4),matrix(0,nrow = 12,dimnames = list(NULL, "iis33"))), quiet = TRUE))
  expect_identical(round(full_isat_result,5),structure(c(222.68878, 220.72002, 222.37005, 224.08472, 224.17809, 
                                                         224.65727, 225.60892, 226.97303, 227.75937, 228.69805, 229.93318, 
                                                         230.49231), index = c(71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 
                                                                               81, 82), class = "zoo"))
  
})
