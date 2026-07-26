test_that("predict.isat",{
  
  ##simulate data:
  set.seed(123)
  dgp.n <- 50
  y <- rnorm(dgp.n) #or: y <- rt(dgp.n, df=4.1)
  mX <- matrix(rnorm(dgp.n*3),dgp.n,3)
  y[1:10] <- y[1:10] + 4 #step-shift
  # plot(as.zoo(y), col="blue") #plot
  
  isatmod <- isat(y, ar=1:2, mxreg=mX, iis=TRUE, sis=TRUE, tis=TRUE, print.searchinfo = FALSE)
  expect_true(all(c("tis7", "tis16") %in% isatmod$ISnames))
  
  mX.pred <- cbind(matrix(rnorm(10*3),10,3))
  mx.pred.tis <- as.matrix(tail(tim(x = dgp.n + 10, which.ones = c(7,16)), 10))
  mX.pred.final <- cbind(mX.pred, mx.pred.tis)
  
  expect_silent(predict(isatmod, newmxreg = mX.pred.final, n.ahead = nrow(mX.pred), quiet = TRUE))
  
  expect_message(predict(isatmod, newmxreg = mX.pred.final, n.ahead = nrow(mX.pred), quiet = FALSE), regexp = "You have provided new data for the following IIS, SIS, or TIS indicators")
  
  # only provide one tim
  expect_message(a <- predict(isatmod, newmxreg = mX.pred.final[,c(1,2,3,4)], n.ahead = nrow(mX.pred), quiet = FALSE), regexp = "You have provided new data for the following IIS, SIS, or TIS indicators")
  # check whether the order in mxreg matters
  expect_message(b <- predict(isatmod, newmxreg = mX.pred.final[,c(1,2,3,5)], n.ahead = nrow(mX.pred), quiet = FALSE), regexp = "You have provided new data for the following IIS, SIS, or TIS indicators")
  expect_equal(a, b)
  
  expect_message(pred_fully_spec <- predict(isatmod, newmxreg = mX.pred.final, n.ahead = nrow(mX.pred)), regexp = "You have provided new data for the following IIS, SIS, or TIS indicators")
  expect_silent(pred_onlymx_spec <- predict(isatmod, newmxreg = mX.pred, n.ahead = nrow(mX.pred)))
  # check that those are all equal
  expect_equal(pred_fully_spec, pred_onlymx_spec)

  
  pred_second_tis <- predict(isatmod,newmxreg = mX.pred.final[, c(1:3, 5), drop = FALSE],n.ahead = nrow(mX.pred),quiet = TRUE)
  expect_equal(pred_fully_spec, pred_second_tis)
  
  
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
  
  

  
  
  # check ordering and partial provision of supplied indicators
  indicators <- mx.pred.tis
  n.ahead <- nrow(mX.pred)
  
  ## All indicators supplied in reverse order
  pred_reverse <- predict(
    isatmod,
    newmxreg = cbind(
      mX.pred,
      indicators[, rev(seq_len(ncol(indicators))), drop = FALSE]
    ),
    n.ahead = n.ahead,
    quiet = TRUE
  )
  
  expect_equal(pred_reverse, pred_fully_spec)
  
  ## Only the first indicator supplied
  pred_first <- predict(
    isatmod,
    newmxreg = cbind(
      mX.pred,
      indicators[, 1L, drop = FALSE]
    ),
    n.ahead = n.ahead,
    quiet = TRUE
  )
  
  expect_equal(pred_first, pred_fully_spec)
  
  ## Only the second indicator supplied
  pred_second <- predict(
    isatmod,
    newmxreg = cbind(
      mX.pred,
      indicators[, 2L, drop = FALSE]
    ),
    n.ahead = n.ahead,
    quiet = TRUE
  )
  
  expect_equal(pred_second, pred_fully_spec)
  
  ## Indicators supplied between the ordinary regressors
  pred_interspersed <- predict(
    isatmod,
    newmxreg = cbind(
      mX.pred[, 1L, drop = FALSE],
      indicators[, 2L, drop = FALSE],
      mX.pred[, 2:3, drop = FALSE],
      indicators[, 1L, drop = FALSE]
    ),
    n.ahead = n.ahead,
    quiet = TRUE
  )
  
  expect_equal(pred_interspersed, pred_fully_spec)
  
})



test_that("predict.isat with newdata",{
  
  ##step indicator saturation:
  set.seed(123)
  y <- rnorm(30)
  isatmod <- isat(y, print.searchinfo = FALSE)
  
  ##generate forecasts of the simplified (specific) model:
  pred.isatmod <- predict(isatmod, newmxreg=matrix(1,12,1), plot=TRUE)
  expect_equal(round(pred.isatmod,5),structure(c(-0.0471, -0.0471, -0.0471, -0.0471, -0.0471, -0.0471, 
                                                 -0.0471, -0.0471, -0.0471, -0.0471, -0.0471, -0.0471), index = c(31, 
                                                                                                                  32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42), class = "zoo"))
  
  
  # test that predict works without covariates (other than mc and ar)
  set.seed(123)
  y2 <- rnorm(50)
  y2[1:20] <- y2[1:20] - 5
  
  mod2 <- isat(y2, iis = FALSE, sis = TRUE, include.gum = FALSE,
               print.searchinfo = FALSE)
  
  expect_gt(length(mod2$ISnames), 0)
  expect_silent(predict(mod2, n.ahead = 4))
  
  
  
  set.seed(123)
  y2 <- rnorm(50)
  y2[1:20] <- y2[1:20] - 5
  
  mod2 <- isat(y2, iis = FALSE, sis = TRUE, include.gum = FALSE,
               print.searchinfo = FALSE)
  
  expect_true(length(mod2$ISnames) > 0)
  expect_silent(predict(mod2, n.ahead = 4))
  
  mod2.predict.old <- predict(mod2, n.ahead = 4, newmxreg = data.frame(sis21 = c(1,1,1,1)), quiet = TRUE)
  expect_equal(predict(mod2, n.ahead = 4), mod2.predict.old)
  
  
  
})
