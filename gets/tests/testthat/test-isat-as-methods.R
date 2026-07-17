test_that("test that as.lm works",{
  
  set.seed(123)
  xregs <- matrix(rnorm(4*70), 70, 4)
  y <- arima.sim(list(ar=0.4), 70) + (xregs[,2]*0.3)
  mx_result <- arx(y, ar=1:2, mxreg=xregs)
  lm_result <- as.lm(mx_result)
  
  expect_identical(coef(mx_result),coef(lm_result))
  
})

set.seed(123)
xregs <- matrix(rnorm(4*70), 70, 4)
y <- arima.sim(list(ar=0.4), 70) + 100 #+ (xregs[,2]*0.3) + 100
y[33] <- y[33] + 30
y[1:20] <- y[1:20] - 5
y[50:70] <- y[50:70] + 1:21

suppressMessages(isat_result <- isat(y, mc = TRUE, mxreg = xregs, iis = TRUE, sis = TRUE, tis = TRUE, print.searchinfo = FALSE))

test_that("test that as.arx.isat works",{
  
  arx_result <- as.arx(isat_result)
  
  # calculate a manual arx
  indic_matrix <- data.frame(iis33 = iim(x = 70, which.ones = 33),
                             sis22 = sim(x = 70, which.ones = 22),
                             tis49 = tim(x = 70, which.ones = 49))
  xregs_manual <- cbind(xregs, indic_matrix)
  colnames(xregs_manual) <- c(paste0("mxreg", 1:4), "iis33", "sis22", "tis49")
  manual_arx_result <- arx(y, mc = TRUE, mxreg = xregs_manual)
  
  expect_identical(coef(isat_result), coef(arx_result))
  
  manual_arx_result$call <- NULL
  arx_result$call <- NULL
  manual_arx_result$date <- NULL
  arx_result$date <- NULL
  expect_identical(manual_arx_result, arx_result, ignore_attr = TRUE)
  
  data(Nile)
  isat_object <- isat(Nile, sis=TRUE, iis=FALSE, plot=FALSE, t.pval=0.005, print.searchinfo = FALSE)
  arx_obj <- as.arx(isat_object)
  expect_identical(coef(isat_result), coef(arx_result))  
  
  isat_object <- isat(Nile, ar = 1:2, sis=TRUE, iis=FALSE, plot=FALSE, t.pval=0.005, print.searchinfo = FALSE)
  arx_obj <- as.arx(isat_object)
  expect_identical(coef(isat_result), coef(arx_result)) 
  
})


test_that("test that as.isat works",{
  
  # ##Simulate from an AR(1):
  # set.seed(123)
  # xregs <- matrix(rnorm(2*80), 80, 2)
  # y <- arima.sim(list(ar=0.7), 80) + (xregs[,2]*0.5)
  # gum01 <- arx(y, mc=TRUE, ar=1:2, mxreg=xregs)
  # meanmod01 <- getsm(gum01, print.searchinfo = FALSE)
  
  
  set.seed(123)
  xregs <- matrix(rnorm(4*70), 70, 4)
  y <- arima.sim(list(ar=0.4), 70) + 100 + (xregs[,2]*0.3) + 100
  y[33] <- y[33] + 30
  y[1:20] <- y[1:20] - 5
  y[50:70] <- y[50:70] + 1:21
  
  indic_matrix <- data.frame(iis33 = iim(x = 70, which.ones = 33),
                             sis22 = sim(x = 70, which.ones = 22),
                             tis49 = tim(x = 70, which.ones = 49),
                             
                             # adding another irrelevant indicator to test that as.isat() can handle this
                             iis60 = iim(x = 70, which.ones = 60)
  )
  xregs_manual <- cbind(xregs, indic_matrix)
  colnames(xregs_manual) <- c(paste0("mxreg", 1:4), "iis33", "sis22", "tis49", "iis60")
  manual_arx_result <- arx(y, mc = TRUE, mxreg = xregs_manual)  
  
  manual_arx_result.isat <- as.isat(manual_arx_result)
  expect_s3_class(manual_arx_result.isat, "isat")
  # because we added the irrelevant indicator, expect the note element
  expect_output(print(manual_arx_result.isat), "Result includes indicators that exceed the target p-value ")
  
  # run a manual isat with the same indicators and compare the results
  # we now expect differences
  suppressMessages(isat_result.full <- isat(y, mc = TRUE, mxreg = xregs_manual[,c("mxreg1","mxreg2", "mxreg3", "mxreg4","iis60")], iis = TRUE, sis = TRUE, tis = TRUE, print.searchinfo = FALSE))
  expect_true(class(isat_result.full) == "isat")
  expect_true(!identical(isat_result.full$ISnames,manual_arx_result.isat$ISnames))
  
  # now lets create a real one
  xregs_manual <- cbind(xregs, indic_matrix[,c("iis33", "sis22", "tis49")])
  colnames(xregs_manual) <- c(paste0("mxreg", 1:4), "iis33", "sis22", "tis49")
  manual_arx_result_2 <- arx(y, mc = TRUE, mxreg = xregs_manual)  
  manual_arx_result.correct.isat <- as.isat(manual_arx_result_2)
  
  suppressMessages(isat_result.full.correct <- isat(y, mc = TRUE, mxreg = xregs_manual[,c("mxreg1","mxreg2", "mxreg3", "mxreg4")], iis = TRUE, sis = TRUE, tis = TRUE, print.searchinfo = FALSE))
  
  # rule out the irrelevant differences
  isat_result.full.correct$call <- NULL
  manual_arx_result.correct.isat$call <- NULL
  isat_result.full.correct$date <- NULL
  manual_arx_result.correct.isat$date <- NULL
  isat_result.full.correct$time.started <- NULL
  manual_arx_result.correct.isat$time.started <- NULL
  isat_result.full.correct$time.finished <- NULL
  manual_arx_result.correct.isat$time.finished <- NULL
  
  
  isat_result.full.correct$no.of.estimations <- NULL
  manual_arx_result.correct.isat$no.of.estimations <- NULL
  isat_result.full.correct$no.of.getsFun.calls <- NULL
  manual_arx_result.correct.isat$no.of.getsFun.calls <- NULL
  
  # paths/terminal based differences
  isat_result.full.correct$paths <- NULL
  manual_arx_result.correct.isat$paths <- NULL
  isat_result.full.correct$terminals <- NULL
  manual_arx_result.correct.isat$terminals <- NULL
  isat_result.full.correct$terminals.results <- NULL
  manual_arx_result.correct.isat$terminals.results <- NULL
  expect_true(is.list(manual_arx_result.correct.isat$ISfinalmodels))
  expect_true(is.list(isat_result.full.correct$ISfinalmodels))
  
  isat_result.full.correct$ISfinalmodels <- NULL
  manual_arx_result.correct.isat$ISfinalmodels <- NULL
  
  isat_result.full.correct$best.terminal <- NULL
  manual_arx_result.correct.isat$best.terminal <- NULL
  
  isat_result.full.correct$specific.spec <- NULL
  manual_arx_result.correct.isat$specific.spec <- NULL
  
  expect_identical(isat_result.full.correct, manual_arx_result.correct.isat)
  
  
  
  # Now let's include the gets element
  gets_result <- getsm(manual_arx_result,ar.LjungB = NULL, arch.LjungB = NULL, print.searchinfo = FALSE, keep = 5)
  gets_result.isat <- as.isat(gets_result)
  
  # further isat results
  suppressMessages(isat_result.sele <- isat(y, mc = TRUE, mxreg = xregs_manual[,c("mxreg1","mxreg2","mxreg4"), drop = FALSE], iis = TRUE, sis = TRUE, tis = TRUE, print.searchinfo = FALSE))
  
  # rule out the irrelevant differences
  isat_result.sele$call <- NULL
  gets_result.isat$call <- NULL
  isat_result.sele$date <- NULL
  gets_result.isat$date <- NULL
  isat_result.sele$time.started <- NULL
  gets_result.isat$time.started <- NULL
  isat_result.sele$time.finished <- NULL
  gets_result.isat$time.finished <- NULL
  
  
  isat_result.sele$no.of.estimations <- NULL
  gets_result.isat$no.of.estimations <- NULL
  isat_result.sele$no.of.getsFun.calls <- NULL
  gets_result.isat$no.of.getsFun.calls <- NULL
  
  # paths/terminal based differences
  isat_result.sele$paths <- NULL
  gets_result.isat$paths <- NULL
  isat_result.sele$terminals <- NULL
  gets_result.isat$terminals <- NULL
  isat_result.sele$terminals.results <- NULL
  gets_result.isat$terminals.results <- NULL
  expect_true(is.list(gets_result.isat$ISfinalmodels))
  expect_true(is.list(isat_result.sele$ISfinalmodels))
  
  isat_result.sele$ISfinalmodels <- NULL
  gets_result.isat$ISfinalmodels <- NULL
  
  isat_result.sele$best.terminal <- NULL
  gets_result.isat$best.terminal <- NULL
  
  isat_result.sele$specific.spec <- NULL
  gets_result.isat$specific.spec <- NULL
  
  expect_equal(isat_result.sele, gets_result.isat)
  
  # check the arguments 
  expect_error(as.isat(manual_arx_result, iis = TRUE, sis = TRUE, tis = TRUE, print.searchinfo = FALSE), 
               "Please use 'indicator_regex' argument to specify indicator patterns instead of 'iis', 'sis', 'tis', or 'uis'")
  expect_error(as.isat(manual_arx_result, indicator_regex = list(iis = "^iis")), regexp = "must be a list with named elements 'iis', 'sis', 'tis', and 'uis'") 
  expect_error(as.isat(manual_arx_result, indicator_regex = list("^iis")), regexp = "must be a list with named elements 'iis', 'sis', 'tis', and 'uis'") 
  
  
})

