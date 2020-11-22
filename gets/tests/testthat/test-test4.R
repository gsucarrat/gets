


test_that("##user-defined Shapiro-Wilks test for normality in the residuals",{
  
  diagnostics <- gets::diagnostics
  
  #environment(diagnostics)
  #environment(SWtest) <- environment(diagnostics)
  
  set.seed(123)
  vY <- rnorm(20)
  mX <- matrix(rnorm(3*20), 20, 3)
  
  x <- ols(vY, mX, method=3)
  
  SWtest <- function(x, ...){
    tmp <- shapiro.test(x$residuals)
    result <- c(tmp$statistic, NA, tmp$p.value)
    return(result)
  }
  
  expect_true(diagnostics(x, user.fun=list(name="SWtest", pval=0.025,envir = environment(SWtest)),verbose=FALSE))
})




# 
# 
# expect_true(is.matrix(diagnostics(x, user.fun=list(name="SWtest", pval=0.025))) & 
#               nrow(diagnostics(x, user.fun=list(name="SWtest", pval=0.025)))==3 & 
#               row.names(diagnostics(x, user.fun=list(name="SWtest", pval=0.025)))[3]=="G-test")
# 
# 
# 
# 
# 
# 
# expect_true(diagnostics(x, user.fun=list(name="SWtest", pval=0.025), verbose=FALSE))
# expect_false(diagnostics(x, user.fun=list(name="SWtest", pval=0.85), verbose=FALSE))
# 
# 
# 
# 
# 
# 




# 
# 
# df <- data.frame(x=rnorm(20),y=rnorm(20))
# 
# larger_function <- function(x,estimator = NULL){
#   #browser()
#   result <- do.call(what = estimator$name,args = list(x = x$x,y = x$y))
#   return(result)
#   
#   
# }
# 
# simple_function <- function(x,y){
#   return(x + y)
# }
# 
# larger_function(df,estimator = list(name="simple_function"))
# 
# 
# test_that("Testing Function",{
#   expect_true(length(larger_function(df,estimator = list(name="simple_function")))==20)
# })
