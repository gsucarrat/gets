# Indicator Saturation (S3 Method)
#
# @param lmobject `lm` object 
# @param ... Further arguments passed to `isat()`
# S3 Method to use an `lm` object in the gets function `isat()`
# @return A `isat` object
# @export
#
# @examples
# Simplest use

# General-to-Specific (GETS) Modelling for an `lm` object (S3 Method)
#
# @param lmobject `lm` object 
# @param ... Further arguments passed to `getsm()`
#
# @return A `getsm` object
# @export
#
# @examples
# Simplest use


# Estimate an AR-X model with log-ARCH-X errors (S3 Method for lm)
#
# @param lmobject `lm` object 
# @param ... Further arguments passed to `arx()`
#
# @return A `arx` object
# @export
#
# @examples
# # Simplest use
mtcars

# arx testing -------------------------------------------------------------

test_that("Check coefficients are equal between lm and arx",{

  # A number of models with varying formula arguments
  lm_model1 <- lm(mpg~ cyl  + gear + carb + hp + drat, mtcars)
  lm_model2 <- lm(log(mpg)~ cyl  + as.factor(gear) + carb + hp + I(hp*hp) + drat + I(drat*carb), mtcars)
  lm_model3 <- lm(cyl ~ log(mpg) + as.factor(gear) + log(carb) + hp + I(hp*hp) + lag(drat) + I(drat*carb) + mpg*carb + hp:drat, mtcars)
  lm_model4 <- lm(cyl ~ log(mpg) + as.factor(gear) + log(carb) + hp + I(hp*hp) + lag(drat) + I(drat*carb) + mpg*carb + hp:drat-1, mtcars)
  
  arx_model1 <- arx(lm_model1)
  arx_model2 <- arx(lm_model2)
  arx_model3 <- arx(lm_model3)
  arx_model4 <- arx(lm_model4)
  
  coef_lm_model1 <- coef(lm_model1)
  coef_lm_model2 <- coef(lm_model2)
  coef_lm_model3 <- coef(lm_model3)
  coef_lm_model4 <- coef(lm_model4)
  
  coef_arx_model1 <- coef(arx_model1)
  coef_arx_model2 <- coef(arx_model2)
  coef_arx_model3 <- coef(arx_model3)
  coef_arx_model4 <- coef(arx_model4)
  
  # modify the intercept name
  names(coef_lm_model1)[which(names(coef_lm_model1)=="(Intercept)")] <- "mconst"
  names(coef_lm_model2)[which(names(coef_lm_model2)=="(Intercept)")] <- "mconst"
  names(coef_lm_model3)[which(names(coef_lm_model3)=="(Intercept)")] <- "mconst"
  names(coef_lm_model4)[which(names(coef_lm_model4)=="(Intercept)")] <- "mconst"
  
  
  expect_equal(coef_lm_model1, coef_arx_model1)
  expect_equal(coef_lm_model2, coef_arx_model2)
  expect_equal(coef_lm_model3, coef_arx_model3)
  expect_equal(coef_lm_model4, coef_arx_model4)
  
  
})


test_that("Check that the coefficients are equal between lm and arx but with vcov.type the S.E. are different",{
  
  # A number of models with varying formula arguments
  lm_model1 <- lm(mpg~ cyl  + gear + carb + hp + drat, mtcars)
  lm_model2 <- lm(log(mpg)~ cyl  + as.factor(gear) + carb + hp + I(hp*hp) + drat + I(drat*carb), mtcars)
  lm_model3 <- lm(cyl ~ log(mpg) + as.factor(gear) + log(carb) + hp + I(hp*hp) + lag(drat) + I(drat*carb) + mpg*carb + hp:drat, mtcars)
  lm_model4 <- lm(cyl ~ log(mpg) + as.factor(gear) + log(carb) + hp + I(hp*hp) + lag(drat) + I(drat*carb) + mpg*carb + hp:drat-1, mtcars)
  
  arx_model1 <- arx(lm_model1, vcov.type = "white")
  arx_model2 <- arx(lm_model2, vcov.type = "white")
  arx_model3 <- arx(lm_model3, vcov.type = "newey-west")
  arx_model4 <- arx(lm_model4, vcov.type = "newey-west")
  
  coef_lm_model1 <- coef(lm_model1)
  coef_lm_model2 <- coef(lm_model2)
  coef_lm_model3 <- coef(lm_model3)
  coef_lm_model4 <- coef(lm_model4)
  
  coef_arx_model1 <- coef(arx_model1)
  coef_arx_model2 <- coef(arx_model2)
  coef_arx_model3 <- coef(arx_model3)
  coef_arx_model4 <- coef(arx_model4)
  
  # modify the intercept name
  names(coef_lm_model1)[which(names(coef_lm_model1)=="(Intercept)")] <- "mconst"
  names(coef_lm_model2)[which(names(coef_lm_model2)=="(Intercept)")] <- "mconst"
  names(coef_lm_model3)[which(names(coef_lm_model3)=="(Intercept)")] <- "mconst"
  names(coef_lm_model4)[which(names(coef_lm_model4)=="(Intercept)")] <- "mconst"
  
  
  expect_equal(coef_lm_model1, coef_arx_model1)
  expect_equal(coef_lm_model2, coef_arx_model2)
  expect_equal(coef_lm_model3, coef_arx_model3)
  expect_equal(coef_lm_model4, coef_arx_model4)
  
  expect_false(all(vcov(lm_model1)==vcov(arx_model1)))
  expect_false(all(vcov(lm_model2)==vcov(arx_model2)))
  expect_false(all(vcov(lm_model3)==vcov(arx_model3)))
  expect_false(all(vcov(lm_model4)==vcov(arx_model4)))

})



# gets testing -----------------------------------------------------------

test_that("Test that gets works with lm models",{
  
  # A number of models with varying formula arguments
  lm_model1 <- lm(mpg~ cyl  + gear + carb + hp + drat, mtcars)
  lm_model2 <- lm(log(mpg)~ cyl  + as.factor(gear) + carb + hp + I(hp*hp) + drat + I(drat*carb), mtcars)
  lm_model3 <- lm(cyl ~ log(mpg) + as.factor(gear) + log(carb) + hp + I(hp*hp) + lag(drat) + I(drat*carb) + mpg*carb + hp:drat, mtcars)
  lm_model4 <- lm(cyl ~ log(mpg) + as.factor(gear) + log(carb) + hp + I(hp*hp) + lag(drat) + I(drat*carb) + mpg*carb + hp:drat-1, mtcars)
  
  
  gets_model1 <- gets(lm_model1, print.searchinfo = FALSE)
  gets_model2 <- gets(lm_model2, print.searchinfo = FALSE)
  gets_model3 <- gets(lm_model3, print.searchinfo = FALSE)
  gets_model4 <- gets(lm_model4, print.searchinfo = FALSE)
  
  expect_output(print(gets_model1))
  expect_output(print(gets_model2))
  expect_output(print(gets_model3))
  expect_output(print(gets_model4))
  
  # Check that all variable names in gets are also in the lm object
  coef_lm_model1 <- coef(lm_model1)
  coef_lm_model2 <- coef(lm_model2)
  coef_lm_model3 <- coef(lm_model3)
  coef_lm_model4 <- coef(lm_model4)
  
  # modify the intercept name
  names(coef_lm_model1)[which(names(coef_lm_model1)=="(Intercept)")] <- "mconst"
  names(coef_lm_model2)[which(names(coef_lm_model2)=="(Intercept)")] <- "mconst"
  names(coef_lm_model3)[which(names(coef_lm_model3)=="(Intercept)")] <- "mconst"
  names(coef_lm_model4)[which(names(coef_lm_model4)=="(Intercept)")] <- "mconst"
  
  expect_true(all(names(coef(gets_model1)) %in% names(coef_lm_model1)))
  expect_true(all(names(coef(gets_model2)) %in% names(coef_lm_model2)))
  expect_true(all(names(coef(gets_model3)) %in% names(coef_lm_model3)))
  expect_true(all(names(coef(gets_model4)) %in% names(coef_lm_model4)))
  
  # Check that the number of variables is always smaller than the original model
  expect_true(length(coef(gets_model1)) <= length(coef(lm_model1)))
  expect_true(length(coef(gets_model2)) <= length(coef(lm_model2)))
  expect_true(length(coef(gets_model3)) <= length(coef(lm_model3)))
  expect_true(length(coef(gets_model4)) <= length(coef(lm_model4)))
  
})


# isat testing ------------------------------------------------------------


test_that("Test that isat works with lm models",{
  data("hpdata", package = "gets")
  lm_model1 <- lm(GD ~ FMBASE + FSDJ + LHC + MU, hpdata)
  lm_model2 <- lm(log(GD) ~ log(FMBASE) + FSDJ + LHC + MU, hpdata)
  lm_model3 <- lm(log(GD) ~ log(FMBASE) + FSDJ*LHC + MU, hpdata)
  lm_model4 <- lm(GD ~ FMBASE + FSDJ + LHC + MU + I(MU*MU)-1, hpdata)
  
  # A number of models with varying formula arguments
  # lm_model1 <- lm(mpg~ cyl  + gear + carb + hp + drat, mtcars)
  # lm_model2 <- lm(log(mpg)~ cyl  + as.factor(gear) + carb + hp + I(hp*hp) + drat + I(drat*carb), mtcars)
  # lm_model3 <- lm(cyl ~ log(mpg) + as.factor(gear) + log(carb) + hp + I(hp*hp) + lag(drat) + I(drat*carb) + mpg*carb + hp:drat, mtcars)
  # lm_model4 <- lm(cyl ~ log(mpg) + as.factor(gear) + log(carb) + hp + I(hp*hp) + lag(drat) + I(drat*carb) + mpg*carb + hp:drat-1, mtcars)
  # 
  # 
  isat_model1 <- isat(lm_model1, print.searchinfo = FALSE, t.pval = 0.001)
  isat_model2 <- isat(lm_model2, print.searchinfo = FALSE, t.pval = 0.001)
  isat_model3 <- isat(lm_model3, print.searchinfo = FALSE, t.pval = 0.001)
  isat_model4 <- isat(lm_model4, print.searchinfo = FALSE, t.pval = 0.001, ar = 1) # adding ar term
  
  expect_output(print(isat_model1))
  expect_output(print(isat_model2))
  expect_output(print(isat_model3))
  expect_output(print(isat_model4))
  
  # Check that all variable names in isat are also in the lm object
  coef_lm_model1 <- coef(lm_model1)
  coef_lm_model2 <- coef(lm_model2)
  coef_lm_model3 <- coef(lm_model3)
  coef_lm_model4 <- coef(lm_model4)
  
  # modify the intercept name
  names(coef_lm_model1)[which(names(coef_lm_model1)=="(Intercept)")] <- "mconst"
  names(coef_lm_model2)[which(names(coef_lm_model2)=="(Intercept)")] <- "mconst"
  names(coef_lm_model3)[which(names(coef_lm_model3)=="(Intercept)")] <- "mconst"
  names(coef_lm_model4)[which(names(coef_lm_model4)=="(Intercept)")] <- "mconst"
  
  
  # Check that the number of variables is always equal or larger than the original model
  expect_true(length(coef(isat_model1)) >= length(coef(lm_model1)))
  expect_true(length(coef(isat_model2)) >= length(coef(lm_model2)))
  expect_true(length(coef(isat_model3)) >= length(coef(lm_model3)))
  expect_true(length(coef(isat_model4)) >= length(coef(lm_model4)))
  
})



test_that("Test that isat works with arx models",{
  data("hpdata", package = "gets")
  arx_model1 <- arx(y = hpdata$GD, mxreg = hpdata[,c("FMBASE","FSDJ","LHC","MU")])
  arx_model2 <- arx(y = log(hpdata$GD), mxreg = hpdata[,c("FMBASE","FSDJ","LHC","MU")])
  arx_model3 <- arx(y = hpdata$GD, mxreg = hpdata[,c("FMBASE","FSDJ","LHC","MU")])
  
  # A number of models with varying formula arguments
  # lm_model1 <- lm(mpg~ cyl  + gear + carb + hp + drat, mtcars)
  # lm_model2 <- lm(log(mpg)~ cyl  + as.factor(gear) + carb + hp + I(hp*hp) + drat + I(drat*carb), mtcars)
  # lm_model3 <- lm(cyl ~ log(mpg) + as.factor(gear) + log(carb) + hp + I(hp*hp) + lag(drat) + I(drat*carb) + mpg*carb + hp:drat, mtcars)
  # lm_model4 <- lm(cyl ~ log(mpg) + as.factor(gear) + log(carb) + hp + I(hp*hp) + lag(drat) + I(drat*carb) + mpg*carb + hp:drat-1, mtcars)
  # 
  # 
  isat_model1 <- isat(arx_model1, print.searchinfo = FALSE, t.pval = 0.001)
  isat_model2 <- isat(arx_model2, print.searchinfo = FALSE, t.pval = 0.001)
  isat_model3 <- isat(arx_model3, print.searchinfo = FALSE, t.pval = 0.001)
  
  expect_output(print(isat_model1))
  expect_output(print(isat_model2))
  expect_output(print(isat_model3))
  
  # Check that all variable names in isat are also in the lm object
  coef_arx_model1 <- coef(arx_model1)
  coef_arx_model2 <- coef(arx_model2)
  coef_arx_model3 <- coef(arx_model3)

  
  # Check that the number of variables is always equal or larger than the original model
  expect_true(length(coef(isat_model1)) >= length(coef(arx_model1)))
  expect_true(length(coef(isat_model2)) >= length(coef(arx_model2)))
  expect_true(length(coef(isat_model3)) >= length(coef(arx_model3)))
  
})



# Special case testing ----------------------------------------------------


test_that("Check that the error message displays when weights are used in lm",{
  lm_model5 <- lm(mpg~ cyl + hp, weights = rnorm(nrow(mtcars),mean = 100), mtcars)
  expect_error(arx(lm_model5))
  
})

# 
# 
# For the examples
# 
# 
# # More complicated model
# model2 <- lm(log(mpg)~ cyl  + as.factor(gear) + carb + hp + I(hp*hp) + drat + I(drat*carb), mtcars)
# arx(model2)
# 
# model3 <- lm(log(mpg)~ cyl  + as.factor(gear) + log(carb) + hp + I(hp*hp) + lag(drat) + I(drat*carb) + cyl*carb + hp:drat, mtcars)
# summary(model3)
# arx(model3)
# 
# getsm(model3)
# 
# isat(model3)
# 
# 
# 
# 
# 
# 
# model1 <- lm(mpg~ cyl  + gear + carb + hp + drat, mtcars)
# isat.model1 <- isat(model1)
# 
# varnames_lm <- names(model1$coefficients)
# varnames_lm <- setdiff(varnames_lm,"(Intercept)")
# names(coef(isat.model1))
# 
# # With further options
# isat(model1, t.pval = 0.01, iis = TRUE, vcov.type = "white")
# 
# # More complicated model
# model2 <- lm(log(mpg)~ cyl  + as.factor(gear) + carb + hp + I(hp*hp) + drat + I(drat*carb), mtcars)
# isat(model2)
# 
# 
# 
# model1 <- lm(mpg~ cyl  + gear + carb + hp + drat, mtcars)
# getsm(model1)
# 
# # With further options
# getsm(model1, vcov.type = "white")
# 
# # More complicated model
# model2 <- lm(log(mpg)~ cyl  + as.factor(gear) + carb + hp + I(hp*hp) + drat + I(drat*carb), mtcars)
# getsm(model2)
# 



