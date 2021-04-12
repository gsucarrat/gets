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
model1 <- lm(mpg~ cyl  + gear + carb + hp + drat, mtcars)
isat(model1)

# With further options
isat(model1, t.pval = 0.01, iis = TRUE, vcov.type = "white")

# More complicated model
model2 <- lm(log(mpg)~ cyl  + as.factor(gear) + carb + hp + I(hp*hp) + drat + I(drat*carb), mtcars)
isat(model2)


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
model1 <- lm(mpg~ cyl  + gear + carb + hp + drat, mtcars)
getsm(model1)

# With further options
getsm(model1, vcov.type = "white")

# More complicated model
model2 <- lm(log(mpg)~ cyl  + as.factor(gear) + carb + hp + I(hp*hp) + drat + I(drat*carb), mtcars)
getsm(model2)

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
model1 <- lm(mpg~ cyl  + gear + carb + hp + drat, mtcars)
arx(model1)

# With further options
arx(model1, t.pval = 0.01, vcov.type = "white")

# More complicated model
model2 <- lm(log(mpg)~ cyl  + as.factor(gear) + carb + hp + I(hp*hp) + drat + I(drat*carb), mtcars)
arx(model2)


