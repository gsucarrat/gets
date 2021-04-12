
isat.lm <- function(lmobject, ...){
  data <- model.frame(lmobject)
  y <- data.frame(y = data[,1])
  names(y) <- names(data)[1] 
  
  
  x_vars <- data[,2:ncol(data)]
  mxreg <- model.matrix(lmobject$terms, data)
  # Deal with intercept
  mc <- ifelse(attr(lmobject$terms, "intercept")==1, TRUE, FALSE)
  mxreg <- mxreg[,!colnames(mxreg) == "(Intercept)"] # remove the intercept
  
  out <- gets::isat(y = y, mxreg = mxreg, mc = mc, ...)
  return(out)
}

getsm.lm <- function(lmobject, ...){
  arx_object <- arx.lm(lmobject)
  getsm(arx_object, ...)
}


arx.lm <- function(lmobject, ...){
  data <- model.frame(lmobject)
  y <- data.frame(y = data[,1])
  names(y) <- names(data)[1] 
  
  
  x_vars <- data[,2:ncol(data)]
  mxreg <- model.matrix(lmobject$terms, data)
  
  # Deal with intercept
  mc <- ifelse(attr(lmobject$terms, "intercept")==1, TRUE, FALSE)
  mxreg <- mxreg[,!colnames(mxreg) == "(Intercept)"] # remove the intercept
  out <- arx(y = y, mxreg = mxreg, mc = mc, ...)
  return(out)
}
