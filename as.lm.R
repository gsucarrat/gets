###=========================
###test that as.lm() works:
#
#set.seed(123)
#y <- rnorm(30)
#arxmod <- arx(y, mc=TRUE, ar=1:3)
#as.lm(arxmod)
#
#getsmod <- getsm(arxmod, keep=1)
#as.lm(getsmod)
#
#isatmod <- isat(y)
#as.lm(isatmod)

##==================================================
##convert to model of class 'lm':
as.lm <- function(object)
{

  ##what kind of class?:
  objectClass <- class(object)
  classOK <-
    ifelse( objectClass %in% c("arx","gets","isat"), TRUE, FALSE)

  ##class not OK:
  if(!classOK){
    stop("'object' must be of class 'arx', 'gets' or 'isat'")
  }

  ##class OK:
  if(classOK){
    y <- object$aux$y
    x <- object$aux$mX
    colnames(x) <- object$aux$mXnames
    result <- lm(y ~ x - 1)
  }
  
  ##return result:
  return(result)
  
}