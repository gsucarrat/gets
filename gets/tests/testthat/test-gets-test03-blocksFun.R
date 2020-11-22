##################################################
## Test file for gets package. First created
## 22 July 2020, Zaragoza.
##
## 1 INITIATE
## 2 TEST MAIN getsFun ARGUMENTS
## 3 TEST SPECIAL CASES
## 7 SIMULATIONS (FOR THE FUTURE)
##
##################################################

##################################################
##1 INITIATE
##################################################

# not needed. Simply click "Run tests" in RStudio in the top right or devtools::test_file("test-filename.R")
# to test this manually (line by line), first load the gets package using:
# devtools::load_all()

##################################################
## 2 TEST MAIN blocksFun ARGUMENTS
##################################################

##dgp1:
set.seed(123)
vY <- rnorm(20)
mX <- matrix(rnorm(length(vY)*10), length(vY), 10)
vY <- mX[,10] + 0.1*rnorm(NROW(mX))

##check if each argument works:

# Visual Inspection
blocksFun(vY, mX)$specific.spec
blocksFun(vY, mX,untransformed.residuals=NA)$specific.spec #should have no effect

test_that("argument untransformed.residuals has no effect", {
  expect_message(blocksFun(vY, mX)$specific.spec)
  expect_message(blocksFun(vY, mX,untransformed.residuals=NA)$specific.spec)
  m1 <- capture_messages(blocksFun(vY, mX)$specific.spec)
  m2 <- capture_messages(blocksFun(vY, mX,untransformed.residuals=NA)$specific.spec)
  expect_identical(m1, m2)
})

mXX <- list(); mXX[[1]] <- mX[,1:5]; mXX[[2]] <- mX[,6:10]
myblocks <- list(); myblocks[[1]] <- list(1:3,4:6,7:10)
blocksFun(vY, mXX)$specific.spec
blocksFun(vY, mX, blocks=myblocks)$blocks
blocksFun(vY, mX, no.of.blocks=4)$blocks
blocksFun(vY, mX, max.block.size=2)$blocks
blocksFun(vY, mX, ratio.threshold=0.4)$blocks

test_that("block variation works", {
  expect_message(blocksFun(vY, mXX)$specific.spec)
  
  expect_message(blocksFun(vY, mX, blocks=myblocks)$blocks)
  a <- blocksFun(vY, mX, blocks=myblocks)$blocks
  expect_type(a, "list")
  expect_length(a, 1)
  expect_length(a$mX, 3)
  expect_equal(a$mX[[1]], 1:3)
  expect_equal(a$mX[[2]], 4:6)
  expect_equal(a$mX[[3]], 7:10)
  
  expect_message(blocksFun(vY, mX, no.of.blocks=4)$blocks)
  b <- blocksFun(vY, mX, no.of.blocks=4)$blocks
  expect_type(b, "list")
  expect_length(b, 1)
  expect_length(b$mX, 4)
  expect_equal(b$mX[[1]], 1:3)
  expect_equal(b$mX[[2]], 4:6)
  expect_equal(b$mX[[3]], 7:9)
  expect_equal(b$mX[[4]], 10)
  
  expect_message(blocksFun(vY, mX, max.block.size=2)$blocks)
  c <- blocksFun(vY, mX, max.block.size=2)$blocks
  expect_type(c, "list")
  expect_length(c, 1)
  expect_length(c$mX, 5)
  expect_equal(c$mX[[1]], 1:2)
  expect_equal(c$mX[[2]], 3:4)
  expect_equal(c$mX[[3]], 5:6)
  expect_equal(c$mX[[4]], 7:8)
  expect_equal(c$mX[[5]], 9:10)
  
  expect_message(blocksFun(vY, mX, ratio.threshold=0.4)$blocks)
  d <- blocksFun(vY, mX, ratio.threshold=0.4)$blocks
  expect_type(d, "list")
  expect_length(d, 1)
  expect_length(d$mX, 3)
  max_ratio <- max(length(d$mX[[1]])/NCOL(mX), length(d$mX[[2]])/NCOL(mX), length(d$mX[[3]])/NCOL(mX))
  expect_lte(max_ratio, 0.4)
  expect_equal(d$mX[[1]], 1:4)
  expect_equal(d$mX[[2]], 5:8)
  expect_equal(d$mX[[3]], 9:10)
})

blocksFun(vY, mX, gets.of.union=FALSE)$specific.spec
blocksFun(vY, mX, force.invertibility=TRUE)$specific.spec
blocksFun(vY, mX,user.estimator=list(name="ols", method=4))$specific.spec
blocksFun(vY, mX,user.estimator=list(name="ols", method=5))$specific.spec
blocksFun(vY, mX, t.pval=0.5)$specific.spec
blocksFun(vY, mX, wald.pval=0.5)$specific.spec
blocksFun(vY, mX, do.pet=FALSE)$specific.spec
blocksFun(vY, mX, ar.LjungB=c(1,0.99))$specific.spec
blocksFun(vY, mX, arch.LjungB=c(1,0.99))$specific.spec
blocksFun(vY, mX, normality.JarqueB=0.99)$specific.spec
blocksFun(vY, mX, keep=1)$specific.spec
blocksFun(vY, mX, include.gum=TRUE)$specific.spec
blocksFun(vY, mX, include.1cut=TRUE)$specific.spec
blocksFun(vY, mX, include.empty=TRUE)$specific.spec
blocksFun(vY, mX,include.gum=TRUE, include.1cut=TRUE, include.empty=TRUE)$specific.spec
blocksFun(vY, mX, max.paths=1)$specific.spec
blocksFun(vY, mX, parallel.options=2)$specific.spec
# blocksFun(vY, mX, parallel.options=5)$specific.spec #should return error
blocksFun(vY, mX, turbo=TRUE)$specific.spec
blocksFun(vY, mX, force.invertibility=FALSE)$specific.spec
blocksFun(vY, mX, tol=1)$specific.spec
blocksFun(vY, mX, LAPACK=TRUE)$specific.spec
# blocksFun(vY, mX, max.regs=1)$specific.spec #should give error
blocksFun(vY, mX, print.searchinfo=FALSE)$specific.spec
blocksFun(vY, mX, alarm=TRUE)$specific.spec

# Unit testing

test_that("TEST MAIN blocksFun ARGUMENTS",{
  expect_message(blocksFun(vY, mX)$specific.spec)
  expect_message(blocksFun(vY, mX,untransformed.residuals=NA)$specific.spec) #should have no effect
  
  mXX <- list()
  mXX[[1]] <- mX[,1:5]
  mXX[[2]] <- mX[,6:10]
  expect_message(blocksFun(vY, mXX)$specific.spec)
  myblocks <- list()
  myblocks[[1]] <- list(1:3,4:6,7:10)
  expect_message(blocksFun(vY, mX, blocks=myblocks)$blocks)
  expect_message(blocksFun(vY, mX, no.of.blocks=4)$blocks)
  expect_message(blocksFun(vY, mX, max.block.size=2)$blocks)
  expect_message(blocksFun(vY, mX, ratio.threshold=0.4)$blocks)
  expect_message(blocksFun(vY, mX, gets.of.union=FALSE)$specific.spec)
  expect_message(blocksFun(vY, mX, force.invertibility=TRUE)$specific.spec)
  expect_message(blocksFun(vY, mX,user.estimator=list(name="ols", method=4))$specific.spec)
  expect_message(blocksFun(vY, mX,user.estimator=list(name="ols", method=5))$specific.spec)
  expect_message(blocksFun(vY, mX, t.pval=0.5)$specific.spec)
  expect_message(blocksFun(vY, mX, wald.pval=0.5)$specific.spec)
  expect_message(blocksFun(vY, mX, do.pet=FALSE)$specific.spec)
  expect_message(blocksFun(vY, mX, ar.LjungB=c(1,0.99))$specific.spec)
  expect_message(blocksFun(vY, mX, arch.LjungB=c(1,0.99))$specific.spec)
  expect_message(blocksFun(vY, mX, normality.JarqueB=0.99)$specific.spec)
  expect_message(blocksFun(vY, mX, keep=1)$specific.spec)
  expect_message(blocksFun(vY, mX, include.gum=TRUE)$specific.spec)
  expect_message(blocksFun(vY, mX, include.1cut=TRUE)$specific.spec)
  expect_message(blocksFun(vY, mX, include.empty=TRUE)$specific.spec)
  expect_message(blocksFun(vY, mX,include.gum=TRUE, include.1cut=TRUE, include.empty=TRUE)$specific.spec)
  expect_message(blocksFun(vY, mX, max.paths=1)$specific.spec)
  expect_message(blocksFun(vY, mX, parallel.options=2)$specific.spec)
  
  expect_error(blocksFun(vY, mX, parallel.options=5)$specific.spec) #should return error
  
  expect_message(blocksFun(vY, mX, turbo=TRUE)$specific.spec)
  expect_message(blocksFun(vY, mX, force.invertibility=FALSE)$specific.spec)
  expect_message(blocksFun(vY, mX, tol=1)$specific.spec)
  expect_message(blocksFun(vY, mX, LAPACK=TRUE)$specific.spec)
  
  expect_error(blocksFun(vY, mX, max.regs=1)$specific.spec) #should give error
  
  expect_silent(blocksFun(vY, mX, print.searchinfo=FALSE)$specific.spec)
  expect_message(blocksFun(vY, mX, alarm=TRUE)$specific.spec)
})


##################################################
## 3 TEST SPECIAL CASES
##################################################

##a dgp:
set.seed(123)
y <- rnorm(20)
x <- matrix(rnorm(length(y)*10), length(y), 10)
mXX <- list()
mXX[[1]] <- x[,1:5]
mXX[[2]] <- x[,6:10]

##list without names:
blocksFun(y, mXX)

test_that("blocksFun Special case - list without names",{
  expect_true(is.list(blocksFun(y, mXX)))
})

##a single block:
myblocks <- list()
myblocks[[1]] <- list(1:10)
blocksFun(y, x, blocks=myblocks)

test_that("blocksFun Special case - single block",{
  expect_message(blocksFun(y, x, blocks=myblocks))
  expect_true(is.list(blocksFun(y, x, blocks=myblocks)))
})

##keep:
blocksFun(y, mXX, keep=1)
mykeep <- list()
mykeep[[1]] <- 2
mykeep[[2]] <- 3
blocksFun(y, mXX, keep=mykeep)

test_that("blocksFun Special Case - keep",{
  expect_true(is.list(blocksFun(y, mXX, keep=1)))
  expect_true(is.list(blocksFun(y, mXX, keep=mykeep)))
  expect_message(blocksFun(y, mXX, keep=1))
  expect_message(blocksFun(y, mXX, keep=mykeep))
})

##compare isat() and blocksFun():
y <- rnorm(50)
mconst <- rep(1,length(y))
x <- list(IIS=cbind(mconst,coredata(iim(y))),
          SIS=cbind(mconst,coredata(sim(y))))
blocksFun(y, x, keep=list(1,1))
isat(y, iis=TRUE, sis=TRUE)
##compare speed (not much difference really...):
system.time(blocksFun(y,x))
system.time(isat(y, iis=TRUE, sis=TRUE))

test_that("compare isat() and blocksFun",{
  blocks <- blocksFun(y, x, keep=list(1,1))
  iis <- isat(y, iis=TRUE, sis=TRUE)
  expect_equal(as.list(blocks$specific.spec),iis$ISfinalmodels,ignore_attr=TRUE)
})

##################################################
## 7 SIMULATIONS (FOR THE FUTURE)
##################################################
