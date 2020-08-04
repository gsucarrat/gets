I AM HERE!!!

#######################################################
## This code replicates the analyses in the JSS-paper
## on the gets package
##
## Note: Some of the examples relies on other R
## packages (all on CRAN): lgarch, MASS, glmnet, and lars
##
#######################################################

#######################################################
##START not in jstatsoft code:

##set working directory:
setwd("C:/Users/sucarrat/Documents/R/gs/gets/devel/")
#setwd(choose.dir())

##load required packages:
require(parallel)
require(zoo)

##remove everything in workspace (.GlobaleEnv:
rm(list=ls())

##load source:
source("gets-base-source.R")
source("gets-isat-source.R")

##END not in jstatsoft code
#######################################################


###=====================================================
### Section 2: An overview, alternatives, and development principles
###=====================================================
#
### 2.1. Variable selection properties of gets, replicating Table 1
###----------------
#
### The following code undertakes the simulations of the Hoover and Perez
### (1999) experiments.
###
### The HP1 experiment: The DGP is y_t = 130*z_t with
### z_t~IIN(0, 1), where the GUM contains 40 (irrelevant)
### variables and an intercept (variable no. 41).
###
### The HP2' experiment: The DGP is y_t = 0.75y_{t-1}
### + 130*(sqrt(7)/4)*z_t, where z_t ~ IIN(0, 1), and
### where the GUM contains a constant and the 40 (irrelevant)
### variables from the HP1 experiment. Note that the DGP in
### HP is written in a different but equivalent form: y_t = 130u_t*,
### where u_t* = 0.75u_{t-1}* + (sqrt(7)/4)*z_t, z_t ~ IIN(0, 1).
###
### The HP7' experiment: The DGP is y_t = 0.75y_{t-1}
### + 1.33x11_t -0.9975x29_t + 6.73*z_t, where z_t ~ IIN(0, 1), and
### where the GUM contains a constant (fixed during simplification)
### and the other variables from the HP1 experiment.
#
#set.seed(123)
#
### required packages
#library("gets")
#library("lgarch")  # needed for the glag function
#
### load and prepare data
#data("hpdata", package = "gets")  # load data
#HPdata <- coredata(as.zoo(hpdata[, -1]))  # convert to matrix, remove date-index
#colnames(HPdata) <- NULL  # remove colnames
#cor(na.trim(HPdata))
#sqrt(diag(var(HPdata, na.rm = TRUE)))
#
### only first 18, because no. 19 is the consumption variable
#vLags <- rep(1, 18)
#vLags[c(2, 3, 5, 9, 10, 18)] <- 2
#dHPdata <- NULL
#for (i in 1:18) {
#  dHPdata <- cbind(dHPdata, gdiff(HPdata[, i], lag = vLags[i]))
#}
#cor(na.trim(dHPdata), use = "everything")
#cor(na.trim(dHPdata), use = "all.obs")
#cor(na.trim(dHPdata), use = "complete.obs")
#cor(na.trim(dHPdata), use = "na.or.complete")
#cor(na.trim(dHPdata), use = "pairwise.complete.obs")
#sqrt(diag(var(dHPdata, na.rm = TRUE)))
#
### add lags to dHPdata
#for (i in 1:18) {
#  dHPdata <- cbind(dHPdata, glag(dHPdata[, i]))
#}
#
### set simulation parameters
#sReps <- 1000
#nobs <- NROW(dHPdata)  # no. of obs. in each replication
#sBurnIn <- 50  # burn-in sample for HP2
#t.pval <- 0.05  # PET sig.level
#
### simulation output
#winners <- list()
#k0.hat <- rep(NA, nobs)
#k1.hat <- rep(NA, nobs)
#
### DGPs
#DGPs <- c("HP1", "HP2", "HP7")
#for (dgp in DGPs) {
#    if (dgp == "HP1") {
#        k1 <- 40  # no. of irrelevant variables
#    }
#    if (dgp == "HP2") {
#        k0 <- 1
#        k1 <- 39  # no. relevant and irrelevant variables
#    }
#    if (dgp == "HP7") {
#        k0 <- 3
#        k1 <- 37  # no. relevant and irrelevant variables
#    }
#
#    ## start time
#    simsStarted <- date()
#
#    ## simulate (replications)
#    for (i in 1:sReps) {
#
#        ## dgp is HP1
#        if (dgp == "HP1") {
#            y <- 130 * rnorm(nobs)
#        }
#
#        ## dgp is HP2
#        if (dgp == "HP2") {
#            u <- rnorm(nobs + sBurnIn)
#            y <- rep(0, nobs + sBurnIn)
#            sigmau <- 130 * (sqrt(7)/4)
#            for (j in 2:length(y)) {
#                y[j] <- 0.75 * y[c(j - 1)] + sigmau * u[j]
#            }
#            y <- y[c(sBurnIn + 1):length(y)]  # delete burn-in values
#        }
#
#        ## dgp is HP7
#        if (dgp == "HP7") {
#            u <- rnorm(nobs)
#            y <- rep(0, nobs)
#            for (j in 3:length(y)) {
#                y[j] <- 0.75 * y[c(j - 1)] + 1.33 * dHPdata[j, 11] -
#                    0.9975 * dHPdata[j, 29] + 6.44 * u[j]
#            }
#        }
#
#        ## create regressors
#        mHP <- dHPdata
#        for (j in 1:4) {
#            mHP <- cbind(mHP, glag(y, k = j))
#        }
#        mHP <- cbind(mHP, 1)
#
#        ## estimate gum
#        HPgum <- arx(y[5:nobs], mxreg = mHP[5:nobs, ], vcov.type = "ordinary")
#
#        ## do gets
#        HPgets <- getsm(HPgum, keep = 41, t.pval = 0.05, do.pet = TRUE, ar.LjungB = NULL,
#                        arch.LjungB = NULL, normality.JarqueB = NULL, include.gum = FALSE, include.empty = FALSE,
#                        info.method = "sc", verbose = FALSE, print.searchinfo = FALSE, estimate.specific = FALSE,
#                        plot = FALSE, alarm = FALSE)
#
#        ## record properties of winner
#        whereWinner <- which.min(HPgets$terminals.results[, 1])
#        winnerSpec <- HPgets$terminals[[whereWinner]]
#        if (dgp == "HP1") {
#            k1.hat[i] <- length(winnerSpec) - 1
#        }
#        if (dgp == "HP2") {
#            k0.hat[i] <- as.numeric(37 %in% winnerSpec)
#            k1.hat[i] <- length(setdiff(winnerSpec, c(37, 41)))
#        }
#        if (dgp == "HP7") {
#            k0.hat[i] <- sum(c(11, 29, 37) %in% winnerSpec)
#            k1.hat[i] <- length(setdiff(winnerSpec, c(11, 29, 37, 41)))
#        }
#
#        ## replication no
#        if (i %in% c(1, 10, 100, 250, 500, 750, 900, 940, 960, 980, 990)) {
#            print(i)
#        }
#
#    }
#
#    cat("DGP = ", dgp, "\n")
#    print(simsStarted)
#    print(date())
#
#    ## print results
#    if (dgp == "HP1") {
#        cat("irrelevance proportion: ", mean(k1.hat/k1), "\n")
#        cat("p(DGP): ", sum(k1.hat == 0)/sReps, "\n")
#    }
#    if (dgp == "HP2") {
#        cat("relevance proportion: ", mean(k0.hat/k0), "\n")
#        cat("irrelevance proportion: ", mean(k1.hat/k1), "\n")
#        cat("p(DGP): ", sum(I(I(k0.hat == 1) + I(k1.hat == 0)) == 2)/sReps, "\n")
#    }
#    if (dgp == "HP7") {
#        cat("relevance proportion: ", mean(k0.hat/k0), "\n")
#        cat("irrelevance proportion: ", mean(k1.hat/k1), "\n")
#        cat("p(DGP): ", sum(I(I(k0.hat == 3) + I(k1.hat == 0)) == 2)/sReps, "\n")
#    }
#}
#
##---------------------------
#### 2.2 Computational Speed Comparisons
#
#k <- 80
#n <- 200
#beta <- rep(0, k)
#alpha = 0.01
#set.seed(123)
#library(gets)
#
####Setting up the DGP:
#x <- as.matrix(replicate(k, rnorm(n)))
#xnames<-paste("x", seq(1:k), sep="")
#colnames(x) <- xnames
#y <- x%*%beta + rnorm(n,0,1)
#
#### Selection using getsm
#m1 <- arx(y, mxreg=x, mc=TRUE)
#start.time <- Sys.time()
#m1.sel <- getsm(m1, t.pval=alpha, print.searchinfo = FALSE)
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#print(paste("Simulation Complete in:"))
#print(time.taken)
#
####### comparison to other packages
#
#### Selection using "glmnet"
#library(glmnet)
#start.time <- Sys.time()
#object1 <- glmnet(x, y, alpha=1)
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#print(paste("Simulation Complete in:"))
#print(time.taken)
#
#### Selection using "lars"
#library(lars)
#start.time <- Sys.time()
#object2 <- lars(x, y, type= c("lasso"))
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#print(paste("Simulation Complete in:"))
#print(time.taken)
#
#### Selection using "backward"
#mod1 <- lm(y ~ x)
#start.time <- Sys.time()
#sel1 <- step(mod1, direction="backward")
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#print(paste("Simulation Complete in:"))
#print(time.taken)
#
#
###---------------------------
### 2.2. Variable selection - comparison of GETS and gets with
### alternatives, replicating Figure 1 and Tables 5-7
### WARNING: Running this code with 1000 replications takes more than 5hrs on a standard 1.8Ghz Processor.
#
#N <- 20  # total number of variables (=20 in paper)
#T <- 100  # total number of observations (denoted n in paper)
#K <- 10  # maximum number of relevant variables (=10 in paper)
#reps <- 1000  # replications in simulation (=1000 in paper)
#alpha <- 0.01  # for selection and testing
#s_fixed <- 0.25  # fixed lasso penalty for non-cross-validated LASSO
#beta_val <- 0.3  # coefficient on variables in DGP
#numcor <- 3  # 3 setups of correlation (uncorrelated, positive correlation, alternating negative/positive correlation)
#rhocorr <- 0.5  # correlation parameter
#
#set.seed(123)
#require("glmnet")
#require("MASS")
#
#start.time <- Sys.time()
#
###---------------------------
### Starting the Simulation Loop
### looping over correlation structures
#
#simulation.data <- list()
#
#for (c in 1:numcor) {
#
#  if (c == 1) {
#    rho <- 0
#    negative <- FALSE
#  }
#  if (c == 2) {
#    rho <- rhocorr
#    negative <- FALSE
#  }
#  if (c == 3) {
#    rho <- rhocorr
#    negative <- TRUE
#  }
#
#  filename <- paste("beta_results_corr_", rho, "_neg_", negative, "_T_", T, "_beta",
#    beta_val, ".csv", sep = "")
#  print(paste("Correlation:", rho, ", negative correlations:", negative, sep = ""))
#  beta <- rep(0, N)
#
#  ## dataframe to store results
#  beta.results <- data.frame(matrix(NA, nrow = (K + 1), ncol = 1))
#  names(beta.results) <- "K"
#
#  for (j in 1:(K + 1)) {
#    ## increase number of relevant variables
#
#    print(paste("Correlation:", rho, ", negative correlations:", negative, sep = ""))
#    print(paste("Variable loop:", j - 1))
#
#    if (j > 1) {
#      beta[1:(j - 1)] <- beta_val
#    }
#    print(beta)
#    beta.results$K[j] <- j - 1
#    lass.results <- data.frame(matrix(NA, reps, 1))
#    names(lass.results) <- "rep"
#
#    ## the simulation loop
#    for (i in 1:reps) {
#
#      print(paste("Variable loop:", j))
#      print(paste("Replication #:", i))
#
#      means <- c(rep(0, N))
#      Sigma <- matrix(rho, ncol = N, nrow = N)
#      diag(Sigma) <- 1
#      x <- mvrnorm(n = T, mu = means, Sigma, empirical = FALSE)
#      if (negative == TRUE) {
#        x[, seq(from = 1, to = N, by = 2)] <- x[, seq(from = 1, to = N, by = 2)] *
#          (-1)
#      }
#      xnames <- paste("x", seq(1:N), sep = "")
#      colnames(x) <- xnames
#
#      y <- x %*% beta + rnorm(T, 0, 1)  # The DGP
#      ################################
#
#      ## crossvalidated LASSO:
#      object1 <- glmnet(x, y, alpha = 1)
#      cv_vals <- cv.glmnet(x, y)
#      coef_cv <- coef(object1, s = cv_vals$lambda.min)
#      ret_cv <- which(coef_cv != 0) - 1
#      ret_cv <- ret_cv[-1]
#      x_lass_cv <- xnames[ret_cv]
#      ######################
#
#      ## lasso with fixed s
#      object2 <- glmnet(x, y, alpha = 1)
#      coef_sfix <- coef(object2, s = s_fixed)
#      ret_sfix <- which(coef_sfix != 0) - 1
#      ret_sfix <- ret_sfix[-1]
#      x_lass_sfix <- xnames[ret_sfix]
#      ######################
#
#      ## getsm
#      mod1 <- arx(y, mxreg = x, mc = TRUE)
#      mod1.sel <- getsm(mod1, t.pval = alpha, keep = 1, arch.LjungB = NULL,
#        ar.LjungB = NULL, print.searchinfo = FALSE)
#      x_getsm <- mod1.sel$aux$mXnames
#      ################
#
#      ## 1-cut (significance in GUM)
#      mod1.p <- mod1$mean.results$`p-value`[-1]
#      x_1cut <- xnames[which(mod1.p < alpha)]
#      ##########################
#
#      ## estimating the DGP (optimal, infeasible)
#      if (any(beta != 0)) {
#        dgp <- arx(y, mxreg = x[, (1:max(which(beta != 0)))], mc = TRUE)
#        dgp.p <- dgp$mean.results$`p-value`[-1]
#        x_dgp <- xnames[which(dgp.p < alpha)]
#      }
#
#      ## computing gauge and potency
#
#      ## Total retained
#      getsm.ret <- length(x_getsm) - 1  #-1 due to constant
#      lass_cv.ret <- length(x_lass_cv)
#      lass_fix.ret <- length(x_lass_sfix)
#      x1cut.ret <- length(x_1cut)
#
#      ## Relevant retained
#      rel_names <- xnames[which(beta != 0)]
#      tot.rel <- length(rel_names)
#
#      getsm.rel.ret <- length(which(x_getsm %in% rel_names == TRUE))
#      lass_cv.rel.ret <- length(which(x_lass_cv %in% rel_names == TRUE))
#      lass_fix.rel.ret <- length(which(x_lass_sfix %in% rel_names == TRUE))
#      x1cut.rel.ret <- length(which(x_1cut %in% rel_names == TRUE))
#
#      ## irrelevant retained
#      getsm.irr.ret <- getsm.ret - getsm.rel.ret
#      lass_cv.irr.ret <- lass_cv.ret - lass_cv.rel.ret
#      lass_fix.irr.ret <- lass_fix.ret - lass_fix.rel.ret
#      x1cut.irr.ret <- x1cut.ret - x1cut.rel.ret
#
#      ## gauges
#      tot.irr <- length(which(beta == 0))
#
#      if (tot.irr != 0) {
#        gauge.gets <- getsm.irr.ret/tot.irr
#        gauge.lass_cv <- lass_cv.irr.ret/tot.irr
#        gauge.lass_fix <- lass_fix.irr.ret/tot.irr
#        gauge.1cut <- x1cut.irr.ret/tot.irr
#      } else {
#        gauge.gets <- NA
#        gauge.lass_cv <- NA
#        gauge.lass_fix <- NA
#        gauge.1cut <- NA
#      }
#
#      ## potencies
#      if (tot.rel != 0) {
#        pot.gets <- getsm.rel.ret/tot.rel
#        pot.lass_cv <- lass_cv.rel.ret/tot.rel
#        pot.lass_fix <- lass_fix.rel.ret/tot.rel
#        pot.1cut <- x1cut.rel.ret/tot.rel
#        pot.dgp <- length(x_dgp)/tot.rel
#
#      } else {
#        pot.gets <- NA
#        pot.lass_cv <- NA
#        pot.lass_fix <- NA
#        pot.1cut <- NA
#        pot.dgp <- NA
#      }
#
#      ## storing results
#      lass.results$rep[i] <- i
#      lass.results$gauge.gets[i] <- gauge.gets
#      lass.results$gauge.lass_cv[i] <- gauge.lass_cv
#      lass.results$gauge.lass_fix[i] <- gauge.lass_fix
#      lass.results$gauge.1cut[i] <- gauge.1cut
#      lass.results$pot.gets[i] <- pot.gets
#      lass.results$pot.lass_cv[i] <- pot.lass_cv
#      lass.results$pot.lass_fix[i] <- pot.lass_fix
#      lass.results$pot.1cut[i] <- pot.1cut
#      lass.results$pot.dgp[i] <- pot.dgp
#
#    } ## Replication Simulation Loop closed
#
#    beta.results$gauge.gets[j] <- mean(lass.results$gauge.gets, na.rm = TRUE)
#    beta.results$gauge.lass_cv[j] <- mean(lass.results$gauge.lass_cv, na.rm = TRUE)
#    beta.results$gauge.lass_fix[j] <- mean(lass.results$gauge.lass_fix, na.rm = TRUE)
#    beta.results$gauge.1cut[j] <- mean(lass.results$gauge.1cut, na.rm = TRUE)
#    beta.results$pot.gets[j] <- mean(lass.results$pot.gets, na.rm = TRUE)
#    beta.results$pot.lass_cv[j] <- mean(lass.results$pot.lass_cv, na.rm = TRUE)
#    beta.results$pot.lass_fix[j] <- mean(lass.results$pot.lass_fix, na.rm = TRUE)
#    beta.results$pot.1cut[j] <- mean(lass.results$pot.1cut, na.rm = TRUE)
#    beta.results$pot.dgp[j] <- mean(lass.results$pot.dgp, na.rm = TRUE)
#  } ## looping over number of variables closed
#
#  simulation.data[[c]] <- beta.results
#
#} ## correlation loop closed
#
### Simulation Results (replicating Tables 5-7)
#print(simulation.data)
#
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken
#print(paste("Simulation Complete in:"))
#print(time.taken)
#
### Replicating Figure 1
#
### load data generated from simulation loop above
#dat1 <- simulation.data[[1]]  # uncorrelated
#dat2 <- simulation.data[[2]]  # positively correlated
#dat3 <- simulation.data[[3]]  # positive/negative correlations
#
### record current graph-values
#def.par <- par(no.readonly = TRUE)
#
#pdf("fig_getsm_performance_v2.pdf", width=11, height=6)
#
#op <- par(mfrow = c(2, 3), oma = c(3, 2, 1, 0) + 0.1, mar = c(3.9, 4, 1, 1) + 0.1)
#
#xl <- 2 # linewidth
#
### Plot of the Gauge
#plot(dat1$K, dat1$gauge.gets, type = "l", ylim = c(-0.01, 1), cex.lab = 1.5, col = "red",
#  ylab = "Gauge", main = "Uncorrelated Regr.", xlab = "", lwd = xl, xaxt = "n",
#  xlim = c(0, 10))
#lines(dat1$K, dat1$gauge.lass_cv, col = "purple", lwd = xl)
#lines(dat1$K, dat1$gauge.lass_fix, col = "purple", lty = 2, lwd = xl)
#lines(dat1$K, dat1$gauge.1cut, col = "blue", lwd = xl)
#abline(h = 0.01, lty = 3)
#axis(1, at = seq(0, 10, by = 1), las = 1, cex.axis = 1.2)
#
#plot(dat2$K, dat2$gauge.gets, type = "l", ylim = c(-0.01, 1), col = "red", cex.lab = 1.5,
#  ylab = "Gauge", main = "Pos. Correlated Regr.", xlab = "", lwd = xl, xaxt = "n",
#  xlim = c(0, 10))
#lines(dat2$K, dat2$gauge.lass_cv, col = "purple", lwd = xl)
#lines(dat2$K, dat2$gauge.lass_fix, col = "purple", lty = 2, lwd = xl)
#lines(dat2$K, dat2$gauge.1cut, col = "blue", lwd = xl)
#abline(h = 0.01, lty = 3)
#legend(1, 0.95, c("getsm", "LassCV", "LassFix", "1-cut", "Nominal 1%"), bg = NA,
#  bty = "n", title.adj = -0.03, lty = c(1, 1, 2, 1, 3), col = c("red", "purple",
#    "purple", "blue", "black"), lwd = 2, cex = 1.5, pt.cex = 1, x.intersp = 0.5,
#  y.intersp = 0.7)
#axis(1, at = seq(0, 10, by = 1), las = 1, cex.axis = 1.2)
#
#plot(dat3$K, dat3$gauge.gets, type = "l", ylim = c(-0.01, 1), col = "red", cex.lab = 1.5,
#  ylab = "Gauge", main = "Pos./Neg. Correlated Regr.", xlab = "", lwd = xl, xaxt = "n",
#  xlim = c(0, 10))
#lines(dat3$K, dat3$gauge.lass_cv, col = "purple", lwd = xl)
#lines(dat3$K, dat3$gauge.lass_fix, col = "purple", lty = 2, lwd = xl)
#lines(dat3$K, dat3$gauge.1cut, col = "blue", lwd = xl)
#abline(h = 0.01, lty = 3)
#axis(1, at = seq(0, 10, by = 1), las = 1, cex.axis = 1.2)
#
### Plot of the Potency
#plot(dat1$K, dat1$pot.gets, type = "l", ylim = c(-0.01, 1), col = "red", cex.lab = 1.5,
#  ylab = "Potency", xlab = "#Relevant Regressors (out of 20)", lwd = xl, xaxt = "n",
#  xlim = c(0, 10))
#lines(dat1$K, dat1$pot.lass_cv, col = "purple", lwd = xl)
#lines(dat1$K, dat1$pot.lass_fix, col = "purple", lty = 2, lwd = xl)
#lines(dat1$K, dat1$pot.1cut, col = "blue", lwd = xl)
#lines(dat1$K, dat1$pot.dgp, col = "gray55", lty = 4, lwd = xl)
#legend(1, 0.4, c("getsm", "LassCV", "LassFix", "1-cut", "DGP"), bg = NA, bty = "n",
#  title.adj = -0.03, lty = c(1, 1, 2, 1, 4), col = c("red", "purple", "purple",
#    "blue", "gray55"), lwd = 2, cex = 1.5, pt.cex = 1, x.intersp = 0.5, y.intersp = 0.7)
#axis(1, at = seq(0, 10, by = 1), las = 1, cex.axis = 1.2)
#
#plot(dat2$K, dat2$pot.gets, type = "l", ylim = c(-0.01, 1), col = "red", cex.lab = 1.5,
#  ylab = "Potency", xlab = "#Relevant Regressors (out of 20)", lwd = xl, xaxt = "n",
#  xlim = c(0, 10))
#lines(dat2$K, dat2$pot.lass_cv, col = "purple", lwd = xl)
#lines(dat2$K, dat2$pot.lass_fix, col = "purple", lty = 2, lwd = xl)
#lines(dat2$K, dat2$pot.1cut, col = "blue", lwd = xl)
#lines(dat2$K, dat2$pot.dgp, col = "gray55", lty = 4, lwd = xl)
#axis(1, at = seq(0, 10, by = 1), las = 1, cex.axis = 1.2)
#
#plot(dat3$K, dat3$pot.gets, type = "l", ylim = c(-0.01, 1), col = "red", cex.lab = 1.5,
#  ylab = "Potency", xlab = "#Relevant Regressors (out of 20)", lwd = xl, xaxt = "n",
#  xlim = c(0, 10))
#lines(dat3$K, dat3$pot.lass_cv, col = "purple", lwd = xl)
#lines(dat3$K, dat3$pot.lass_fix, col = "purple", lty = 2, lwd = xl)
#lines(dat3$K, dat3$pot.1cut, col = "blue", lwd = xl)
#lines(dat3$K, dat3$pot.dgp, col = "gray55", lty = 4, lwd = xl)
#axis(1, at = seq(0, 10, by = 1), las = 1, cex.axis = 1.2)
#
#dev.off()
#
### return to default par options
#par(def.par, no.readonly = FALSE)
#

##=====================================================
## Section 4: The ar-x model w/log-arch-x errors
##=====================================================

## 4.1. Simulation
##----------------

set.seed(123)  # for reproducability
y <- arima.sim(list(ar = 0.4), 100)
library("lgarch")
eps <- lgarchSim(100, arch = 0.3, garch = 0)
yy <- arima.sim(list(ar = 0.4), 100, innov = eps)
plot(as.zoo(cbind(y, yy, eps)))

## 4.2. arx: Estimation
##---------------------

##library("gets") #this line is not commented out in the jstatsoftcode
mod01 <- arx(y, ar = 1)
mod01
mX <- matrix(rnorm(100 * 5), 100, 5)
mod02 <- arx(y, mc = TRUE, ar = 1:2, mxreg = mX, vcov.type = "white")
mod02
mod03 <- arx(eps, arch = 1)
mod03
mod04 <- arx(eps, arch = 1:3, asym = 2, vxreg = log(mX^2))
mod04
mod05 <- arx(yy, mc = TRUE, ar = 1:2, mxreg = mX, arch = 1:3, asym = 2, 
  vxreg = log(mX^2), vcov.type = "white")
mod05

## 4.4. Example: A model of quarterly inflation with time-varying conditional
## variance
##-----------------------------------------------------

data("infldata", package = "gets")
infldata <- zooreg(infldata[, -1], frequency = 4, start = c(1989, 1))
inflMod01 <- inflMod01 <- arx(infldata[, "infl"], mc = TRUE, ar = 1:4, mxreg = infldata[, 
  2:4], vcov.type = "white")
inflMod01
inflMod02 <- inflMod02 <- arx(infldata[, "infl"], mc = TRUE, ar = 1:4, mxreg = infldata[, 
  2:4], arch = 1:4, vxreg = infldata[, 2:4], vcov.type = "white")
inflMod02
info.criterion(as.numeric(logLik(inflMod01)), n = 104, k = 8 + 1)
info.criterion(as.numeric(logLik(inflMod02)), n = 100, k = 8 + 8)

## 4.5. Example: A model of daily SP500 volatility
##------------------------------------------------

data("sp500data", package = "gets")
sp500data <- zoo(sp500data[, -1], order.by = as.Date(sp500data[, "Date"]))
sp500data <- window(sp500data, start = as.Date("1983-07-01"))
head(sp500data)
plot(sp500data)

## make log-returns in %
sp500Ret <- diff(log(sp500data[, "Adj.Close"])) * 100

## make lagged volatility proxy (range-based)
relrange <- (log(sp500data[, "High"]) - log(sp500data[, "Low"])) * 100
volproxy <- log(relrange^2)
volproxylag <- lag(volproxy, k = -1)

## make volume variable
volume <- log(sp500data[, "Volume"])
volumediff <- diff(volume) * 100
volumedifflag <- lag(volumediff, k = -1)

## make day-of-the-week dummies
sp500Index <- index(sp500Ret)
days <- weekdays(sp500Index)
days <- union(days, days)
dTue <- zoo(as.numeric(weekdays(sp500Index) == days[1]), order.by = sp500Index)
dWed <- zoo(as.numeric(weekdays(sp500Index) == days[2]), order.by = sp500Index)
dThu <- zoo(as.numeric(weekdays(sp500Index) == days[3]), order.by = sp500Index)
dFri <- zoo(as.numeric(weekdays(sp500Index) == days[4]), order.by = sp500Index)

## estimate log-arch(5)-x
sp500Mod01 <- arx(sp500Ret, arch = 1:5, log.ewma = c(5, 20, 60, 120), asym = 1, vxreg = cbind(volproxylag, 
  volumedifflag, dTue, dWed, dThu, dFri))
sp500Mod01

## comparison with log-garch(1, 1)
library("lgarch")
sp500Mod02 <- lgarch(sp500Ret)
sp500Mod02
logLik(sp500Mod02)
info.criterion(as.numeric(logLik(sp500Mod01)), n = 8235, k = 17)
info.criterion(as.numeric(logLik(sp500Mod02)), n = 8240, k = 3)


##=====================================================
## Section 5: gets-modeling
##=====================================================

## 5.1 getsm: Modeling the mean
##------------------------------

getsm05 <- getsm(mod05)
getsm05
getsm05a <- getsm(mod05, t.pval = 0.01, arch.LjungB = NULL)
getsm05a
getsm05b <- getsm(mod05, keep = 1)
getsm05b

## 5.2. getsv: Modeling the log-variance
##---------------------------------------

getsv05 <- getsv(mod05)
getsv05
mod06 <- arx(residuals(getsm05), arch = 1:3, asym = 2, vxreg = log(mX^2))
getsv06 <- getsv(mod06)
getsv06
getsv06b <- getsv(mod06, keep = 1:4, ar.LjungB = NULL)
getsv06b

## 5.4. Example: A parsimonious model of quarterly inflation
##------------------------------------------------

## gets modeling
inflMod03 <- getsm(inflMod02)
inflMod03
inflMod04 <- arx(residuals(inflMod03), arch = 1:4, vxreg = infldata[, 2:4])
inflMod05 <- getsv(inflMod04, ar.LjungB = list(lag = 5, pval = 0.025))
inflMod05

## out-of-sample forecasts of sd
inflMod06 <- inflMod06 <- arx(infldata[, "infl"], mc = TRUE, ar = c(1, 4), arch = 1:2, 
  vxreg = infldata[, 2:4], vcov.type = "white")
inflMod06
newvxreg <- matrix(0, 4, 3)
colnames(newvxreg) <- c("q2dum", "q3dum", "q4dum")
newvxreg[2, "q2dum"] <- 1
newvxreg[3, "q3dum"] <- 1
newvxreg[4, "q4dum"] <- 1
set.seed(123)  # for reproducability
predict(inflMod06, n.ahead = 4, spec = "variance", newvxreg = newvxreg)
2 - 2 * sqrt(1.0448239)
2 + 2 * sqrt(1.0448239)
2 - 2 * sqrt(0.2101471)
2 + 2 * sqrt(0.2101471)

## 5.5. Example: A parsimonious model of daily SP500 volatility
##--------------------------------------------------

sp500Mod03 <- getsv(sp500Mod01, t.pval = 0.001, arch.LjungB = NULL)
sp500Mod03
info.criterion(as.numeric(logLik(sp500Mod03)), n = 8235, k = 7)


##=====================================================
## Section 6: indicator saturation
##=====================================================

options(plot = TRUE)  # turn on plotting
so2 <- data("so2data", package = "gets")
yso2 <- zoo(so2data[, "DLuk_tot_so2"], order.by = so2data[, "year"])
sis <- isat(yso2, t.pval = 0.01)
sis
##check:
coefpaper <- c(0.01465385, -0.04332051, -0.11693333, 0.12860000,
  -0.28400000, 0.24550000, -0.11550000)
coefsis <- coef(sis)
all( abs(coefpaper-coefsis) < 1e-8 )

x1972 <- zoo(sim(so2data[, "year"])[, 26], order.by = so2data[, "year"])
sis2 <- isat(yso2, t.pval = 0.01, mxreg = x1972)
##check:
coefpaper <- c(0.01465385, -0.04332051, -0.11693333, 0.12860000,
  -0.28400000, 0.24550000, -0.11550000)
coefsis2 <- coef(sis2)
all( abs(coefpaper-coefsis2) < 1e-8 )

sisvar <- isatvar(sis)
sisvar
##check:
sisvarpaper <- #first two years, i.e. 1946 and 1947, as in paper
  c(0.00000000, 0.01465385, 6.291637e-05, 0.007931984,
  0.00000000, 0.01465385, 6.291637e-05, 0.007931984)
sisvartmp <- c(as.vector(sisvar[1,]), as.vector(sisvar[2,]))
all( abs(sisvarpaper-sisvartmp) < 1e-8 )

iis <- isat(yso2, ar = 1, sis = FALSE, iis = TRUE, t.pval = 0.05)
## With Correction (conscorr = TRUE, effcorr = TRUE)
iisvar <- isatvar(iis, conscorr = TRUE, effcorr = TRUE)
##check:
iisvarpaper <- #first two years, i.e. 1946 and 1947, as in paper
  c(0.00000000, -0.006210179, 7.225479e-05, 0.008500282,
  0.00000000, -0.006210179, 7.225479e-05, 0.008500282)
iisvartmp <- c(as.vector(iisvar[1,]), as.vector(iisvar[2,]))
##not all are true (isatvar has changed?)
abs(iisvarpaper-iisvartmp) < 1e-8

## Without Correction (conscorr = FALSE, effcorr = FALSE)
isatvar(iis, conscorr = FALSE, effcorr = FALSE)

bcorr <- biascorr(b = sisvar[, "const.path"], b.se = sisvar[, "const.se"], p.alpha = 0.01, 
  T = length(sisvar[, "const.path"]))
bcorr

isattest(sis, hnull = 0, lr = FALSE, ci.pval = 0.99, plot.turn = TRUE, biascorr = TRUE)


### Replicating Figure 2 and Table 8: False Detection Rate of isat
###------------------------------------------------
#### Warning: This code will take a long time to run when using 1000 replications.
#
#library("gets")
#set.seed(123)
#T_sim <- c(30, 50, 100, 200, 300)  # sample sizes
#reps <- 1000  # number of replications in the simulation, set to 1000 to match paper.
#alpha <- 0.01  # target-level of significance of selection
#
#sample.res <- data.frame(matrix(NA, length(T_sim), 1))
#names(sample.res) <- "T"
#
#start.time <- Sys.time()
#for (j in 1:length(T_sim)) {
#
#  T <- T_sim[j]
#  print(T)
#
#  print(paste("New sample: ", j, ":", T))
#
#  sim.res <- data.frame(matrix(NA, reps, 1))
#  names(sim.res) <- "rep"
#
#  for (i in 1:reps) {
#
#    print(paste("sample: ", j, ":", T))
#    print(paste("rep:", i))
#
#    y <- rnorm(T, 0, 1)
#
#    ## iis
#    iis.m <- isat(y, t.pval = alpha, iis = TRUE, sis = FALSE, tis = FALSE, print.searchinfo = FALSE)
#
#    ## iis gauge
#    iis.ret <- length(iis.m$ISnames)
#    iis.gauge <- iis.ret/T
#
#    ## sis
#    sis.m <- isat(y, t.pval = alpha, iis = FALSE, sis = TRUE, tis = FALSE, print.searchinfo = FALSE)
#    sis.ret <- length(sis.m$ISnames)
#    sis.gauge <- sis.ret/(T - 1)
#
#    ## tis
#    tis.m <- isat(y, t.pval = alpha, iis = FALSE, sis = FALSE, tis = TRUE, print.searchinfo = FALSE)
#    tis.ret <- length(tis.m$ISnames)
#    tis.gauge <- tis.ret/(T)
#
#    sim.res$rep[i] <- i
#    sim.res$iis.gauge[i] <- iis.gauge
#    sim.res$sis.gauge[i] <- sis.gauge
#    sim.res$tis.gauge[i] <- tis.gauge
#
#  } ## i loop closed
#
#  sample.res$T[j] <- T
#  sample.res$iis.gauge[j] <- mean(sim.res$iis.gauge)
#  sample.res$sis.gauge[j] <- mean(sim.res$sis.gauge)
#  sample.res$tis.gauge[j] <- mean(sim.res$tis.gauge)
#
#} ## j sample loop closed
#
#
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken
#
### Replicating Figure 2
###----------------------
#
### record current graph-values
#def.par <- par(no.readonly = TRUE)
#
#pdf("fig_isat_performance_v2.pdf", width=11, height=6)
#
#op <- par(mfrow = c(1, 1), oma = c(3, 2, 1, 0) + 0.1, mar = c(3.9, 4, 1, 1) + 0.1)
#
#xl <- 2
#plot(sample.res$T, sample.res$iis.gauge, type = "l", ylim = c(-0.01, 0.1), cex.lab = 1.5,
#  col = "red", ylab = "Gauge", main = "isat: False Detection rate", xlab = "Sample Size n",
#  lwd = xl)
#lines(sample.res$T, sample.res$sis.gauge, col = "blue", lwd = xl)
#lines(sample.res$T, sample.res$tis.gauge, col = "darkgreen", lty = 2, lwd = xl)
#abline(h = 0.01, lty = 3, lwd = 2)
#legend(100, 0.1, c("IIS", "SIS", "TIS", "Nominal 1%"), bg = NA, bty = "n", title.adj = -0.03,
#  lty = c(1, 1, 2, 3), col = c("red", "blue", "darkgreen", "black"), lwd = 2, cex = 1.5,
#  pt.cex = 1, x.intersp = 0.5, y.intersp = 0.7)
#
#dev.off()
#
### Replicating Table 8
#
#print(sample.res)
#
### return to default par options
#par(def.par, no.readonly = FALSE)


##=====================================================
## Section 7: export to eviews and stata, and
## generation of latex-code
##=====================================================

eviews(getsm05)
stata(getsm05)
