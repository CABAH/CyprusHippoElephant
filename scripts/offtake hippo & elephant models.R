#############################################################################################################################
## dwarf hippo & elephant model (Cyprus)
## species: Phanourios minor (PM)
##          Palaeoloxodon cypriotes (PC)
## Corey Bradshaw
## October 2023
## Flinders University
#############################################################################################################################

# libraries
library(plotly)
library(ggpubr)

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

## source
source("matrixOperators.r")


################################################################
## offtake models (requires running base models script first) ##
################################################################


itdiv <- iter/10
Q.ext.thresh <- 100/2 # quasi-extinction threshold


## PROGRESSIVELY REMOVE INCREASING NUMBER OF INDIVIDUALS RANDOMLY FROM N VECTOR ##
iter <- 1000
ind.rem.vec <- seq(0,1600,50)

PM.Q.ext.pr <- PC.Q.ext.pr <- rep(NA,length(ind.rem.vec))

for (k in 1:length(ind.rem.vec)) {
  
  ## PHANOURIOS (PM)
  ## set storage matrices & vectors
  PM.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
  
  for (e in 1:iter) {
    PM.popmat <- PM.popmat.orig
    
    PM.n.mat <- matrix(0, nrow=PM.age.max+1,ncol=(t+1))
    PM.n.mat[,1] <- PM.init.vec
    
    for (i in 1:t) {
      # stochastic survival values
      PM.s.alpha <- estBetaParams(PM.Sx, PM.s.sd.vec^2)$alpha
      PM.s.beta <- estBetaParams(PM.Sx, PM.s.sd.vec^2)$beta
      PM.s.stoch <- rbeta(length(PM.s.alpha), PM.s.alpha, PM.s.beta)
      
      if (rbinom(1, 1, 0.14/PM.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        PM.s.stoch <- PM.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      PM.fert.stch <- rnorm(length(PM.popmat[,1]), PM.pred.p.mm, PM.m.sd.vec)
      
      PM.totN.i <- sum(PM.n.mat[,i], na.rm=T)
      PM.pred.red <- PM.a.lp/(1+(PM.totN.i/PM.b.lp)^PM.c.lp)
      
      diag(PM.popmat[2:(PM.age.max+1),]) <- (PM.s.stoch[-(PM.age.max+1)])*PM.pred.red
      PM.popmat[PM.age.max+1,PM.age.max+1] <- (PM.s.stoch[PM.age.max+1])*PM.pred.red
      PM.popmat[1,] <- ifelse(PM.fert.stch < 0, 0, PM.fert.stch)
      PM.n.mat[,i+1] <- PM.popmat %*% PM.n.mat[,i]
      
      PM.n.mat[,i+1] <- ifelse((PM.n.mat[,i+1] - round((ind.rem.vec[k]*PM.ssd),0)) < 0, 0, (PM.n.mat[,i+1] - round((ind.rem.vec[k]*PM.ssd),0))) # remove animals following ssd
      
    } # end i loop
    
    PM.n.sums.mat[e,] <- (as.vector(colSums(PM.n.mat)))
    
    #if (e %% itdiv==0) print(e) 
    
  } # end e loop
  
  # total N
  PM.n.md <- apply(PM.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  PM.n.up <- apply(PM.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  PM.n.lo <- apply(PM.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  PM.Q.ext.mat <- ifelse(PM.n.sums.mat < Q.ext.thresh, 1, 0)
  PM.Q.ext.sum <- apply(PM.Q.ext.mat[,ceiling(PM.gen.l):dim(PM.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  PM.Q.ext.pr[k] <- length(which(PM.Q.ext.sum > 0)) / iter
  print("PHANOURIOS")
  
  
  ## PALAEOLOXODON (PC)
  ## set storage matrices & vectors
  PC.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
  
  for (e in 1:iter) {
    PC.popmat <- PC.popmat.orig
    
    PC.n.mat <- matrix(0, nrow=PC.age.max+1,ncol=(t+1))
    PC.n.mat[,1] <- PC.init.vec
    
    for (i in 1:t) {
      # stochastic survival values
      PC.s.alpha <- estBetaParams(PC.Sx, PC.s.sd.vec^2)$alpha
      PC.s.beta <- estBetaParams(PC.Sx, PC.s.sd.vec^2)$beta
      PC.s.stoch <- rbeta(length(PC.s.alpha), PC.s.alpha, PC.s.beta)
      
      if (rbinom(1, 1, 0.14/PC.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        PC.s.stoch <- PC.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertilty sampler (gaussian)
      PC.fert.stch <- rnorm(length(PC.popmat[,1]), PC.pred.p.mm, PC.m.sd.vec)
      
      PC.totN.i <- sum(PC.n.mat[,i], na.rm=T)
      PC.pred.red <- PC.a.lp/(1+(PC.totN.i/PC.b.lp)^PC.c.lp)
      
      diag(PC.popmat[2:(PC.age.max+1),]) <- (PC.s.stoch[-(PC.age.max+1)])*PC.pred.red
      PC.popmat[PC.age.max+1,PC.age.max+1] <- (PC.s.stoch[PC.age.max+1])*PC.pred.red
      PC.popmat[1,] <- ifelse(PC.fert.stch < 0, 0, PC.fert.stch)
      PC.n.mat[,i+1] <- PC.popmat %*% PC.n.mat[,i]
      
      PC.n.mat[,i+1] <- ifelse((PC.n.mat[,i+1] - round((ind.rem.vec[k]*PC.ssd),0)) < 0, 0, (PC.n.mat[,i+1] - round((ind.rem.vec[k]*PC.ssd),0))) # remove animals following ssd
      
    } # end i loop
    
    PC.n.sums.mat[e,] <- (as.vector(colSums(PC.n.mat)))
    
    #if (e %% itdiv==0) print(e) 
    
  } # end e loop
  
  # total N
  PC.n.md <- apply(PC.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  PC.n.up <- apply(PC.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  PC.n.lo <- apply(PC.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # quasi-extinction probability
  PC.Q.ext.mat <- ifelse(PC.n.sums.mat < Q.ext.thresh, 1, 0)
  PC.Q.ext.sum <- apply(PC.Q.ext.mat[,ceiling(PC.gen.l):dim(PC.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  PC.Q.ext.pr[k] <- length(which(PC.Q.ext.sum > 0)) / iter
  print("PALAEOLOXODON")
  
  print(ind.rem.vec[k]) 
} # end k loop

plot(ind.rem.vec, PM.Q.ext.pr, type="l", xlab="individuals removed/year", ylab="Pr(quasi-extinction)", ylim=c(0,1))
lines(ind.rem.vec, PC.Q.ext.pr, lty=2)

## compile species outputs
spp.mass.vec <- c(PM.mass,PC.mass)

spp.gen.l.vec <- c(PM.gen.l,PC.gen.l)

PM.auc <- sum(sum(PM.Q.ext.pr)/length(ind.rem.vec))
PC.auc <- sum(sum(PC.Q.ext.pr)/length(ind.rem.vec))

Qext.auc <- c(PM.auc,PC.auc)

labs.vec <- c("PM","PC")
Qextpr.dat <- data.frame(ind.rem.vec,PM.Q.ext.pr,PC.Q.ext.pr)
colnames(Qextpr.dat) <- c("IndRem",labs.vec)
dat.out <- data.frame(spp.mass.vec, Qext.auc)
colnames(dat.out) <- c("M","AUC")
rownames(dat.out) <- labs.vec
