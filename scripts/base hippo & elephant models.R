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


######################################################
## BASE MODELS
######################################################

############################
## Phanourios minor (PM)
## body mass estimates
PM.mass <- 132 # Phanourios minor (Lomolino et al. 2013 https://doi.org/10.1111/jbi.12096)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
PM.rm.pred <- 10^(0.6914 - (0.2622*log10(PM.mass*1000)))
PM.rm.pred
PM.lm.pred <- exp(PM.rm.pred)
PM.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
PM.D.pred <- (10^(4.196 - (0.74*log10(PM.mass*1000))))/2 # divided by 2 for females only
PM.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
PM.age.max <- round(10^(0.89 + (0.13*log10(PM.mass*1000))), 0)
PM.age.max

## age vector
PM.age.vec <- 0:PM.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
PM.F.pred <- exp(2.719 - (0.211*log(PM.mass*1000)))/2 # divided by 2 for females
PM.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
PM.alpha <- ceiling(exp(-1.34 + (0.214*log(PM.mass*1000))))
PM.alpha

## define m function with age
PM.m.vec <- c(rep(0, PM.alpha-1), rep(0.75*PM.F.pred, round(PM.alpha/2,0)), rep(PM.F.pred, (PM.age.max+1-((PM.alpha-1+round(PM.alpha/2,0))))))
PM.m.sd.vec <- 0.05*PM.m.vec
plot(PM.age.vec, PM.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
PM.m.dat <- data.frame(PM.age.vec, PM.m.vec)
param.init <- c(0.6, 4, -5)
PM.fit.logp <- nls(PM.m.vec ~ a / (1+(PM.age.vec/b)^c), 
                data = PM.m.dat,
                algorithm = "port",
                start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                trace = TRUE,      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PM.fit.logp.summ <- summary(PM.fit.logp)
plot(PM.age.vec, PM.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
PM.age.vec.cont <- seq(0,max(PM.age.vec),1)
PM.pred.p.m <- coef(PM.fit.logp)[1] / (1+(PM.age.vec.cont/coef(PM.fit.logp)[2])^coef(PM.fit.logp)[3])
PM.pred.p.mm <- ifelse(PM.pred.p.m > 1, 1, PM.pred.p.m)
lines(PM.age.vec.cont, PM.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
PM.s.tran <- ln.a.s + b.s*log(PM.mass*1000) + log(1)
PM.s.ad.yr <- exp(-exp(PM.s.tran))
PM.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.965*PM.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 1.9 # rate of mortality decline (also known as bt)
a2 <- 1 - PM.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.05 # rate of mortality increase
longev <- PM.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
#plot(x,l.x,type="l")
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
PM.Sx <- c(0.935*PM.s.ad.yr, 1 - qx)
plot(x, PM.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
PM.s.sd.vec <- 0.05*PM.Sx

# qx from Nicolaou et al. 2020 Quat Int (doi:10.1016/j.quaint.2020.09.016)
qx.Nic <- read.table("qx-Nicolaou.csv", header=T, sep=",")
PM.Sx.Nic <- c(1 - qx.Nic$qx)
plot(qx.Nic$x, PM.Sx.Nic, type="l")
lines(x, PM.Sx, lty=2, col="red")
PM.popmat.Nic <- matrix(data = 0, nrow=max(qx.Nic$x)+1, ncol=max(qx.Nic$x)+1)
diag(PM.popmat.Nic[2:(max(qx.Nic$x)+1),]) <- c(PM.Sx[1],PM.Sx.Nic[1:(max(qx.Nic$x)-1)])
PM.popmat.Nic[max(qx.Nic$x)+1,max(qx.Nic$x)+1] <- PM.Sx.Nic[max(qx.Nic$x)]
PM.popmat.Nic[1,] <- PM.pred.p.mm[1:21]
colnames(PM.popmat.Nic) <- c(0:max(qx.Nic$x))
rownames(PM.popmat.Nic) <- c(0:max(qx.Nic$x))
max.lambda(PM.popmat.Nic) ## 1-yr lambda
PM.lm.pred
max.r(PM.popmat.Nic) # rate of population change, 1-yr

## create matrix
PM.popmat <- matrix(data = 0, nrow=PM.age.max+1, ncol=PM.age.max+1)
diag(PM.popmat[2:(PM.age.max+1),]) <- PM.Sx[-(PM.age.max+1)]
PM.popmat[PM.age.max+1,PM.age.max+1] <- PM.Sx[PM.age.max+1]
PM.popmat[1,] <- PM.pred.p.mm
colnames(PM.popmat) <- c(0:PM.age.max)
rownames(PM.popmat) <- c(0:PM.age.max)
PM.popmat.orig <- PM.popmat ## save original matrix

## matrix properties
max.lambda(PM.popmat.orig) ## 1-yr lambda
PM.lm.pred
max.r(PM.popmat.orig) # rate of population change, 1-yr
PM.ssd <- stable.stage.dist(PM.popmat.orig) ## stable stage distribution
plot(PM.age.vec, PM.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(PM.popmat.orig, PM.age.max) # reproductive value
PM.gen.l <- G.val(PM.popmat.orig, PM.age.max) # mean generation length
PM.gen.l

## initial population vector
area <- 11193.62388 # km2 @14.125 ka
PM.pop.found <- round(area*PM.D.pred, 0) # population size
PM.pop.found
PM.init.vec <- PM.ssd * PM.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*PM.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

PM.tot.F <- sum(PM.popmat.orig[1,])
PM.popmat <- PM.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
PM.n.mat <- matrix(0, nrow=PM.age.max+1,ncol=(t+1))
PM.n.mat[,1] <- PM.init.vec

## set up projection loop
for (i in 1:t) {
  PM.n.mat[,i+1] <- PM.popmat %*% PM.n.mat[,i]
}

PM.n.pred <- colSums(PM.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(PM.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
PM.K.max <- 1*PM.pop.found
PM.K.vec <- c(1, PM.K.max/2, 0.75*PM.K.max, PM.K.max) 
PM.red.vec <- c(1,0.935,0.855,0.79)
plot(PM.K.vec, PM.red.vec,pch=19,type="b")
PM.Kred.dat <- data.frame(PM.K.vec, PM.red.vec)

# logistic power function a/(1+(x/b)^c)
PM.param.init <- c(1, 2*PM.K.max, 2)
PM.fit.lp <- nls(PM.red.vec ~ a/(1+(PM.K.vec/b)^c), 
                 data = PM.Kred.dat,
                 algorithm = "port",
                 start = c(a = PM.param.init[1], b = PM.param.init[2], c = PM.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PM.fit.lp.summ <- summary(PM.fit.lp)
plot(PM.K.vec, PM.red.vec, pch=19,xlab="N",ylab="reduction factor")
PM.K.vec.cont <- seq(1,2*PM.pop.found,1)
PM.pred.lp.fx <- coef(PM.fit.lp)[1]/(1+(PM.K.vec.cont/coef(PM.fit.lp)[2])^coef(PM.fit.lp)[3])
lines(PM.K.vec.cont, PM.pred.lp.fx, lty=3,lwd=3,col="red")

PM.a.lp <- coef(PM.fit.lp)[1]
PM.b.lp <- coef(PM.fit.lp)[2]
PM.c.lp <- coef(PM.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
PM.n.mat <- matrix(0, nrow=PM.age.max+1, ncol=(t+1))
PM.n.mat[,1] <- PM.init.vec
PM.popmat <- PM.popmat.orig

## set up projection loop
for (i in 1:t) {
  PM.totN.i <- sum(PM.n.mat[,i])
  PM.pred.red <- as.numeric(PM.a.lp/(1+(PM.totN.i/PM.b.lp)^PM.c.lp))
  diag(PM.popmat[2:(PM.age.max+1),]) <- (PM.Sx[-(PM.age.max+1)])*PM.pred.red
  PM.popmat[PM.age.max+1,PM.age.max+1] <- (PM.Sx[PM.age.max+1])*PM.pred.red
  PM.popmat[1,] <- PM.pred.p.mm
  PM.n.mat[,i+1] <- PM.popmat %*% PM.n.mat[,i]
}

PM.n.pred <- colSums(PM.n.mat)
plot(yrs, PM.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=PM.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 100
itdiv <- iter/10

PM.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
PM.s.arr <- PM.m.arr <- array(data=NA, dim=c(t+1, PM.age.max+1, iter))

for (e in 1:iter) {
  PM.popmat <- PM.popmat.orig
  
  PM.n.mat <- matrix(0, nrow=PM.age.max+1,ncol=(t+1))
  PM.n.mat[,1] <- PM.init.vec
  
  for (i in 1:t) {
    # stochastic survival values
    PM.s.alpha <- estBetaParams(PM.Sx, PM.s.sd.vec^2)$alpha
    PM.s.beta <- estBetaParams(PM.Sx, PM.s.sd.vec^2)$beta
    PM.s.stoch <- rbeta(length(PM.s.alpha), PM.s.alpha, PM.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    PM.fert.stch <- rnorm(length(PM.popmat[,1]), PM.pred.p.mm, PM.m.sd.vec)
    PM.m.arr[i,,e] <- ifelse(PM.fert.stch < 0, 0, PM.fert.stch)
    
    PM.totN.i <- sum(PM.n.mat[,i], na.rm=T)
    PM.pred.red <- PM.a.lp/(1+(PM.totN.i/PM.b.lp)^PM.c.lp)
    
    diag(PM.popmat[2:(PM.age.max+1),]) <- (PM.s.stoch[-(PM.age.max+1)])*PM.pred.red
    PM.popmat[PM.age.max+1,PM.age.max+1] <- (PM.s.stoch[PM.age.max+1])*PM.pred.red
    PM.popmat[1,] <- PM.m.arr[i,,e]
    PM.n.mat[,i+1] <- PM.popmat %*% PM.n.mat[,i]
    
    PM.s.arr[i,,e] <- PM.s.stoch * PM.pred.red
    
  } # end i loop
  
  PM.n.sums.mat[e,] <- ((as.vector(colSums(PM.n.mat))/PM.pop.found))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

PM.n.md <- apply(PM.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
PM.n.up <- apply(PM.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PM.n.lo <- apply(PM.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

par(mfrow=c(1,3))
plot(yrs,PM.n.md,type="l", main = "", xlab="year", ylab="pN1", lwd=2, ylim=c(0.95*min(PM.n.lo),1.05*max(PM.n.up)))
lines(yrs,PM.n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,PM.n.up,lty=2,col="red",lwd=1.5)

PM.s.add <- PM.m.add  <- rep(0, PM.age.max+1)
for (m in 1:iter) {
  PM.s.add <- rbind(PM.s.add, PM.s.arr[ceiling(PM.gen.l):(t+1),,m])
  PM.m.add <- rbind(PM.m.add, PM.m.arr[ceiling(PM.gen.l):(t+1),,m])
}
PM.s.add <- PM.s.add[-1,]
PM.m.add <- PM.m.add[-1,]

PM.s.md <- apply(PM.s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
PM.s.up <- apply(PM.s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PM.s.lo <- apply(PM.s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(PM.age.vec,PM.s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(PM.s.lo),1.05*max(PM.s.up)))
lines(PM.age.vec,PM.s.lo,lty=2,col="red",lwd=1.5)
lines(PM.age.vec,PM.s.up,lty=2,col="red",lwd=1.5)

PM.m.md <- apply(PM.m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
PM.m.up <- apply(PM.m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PM.m.lo <- apply(PM.m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(PM.age.vec,PM.m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(PM.m.lo),1.05*max(PM.m.up)))
lines(PM.age.vec,PM.m.lo,lty=2,col="red",lwd=1.5)
lines(PM.age.vec,PM.m.up,lty=2,col="red",lwd=1.5)
par(mfrow=c(1,1))




############################
## Palaeoloxodon cypriotes (PC)

# mass
PC.mass <- 531 # Palaeoloxodon cypriotes (Lomolino et al. 2013 https://doi.org/10.1111/jbi.12096)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
PC.rm.pred <- 10^(0.6914 - (0.2622*log10(PC.mass*1000)))
PC.rm.pred
PC.lm.pred <- exp(PC.rm.pred)
PC.lm.pred

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
PC.D.pred <- (10^(4.196 - (0.74*log10(PC.mass*1000))))/2 # divided by 2 for females only
PC.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
PC.age.max <- round(10^(0.89 + (0.13*log10(PC.mass*1000))), 0)
PC.age.max
# predicted max longevity for ungulates
# Y=a+bx, a=1.21941, b=0.135772
10^(1.270706+(log10(PC.mass)*0.124650))  # https://doi.org/10.1038%2Fs41598-021-02192-4 SI

## age vector
PC.age.vec <- 0:PC.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
PC.F.pred <- exp(2.719 - (0.211*log(PC.mass*1000)))/2 # divided by 2 for females
PC.F.pred

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
PC.alpha <- ceiling(exp(-1.34 + (0.214*log(PC.mass*1000))))
PC.alpha
round(10^(0.052366+(log10(531)*0.200756)),0) # https://doi.org/10.1038%2Fs41598-021-02192-4 SI
 
## define m function with age
PC.m.vec <- c(rep(0, PC.alpha-1), rep(0.75*PC.F.pred, round(PC.alpha/2,0)), rep(PC.F.pred, (PC.age.max+1-((PC.alpha-1+round(PC.alpha/2,0))))))
PC.m.sd.vec <- 0.05*PC.m.vec
plot(PC.age.vec, PC.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
PC.m.dat <- data.frame(PC.age.vec, PC.m.vec)
param.init <- c(0.1, 2, -5)
PC.fit.logp <- nls(PC.m.vec ~ a / (1+(PC.age.vec/b)^c), 
                   data = PC.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PC.fit.logp.summ <- summary(PC.fit.logp)
plot(PC.age.vec, PC.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
PC.age.vec.cont <- seq(0,max(PC.age.vec),1)
PC.pred.p.m <- coef(PC.fit.logp)[1] / (1+(PC.age.vec.cont/coef(PC.fit.logp)[2])^coef(PC.fit.logp)[3])
PC.pred.p.mm <- ifelse(PC.pred.p.m > 1, 1, PC.pred.p.m)
lines(PC.age.vec.cont, PC.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
PC.s.tran <- ln.a.s + b.s*log(PC.mass*1000) + log(1)
PC.s.ad.yr <- exp(-exp(PC.s.tran))
PC.s.ad.yr

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (1.01*PC.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 2.7 # rate of mortality decline (also known as bt)
a2 <- 1 - PC.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.02 # rate of mortality increase
longev <- PC.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
#plot(x,l.x,type="l")
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
PC.Sx <- c(0.995*PC.s.ad.yr, 1 - qx)
plot(x, PC.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
PC.s.sd.vec <- 0.05*PC.Sx

## create matrix
PC.popmat <- matrix(data = 0, nrow=PC.age.max+1, ncol=PC.age.max+1)
diag(PC.popmat[2:(PC.age.max+1),]) <- PC.Sx[-(PC.age.max+1)]
PC.popmat[PC.age.max+1,PC.age.max+1] <- PC.Sx[PC.age.max+1]
PC.popmat[1,] <- PC.pred.p.mm
colnames(PC.popmat) <- c(0:PC.age.max)
rownames(PC.popmat) <- c(0:PC.age.max)
PC.popmat.orig <- PC.popmat ## save original matrix

## matrix properties
max.lambda(PC.popmat.orig) ## 1-yr lambda
PC.lm.pred
max.r(PC.popmat.orig) # rate of population change, 1-yr
PC.ssd <- stable.stage.dist(PC.popmat.orig) ## stable stage distribution
plot(PC.age.vec, PC.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(PC.popmat.orig, PC.age.max) # reproductive value
PC.gen.l <- G.val(PC.popmat.orig, PC.age.max) # mean generation length
PC.gen.l

## initial population vector
PC.pop.found <- round(area*PC.D.pred, 0) # founding population size
PC.pop.found
PC.init.vec <- PC.ssd * PC.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*PC.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

PC.tot.F <- sum(PC.popmat.orig[1,])
PC.popmat <- PC.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
PC.n.mat <- matrix(0, nrow=PC.age.max+1,ncol=(t+1))
PC.n.mat[,1] <- PC.init.vec

## set up projection loop
for (i in 1:t) {
  PC.n.mat[,i+1] <- PC.popmat %*% PC.n.mat[,i]
}

PC.n.pred <- colSums(PC.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(PC.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
PC.K.max <- 1*PC.pop.found
PC.K.vec <- c(1, PC.K.max/2, 0.75*PC.K.max, PC.K.max) 
PC.red.vec <- c(1,0.96,0.90,0.829)
plot(PC.K.vec, PC.red.vec,pch=19,type="b")
PC.Kred.dat <- data.frame(PC.K.vec, PC.red.vec)

# logistic power function a/(1+(x/b)^c)
PC.param.init <- c(1, 2*PC.K.max, 2)
PC.fit.lp <- nls(PC.red.vec ~ a/(1+(PC.K.vec/b)^c), 
                 data = PC.Kred.dat,
                 algorithm = "port",
                 start = c(a = PC.param.init[1], b = PC.param.init[2], c = PC.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
PC.fit.lp.summ <- summary(PC.fit.lp)
plot(PC.K.vec, PC.red.vec, pch=19,xlab="N",ylab="reduction factor")
PC.K.vec.cont <- seq(1,2*PC.pop.found,1)
PC.pred.lp.fx <- coef(PC.fit.lp)[1]/(1+(PC.K.vec.cont/coef(PC.fit.lp)[2])^coef(PC.fit.lp)[3])
lines(PC.K.vec.cont, PC.pred.lp.fx, lty=3,lwd=3,col="red")

PC.a.lp <- coef(PC.fit.lp)[1]
PC.b.lp <- coef(PC.fit.lp)[2]
PC.c.lp <- coef(PC.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
PC.n.mat <- matrix(0, nrow=PC.age.max+1, ncol=(t+1))
PC.n.mat[,1] <- PC.init.vec
PC.popmat <- PC.popmat.orig

## set up projection loop
for (i in 1:t) {
  PC.totN.i <- sum(PC.n.mat[,i])
  PC.pred.red <- as.numeric(PC.a.lp/(1+(PC.totN.i/PC.b.lp)^PC.c.lp))
  diag(PC.popmat[2:(PC.age.max+1),]) <- (PC.Sx[-(PC.age.max+1)])*PC.pred.red
  PC.popmat[PC.age.max+1,PC.age.max+1] <- (PC.Sx[PC.age.max+1])*PC.pred.red
  PC.popmat[1,] <- PC.pred.p.mm
  PC.n.mat[,i+1] <- PC.popmat %*% PC.n.mat[,i]
}

PC.n.pred <- colSums(PC.n.mat)
plot(yrs, PC.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=PC.pop.found, lty=2, col="red", lwd=2)


## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 100
itdiv <- iter/10

PC.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
PC.s.arr <- PC.m.arr <- array(data=NA, dim=c(t+1, PC.age.max+1, iter))

for (e in 1:iter) {
  PC.popmat <- PC.popmat.orig
  
  PC.n.mat <- matrix(0, nrow=PC.age.max+1,ncol=(t+1))
  PC.n.mat[,1] <- PC.init.vec

  for (i in 1:t) {
    # stochastic survival values
    PC.s.alpha <- estBetaParams(PC.Sx, PC.s.sd.vec^2)$alpha
    PC.s.beta <- estBetaParams(PC.Sx, PC.s.sd.vec^2)$beta
    PC.s.stoch <- rbeta(length(PC.s.alpha), PC.s.alpha, PC.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    PC.fert.stch <- rnorm(length(PC.popmat[,1]), PC.pred.p.mm, PC.m.sd.vec)
    PC.m.arr[i,,e] <- ifelse(PC.fert.stch < 0, 0, PC.fert.stch)
    
    PC.totN.i <- sum(PC.n.mat[,i], na.rm=T)
    PC.pred.red <- PC.a.lp/(1+(PC.totN.i/PC.b.lp)^PC.c.lp)
    
    diag(PC.popmat[2:(PC.age.max+1),]) <- (PC.s.stoch[-(PC.age.max+1)])*PC.pred.red
    PC.popmat[PC.age.max+1,PC.age.max+1] <- (PC.s.stoch[PC.age.max+1])*PC.pred.red
    PC.popmat[1,] <- PC.m.arr[i,,e]
    PC.n.mat[,i+1] <- PC.popmat %*% PC.n.mat[,i]

    PC.s.arr[i,,e] <- PC.s.stoch * PC.pred.red
    
  } # end i loop
  
  PC.n.sums.mat[e,] <- ((as.vector(colSums(PC.n.mat))/PC.pop.found))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

PC.n.md <- apply(PC.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
PC.n.up <- apply(PC.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PC.n.lo <- apply(PC.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

par(mfrow=c(1,3))
plot(yrs,PC.n.md,type="l", main = "", xlab="year", ylab="pN1", lwd=2, ylim=c(0.95*min(PC.n.lo),1.05*max(PC.n.up)))
lines(yrs,PC.n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,PC.n.up,lty=2,col="red",lwd=1.5)

PC.s.add <- PC.m.add  <- rep(0, PC.age.max+1)
for (m in 1:iter) {
  PC.s.add <- rbind(PC.s.add, PC.s.arr[ceiling(PC.gen.l):(t+1),,m])
  PC.m.add <- rbind(PC.m.add, PC.m.arr[ceiling(PC.gen.l):(t+1),,m])
}
PC.s.add <- PC.s.add[-1,]
PC.m.add <- PC.m.add[-1,]

PC.s.md <- apply(PC.s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
PC.s.up <- apply(PC.s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PC.s.lo <- apply(PC.s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(PC.age.vec,PC.s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(PC.s.lo),1.05*max(PC.s.up)))
lines(PC.age.vec,PC.s.lo,lty=2,col="red",lwd=1.5)
lines(PC.age.vec,PC.s.up,lty=2,col="red",lwd=1.5)

PC.m.md <- apply(PC.m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
PC.m.up <- apply(PC.m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
PC.m.lo <- apply(PC.m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(PC.age.vec,PC.m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(PC.m.lo),1.05*max(PC.m.up)))
lines(PC.age.vec,PC.m.lo,lty=2,col="red",lwd=1.5)
lines(PC.age.vec,PC.m.up,lty=2,col="red",lwd=1.5)
par(mfrow=c(1,1))


