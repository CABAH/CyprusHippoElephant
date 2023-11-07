#############################################################################################################################
## dwarf hippo & elephant model (Cyprus)
## species: Phanourios minor (PM)
##          Palaeoloxodon cypriotes (PC)
## global sensitivity analysis
## Corey Bradshaw
## November 2023
## Flinders University
#############################################################################################################################

##########################################################
## SET PARAMETER LIMITS FOR GLOBAL SENSITIVITY ANALYSIS ##
##########################################################
#
## testing effect on extinction risk (Phanourios and Palaeoloxodon, respectively)
## of modifying the following variables:
##
###########################################################################################################
###########################################################################################################
## relative prey choice probability (PM vs. PC):
# 1. PM.prob3 = 0.7765041 set in model; ± 10%; 0.7765041-(0.1*0.7765041) to 0.7765041+(0.1*0.7765041)
##
## standard deviation of relative prey choice probability (PM vs. PC):
# 2. PM.prob3.sd.pc = 0.05 set in model; 0.025 to 0.1
##
## proportion of meat in human diet
# 3. propdietmeat = 0.65 set in model; 0.5 to 0.7
##
## energy intake by hunter-gatherers (female)
# 4. f.E.int.m = 1877 in model; ± 10%; 1877-(0.1*1877) to 1877+(0.1*1877)
##
## energy intake by hunter-gatherers (male)
# 5. m.E.int.m <- 2649 set in model; ± 10%; 2649-(0.1*2649) to 2649+(0.1*2649)
##
## founding population size (PM)
# 6. PM.pop.found = 14280 set in model; ± 20%; round(14280-(0.2*14280),0) to round(14280+(0.2*14280),0)
##
## founding population size (PC)
# 7. PC.pop.found = 5098 set in model; ± 20%; round(5098-(0.2*5098),0) to round(5098+(0.2*5098),0)
##
## relative proportion of female vs. male prey harvested
# 8. FvsMhunt = 0.5 set in model; 0.4 to 0.6
##
## kCal per 100 grams of meat
# 9. kCalmeat = 130 set in model; ± 10%; 130-(0.1*130) to 130+(0.1*130)
##
## median proportion meat edible per carcase
# 10. mpsamp = 0.3156746 set in model; 0.25 to 0.4
###########################################################################################################
###########################################################################################################

## remove everything
rm(list = ls())

# libraries
library(plotly)
library(ggpubr)
library(truncnorm)
library(doSNOW)
library(iterators)
library(snow)
library(foreach)
library(lhs)
library(data.table)
library(dismo)
library(gbm)

  ## functions
  # beta distribution shape parameter estimator function
  estBetaParams <- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }
  estBetaParams2 <- function(mu, var) {
    alpha <- -((mu * (var + mu^2 - mu))/var)
    beta <- ((var + mu^2 - mu)*(mu - 1))/var
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
  source("~/Documents/Papers/Other/Global human population/ms/PNAS/R1/matrixOperators.r")
  
  ## Set up parallel processing (nproc is the number of processing cores to use)
  nproc <- 6
  cl.tmp = makeCluster(rep('localhost',nproc), type='SOCK')
  registerDoSNOW(cl.tmp)
  getDoParWorkers()
  
    
      ######################################################
      ## BASE MODELS
      ######################################################
      
      ############################
      ## Phanourios minor (PM)
      ## body mass estimates
      PM.mass <- 132 # Phanourios minor (Lomolino et al. 2013 https://doi.org/10.1111/jbi.12096)
      
      ## predicted rm (from Henneman 1983 Oecologia 56:104-108)
      PM.rm.pred <- 10^(0.6914 - (0.2622*log10(PM.mass*1000)))
      PM.lm.pred <- exp(PM.rm.pred)
    
      ## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
      PM.D.pred <- (10^(4.196 - (0.74*log10(PM.mass*1000))))/2 # divided by 2 for females only
    
      ## max age
      PM.age.max <- round(10^(0.89 + (0.13*log10(PM.mass*1000))), 0)
    
      ## age vector
      PM.age.vec <- 0:PM.age.max
      
      ## fertility
      PM.F.pred <- exp(2.719 - (0.211*log(PM.mass*1000)))/2 # divided by 2 for females
    
      ## age at primiparity
      PM.alpha <- ceiling(exp(-1.34 + (0.214*log(PM.mass*1000))))
    
      ## define m function with age
      PM.m.vec <- c(rep(0, PM.alpha-1), rep(0.75*PM.F.pred, round(PM.alpha/2,0)), rep(PM.F.pred, (PM.age.max+1-((PM.alpha-1+round(PM.alpha/2,0))))))
      PM.m.sd.vec <- 0.05*PM.m.vec
    
      # logistic power function y = a / (1+(x/b)^c)
      PM.m.dat <- data.frame(PM.age.vec, PM.m.vec)
      param.init <- c(0.6, 4, -5)
      PM.fit.logp <- nls(PM.m.vec ~ a / (1+(PM.age.vec/b)^c), 
                      data = PM.m.dat,
                      algorithm = "port",
                      start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                      trace = F,      
                      nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
      PM.fit.logp.summ <- summary(PM.fit.logp)
      PM.age.vec.cont <- seq(0,max(PM.age.vec),1)
      PM.pred.p.m <- coef(PM.fit.logp)[1] / (1+(PM.age.vec.cont/coef(PM.fit.logp)[2])^coef(PM.fit.logp)[3])
      PM.pred.p.mm <- ifelse(PM.pred.p.m > 1, 1, PM.pred.p.m)
      
      ## survival
      ln.a.s <- -0.5; b.s <- -0.25
      PM.s.tran <- ln.a.s + b.s*log(PM.mass*1000) + log(1)
      PM.s.ad.yr <- exp(-exp(PM.s.tran))
      
      # Siler hazard h(x) (Gurven et al. 2007)
      a1 <- 1 - (0.965*PM.s.ad.yr) # initial infant mortality rate (also known as αt)
      b1 <- 1.9 # rate of mortality decline (also known as bt)
      a2 <- 1 - PM.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
      a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
      b3 <- 0.05 # rate of mortality increase
      longev <- PM.age.max
      x <- seq(0,longev,1) # age vector
      h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
      l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
      init.pop <- 10000
      lx <- round(init.pop*l.x,0)
      len.lx <- length(lx)
      dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
      qx <- dx/lx[1:(length(lx)-1)]
      PM.Sx <- c(0.935*PM.s.ad.yr, 1 - qx)
      PM.s.sd.vec <- 0.05*PM.Sx
      
      # qx from Nicolaou et al. 2020 QI
      qx.Nic <- read.table("qx-Nicolaou.csv", header=T, sep=",")
      PM.Sx.Nic <- c(1 - qx.Nic$qx)
      PM.popmat.Nic <- matrix(data = 0, nrow=max(qx.Nic$x)+1, ncol=max(qx.Nic$x)+1)
      diag(PM.popmat.Nic[2:(max(qx.Nic$x)+1),]) <- c(PM.Sx[1],PM.Sx.Nic[1:(max(qx.Nic$x)-1)])
      PM.popmat.Nic[max(qx.Nic$x)+1,max(qx.Nic$x)+1] <- PM.Sx.Nic[max(qx.Nic$x)]
      PM.popmat.Nic[1,] <- PM.pred.p.mm[1:21]
      colnames(PM.popmat.Nic) <- c(0:max(qx.Nic$x))
      rownames(PM.popmat.Nic) <- c(0:max(qx.Nic$x))
      
      ## create matrix
      PM.popmat <- matrix(data = 0, nrow=PM.age.max+1, ncol=PM.age.max+1)
      diag(PM.popmat[2:(PM.age.max+1),]) <- PM.Sx[-(PM.age.max+1)]
      PM.popmat[PM.age.max+1,PM.age.max+1] <- PM.Sx[PM.age.max+1]
      PM.popmat[1,] <- PM.pred.p.mm
      colnames(PM.popmat) <- c(0:PM.age.max)
      rownames(PM.popmat) <- c(0:PM.age.max)
      PM.popmat.orig <- PM.popmat ## save original matrix
      
      ## matrix properties
      PM.ssd <- stable.stage.dist(PM.popmat.orig) ## stable stage distribution
      PM.gen.l <- G.val(PM.popmat.orig, PM.age.max) # mean generation length
      
      ## initial population vector
      area <- 11193.62388 # km2 @14.125 ka
      PM.init.vec <- PM.ssd * round(area*PM.D.pred, 0)
      
      #################
      ## project
      ## set time limit for projection in 1-yr increments
      yr.st <- 1
      yr.end <- round(40*PM.gen.l, 0) # set projection end date
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
      } # end i loop
      
      PM.n.pred <- colSums(PM.n.mat)
      yrs <- seq(yr.st, yr.end, 1)
    
      # compensatory density feedback
      PM.K.max <- 1*round(area*PM.D.pred, 0)
      PM.K.vec <- c(1, PM.K.max/2, 0.75*PM.K.max, PM.K.max) 
      PM.red.vec <- c(1,0.935,0.855,0.79)
      PM.Kred.dat <- data.frame(PM.K.vec, PM.red.vec)
      
      # logistic power function a/(1+(x/b)^c)
      PM.param.init <- c(1, 2*PM.K.max, 2)
      PM.fit.lp <- nls(PM.red.vec ~ a/(1+(PM.K.vec/b)^c), 
                       data = PM.Kred.dat,
                       algorithm = "port",
                       start = c(a = PM.param.init[1], b = PM.param.init[2], c = PM.param.init[3]),
                       trace = F,      
                       nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
      PM.fit.lp.summ <- summary(PM.fit.lp)
      PM.K.vec.cont <- seq(1,2*round(area*PM.D.pred, 0),1)
      PM.pred.lp.fx <- coef(PM.fit.lp)[1]/(1+(PM.K.vec.cont/coef(PM.fit.lp)[2])^coef(PM.fit.lp)[3])
    
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
      } # end i loop
      
      PM.n.pred <- colSums(PM.n.mat)
      
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
        
        PM.n.sums.mat[e,] <- ((as.vector(colSums(PM.n.mat))/round(area*PM.D.pred, 0)))
        
      } # end e loop
      
      PM.n.md <- apply(PM.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
      PM.n.up <- apply(PM.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      PM.n.lo <- apply(PM.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
      
      PM.s.add <- PM.m.add  <- rep(0, PM.age.max+1)
      for (m in 1:iter) {
        PM.s.add <- rbind(PM.s.add, PM.s.arr[ceiling(PM.gen.l):(t+1),,m])
        PM.m.add <- rbind(PM.m.add, PM.m.arr[ceiling(PM.gen.l):(t+1),,m])
      } # end m loop
      
      PM.s.add <- PM.s.add[-1,]
      PM.m.add <- PM.m.add[-1,]
      
      PM.s.md <- apply(PM.s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
      PM.s.up <- apply(PM.s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      PM.s.lo <- apply(PM.s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
      
      PM.m.md <- apply(PM.m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
      PM.m.up <- apply(PM.m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      PM.m.lo <- apply(PM.m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
      
      
      
      ############################
      ## Palaeoloxodon cypriotes (PC)
      
      # mass
      PC.mass <- 531 # Palaeoloxodon cypriotes (Lomolino et al. 2013 https://doi.org/10.1111/jbi.12096)
      
      ## predicted rm (from Henneman 1983 Oecologia 56:104-108)
      PC.rm.pred <- 10^(0.6914 - (0.2622*log10(PC.mass*1000)))
      PC.lm.pred <- exp(PC.rm.pred)
    
      ## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
      ## log10D = 4.196 − 0.74*(log10m)
      PC.D.pred <- (10^(4.196 - (0.74*log10(PC.mass*1000))))/2 # divided by 2 for females only
    
      ## max age
      ## non-volant birds & mammals (Healy K et al. 2014 PRSB)
      PC.age.max <- round(10^(0.89 + (0.13*log10(PC.mass*1000))), 0)

      ## age vector
      PC.age.vec <- 0:PC.age.max
      
      ## fertility
      ## total fecundity from Allainé et al. 1987 (Oecologia)
      PC.F.pred <- exp(2.719 - (0.211*log(PC.mass*1000)))/2 # divided by 2 for females
    
      ## age at primiparity
      ## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
      PC.alpha <- ceiling(exp(-1.34 + (0.214*log(PC.mass*1000))))
    
      ## define m function with age
      PC.m.vec <- c(rep(0, PC.alpha-1), rep(0.75*PC.F.pred, round(PC.alpha/2,0)), rep(PC.F.pred, (PC.age.max+1-((PC.alpha-1+round(PC.alpha/2,0))))))
      PC.m.sd.vec <- 0.05*PC.m.vec
    
      # fit sigmoidal function
      # logistic power function y = a / (1+(x/b)^c)
      PC.m.dat <- data.frame(PC.age.vec, PC.m.vec)
      param.init <- c(0.1, 2, -5)
      PC.fit.logp <- nls(PC.m.vec ~ a / (1+(PC.age.vec/b)^c), 
                         data = PC.m.dat,
                         algorithm = "port",
                         start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                         trace = F,      
                         nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
      PC.fit.logp.summ <- summary(PC.fit.logp)
      PC.age.vec.cont <- seq(0,max(PC.age.vec),1)
      PC.pred.p.m <- coef(PC.fit.logp)[1] / (1+(PC.age.vec.cont/coef(PC.fit.logp)[2])^coef(PC.fit.logp)[3])
      PC.pred.p.mm <- ifelse(PC.pred.p.m > 1, 1, PC.pred.p.m)
    
      ## survival
      ln.a.s <- -0.5; b.s <- -0.25
      PC.s.tran <- ln.a.s + b.s*log(PC.mass*1000) + log(1)
      PC.s.ad.yr <- exp(-exp(PC.s.tran))
    
      # Siler hazard h(x) (Gurven et al. 2007)
      a1 <- 1 - (1.01*PC.s.ad.yr) # initial infant mortality rate (also known as αt)
      b1 <- 2.7 # rate of mortality decline (also known as bt)
      a2 <- 1 - PC.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
      a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
      b3 <- 0.02 # rate of mortality increase
      longev <- PC.age.max
      x <- seq(0,longev,1) # age vector
      h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
      l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
      init.pop <- 10000
      lx <- round(init.pop*l.x,0)
      len.lx <- length(lx)
      dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
      qx <- dx/lx[1:(length(lx)-1)]
      PC.Sx <- c(0.995*PC.s.ad.yr, 1 - qx)
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
      PC.ssd <- stable.stage.dist(PC.popmat.orig) ## stable stage distribution
      PC.gen.l <- G.val(PC.popmat.orig, PC.age.max) # mean generation length
      
      ## initial population vector
      PC.init.vec <- PC.ssd * round(area*PC.D.pred, 0)
      
      #################
      ## project
      ## set time limit for projection in 1-yr increments
      yr.st <- 1
      yr.end <- round(40*PC.gen.l, 0) # set projection end date
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
      } # end i loop
      
      PC.n.pred <- colSums(PC.n.mat)
      yrs <- seq(yr.st, yr.end, 1)
    
      # compensatory density feedback
      PC.K.max <- 1*round(area*PC.D.pred, 0)
      PC.K.vec <- c(1, PC.K.max/2, 0.75*PC.K.max, PC.K.max) 
      PC.red.vec <- c(1,0.96,0.90,0.829)
      PC.Kred.dat <- data.frame(PC.K.vec, PC.red.vec)
      
      # logistic power function a/(1+(x/b)^c)
      PC.param.init <- c(1, 2*PC.K.max, 2)
      PC.fit.lp <- nls(PC.red.vec ~ a/(1+(PC.K.vec/b)^c), 
                       data = PC.Kred.dat,
                       algorithm = "port",
                       start = c(a = PC.param.init[1], b = PC.param.init[2], c = PC.param.init[3]),
                       trace = F,      
                       nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
      PC.fit.lp.summ <- summary(PC.fit.lp)
      PC.K.vec.cont <- seq(1,2*round(area*PC.D.pred,0),1)
      PC.pred.lp.fx <- coef(PC.fit.lp)[1]/(1+(PC.K.vec.cont/coef(PC.fit.lp)[2])^coef(PC.fit.lp)[3])
    
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
      } # end i loop
      
      PC.n.pred <- colSums(PC.n.mat)
      
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
        
        PC.n.sums.mat[e,] <- ((as.vector(colSums(PC.n.mat))/round(area*PC.D.pred, 0)))
        
      } # end e loop
      
      PC.n.md <- apply(PC.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
      PC.n.up <- apply(PC.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      PC.n.lo <- apply(PC.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
      
      PC.s.add <- PC.m.add <- rep(0, PC.age.max+1)
      for (m in 1:iter) {
        PC.s.add <- rbind(PC.s.add, PC.s.arr[ceiling(PC.gen.l):(t+1),,m])
        PC.m.add <- rbind(PC.m.add, PC.m.arr[ceiling(PC.gen.l):(t+1),,m])
      } # end m loop
      
      PC.s.add <- PC.s.add[-1,]
      PC.m.add <- PC.m.add[-1,]
      
      PC.s.md <- apply(PC.s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
      PC.s.up <- apply(PC.s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      PC.s.lo <- apply(PC.s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
      
      PC.m.md <- apply(PC.m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
      PC.m.up <- apply(PC.m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
      PC.m.lo <- apply(PC.m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
      
    
      
      #####################
      ## meat equivalents
      #####################
      # sexual maturity
      # hippos: in the wild: males 6-13 yrs.; females 7-15 yrs https://ielc.libguides.com/sdzg/factsheets/hippopotamus/reproduction
      PM.mat.ratio <- mean(c(6,13))/mean(c(7,15))
      
      ## growth functions
      ## PM hippos
      PM.f.age.vec <- c(0,PM.alpha,PM.age.max); PM.m.age.vec <- c(0,PM.alpha*PM.mat.ratio,PM.age.max)
      PM2CLmassratio <- PM.mass/mean(c(179,273))
      PM.f.wt.vec <- c(4.5,0.78*179,179)*PM2CLmassratio; PM.m.wt.vec <- c(6.2,0.65*273,273)*PM2CLmassratio
      
      ## Create data frame
      PM.gr.data <- data.frame(PM.f.wt.vec,PM.f.age.vec,PM.m.wt.vec,PM.m.age.vec)
    
      ## Von Bertalanffy growth function: M(t) = Mmax - (Mmax - M0) exp(-kt)
      PM.fit.gr.f <- nls(PM.gr.data$PM.f.wt.vec ~ PM.gr.data$PM.f.wt.vec[3] - (PM.gr.data$PM.f.wt.vec[3] - PM.gr.data$PM.f.wt.vec[1]) * exp(-k.coeff*PM.gr.data$PM.f.age.vec),
                      data = PM.gr.data,
                      start = list(k.coeff = 0.2),
                      trace = F)
      PM.sum.fit.gr.f <- summary(PM.fit.gr.f)
      
      PM.fit.gr.m <- nls(PM.gr.data$PM.m.wt.vec ~ PM.gr.data$PM.m.wt.vec[3] - (PM.gr.data$PM.m.wt.vec[3] - PM.gr.data$PM.m.wt.vec[1]) * exp(-k.coeff*PM.gr.data$PM.m.age.vec),
                      data = PM.gr.data,
                      start = list(k.coeff = 0.24),
                      trace = F)
      PM.sum.fit.gr.m <- summary(PM.fit.gr.m)
      
      ## Coefficients from fit
      PM.coeff.fit.gr.f <- as.numeric(PM.sum.fit.gr.f$parameters)
      PM.coeff.fit.gr.m <- as.numeric(PM.sum.fit.gr.m$parameters)
      
      ## Predict new q vector with integer age inputs
      PM.grf <- PM.gr.data$PM.f.wt.vec[3] - (PM.gr.data$PM.f.wt.vec[3] - PM.gr.data$PM.f.wt.vec[1]) * exp(-PM.coeff.fit.gr.f[1]*PM.age.vec)
      PM.grm <- PM.gr.data$PM.m.wt.vec[3] - (PM.gr.data$PM.m.wt.vec[3] - PM.gr.data$PM.m.wt.vec[1]) * exp(-PM.coeff.fit.gr.m[1]*PM.age.vec)
      
      PM.VB.out <- data.frame(PM.age.vec,PM.grf,PM.grm)
      
      ## PC elephant
      mean.m.La <- mean(c(4000,6300)) # https://ielc.libguides.com/sdzg/factsheets/african_elephant/characteristics
      mean.f.La <- mean(c(2400,3500)) # https://ielc.libguides.com/sdzg/factsheets/african_elephant/characteristics
      PC.Wmax.ratio <- mean.m.La/mean.f.La
      f2avg <- mean.f.La/mean(c(mean.m.La,mean.f.La))
      PC.mass.f <- PC.mass*f2avg
      PC.mass.m <- PC.mass.f*PC.Wmax.ratio
      PC.age.vec <- 0:PC.age.max
      PC.grf <- PC.mass.f*(1 - exp(-0.092*(PC.age.vec+6.15)))^3 # Sukumar, R., Joshi, N.V. & Krishnamurthy, V. (1988). Growth in the Asian elephant. Proceedings of the Indian Academy of Sciences Animal Sciences, 97, 561-571
      PC.grm <- PC.mass.m*(1 - exp(-0.149*(PC.age.vec+3.16)))^3 # Sukumar, R., Joshi, N.V. & Krishnamurthy, V. (1988). Growth in the Asian elephant. Proceedings of the Indian Academy of Sciences Animal Sciences, 97, 561-571
    
      PC.VB.out <- data.frame(0:PC.age.max,PC.grf,PC.grm)
    
      # tissue mass vs. total mass land mammals (Anderson et al. 1979 Scaling of supportive tissue mass
      PC.tiss.m <- (0.033*((PC.grm*1000)^1.090))/1000
      PC.tiss.m.prop <- PC.tiss.m/PC.grm
      
      sdsamp <- 0.09528481 # standard error on edible meat proportion per animal

      # single-sampling test of parameter ranges
      # PM.prob3 <- runif(1, min=min(c(0.7765041-(0.1*0.7765041), 0.7765041+(0.1*0.7765041))),
      #                   max=max(c(0.7765041-(0.1*0.7765041), 0.7765041+(0.1*0.7765041)))) # 1
      # PM.prob3.sd.pc <- runif(1, min=0.025, max=0.1) # 2
      # propdietmeat <- runif(1, min=0.5, max=0.7) # 3
      # f.E.int.m <- runif(1, min=min(c(1877-(0.1*1877), 1877+(0.1*1877))),
      #                    max=max(c(1877-(0.1*1877), 1877+(0.1*1877)))) # 4
      # m.E.int.m <- runif(1, min=min(c(2649-(0.1*2649), 2649+(0.1*2649))),
      #                    max=max(c(2649-(0.1*2649), 2649+(0.1*2649)))) # 5
      # PM.pop.found <- runif(1, min=min(c(round(14280-(0.2*14280),0), round(14280+(0.2*14280),0))),
      #                       max=max(c(round(14280-(0.2*14280),0), round(14280+(0.2*14280),0)))) # 6
      # PC.pop.found <- runif(1, min=min(c(round(5098-(0.2*5098),0), round(5098+(0.2*5098),0))),
      #                       max=max(c(round(5098-(0.2*5098),0), round(5098+(0.2*5098),0)))) # 7
      # FvsMhunt <- runif(1, min=0.4, max=0.6) # 8
      # kCalmeat <- runif(1, min=min(c(130-(0.1*130), 130+(0.1*130))), max=max(c(130-(0.1*130), 130+(0.1*130)))) # 9
      # mpsamp <- runif(1, min=0.25, max=0.4) # 10
      
      
      ## ancient human founding population simulation function
      hipelemeat_sim <- function(input, dir.nm, rowNum) {
        
        ## assign all parameter values
        for (d in 1:ncol(input)) {assign(names(input)[d], input[,d])}
        
        ## energy intake by hunter-gatherers 10.1371/journal.pone.0040503
        #############################################################################################
        ######################## hypercube-sampled variable: f.E.int.m ##############################
        f.E.int.sd.pc <- 0.1939265; f.E.int.sd <- f.E.int.sd.pc*f.E.int.m
        #############################################################################################
        #############################################################################################
        
        #############################################################################################
        ######################## hypercube-sampled variable: m.E.int.m ##############################
        m.E.int.sd.pc <- 0.1491129; m.E.int.sd <- m.E.int.sd.pc*m.E.int.m
        #############################################################################################
        #############################################################################################
        
        # elephant meat
        # average 130 kCal / 100 g (Lupo & Schmitt 2016 J Anthropol Archaeol 44:185–197)
        #############################################################################################
        ######################## hypercube-sampled variable: kCalmeat ###############################
        kCalpgmeat <- kCalmeat/100
        #############################################################################################
        #############################################################################################
        
        gmeatpkCal <- 1/kCalpgmeat
        
        #############################################################################################
        ######################## hypercube-sampled variable: propdietmeat ###########################
        ######################## hypercube-sampled variable: f.E.int.m ##############################
        mkgmeatpday.f <- ((f.E.int.m*propdietmeat)*gmeatpkCal)/1000
        #############################################################################################
        #############################################################################################
        
        mkgmeatpyr.f <- 365*mkgmeatpday.f
        
        #############################################################################################
        ######################## hypercube-sampled variable: f.E.int.m ##############################
        ######################## hypercube-sampled variable: m.E.int.m ##############################
        mkgmeatpyr.m <- m.E.int.m/f.E.int.m * mkgmeatpyr.f
        #############################################################################################
        #############################################################################################
        
        ## translate annual removals into meat
        PM.f.harvWt <- sum(0.5 * 900 * PM.ssd * PM.grf)
        PM.m.harvWt <- sum(0.5 * 900 * PM.ssd * PM.grm)
        
        #############################################################################################
        ######################## hypercube-sampled variable: mpsamp #################################
        PM.tot.meat <- mpsamp*(PM.f.harvWt+PM.m.harvWt)
        #############################################################################################
        #############################################################################################
        
        PM.N.fem.meatEquiv <- PM.tot.meat/mkgmeatpyr.f
        PM.N.mal.meatEquiv <- PM.tot.meat/mkgmeatpyr.m
        PM.N.tot.meatEquiv <- sum(c(PM.N.fem.meatEquiv,PM.N.mal.meatEquiv))
        
        PC.f.harvWt <- sum(0.5 * 350 * PC.ssd * PC.grf)
        PC.m.harvWt <- sum(0.5 * 350 * PC.ssd * PC.grm)
        
        #############################################################################################
        ######################## hypercube-sampled variable: mpsamp #################################
        PC.tot.meat <- mpsamp*(PC.f.harvWt+PC.m.harvWt)
        #############################################################################################
        #############################################################################################
    
        PC.N.fem.meatEquiv <- PC.tot.meat/mkgmeatpyr.f
        PC.N.mal.meatEquiv <- PC.tot.meat/mkgmeatpyr.m
        PC.N.tot.meatEquiv <- sum(c(PC.N.fem.meatEquiv,PC.N.mal.meatEquiv))
        
        
        ##############################################################################
        ## incorporate meat and human diet components to estimate equivalent number ##
        ## of humans required to achieve relative probabilities of extinction       ##
        ##############################################################################
        age.req <- c(2,6,11,16,21)
        rel.prot <- c(9.2,13.5,24,37,40)
        prop.intake.age <- rel.prot/max(rel.prot)
    
        # fit sigmoidal function
        # logistic power function Y=YM*Y0/((YM-Y0)*exp(-k*x) +Y0)
        relintake.dat <- data.frame(age.req, prop.intake.age)
        
        param.init <- c(1.138, 0.1393, 0.1983)
        relintake.fit <- nls(prop.intake.age ~ a*b/((a-b)*exp(-c*age.req) + b), 
                           data = relintake.dat,
                           algorithm = "port",
                           start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                           trace = F,      
                           nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
        
        relintake.fit.logp.summ <- summary(relintake.fit)
        relintake.age.vec.cont <- seq(0,max(age.req),1)
        relintake.pred1 <- coef(relintake.fit.logp.summ)[1]*coef(relintake.fit.logp.summ)[2]/((coef(relintake.fit.logp.summ)[1]-coef(relintake.fit.logp.summ)[2])*exp(-coef(relintake.fit.logp.summ)[3]*relintake.age.vec.cont) + coef(relintake.fit.logp.summ)[2])
        relintake.pred <- ifelse(relintake.pred1 > 1, 1, relintake.pred1)
        relintake.out <- data.frame(relintake.age.vec.cont,relintake.pred)
        
        # starting proportion of total meat intake coming from sources other than hippos/elephants (e.g., seafood, genets)
        meat.other.prop.st <- 0.33
        meat.other.prop.en <- 1
        meat.other.vec <- c(meat.other.prop.st, meat.other.prop.st*1.2*(0.5/meat.other.prop.st),
                            meat.other.prop.st*1.4*(0.5/meat.other.prop.st),meat.other.prop.st*1.5*(0.5/meat.other.prop.st),
                            meat.other.prop.st*1.7*(0.5/meat.other.prop.st),meat.other.prop.st*1.76*(0.5/meat.other.prop.st),
                            meat.other.prop.st*1.9*(0.5/meat.other.prop.st),meat.other.prop.st*1.96*(0.5/meat.other.prop.st),
                            meat.other.prop.en)
        propRmeatOther <- data.frame(pRem=c(1,.9,.8,.7,.5,.4,.2,.1,0.01), pOth=meat.other.vec)
        Prem.vec.cont <- seq(1, 0, -0.01)
        Bcoeffs <- c(-2832,1,3.397,-0.2868)
        Poth.pred <- Bcoeffs[1]+((Bcoeffs[2]-Bcoeffs[1])/(1+exp((Bcoeffs[3]-Prem.vec.cont)/Bcoeffs[4])))
        
        human.ssd <- read.table("ssdHuman.csv", sep=",", header=T)
        relintake.full.prop <- c(relintake.pred,rep(1,59))
        relmeatpyr.age.f <- mkgmeatpyr.f*relintake.full.prop
        relmeatpyr.age.m <- mkgmeatpyr.m*relintake.full.prop
        relintakeMeat.out <- data.frame(relmeatpyr.age.f,relmeatpyr.age.m)
        
        ## prey selection
        # energy / total handling time; Lupo, & Schmitt J Anthropol Archaeol 2016 44:185–197
        #############################################################################################
        ######################## hypercube-sampled variable: PM.prob3 ###############################
        ######################## hypercube-sampled variable: PM.prob3.sd.pc #########################
        PM.prob3.sd <- PM.prob3.sd.pc*PM.prob3
        #############################################################################################
        #############################################################################################
        
        iter <- 1000
        itdiv <- iter/10
        Q.ext.thresh <- 100/2 # quasi-extinction threshold
        
        #############################################################################################
        ######################## hypercube-sampled variable: PM.pop.found ###########################
        PM.init.vec <- PM.ssd * PM.pop.found
        #############################################################################################
        #############################################################################################
    
        #############################################################################################
        ######################## hypercube-sampled variable: PC.pop.found ###########################
        PC.init.vec <- PC.ssd * PC.pop.found
        #############################################################################################
        #############################################################################################
        
        yr.st <- 1
        yr.end <- round(80*PM.gen.l, 0) # set projection end date
        t <- (yr.end - yr.st)
        
        # set mid-range human population and test effect of varying parameter values on extinction risk for both species
        peopleN <- 5000 
        
          PM.n.sums.mat <- PC.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
          
          for (e in 1:iter) {
            PM.popmat <- PM.popmat.orig
            PC.popmat <- PC.popmat.orig
            
            PM.n.mat <- matrix(0, nrow=PM.age.max+1,ncol=(t+1))
            PC.n.mat <- matrix(0, nrow=PC.age.max+1,ncol=(t+1))
            PM.n.mat[,1] <- PM.init.vec
            PC.n.mat[,1] <- PC.init.vec
            
            for (i in 1:t) {
              
              # PHANOURIOS
              # survival
              PM.s.alpha <- estBetaParams(PM.Sx, PM.s.sd.vec^2)$alpha
              PM.s.beta <- estBetaParams(PM.Sx, PM.s.sd.vec^2)$beta
              PM.s.stoch <- rbeta(length(PM.s.alpha), PM.s.alpha, PM.s.beta)
              
              if (rbinom(1, 1, 0.14/PM.gen.l) == 1) { # catastrophe
                cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
                cat.beta <- estBetaParams(0.5, 0.05^2)$beta
                PM.s.stoch <- PM.s.stoch * (rbeta(1, cat.alpha, cat.beta))
              } # end if
              
              # stochastic fertility sampler (Gaussian)
              PM.fert.stch <- rnorm(length(PM.popmat[,1]), PM.pred.p.mm, PM.m.sd.vec)
              
              PM.totN.i <- sum(PM.n.mat[,i], na.rm=T)
              PM.pred.red <- PM.a.lp/(1+(PM.totN.i/PM.b.lp)^PM.c.lp)
              
              diag(PM.popmat[2:(PM.age.max+1),]) <- (PM.s.stoch[-(PM.age.max+1)])*PM.pred.red
              PM.popmat[PM.age.max+1,PM.age.max+1] <- (PM.s.stoch[PM.age.max+1])*PM.pred.red
              PM.popmat[1,] <- ifelse(PM.fert.stch < 0, 0, PM.fert.stch)
              PM.currN <- PM.popmat %*% PM.n.mat[,i]
              PM.n.mat[,i+1] <- PM.currN
              
              # PALAEOLOXODON
              # survival
              PC.s.alpha <- estBetaParams(PC.Sx, PC.s.sd.vec^2)$alpha
              PC.s.beta <- estBetaParams(PC.Sx, PC.s.sd.vec^2)$beta
              PC.s.stoch <- rbeta(length(PC.s.alpha), PC.s.alpha, PC.s.beta)
              
              if (rbinom(1, 1, 0.14/PC.gen.l) == 1) { # catastrophe
                cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
                cat.beta <- estBetaParams(0.5, 0.05^2)$beta
                PC.s.stoch <- PC.s.stoch * (rbeta(1, cat.alpha, cat.beta))
              } # end if
              
              # stochastic fertility sampler (Gaussian)
              PC.fert.stch <- rnorm(length(PC.popmat[,1]), PC.pred.p.mm, PC.m.sd.vec)
              
              PC.totN.i <- sum(PC.n.mat[,i], na.rm=T)
              PC.pred.red <- PC.a.lp/(1+(PC.totN.i/PC.b.lp)^PC.c.lp)
              
              diag(PC.popmat[2:(PC.age.max+1),]) <- (PC.s.stoch[-(PC.age.max+1)])*PC.pred.red
              PC.popmat[PC.age.max+1,PC.age.max+1] <- (PC.s.stoch[PC.age.max+1])*PC.pred.red
              PC.popmat[1,] <- ifelse(PC.fert.stch < 0, 0, PC.fert.stch)
              PC.n.mat[,i+1] <- PC.popmat %*% PC.n.mat[,i]
              
              # HUMANS
              ppl.real <- as.vector(round(peopleN*human.ssd, 0))$x
              f.real <- ifelse(ppl.real < 2, rbinom(length(ppl.real),size=1,prob=0.5), round(ppl.real/2, 0))
              m.real <- ifelse(ppl.real-f.real < 0, 0, ppl.real-f.real)
              
              # meat requirements
              meatpyr.f <- f.real*relmeatpyr.age.f
              meatpyr.m <- m.real*relmeatpyr.age.m
              totmeatpyr <- sum(meatpyr.f + meatpyr.m)
              
              # equivalent meat in PC and PM (if only source of meat)
              #############################################################################################
              ######################## hypercube-sampled variable: PM.prob3 ###############################
              PM.prob3.alpha <- estBetaParams2(mu=round(PM.prob3,3), var=round(PM.prob3.sd^2, 3))$alpha
              PM.prob3.beta <- estBetaParams2(mu=round(PM.prob3,3), var=round(PM.prob3.sd^2, 3))$beta
              #############################################################################################
              #############################################################################################
              
              PM.prob3.stoch <- rbeta(length(PM.prob3.alpha), PM.prob3.alpha, PM.prob3.beta)
              
              PM.meat.tot <- PM.prob3.stoch*totmeatpyr
              PC.meat.tot <- totmeatpyr - PM.meat.tot
              
              PM.meat.f <- PM.prob3.stoch*sum(meatpyr.f)
              PC.meat.f <- sum(meatpyr.f) - PM.meat.f
              PM.meat.m <- PM.prob3.stoch*sum(meatpyr.m)
              PC.meat.m <- sum(meatpyr.m) - PM.meat.m
              
              # meat equivalents in n PM and PC
              # stochastic expression of edible meat return
              #############################################################################################
              ######################## hypercube-sampled variable: mpsamp #################################
              psamp.alpha <- estBetaParams(mpsamp, sdsamp^2)$alpha
              psamp.beta <- estBetaParams(mpsamp, sdsamp^2)$beta
              #############################################################################################
              #############################################################################################
              
              psamp.stoch <- rbeta(length(psamp.alpha), psamp.alpha, psamp.beta)
              
              # relative increase in 'meat.other.prop' with declining animal densities
              PM.pInit.remain <- sum(PM.n.mat[,i])/sum(PM.n.mat[,1])
              PC.pInit.remain <- sum(PC.n.mat[,i])/sum(PC.n.mat[,1])
              
              PM.meat.oth.t1 <- Bcoeffs[1]+((Bcoeffs[2]-Bcoeffs[1])/(1+exp((Bcoeffs[3]-PM.pInit.remain)/Bcoeffs[4])))
              PC.meat.oth.t1 <- Bcoeffs[1]+((Bcoeffs[2]-Bcoeffs[1])/(1+exp((Bcoeffs[3]-PC.pInit.remain)/Bcoeffs[4])))
              
              PM.meat.oth.t2 <- ifelse(PM.meat.oth.t1 > meat.other.prop.en, 1, PM.meat.oth.t1)
              PM.meat.oth.t <- ifelse(PM.meat.oth.t2 < meat.other.prop.st, 1, PM.meat.oth.t2)
              PC.meat.oth.t2 <- ifelse(PC.meat.oth.t1 > meat.other.prop.en, 1, PC.meat.oth.t1)
              PC.meat.oth.t <- ifelse(PC.meat.oth.t2 < meat.other.prop.st, 1, PC.meat.oth.t2)
              
              # correct for relative 'meat other' proportion
              PM.totf.wt <- PM.meat.f/psamp.stoch*(1-PM.meat.oth.t)
              PC.totf.wt <- PC.meat.f/psamp.stoch*(1-PC.meat.oth.t)
        
              PM.totm.wt <- PM.meat.m/psamp.stoch*(1-PM.meat.oth.t)
              PC.totm.wt <- PC.meat.m/psamp.stoch*(1-PC.meat.oth.t)
              
              # assume FvsMhunt females harvested vs. males
              # PM - female
              PMwtvecf <- (1000*PM.ssd*PM.grf)
              PMwtvecfwt <- PMwtvecf/sum(PMwtvecf)
              
              #############################################################################################
              ######################## hypercube-sampled variable: FvsMhunt ###############################
              PMf.wtXage <- FvsMhunt*(PM.totf.wt)*PMwtvecfwt
              #############################################################################################
              #############################################################################################
              
              PMf.NxAge <- as.vector(round(PMf.wtXage/PM.grf, 0))
              
              # PM - male
              PMwtvecm <- (1000*PM.ssd*PM.grm)
              PMwtvecmwt <- PMwtvecm/sum(PMwtvecm)
              
              #############################################################################################
              ######################## hypercube-sampled variable: FvsMhunt ###############################
              PMm.wtXage <- (1-FvsMhunt)*(PM.totm.wt)*PMwtvecmwt
              #############################################################################################
              #############################################################################################
              
              PMm.NxAge <- as.vector(round(PMm.wtXage/PM.grm, 0))
              
              # PC - female
              PCwtvecf <- (1000*PC.ssd*PC.grf)
              PCwtvecfwt <- PCwtvecf/sum(PCwtvecf)
              
              #############################################################################################
              ######################## hypercube-sampled variable: FvsMhunt ###############################
              PCf.wtXage <- FvsMhunt*(PC.totf.wt)*PCwtvecfwt
              #############################################################################################
              #############################################################################################
              
              PCf.NxAge <- as.vector(round(PCf.wtXage/PC.grf, 0))
              
              # PC - male
              PCwtvecm <- (1000*PC.ssd*PC.grm)
              PCwtvecmwt <- PCwtvecm/sum(PCwtvecm)
              
              #############################################################################################
              ######################## hypercube-sampled variable: FvsMhunt ###############################
              PCm.wtXage <- (1-FvsMhunt)*(PC.totm.wt)*PCwtvecmwt
              #############################################################################################
              #############################################################################################
              
              PCm.NxAge <- as.vector(round(PCm.wtXage/PC.grm, 0))
              
              # if not enough PM meat to remove relative to size of extant PM population,
              # transfer meat equivalent to PC population
              PCf.NxAgeADD <- rep(0,length(PCwtvecfwt))
              if ((sum(PM.currN) - sum(PMf.NxAge)) < 0) {
                PMfn2twt <- sum(PMf.NxAge*PM.grf) # transform extra female PMs to total female weight
                PMfn2twtCor <- PMfn2twt*(1-PC.meat.oth.t) # correct for meat 'other' proportion
                PCf.wtXageADD <- (PMfn2twtCor)*PCwtvecfwt # weight by PC age
                PCf.NxAgeADD <- as.vector(round(PCf.wtXageADD/PC.grf, 0)) # n PC by age
               } # end if
              
              # REMOVALS FROM n MAT - PM
              PM.n.mat[,i+1] <- ifelse((PM.n.mat[,i+1] - PMf.NxAge) < 0, 0, (PM.n.mat[,i+1] - PMf.NxAge))
        
              # REMOVALS FROM n MAT - PC
              PC.n.mat[,i+1] <- ifelse((PC.n.mat[,i+1] - (PCf.NxAge+PCf.NxAgeADD)) < 0, 0, (PC.n.mat[,i+1] - (PCf.NxAge+PCf.NxAgeADD)))
              
            } # end i loop
            
            PM.n.sums.mat[e,] <- (as.vector(colSums(PM.n.mat)))
            PC.n.sums.mat[e,] <- (as.vector(colSums(PC.n.mat)))
            
          } # end e loop
          
          # quasi-extinction probability
          PM.Q.ext.mat <- ifelse(PM.n.sums.mat < Q.ext.thresh, 1, 0)
          PM.Q.ext.sum <- apply(PM.Q.ext.mat[,ceiling(PM.gen.l):dim(PM.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
    
          PC.Q.ext.mat <- ifelse(PC.n.sums.mat < Q.ext.thresh, 1, 0)
          PC.Q.ext.sum <- apply(PC.Q.ext.mat[,ceiling(PC.gen.l):dim(PC.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
    
          # save
          input$PM.Qext.pr <- length(which(PM.Q.ext.sum > 0)) / iter
          input$PC.Qext.pr <- length(which(PC.Q.ext.sum > 0)) / iter
          save.nm <- paste0('res',sprintf("%09.0f", rowNum))
          assign(save.nm, input)
          save(list=save.nm,file=paste(dir.nm,save.nm,sep='/'))
          
} # end d lh loop


## parameter ranges
ranges <- list()
ranges$PM.prob3 <- c(0.7765041-(0.1*0.7765041), 0.7765041+(0.1*0.7765041)) # 1
ranges$PM.prob3.sd.pc <- c(0.025, 0.1) # 2
ranges$propdietmeat <- c(0.5, 0.7) # 3
ranges$f.E.int.m <- c(1877-(0.1*1877), 1877+(0.1*1877)) # 4
ranges$m.E.int.m <- c(2649-(0.1*2649), 2649+(0.1*2649)) # 5
ranges$PM.pop.found <- c(round(14280-(0.2*14280),0), round(14280+(0.2*14280),0)) # 6
ranges$PC.pop.found <- c(round(5098-(0.2*5098),0), round(5098+(0.2*5098),0)) # 7
ranges$FvsMhunt <- c(0.4, 0.6) # 8
ranges$kCalmeat <- c(130-(0.1*130), 130+(0.1*130)) # 9
ranges$mpsamp <- c(0.25, 0.4) # 10

## create hypercube
nSamples <- 1000
lh <- data.frame(randomLHS(n=nSamples, k=length(ranges)))
names(lh) <- names(ranges)

## convert parameters to required scale
for (j in 1:ncol(lh)) {
  par <- names(lh)[j]
  lh[,par] <- qunif(lh[,j], min=ranges[[par]][1], max=ranges[[par]][2]) ## continuous
}

## number of iterations for each parameter set
lh$iter <- 1

## folder for saving the results of each row
## we could just store in memory, but then if something breaks we will lose the lot
dir.nm <- 'GSAhipeleharv'

if (dir.exists(dir.nm) == F) {
  dir.create(dir.nm)
}

## run in parallel
res <- foreach(rowNum=1:nrow(lh),.verbose=T) %do% {hipelemeat_sim(input=lh[rowNum,],dir.nm=dir.nm,rowNum=rowNum)}

## retrieve results
res.nms <- list.files(dir.nm)
res.list <- lapply(res.nms, function(x) {load(paste(dir.nm,x,sep='/')) ; print(x) ; return(eval(as.name(x)))})
dat <- rbindlist(res.list)
head(dat)
dim(dat)[1]
tail(dat)
sum(is.na(dat$PM.Qext.pr))
sum(is.na(dat$PC.Qext.pr))



######################################
## Boosted Regression Tree Emulator ##
######################################
## PM
dat.nona <- data.frame(na.omit(dat[!is.infinite(rowSums(dat)),]))
dim(dat.nona)[1]
brt.fit <- gbm.step(dat.nona, gbm.x = attr(dat.nona, "names")[1:10], gbm.y = attr(dat.nona, "names")[12], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.008, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit)
dim(dat.nona)[1]
D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)

CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
CV.cor
CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
CV.cor.se
print(c(CV.cor, CV.cor.se))

## PC
dat.nona <- data.frame(na.omit(dat[!is.infinite(rowSums(dat)),]))
dim(dat.nona)[1]
brt.fit <- gbm.step(dat.nona, gbm.x = attr(dat.nona, "names")[1:10], gbm.y = attr(dat.nona, "names")[13], family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.008, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit)
dim(dat.nona)[1]
D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)

CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
CV.cor
CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
CV.cor.se
print(c(CV.cor, CV.cor.se))

  
save.image("~/Documents/Papers/Palaeo/Cyprus/data/palaeo/pplExtGSA.RData")

