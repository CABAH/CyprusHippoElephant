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

################################################################################
## translating human dietary requirements to meat equivalents to estimate the ##
## probability of driving both megafauna species to extinction                ##
## (requires running 'base' models first')
################################################################################

# pygmy hippo weight: weighs 180–275 kg
# male pygmy hippos have a range of 1.5 km2, while a female ranges from 0.4 to 0.6 km2.
# at birth, pygmy hippos weigh 4.5-6.2 kg
# During the Pleistocene, Cyprus was inhabited by Phanourios minor, the smallest Mediterranean dwarf
# hippopotamus known, featuring a body mass of approximately 130 kg (Lomolino et al., 2013). The
# strongly endemic and impoverished environment of P. minor was characterised by few competition
# and a lack of predation, since larger carnivores were absent (Marra, 2005a)
# https://www.hippohaven.com/pygmy-hippopotamus/
#
# sexual maturity
# hippos: in the wild: males 6-13 yrs.; females 7-15 yrs https://ielc.libguides.com/sdzg/factsheets/hippopotamus/reproduction
PM.mat.ratio <- mean(c(6,13))/mean(c(7,15))

## growth functions
## PM hippos
PM.f.age.vec <- c(0,PM.alpha,PM.age.max); PM.m.age.vec <- c(0,PM.alpha*PM.mat.ratio,PM.age.max)
# sexual maturity 0.78 of mass (F); 0.65 of mass (M) https://doi.org/10.3957/056.053.0103
# ratio of P. minor to C. liberiensis (extant pygmy hippo)
PM2CLmassratio <- PM.mass/mean(c(179,273))

PM.f.wt.vec <- c(4.5,0.78*179,179)*PM2CLmassratio; PM.m.wt.vec <- c(6.2,0.65*273,273)*PM2CLmassratio

## Create data frame
PM.gr.data <- data.frame(PM.f.wt.vec,PM.f.age.vec,PM.m.wt.vec,PM.m.age.vec)
PM.gr.data

## model formula
## Von Bertalanffy growth function: M(t) = Mmax - (Mmax - M0) exp(-kt)
## M0 = mean birth mass
## Mmax = mean maximum mass
## k = rate constant per year

PM.fit.gr.f <- nls(PM.gr.data$PM.f.wt.vec ~ PM.gr.data$PM.f.wt.vec[3] - (PM.gr.data$PM.f.wt.vec[3] - PM.gr.data$PM.f.wt.vec[1]) * exp(-k.coeff*PM.gr.data$PM.f.age.vec),
                data = PM.gr.data,
                start = list(k.coeff = 0.2),
                trace = TRUE)
PM.sum.fit.gr.f <- summary(PM.fit.gr.f)

PM.fit.gr.m <- nls(PM.gr.data$PM.m.wt.vec ~ PM.gr.data$PM.m.wt.vec[3] - (PM.gr.data$PM.m.wt.vec[3] - PM.gr.data$PM.m.wt.vec[1]) * exp(-k.coeff*PM.gr.data$PM.m.age.vec),
                data = PM.gr.data,
                start = list(k.coeff = 0.24),
                trace = TRUE)
PM.sum.fit.gr.m <- summary(PM.fit.gr.m)

## Coefficients from fit
PM.coeff.fit.gr.f <- as.numeric(PM.sum.fit.gr.f$parameters)
PM.coeff.fit.gr.m <- as.numeric(PM.sum.fit.gr.m$parameters)


## Predict new q vector with integer age inputs
PM.grf <- PM.gr.data$PM.f.wt.vec[3] - (PM.gr.data$PM.f.wt.vec[3] - PM.gr.data$PM.f.wt.vec[1]) * exp(-PM.coeff.fit.gr.f[1]*PM.age.vec)
PM.grm <- PM.gr.data$PM.m.wt.vec[3] - (PM.gr.data$PM.m.wt.vec[3] - PM.gr.data$PM.m.wt.vec[1]) * exp(-PM.coeff.fit.gr.m[1]*PM.age.vec)

row <- 1
col <- 2
par(mfrow=c(row,col))
plot(PM.age.vec,PM.grf,xlab="age (years)",ylab="mass (kg)",ylim=range(c(0,max(m.wt.vec))),type="l")
lines(PM.age.vec,PM.grm,col="red")
title(main="hippo")

PM.VB.out <- data.frame(PM.age.vec,PM.grf,PM.grm)

## PC elephant
# growth constant K in P. falconeri is far lower (0.055–0.079) than that of
# L. africana (0.1666) https://doi.org/10.1038%2Fs41598-021-02192-4
# female growth W = Wmax*(1 - exp(-0.092*(x+6.15)))^3
# male growth  W = Wmax*(1 - exp(-0.149*(x+3.16)))^3
mean.m.La <- mean(c(4000,6300)) # https://ielc.libguides.com/sdzg/factsheets/african_elephant/characteristics
mean.f.La <- mean(c(2400,3500)) # https://ielc.libguides.com/sdzg/factsheets/african_elephant/characteristics
PC.Wmax.ratio <- mean.m.La/mean.f.La
f2avg <- mean.f.La/mean(c(mean.m.La,mean.f.La))
PC.mass.f <- PC.mass*f2avg
PC.mass.m <- PC.mass.f*PC.Wmax.ratio
PC.age.vec <- 0:PC.age.max
PC.grf <- PC.mass.f*(1 - exp(-0.092*(PC.age.vec+6.15)))^3 # Sukumar, R., Joshi, N.V. & Krishnamurthy, V. (1988). Growth in the Asian elephant. Proceedings of the Indian Academy of Sciences Animal Sciences, 97, 561-571
plot(0:PC.age.max, PC.grf, type="l",xlab="age (years)",ylab="mass (kg)",ylim=range(c(0,max(PC.mass.m))))
PC.grm <- PC.mass.m*(1 - exp(-0.149*(PC.age.vec+3.16)))^3 # Sukumar, R., Joshi, N.V. & Krishnamurthy, V. (1988). Growth in the Asian elephant. Proceedings of the Indian Academy of Sciences Animal Sciences, 97, 561-571
lines(PC.age.vec,PC.grm,col="red")
title(main="elephant")
par(mfrow=c(1,1))
PC.VB.out <- data.frame(0:PC.age.max,PC.grf,PC.grm)

## ungulate edible weights https://www.gov.nt.ca/sites/ecc/files/weights_of_wildlife.pdf
edw.bgcaribou <- c(36.4,45,45.45,58.2,36,55,48,41,37,37,48,48,48,45,45,90,36,50,45,37,36,33,29)
edw.wcaribou <- c(77,61.8,68,95,61.8,50)
edw.moose <- c(159.1,160,199,204.5,199,160,199,199,205,227,199,180,140)
edw.wbison <- c(272.7,273,250,272.7,409)
edw.muskox <- c(137.5,136,110,110,110,95,100,69)

m.bgcaribou <- c(90,135,159,182) # https://www.registrelep-sararegistry.gc.ca/virtual_sara/files/cosewic/sr_Caribou%20Barren-ground_2016_e.pdf; https://www.adfg.alaska.gov/index.cfm?adfg=caribou.main
mm.bgcaribou <- mean(m.bgcaribou)
msd.bgcaribou <- sd(m.bgcaribou)
edpropm.bgcaribou <- mean(edw.bgcaribou)/mm.bgcaribou
edpropsd.bgcaribou <- sd(edw.bgcaribou)/mm.bgcaribou

m.wcaribou <- c(110,150,160,210) # https://www.fws.gov/species/woodland-caribou-rangifer-tarandus-caribou#:~:text=Height%20of%20the%20woodland%20caribou,(160%20to%20210%20kg).
mm.wcaribou <- mean(m.wcaribou)
msd.wcaribou <- sd(m.wcaribou)
edpropm.wcaribou <- mean(edw.wcaribou)/mm.wcaribou
edpropsd.wcaribou <- sd(edw.wcaribou)/mm.wcaribou

m.moose <- c(2000*0.453592, 1200*0.453592) # https://www.adfg.alaska.gov/index.cfm?adfg=moose.main
mm.moose <- mean(m.moose)
msd.moose <- sd(m.moose)
edpropm.moose <- mean(edw.moose)/mm.moose
edpropsd.moose <- sd(edw.moose)/mm.moose

m.muskox <- c(600*0.453592, 800*0.453592, 400*0.453592, 500*0.453592) # https://www.adfg.alaska.gov/index.cfm?adfg=muskox.printerfriendly
mm.muskox <- mean(m.muskox)
msd.muskox <- sd(m.muskox)
edpropm.muskox <- mean(edw.muskox)/mm.muskox
edpropsd.muskox <- sd(edw.muskox)/mm.muskox

m.mass.vec <- c(mm.bgcaribou,mm.wcaribou,mm.moose,mm.muskox)
edpropm.vec <- c(edpropm.bgcaribou,edpropm.wcaribou,edpropm.moose,edpropm.muskox)
plot(m.mass.vec,edpropm.vec)

editer <- 10000
psamp.bgcaribou <- psamp.wcaribou <- psamp.moose <- psamp.muskox <- rep(NA,editer)
for (e in 1:editer) {
  psamp.bgcaribou[e] <- rnorm(1,mean(edw.bgcaribou), sd(edw.bgcaribou)) / rnorm(1,mm.bgcaribou, msd.bgcaribou)
  psamp.wcaribou[e] <- rnorm(1,mean(edw.wcaribou), sd(edw.wcaribou)) / rnorm(1,mm.wcaribou, msd.wcaribou)
  psamp.moose[e] <- rnorm(1,mean(edw.moose), sd(edw.moose)) / rnorm(1,mm.moose, msd.moose)
  psamp.muskox[e] <- rnorm(1,mean(edw.muskox), sd(edw.muskox)) / rnorm(1,mm.muskox, msd.muskox)
}
psamp1.vec <- c(psamp.bgcaribou,psamp.wcaribou,psamp.moose,psamp.muskox)
psamp.vec <- psamp1.vec[which(psamp1.vec < 0.5 & psamp1.vec > 0)]
range(psamp.vec)
mpsamp <- median(psamp.vec)
sdsamp <- sd(psamp.vec)
uppsamp <- quantile(psamp.vec, probs=0.975)
lopsamp <- quantile(psamp.vec, probs=0.025)
mpsamp
sdsamp
uppsamp
lopsamp

## energy intake by hunter-gatherers 10.1371/journal.pone.0040503
# female: 1877 +/- 364 kcal/day; male: 2649 +/- 395 kcal/day
f.E.int.m <- 1877
f.E.int.sd <- 364
m.E.int.m <- 2649
m.E.int.sd <- 395

# proportion diet meat 10.1038/sj.ejcn=1601353
propdietmeat <- 0.65

# elephant meat
# average 130 kCal / 100 g (Lupo & Schmitt 2016 J Anthropol Archaeol 44:185–197; )
kCalpgmeat <- 130/100
gmeatpkCal <- 1/kCalpgmeat
mkgmeatpday.f <- ((f.E.int.m*propdietmeat)*gmeatpkCal)/1000
mkgmeatpyr.f <- 365*mkgmeatpday.f
mkgmeatpyr.m <- m.E.int.m/f.E.int.m * mkgmeatpyr.f

## translate annual removals into meat
PM.f.harvWt <- sum(0.5 * ind.rem.vec[19] * PM.ssd * PM.grf)
PM.m.harvWt <- sum(0.5 * ind.rem.vec[19] * PM.ssd * PM.grm)
PM.tot.meat <- mpsamp*(PM.f.harvWt+PM.m.harvWt)
PM.N.fem.meatEquiv <- PM.tot.meat/mkgmeatpyr.f
PM.N.mal.meatEquiv <- PM.tot.meat/mkgmeatpyr.m
PM.N.tot.meatEquiv <- sum(c(PM.N.fem.meatEquiv,PM.N.mal.meatEquiv))

PC.f.harvWt <- sum(0.5 * ind.rem.vec[8] * PC.ssd * PC.grf)
PC.m.harvWt <- sum(0.5 * ind.rem.vec[8] * PC.ssd * PC.grm)
PC.tot.meat <- mpsamp*(PC.f.harvWt+PC.m.harvWt)
PC.N.fem.meatEquiv <- PC.tot.meat/mkgmeatpyr.f
PC.N.mal.meatEquiv <- PC.tot.meat/mkgmeatpyr.m
PC.N.tot.meatEquiv <- sum(c(PC.N.fem.meatEquiv,PC.N.mal.meatEquiv))


##############################################################################
## incorporate meat and human diet components to estimate equivalent number ##
## of humans required to achieve relative probabilities of extinction       ##
##############################################################################

# Kyriakou et al. 2014
# relative protein requirements for people 
# adult 50 kg: 40 g protein/day
# child 1-3 yrs: 9.2 g/day
# child 4-8 yrs: 13.5 g/day
# child 9-13yrs: 24 g/day
# child 14-18 yrs: 37 g/day

age.req <- c(2,6,11,16,21)
rel.prot <- c(9.2,13.5,24,37,40)
prop.intake.age <- rel.prot/max(rel.prot)
plot(age.req, prop.intake.age, type="b", pch=19, xlab="age", ylab="proportion total adult intake")

# fit sigmoidal function
# logistic power function Y=YM*Y0/((YM-Y0)*exp(-k*x) +Y0)
relintake.dat <- data.frame(age.req, prop.intake.age)

param.init <- c(1.138, 0.1393, 0.1983)
relintake.fit <- nls(prop.intake.age ~ a*b/((a-b)*exp(-c*age.req) + b), 
                   data = relintake.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))

relintake.fit.logp.summ <- summary(relintake.fit)
plot(age.req, prop.intake.age, type="b", pch=19, xlab="age (yrs)", ylab="prop tot adult intake", xlim=c(0,max(age.req)), ylim=c(0,1))
relintake.age.vec.cont <- seq(0,max(age.req),1)
relintake.pred1 <- coef(relintake.fit.logp.summ)[1]*coef(relintake.fit.logp.summ)[2]/((coef(relintake.fit.logp.summ)[1]-coef(relintake.fit.logp.summ)[2])*exp(-coef(relintake.fit.logp.summ)[3]*relintake.age.vec.cont) + coef(relintake.fit.logp.summ)[2])
relintake.pred <- ifelse(relintake.pred1 > 1, 1, relintake.pred1)
lines(relintake.age.vec.cont, relintake.pred,lty=2,lwd=3,col="red")

relintake.out <- data.frame(relintake.age.vec.cont,relintake.pred)


##################################################################################
## PROGRESSIVELY INCREASE NUMBER OF PEOPLE & TRANSLATE TO ANIMALS KILLED        ##
##################################################################################
yr.st <- 1
#************************
yr.end <- round(80*PM.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

iter <- 1000
people.vec <- seq(3000,12000,250)

# starting proportion of total meat intake coming from sources other than hippos/elephants (e.g., seafood, genets)
meat.other.prop.st <- 0.33
meat.other.prop.en <- 1
meat.other.vec <- c(meat.other.prop.st, meat.other.prop.st*1.2*(0.5/meat.other.prop.st),
                    meat.other.prop.st*1.4*(0.5/meat.other.prop.st),meat.other.prop.st*1.5*(0.5/meat.other.prop.st),
                    meat.other.prop.st*1.7*(0.5/meat.other.prop.st),meat.other.prop.st*1.76*(0.5/meat.other.prop.st),
                    meat.other.prop.st*1.9*(0.5/meat.other.prop.st),meat.other.prop.st*1.96*(0.5/meat.other.prop.st),
                    meat.other.prop.en)
propRmeatOther <- data.frame(pRem=c(1,.9,.8,.7,.5,.4,.2,.1,0.01), pOth=meat.other.vec)
plot(propRmeatOther$pRem, propRmeatOther$pOth, pch=19, ylab="prop 'other' meat source", xlab="prop pop remaining")
Prem.vec.cont <- seq(1, 0, -0.01)
Bcoeffs <- c(-2832,1,3.397,-0.2868)
Poth.pred <- Bcoeffs[1]+((Bcoeffs[2]-Bcoeffs[1])/(1+exp((Bcoeffs[3]-Prem.vec.cont)/Bcoeffs[4])))
lines(Prem.vec.cont, Poth.pred, lty=3,lwd=3,col="red")

human.ssd <- read.table("ssdHuman.csv", sep=",", header=T)
relintake.full.prop <- c(relintake.pred,rep(1,59))
relmeatpyr.age.f <- mkgmeatpyr.f*relintake.full.prop
relmeatpyr.age.m <- mkgmeatpyr.m*relintake.full.prop
relintakeMeat.out <- data.frame(relmeatpyr.age.f,relmeatpyr.age.m)

## prey selection
# energetic payoff (e) (equations from Yaworsky et al. 2023-Sci Rep)
PM.e.md <- 43551.8 + PM.mass*546.5
PM.e.lo <- 43551.8-(1.96*1856.3) + PM.mass*(546.5-(1.96*15.6))
PM.e.up <- 43551.8+(1.96*1856.3) + PM.mass*(546.5+(1.96*15.6))
PM.e.sd <- mean(1/1.96*(PM.e.up-PM.e.md), 1/1.96*(PM.e.md-PM.e.lo))

PC.e.md <- 43551.8 + PC.mass*546.5
PC.e.lo <- 43551.8-(1.96*1856.3) + PC.mass*(546.5-(1.96*15.6))
PC.e.up <- 43551.8+(1.96*1856.3) + PC.mass*(546.5+(1.96*15.6))
PC.e.sd <- mean(1/1.96*(PC.e.up-PC.e.md), 1/1.96*(PC.e.md-PC.e.lo))

# pre-acquisition handling costs (c) (equations from Yaworsky et al. 2023-Sci Rep)
PM.c.md <- 514.1590 + PM.mass*0.6252
PM.c.lo <- 514.1590-(1.96*393.8621) + PM.mass*(0.6252-(1.96*0.3309))
PM.c.up <- 514.1590+(1.96*393.8621) + PM.mass*(0.6252+(1.96*0.3309))
PM.c.lo <- ifelse(PM.c.lo < 0, 0, PM.c.lo)
PM.c.sd <- mean(1/1.96*(PM.c.up-PM.c.md), 1/1.96*(PM.c.md-PM.c.lo))

PC.c.md <- 514.1590 + PC.mass*0.6252
PC.c.lo <- 514.1590-(1.96*393.8621) + PC.mass*(0.6252-(1.96*0.3309))
PC.c.up <- 514.1590+(1.96*393.8621) + PC.mass*(0.6252+(1.96*0.3309))
PC.c.lo <- ifelse(PC.c.lo < 0, 0, PC.c.lo)
PC.c.sd <- mean(1/1.96*(PC.c.up-PC.c.md), 1/1.96*(PC.c.md-PC.c.lo))

# post-acquisition handling costs (h) (equations from Yaworsky et al. 2023-Sci Rep)
PM.h.md <- -897.6988 + PM.mass*10.5013
PM.h.lo <- -897.6988-(1.96*297.4462) + PM.mass*(10.5013-(1.96*0.2499))
PM.h.up <- -897.6988+(1.96*297.4462) + PM.mass*(10.5013+(1.96*0.2499))
PM.h.lo <- ifelse(PM.h.lo < 0, 0, PM.h.lo)
PM.h.sd <- mean(1/1.96*(PM.h.up-PM.h.md), 1/1.96*(PM.h.md-PM.h.lo))

PC.h.md <- -897.6988 + PC.mass*10.5013
PC.h.lo <- -897.6988-(1.96*297.4462) + PC.mass*(10.5013-(1.96*0.2499))
PC.h.up <- -897.6988+(1.96*297.4462) + PC.mass*(10.5013+(1.96*0.2499))
PC.h.lo <- ifelse(PC.h.lo < 0, 0, PC.h.lo)
PC.h.sd <- mean(1/1.96*(PC.h.up-PC.h.md), 1/1.96*(PC.h.md-PC.h.lo))

# proportion failed pursuits (p) (equations from Yaworsky et al. 2023-Sci Rep)
PM.p.md <- 0.1753815 + PM.mass*0.0002839
PM.p.lo <- 0.1753815-(1.96*0.2221416) + PC.mass*(0.0002839-(1.96*0.0002211))
PM.p.up <- 0.1753815+(1.96*0.2221416) + PC.mass*(0.0002839+(1.96*0.0002211))
PM.p.lo <- ifelse(PM.p.lo < 0, 0, PM.p.lo)
PM.p.sd <- mean(1/1.96*(PM.p.up-PM.p.md), 1/1.96*(PM.p.md-PM.p.lo))

PC.p.md <- 0.1753815 + PC.mass*0.0002839
PC.p.lo <- 0.1753815-(1.96*0.2221416) + PC.mass*(0.0002839-(1.96*0.0002211))
PC.p.up <- 0.1753815+(1.96*0.2221416) + PC.mass*(0.0002839+(1.96*0.0002211))
PC.p.lo <- ifelse(PC.p.lo < 0, 0, PC.p.lo)
PC.p.sd <- mean(1/1.96*(PC.p.up-PC.p.md), 1/1.96*(PC.p.md-PC.p.lo))

# post-encounter return rate (π) (equations from Yaworsky et al. 2023-Sci Rep)
PM.pi.md <- ((PM.e.md*(1-PM.p.md))/(((1-PM.p.md)*PM.h.md)+PM.c.md))*60
PC.pi.md <- ((PC.e.md*(1-PC.p.md))/(((1-PC.p.md)*PC.h.md)+PC.c.md))*60

perriter <- 100000
PM.prob.vec <- choice.vec <- rep(NA,perriter)
for (p in 1:perriter) {
  PM.e.it <- rtruncnorm(1, a=0, b=Inf, mean=PM.e.md, sd=PM.e.sd)
    PM.p.alpha <- estBetaParams(PM.p.md, PM.p.sd^2)$alpha
    PM.p.beta <- estBetaParams(PM.p.md, PM.p.sd^2)$beta
  PM.p.it <- rbeta(1, PM.p.alpha, PM.p.beta)
  PM.h.it <- rtruncnorm(1, a=0, b=Inf, mean=PM.h.md, sd=PM.h.sd)
  PM.c.it <- rtruncnorm(1, a=0, b=Inf, mean=PM.c.md, sd=PM.c.sd)
  
  PM.pi.it <- ((PM.e.it*(1-PM.p.it))/(((1-PM.p.it)*PM.h.it)+PM.c.it))*60
  PM.pi.it <- ifelse(PM.pi.it == 0, NA, PM.pi.it)
  
  PC.e.it <- rtruncnorm(1, a=0, b=Inf, mean=PC.e.md, sd=PC.e.sd)
  PC.p.alpha <- estBetaParams(PC.p.md, PC.p.sd^2)$alpha
  PC.p.beta <- estBetaParams(PC.p.md, PC.p.sd^2)$beta
  PC.p.it <- rbeta(1, PC.p.alpha, PC.p.beta)
  PC.h.it <- rtruncnorm(1, a=0, b=Inf, mean=PC.h.md, sd=PC.h.sd)
  PC.c.it <- rtruncnorm(1, a=0, b=Inf, mean=PC.c.md, sd=PC.c.sd)
  
  PC.pi.it <- ((PC.e.it*(1-PC.p.it))/(((1-PC.p.it)*PC.h.it)+PC.c.it))*60
  PC.pi.it <- ifelse(PC.pi.it == 0, NA, PC.pi.it)

  choice.vec[p] <- ifelse(PM.pi.it > PC.pi.it, 1, 0)
}

PM.prob3 <- as.vector(sum(choice.vec, na.rm=T)/length(which(is.na(choice.vec)==F)))
PM.prob3.sd <- 0.05*PM.prob3

PM.Q.ext.pr <- PC.Q.ext.pr <- PM.yrsRem.md <- PM.yrsRem.up <- PM.yrsRem.lo <- 
  PC.yrsRem.md <- PC.yrsRem.up <- PC.yrsRem.lo <- rep(NA,length(people.vec))

for (k in 1:length(people.vec)) {
  
  PM.n.sums.mat <- PM.n.remsums.mat <- PC.n.sums.mat <- PC.n.remsums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
  
  for (e in 1:iter) {
    PM.popmat <- PM.popmat.orig
    PC.popmat <- PC.popmat.orig
    
    PM.n.mat <- PM.n.rem <- matrix(0, nrow=PM.age.max+1,ncol=(t+1))
    PC.n.mat <- PC.n.rem <- matrix(0, nrow=PC.age.max+1,ncol=(t+1))
    PM.n.mat[,1] <- PM.init.vec
    PC.n.mat[,1] <- PC.init.vec
    PM.n.rem[,1] <- 0
    PC.n.rem[,1] <- 0
    
    for (i in 1:t) {
      
      # PHANOURIOS
      # stochastic survival values
      PM.s.alpha <- estBetaParams(PM.Sx, PM.s.sd.vec^2)$alpha
      PM.s.beta <- estBetaParams(PM.Sx, PM.s.sd.vec^2)$beta
      PM.s.stoch <- rbeta(length(PM.s.alpha), PM.s.alpha, PM.s.beta)
      
      if (rbinom(1, 1, 0.14/PM.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        PM.s.stoch <- PM.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertility sampler (gaussian)
      PM.fert.stch <- rnorm(length(PM.popmat[,1]), PM.pred.p.mm, PM.m.sd.vec)
      
      PM.totN.i <- sum(PM.n.mat[,i], na.rm=T)
      PM.pred.red <- PM.a.lp/(1+(PM.totN.i/PM.b.lp)^PM.c.lp)
      
      diag(PM.popmat[2:(PM.age.max+1),]) <- (PM.s.stoch[-(PM.age.max+1)])*PM.pred.red
      PM.popmat[PM.age.max+1,PM.age.max+1] <- (PM.s.stoch[PM.age.max+1])*PM.pred.red
      PM.popmat[1,] <- ifelse(PM.fert.stch < 0, 0, PM.fert.stch)
      PM.n.mat[,i+1] <- PM.popmat %*% PM.n.mat[,i]
      
      # PALAEOLOXODON
      PC.s.alpha <- estBetaParams(PC.Sx, PC.s.sd.vec^2)$alpha
      PC.s.beta <- estBetaParams(PC.Sx, PC.s.sd.vec^2)$beta
      PC.s.stoch <- rbeta(length(PC.s.alpha), PC.s.alpha, PC.s.beta)
      
      if (rbinom(1, 1, 0.14/PC.gen.l) == 1) { # catastrophe
        cat.alpha <- estBetaParams(0.5, 0.05^2)$alpha
        cat.beta <- estBetaParams(0.5, 0.05^2)$beta
        PC.s.stoch <- PC.s.stoch * (rbeta(1, cat.alpha, cat.beta)) }
      
      # stochastic fertility sampler (gaussian)
      PC.fert.stch <- rnorm(length(PC.popmat[,1]), PC.pred.p.mm, PC.m.sd.vec)
      
      PC.totN.i <- sum(PC.n.mat[,i], na.rm=T)
      PC.pred.red <- PC.a.lp/(1+(PC.totN.i/PC.b.lp)^PC.c.lp)
      
      diag(PC.popmat[2:(PC.age.max+1),]) <- (PC.s.stoch[-(PC.age.max+1)])*PC.pred.red
      PC.popmat[PC.age.max+1,PC.age.max+1] <- (PC.s.stoch[PC.age.max+1])*PC.pred.red
      PC.popmat[1,] <- ifelse(PC.fert.stch < 0, 0, PC.fert.stch)
      PC.n.mat[,i+1] <- PC.popmat %*% PC.n.mat[,i]
      
      # HUMANS
      ppl.real <- as.vector(round(people.vec[k]*human.ssd, 0))$x
      f.real <- ifelse(ppl.real < 2, rbinom(length(ppl.real),size=1,prob=0.5), round(ppl.real/2, 0))
      m.real <- ifelse(ppl.real-f.real < 0, 0, ppl.real-f.real)
      
      # meat requirements
      meatpyr.f <- f.real*relmeatpyr.age.f
      meatpyr.m <- m.real*relmeatpyr.age.m
      totmeatpyr <- sum(meatpyr.f + meatpyr.m)
      
      # equivalent meat in PC and PM (if only source of meat)
      PM.prob3.alpha <- estBetaParams(PM.prob3, (0.05*PM.prob3)^2)$alpha
      PM.prob3.beta <- estBetaParams(PM.prob3, (0.05*PM.prob3)^2)$beta
      PM.prob3.stoch <- rbeta(length(PM.prob3.alpha), PM.prob3.alpha, PM.prob3.beta)
      
      PM.meat.tot <- PM.prob3.stoch*totmeatpyr
      PC.meat.tot <- totmeatpyr - PM.meat.tot
      
      PM.meat.f <- PM.prob2.stoch*sum(meatpyr.f)
      PC.meat.f <- sum(meatpyr.f) - PM.meat.f
      PM.meat.m <- PM.prob2.stoch*sum(meatpyr.m)
      PC.meat.m <- sum(meatpyr.m) - PM.meat.m
      
      # meat equivalents in n PM and PC
      # stochastic expression of edible meat return
      psamp.alpha <- estBetaParams(mpsamp, sdsamp^2)$alpha
      psamp.beta <- estBetaParams(mpsamp, sdsamp^2)$beta
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
      
      # assume half females & half males harvested
      # PM - female
      PMwtvecf <- (1000*PM.ssd*PM.grf)
      PMwtvecfwt <- PMwtvecf/sum(PMwtvecf)
      PMf.wtXage <- 0.5*(PM.totf.wt)*PMwtvecfwt
      PMf.NxAge <- as.vector(round(PMf.wtXage/PM.grf, 0))
      
      # PM - male
      PMwtvecm <- (1000*PM.ssd*PM.grm)
      PMwtvecmwt <- PMwtvecm/sum(PMwtvecm)
      PMm.wtXage <- 0.5*(PM.totm.wt)*PMwtvecmwt
      PMm.NxAge <- as.vector(round(PMm.wtXage/PM.grm, 0))
      
      # PC - female
      PCwtvecf <- (1000*PC.ssd*PC.grf)
      PCwtvecfwt <- PCwtvecf/sum(PCwtvecf)
      PCf.wtXage <- 0.5*(PC.totf.wt)*PCwtvecfwt
      PCf.NxAge <- as.vector(round(PCf.wtXage/PC.grf, 0))
      
      # PC - male
      PCwtvecm <- (1000*PC.ssd*PC.grm)
      PCwtvecmwt <- PCwtvecm/sum(PCwtvecm)
      PCm.wtXage <- 0.5*(PC.totm.wt)*PCwtvecmwt
      PCm.NxAge <- as.vector(round(PCm.wtXage/PC.grm, 0))
      
     # if not enough PM meat to remove relative to size of extant PM population,
      # transfer meat equivalent to PC population
      PCf.NxAgeADD <- rep(0,length(PCwtvecfwt))
      if ((sum(PM.n.mat[,i+1]) - sum(PMf.NxAge)) < 0) {
        PMfn2twt <- sum(PMf.NxAge*PM.grf) # transform extra female PMs to total female weight
        PMfn2twtCor <- PMfn2twt*(1-PC.meat.oth.t) # correct for meat 'other' proportion
        PCf.wtXageADD <- (PMfn2twtCor)*PCwtvecfwt # weight by PC age
        PCf.NxAgeADD <- as.vector(round(PCf.wtXageADD/PC.grf, 0)) # n PC by age
       }
      
      # REMOVALS FROM n MAT - PM
      PM.n.mat[,i+1] <- ifelse((PM.n.mat[,i+1] - PMf.NxAge) < 0, 0, (PM.n.mat[,i+1] - PMf.NxAge))
      
      # REMOVALS FROM n MAT - PC
      PC.n.mat[,i+1] <- ifelse((PC.n.mat[,i+1] - (PCf.NxAge+PCf.NxAgeADD)) < 0, 0, (PC.n.mat[,i+1] - (PCf.NxAge+PCf.NxAgeADD)))
      
      # animals removed
      PM.n.rem[,i+1] <- ifelse(PM.n.mat[,i+1] == 0, 0, PMf.NxAge)
      PC.n.rem[,i+1] <- ifelse(PC.n.mat[,i+1] == 0, 0, PCf.NxAge)
      
    } # end i loop
    
    PM.n.sums.mat[e,] <- (as.vector(colSums(PM.n.mat)))
    PM.n.remsums.mat[e,] <- (as.vector(colSums(PM.n.rem)))
    
    PC.n.sums.mat[e,] <- (as.vector(colSums(PC.n.mat)))
    PC.n.remsums.mat[e,] <- (as.vector(colSums(PC.n.rem)))
    
    #if (e %% itdiv==0) print(e) 
    
  } # end e loop
  
  # total N
  PM.n.md <- apply(PM.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  PM.n.up <- apply(PM.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  PM.n.lo <- apply(PM.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  PC.n.md <- apply(PC.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  PC.n.up <- apply(PC.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  PC.n.lo <- apply(PC.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # total removed
  PM.rem.md <- apply(PM.n.remsums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  PM.rem.up <- apply(PM.n.remsums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  PM.rem.lo <- apply(PM.n.remsums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  PC.rem.md <- apply(PC.n.remsums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
  PC.rem.up <- apply(PC.n.remsums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  PC.rem.lo <- apply(PC.n.remsums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  # years when > 1 animals removed
  PM.yrsRem.md[k] <- ifelse(length(which(PM.rem.md >= 1)) + 1 >= t, NA, length(which(PM.rem.md >= 1)) + 1)
  PM.yrsRem.up[k] <- ifelse(length(which(PM.rem.up >= 1)) + 1 >= t, NA, length(which(PM.rem.up >= 1)) + 1)
  PM.yrsRem.lo[k] <- ifelse(length(which(PM.rem.lo >= 1)) + 1 >= t, NA, length(which(PM.rem.lo >= 1)) + 1)
  
  PC.yrsRem.md[k] <- ifelse(length(which(PC.rem.md >= 1)) + 1 >= t, NA, length(which(PC.rem.md >= 1)) + 1)
  PC.yrsRem.up[k] <- ifelse(length(which(PC.rem.up >= 1)) + 1 >= t, NA, length(which(PC.rem.up >= 1)) + 1)
  PC.yrsRem.lo[k] <- ifelse(length(which(PC.rem.lo >= 1)) + 1 >= t, NA, length(which(PC.rem.lo >= 1)) + 1)
  
  # quasi-extinction probability
  PM.Q.ext.mat <- ifelse(PM.n.sums.mat < Q.ext.thresh, 1, 0)
  PM.Q.ext.sum <- apply(PM.Q.ext.mat[,ceiling(PM.gen.l):dim(PM.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  PM.Q.ext.pr[k] <- length(which(PM.Q.ext.sum > 0)) / iter
  
  PC.Q.ext.mat <- ifelse(PC.n.sums.mat < Q.ext.thresh, 1, 0)
  PC.Q.ext.sum <- apply(PC.Q.ext.mat[,ceiling(PC.gen.l):dim(PC.Q.ext.mat)[2]], MARGIN=1, sum, na.rm=T)
  PC.Q.ext.pr[k] <- length(which(PC.Q.ext.sum > 0)) / iter

  print(people.vec[k]) 
  if (k > 1) {
    plot(people.vec, PM.Q.ext.pr, type="l", xlab="people on the island", ylab="Pr(quasi-extinction)", ylim=c(0,1))
    lines(people.vec, PC.Q.ext.pr, lty=2, col="red")
    
  } # end if
  
} # end k loop

plot(people.vec, PM.Q.ext.pr, type="l", xlab="people on the island", ylab="Pr(quasi-extinction)", ylim=c(0,1))
lines(people.vec, PC.Q.ext.pr, lty=2)

plot(people.vec, PM.yrsRem.md, type="l", xlab="people on the island", ylab="years to extinction [PM]", ylim=c(min(PM.yrsRem.lo, na.rm=T),max(PM.yrsRem.up, na.rm=T)))
lines(people.vec, PM.yrsRem.up, lty=2, col="red")
lines(people.vec, PM.yrsRem.lo, lty=2, col="red")

plot(people.vec, PC.yrsRem.md, type="l", xlab="people on the island", ylab="years to extinction [PC]", ylim=c(min(PC.yrsRem.lo, na.rm=T),max(PC.yrsRem.up, na.rm=T)))
lines(people.vec, PC.yrsRem.up, lty=2, col="red")
lines(people.vec, PC.yrsRem.lo, lty=2, col="red")

out.dat <- data.frame(people.vec, PM.Q.ext.pr, PC.Q.ext.pr, PM.yrsRem.md, PM.yrsRem.up, PM.yrsRem.lo,
                      PC.yrsRem.md, PC.yrsRem.up, PC.yrsRem.lo)

