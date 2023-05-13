rm(list = ls())
setwd("C:/Users/joaquin/Desktop/new_models")

library(pacman)
pacman::p_load(RandomFields, geoR, fields, prodlim, TMB, INLA, dplyr, tmbstan, rstan, parallel)

options(scipen=999)

# Calculate the number of cores
no_cores <- detectCores() - 5


## From package 'geoR':
# Matern correlation
matern_fun <- function(h, nu, kappa) {
  ifelse(h > 0, besselK(h * kappa, nu) * (h * kappa)^nu / 
           (gamma(nu) * 2^(nu - 1)), 1)
}

## Test:
set.seed(123)
n <- 50
loc <- matrix(runif(2*n, 0, 10), n)
dim(loc)

 
dist_mat <- as.matrix(dist(loc))

sigma_e <- 1.0
sigma_u <- 1.0
kappa <- 5
rho <- sqrt(8)/kappa
nu <- 1

mcor <- matern_fun(dist_mat, nu, kappa)
#mcov <- sigma_e * diag(nrow(mcor)) + sigma_u * mcor

u1 <- rnorm(n,0,1) %*% chol(mcor) #Spatial random effect
u2 <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = mcor)


plot(loc, cex = u1 - min(u1))
plot(loc, cex = u2 - min(u2))

#sim_grf <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = K)

beta0 <- 2.0
beta1 <- 1.0 
x1 <- rnorm(n, 0 ,1)


y1 <- beta0 + beta1*x1 + u1 + rnorm(n, mean = 0, sd = sigma_e)
hist(y1)
# 
y2 <- beta0 + beta1*x1 + u2 + rnorm(n, mean = 0, sd = sigma_e)
hist(y2)
# 
# 
y3 <- rnorm(n, mean = beta0 + beta1*x1 + u1, sd = sigma_e)
hist(y3)

y4 <- rnorm(n, mean = beta0 + beta1*x1 + u2, sd = sigma_e)
hist(y4)








#========================================
# Fit full Gaussian random field in TMB 
#========================================
tmb_data <- list(y = y1, 
                 X = model.matrix(~ 1 + x1),
                 dist_mat=as.matrix(dist(loc)))

tmb_par <- list(beta = c(0.1, 0.1),
                logsigma = 0.1,
                logphi = 0.1,
                logkappa = 0.1,
                u = rep(0, n))



compile("tmbstan_matern.cpp")
dyn.load(dynlib("tmbstan_matern"))


## Create the TMB object
obj <- MakeADFun(data = tmb_data, 
                 parameters = tmb_par, 
                 random = 'u',
                 DLL = 'tmbstan_matern', 
                 hessian = T)

obj$fn() #logLik evaluates
# 
#opt <- nlminb(obj$par, obj$fn, obj$gr, lower=rep(-5,3),upper=rep(3,3))
opt <- nlminb(obj$par, obj$fn, obj$gr)

opt$par
exp(opt$par[3:5])



startTime <- Sys.time()
fit_tmbstan = tmbstan(obj, chains = 2, open_progress = FALSE, 
                      init=list("last.par.best"), 
                      control = list(max_treedepth = 5, adapt_delta = 0.6), 
                      iter = 2000, warmup=1000, cores = no_cores, seed=483892929)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)


s <- summary(fit_tmbstan, probs = c(0.25, 0.75))
s$summary  # all chaines merged


c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")


library(gridExtra)
library(tidyverse)
traceplot(fit_tmbstan, pars=names(obj$par), inc_warmup=TRUE) + 
  theme(strip.text = ggplot2::element_text(size = 16, color = "black"),
        axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14))

#load package
require(MCMCvis)
MCMCsummary(fit_tmbstan, round = 2)
pairs(fit_tmbstan, pars = c("logsigma", "logphi", "logkappa")) # below the diagonal
