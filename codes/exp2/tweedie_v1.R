rm(list = ls())
setwd("C:/Users/joaquin/Desktop/stan_tmb")
library(TMB)
library(mgcv)
library(tmbstan)
library(MASS)
library(patchwork)
library(tidyverse)
library(Matrix)
library(gridExtra)
library(grid)
library(lattice)
library(sp)
library(rgdal)
library(raster)
library(leaflet)
library(mapview)
library(tictoc)
library(parallel)
library(tweedie)

options(scipen=999)

# Calculate the number of cores
no_cores <- detectCores() - 1

#=================================================
# Compilamos el modelo y lo cargamos en R
compile("tweedie_v1.cpp")
dyn.load(dynlib("tweedie_v1"))
#=================================================

#=============

#===========================================================================
#                          Simulated data
#===========================================================================
# Simule new data
n <- 1000
mu_ini  <- 2
phi_ini<- 2                          # Desvest
p_ini  <- 1.5                          # slant parameter of the Skew Normal


#y_sim = rtweedie(n, mu = mu_ini, phi = phi_ini, power = p_ini)
y_sim = tweedie:::rtweedie(n, mu= mu_ini, phi= phi_ini, p= p_ini) 
hist(y_sim)



#=================================================================
#                           TMB modelling
#=================================================================

#======================================
#                TMB data
#======================================
tmb_data = list(y = y_sim)           



#=====================================
#            TMB parameters
#=====================================
tmb_par = list(mu = 1.1, 
               phi = 1.1,
               p = 1.1)



#=====================================
#             Run the model
#=====================================
obj = MakeADFun(data = tmb_data, parameters = tmb_par, DLL = "tweedie_v1", hessian = TRUE)

# lwr <- c(-Inf, -Inf, 0, -Inf)
# upr <- c(Inf, Inf, 0.69, Inf)

tic("TMB Tweddie model fitting")
# opt = nlminb(obj$par,obj$fn,obj$gr, lower=lwr, upper=upr)
opt = nlminb(obj$par,obj$fn,obj$gr)
toc()

rep = sdreport(obj)

param = cbind(opt$par[1], exp(opt$par[2]), exp(opt$par[3]))
colnames(param) = c("mu", "phi", "p")
param












#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================
init.fn <- function()
  split(unname(obj$env$last.par.best),names(obj$env$last.par.best))

tic("Time of estimation")
fit = tmbstan(obj, chains = 2, open_progress = FALSE,# lower = lwr, upper = upr,
              init = init.fn, control = list(max_treedepth = 10, adapt_delta = 0.8), 
              iter = 3000, warmup=1000, cores=no_cores, seed=483892929)
toc()

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")


traceplot(fit, pars=names(obj$par), inc_warmup=TRUE)

#load package
require(MCMCvis)
MCMCsummary(fit, round = 2)

#pairs(fit, pars=names(fit)[-grep('beta0',names(fit))][1:10])
#pairs(fit, pars = c("logphi", "logalpha", "lp__"), las = 1) # below the diagonal
pairs(fit, pars = c("phi", "p")) # below the diagonal


pairs(fit, pars = c("logsigma", "logalpha", "lp__"), las = 1) # below the diagonal


params_cp <- as.data.frame(fit)
names(params_cp) <- gsub("chain:1.", "", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("[", ".", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("]", "", names(params_cp), fixed = TRUE)
params_cp$iter <- 1:4000


# logalpha
par(mfrow=c(3,2),mar=c(4,4,0.5,0.5), oma=c(2,3,1,1))
plot(params_cp$iter, params_cp$logalpha, col=c_dark, pch=16, cex=0.8, type = "l",
     xlab="Iteration", ylab="logalpha", cex.lab=1.3, cex.axis=1.3)

running_means_alpha = sapply(params_cp$iter, function(n) mean(params_cp$logalpha[1:n]))
plot(params_cp$iter, running_means_alpha, col=c_dark, pch=16, cex=0.8,  cex.lab=1.3, cex.axis=1.3,
     xlab="Iteration", ylab="MCMC mean of logalpha")
abline(h=mean(running_means_alpha), col="grey", lty="dashed", lwd=3)



# logsigma
plot(params_cp$iter, params_cp$logsigma, col=c_dark, pch=16, cex=0.8, type = "l",
     xlab="Iteration", ylab="logsigma",  cex.lab=1.3, cex.axis=1.3)

running_means_sigma = sapply(params_cp$iter, function(n) mean(params_cp$logsigma[1:n]))
plot(params_cp$iter, running_means_sigma, col=c_dark, pch=16, cex=0.8,  cex.lab=1.3, cex.axis=1.3,
     xlab="Iteration", ylab="MCMC mean of logsigma")
abline(h=mean(running_means_sigma), col="grey", lty="dashed", lwd=3)


# loglambda
plot(params_cp$iter, params_cp$loglambda, col=c_dark, pch=16, cex=0.8, type = "l",
     xlab="Iteration", ylab="loglambda",  cex.lab=1.3, cex.axis=1.3)

running_means_lambda = sapply(params_cp$iter, function(n) mean(params_cp$loglambda[1:n]))
plot(params_cp$iter, running_means_lambda, col=c_dark, pch=16, cex=0.8,  cex.lab=1.3, cex.axis=1.3,
     xlab="Iteration", ylab="MCMC mean of loglambda")
abline(h=mean(running_means_lambda), col="grey", lty="dashed", lwd=3)
mtext("Convergence of the parameters alpha, sigma and lambda", outer=TRUE,  cex=1, line=-0.5)




# Divergences
divergent = get_sampler_params(fit, inc_warmup=FALSE)[[1]][,'divergent__']
sum(divergent)

## Methods provided by 'rstan'
class(fit)
methods(class ="stanfit")

get_cppo_mode(fit)
get_stancode(fit)
get_stanmodel(fit)
log_prob(fit)
loo(fit)
traceplot(fit)

#launch_shinystan(fit)

## ESS and Rhat from rstan::monitor
mon = monitor(fit)
max(mon$Rhat)
min(mon$Tail_ESS)

# evalaute problem of convergence
sum(mon$Rhat > 1.01)
sum(mon$Tail_ESS < 400)

source('monitornew.R')
source('monitorplot.R')
source('stan_utility.R')

which_min_ess = which.min(mon[1:200, 'Tail_ESS'])
plot_local_ess(fit = fit, par = which_min_ess, nalpha = 10)

plot_quantile_ess(fit = fit, par = which_min_ess, nalpha = 50)

plot_change_ess(fit = fit, par = which_min_ess)

check_rhat(fit)
check_treedepth(fit, 10)
check_energy(fit)  #
check_div(fit)


# Variance
library(MCMCvis)
color_scheme_set("viridisE")
#mcmc_rank_hist(fit, pars = c("sigma_beta_year", "sigma_beta_depth", "sigma_beta_trim", "sigma_beta_destine"), ref_line = TRUE)
mcmc_rank_hist(fit, pars = c("logalpha", "logsigma", "loglambda"), ref_line = TRUE)   # spatial random field


## Extract marginal posteriors
posterior <- as.matrix(fit)

exp(mean(posterior[, "logalpha"]))
exp(mean(posterior[, "logsigma"]))
exp(mean(posterior[, "loglambda"]))

exp(opt$par[3:5])




dyn.unload(dynlib("tweedie_v1"))
