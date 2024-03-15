rm(list = ls())
setwd("C:/Users/cavieresgaet/Desktop/stan_tmb/new_models/output_real_application/gamma_tps")

library(pacman)
pacman::p_load(geoR, fields, prodlim, TMB, INLA, dplyr, tmbstan, rstan, parallel, mgcv, ggplot2, bayesplot, TMBhelper)

options(scipen=999)

# Calculate the number of cores
no_cores <- parallelly::availableCores()  - 1



options(scipen=999)

#=================================================
# Compilamos el modelo y lo cargamos en R
TMB::compile('gamma_tps.cpp')
dyn.load(dynlib("gamma_tps"))
# #=================================================
# 


#=============================================================================
#                              Applied analysis
#=============================================================================
data = read.csv("north2.csv", header = T)
head(data, 3)


#====================================================
#      Obtenemos matriz penalizada desde mgcv
#====================================================
# Creamos una variable artificial
z <- rep(0, length(data$cpue))
data$z <- z
head(data)

#===================
#     mgcv setup
#===================
m0 <- gam(z ~ year + trim + destine + depth  + s(latitude, longitude, bs = "tp", k = length(table(data$site))), 
                                                 data = data, fit = FALSE)


Xtps <- m0$X[, c(-1, -2, -3, -4, -5)]        # Matricial form without intercept and the parameters asociated with the covariates
Stps <- m0$smooth[[1]]$S[[1]]                # Extrtact penelization matrices
dim(Stps)

S_list = list(Stps)
S_combined = .bdiag(S_list)                               # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow))                       # Find dimension of each S (in this case only one)
#----------------------------------------

#===========================================================================
# For report, used for constructing plots
latitude = seq(min(data$latitude), max(data$latitude), by = 1)
longitude = seq(min(data$longitude), max(data$longitude), by = 1)

# by cambia para hacer coincidir los largos de los vectores
tpsReport = PredictMat(m0$smooth[[1]], data = data.frame(latitude, longitude))
tpsReport = list(tpsReport)



#=================================================================
#                           TMB modelling
#=================================================================
data$year <- as.factor(data$year)
data$trim <- as.factor(data$trim)
data$destine <- as.factor(data$destine)


#======================================
#                TMB data
#======================================
tmb_data = list(likelihood = 2, # 1 = lognormal, 2 = Gamma, 3 = Slew normal
                cpue     = as.vector(data$cpue),          # Gamma response
                year = as.numeric(data$year) - 1,
                trim = as.numeric(data$trim) - 1,
                destine = as.numeric(data$destine) - 1,
                depth = as.numeric(data$depth),
                TPS     = Xtps,                               # Design matrix, without intercept and betas
                S     = as(S_combined, "dgTMatrix"),          # Combined penalty matrix
                Sdims = Sdims,
                tpsReport = .bdiag(tpsReport))



#=====================================
#            TMB parameters
#=====================================
tmb_par = list(beta0 = 0.1,
               beta_year = rep(0.1, length(levels(data$year))),
               beta_trim = rep(0.1, length(levels(data$trim))),
               beta_destine = rep(0.1, length(levels(data$destine))),
               beta_depth = 0.1,
               smoothCoefs = rep(rnorm(ncol(Xtps), 0, 1)),              # Spline coefficients
               loglambda = exp(0.1),
               logsigma = exp(0.1))
               

#=====================================
#             GAMMA MODEL
#=====================================
obj <- MakeADFun(tmb_data, random = c("smoothCoefs"), tmb_par, DLL="gamma_tps")


# Gamma
startTime <- Sys.time()
opt <- nlminb(obj$par,obj$fn,obj$gr)
rep <- sdreport(obj)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)

obj$report()$log_lik

#-------------------------------------------

TMBAIC(opt)
TMBAIC(opt2)
TMBAIC(opt3)


saveRDS(obj, file='TMB_tps_gamma.RDS')

#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================
# Calculate the number of cores

init.fn <- function()
  split(unname(obj$env$last.par.best),names(obj$env$last.par.best))


startTime <- Sys.time()
fit_tps_gamma = tmbstan(obj, chains = 3, open_progress = FALSE, 
              init = init.fn, control = list(max_treedepth = 13, adapt_delta = 0.95), 
              iter = 3500, warmup= 700, cores = no_cores, seed=483892929)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)

##==========================================================
##        Time difference of 6.750947 hours !!!!!!
##==========================================================
s <- summary(fit_tps_gamma, probs = c(0.25, 0.75))
s$summary  # all chaines merged

posterior_tps_gamma <- as.matrix(fit_tps_gamma)
saveRDS(posterior_tps_gamma, file='posterior_tps_gamma.RDS')
saveRDS(fit_tps_gamma, file='fit_tps_gamma.RDS')

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

#load packages
require(MCMCvis)
MCMCsummary(fit_tps, round = 2)

#pairs(fit, pars=names(fit)[-grep('beta0',names(fit))][1:10])
pairs(fit_tps_gamma, pars = c("logsigma", "loglambda", "lp__"), las = 1) # below the diagonal
pairs(fit_tps_gamma, pars = c("logsigma", "logomega", "loglambda")) # below the diagonal


params_cp <- as.data.frame(fit_tps)
names(params_cp) <- gsub("chain:1.", "", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("[", ".", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("]", "", names(params_cp), fixed = TRUE)
params_cp$iter <- 1:6900


# logalpha
par(mfrow=c(3,2),mar=c(4, 5,0.5,0.5), oma=c(2,3,1,1))
plot(params_cp$iter, params_cp$loglambda, col=c_dark, pch=16, cex=2.2, type = "l",
     xlab="Iteration", ylab="logalpha", cex.lab=2.0, cex.axis=1.5)

running_means_lambda = sapply(params_cp$iter, function(n) mean(params_cp$loglambda[1:n]))
plot(params_cp$iter, running_means_lambda, col=c_dark, pch=16, cex=1.2,  cex.lab=2.0, cex.axis=1.5,
     xlab="Iteration", ylab="MCMC mean of loglambda")
abline(h=mean(running_means_lambda), col="grey", lty="dashed", lwd=3)



# logsigma
plot(params_cp$iter, params_cp$logsigma, col=c_dark, pch=16, cex=0.8, type = "l",
     xlab="Iteration", ylab="logsigma",  cex.lab=1.3, cex.axis=1.3)

running_means_sigma = sapply(params_cp$iter, function(n) mean(params_cp$logsigma[1:n]))
plot(params_cp$iter, running_means_sigma, col=c_dark, pch=16, cex=0.8,  cex.lab=1.3, cex.axis=1.3,
     xlab="Iteration", ylab="MCMC mean of logsigma")
abline(h=mean(running_means_sigma), col="grey", lty="dashed", lwd=3)


# loglambda
plot(params_cp$iter, params_cp$logomega, col=c_dark, pch=16, cex=0.8, type = "l",
     xlab="Iteration", ylab="logomega",  cex.lab=1.3, cex.axis=1.3)

running_means_omega = sapply(params_cp$iter, function(n) mean(params_cp$logomega[1:n]))
plot(params_cp$iter, running_means_omega, col=c_dark, pch=16, cex=0.8,  cex.lab=1.3, cex.axis=1.3,
     xlab="Iteration", ylab="MCMC mean of logomega")
abline(h=mean(running_means_omega), col="grey", lty="dashed", lwd=3)
mtext("Convergence of the parameters loglambda, logsigma and logomega", outer=TRUE,  cex=1, line=-0.5)




# Divergences
divergent = get_sampler_params(fit_tps, inc_warmup=FALSE)[[1]][,'divergent__']
sum(divergent)

## Methods provided by 'rstan'
class(fit_tps)
methods(class ="stanfit")

get_cppo_mode(fit_tps)
get_stancode(fit_tps)
get_stanmodel(fit_tps)
log_prob(fit_tps)
loo(fit)

## ESS and Rhat from rstan::monitor
mon = monitor(fit_tps)
max(mon$Rhat)
min(mon$Tail_ESS)

# evalaute problem of convergence
sum(mon$Rhat > 1.01)
sum(mon$Tail_ESS < 400)

