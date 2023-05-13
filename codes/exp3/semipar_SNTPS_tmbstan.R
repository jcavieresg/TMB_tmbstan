rm(list = ls())
setwd("C:/Users/joaquin/Desktop/new_models")
library(TMB)
library(mgcv)
library(sn)
library(MASS)
library(tidyverse)
library(gridExtra)
library(grid)
library(lattice)
library(tmbstan)
library(rstan)
library(tmbstan)
library(rstan)
library(bayesplot)
library(geoR)
library(raster)
library(Matrix)


options(scipen=999)

#=================================================
# Compilamos el modelo y lo cargamos en R
TMB::compile('semipar_SNTPS_tmbstan.cpp')
dyn.load(dynlib("semipar_SNTPS_tmbstan"))
# #=================================================
# 


#=============================================================================
#                              Applied analysis
#=============================================================================
data = read.csv("north2.csv", header = T)
head(data, 3)

# data <- data[sample(nrow(data), 500), ]
# dim(data)


#====================================================
#      Obtenemos matriz penalizada desde mgcv
#====================================================
# Creamos una variable artificial
z <- rep(0, length(data$cpue))
data$z <- z
head(data)


# data <- sample_n(data, 2500)
# dim(data)

#===================
#     mgcv setup
#===================
m0 <- gam(z ~ year + trim + destine + depth  + s(latitude, longitude, bs = "tp", k = 13), data = data, fit = FALSE)
#m0 <- gam(z ~ year + depth  + s(latitude, longitude, bs = "tp", k = 13), data = data, fit = FALSE) 


Xtps <- m0$X[, c(-1, -2, -3, -4, -5)]                    # Matricial form without intercept and the parameters asociated with the covariates
#Xtps <- m0$X[, c(-1, -2. -3)]                    # Matricial form without intercept and the parameters asociated with the covariates
Stps <- m0$smooth[[1]]$S[[1]]                # Extrtact penelization matrices
dim(Stps)

S_list = list(Stps)
S_combined = .bdiag(S_list)                               # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow))                       # Find dimension of each S (in this case only one)
#----------------------------------------

#===========================================================================
# For report, used for constructing plots
latitude = seq(min(data$latitude), max(data$latitude), by = 0.1)
longitude = seq(min(data$longitude), max(data$longitude), by = 0.1)

# by cambia para hacer coincidir los largos de los vectores
tpsReport = PredictMat(m0$smooth[[1]], data = data.frame(latitude, longitude))
tpsReport = list(tpsReport)



#-------------------------------------------



# # Plot
# to_plot = TRUE   # generate plots
# 
# class(data)
# 
# if (to_plot){
#   ggplot(data) +
#     geom_point(aes(x = longitude, y = latitude, fill = cpue)) +
#     coord_equal() +
#     scale_fill_viridis_c()
# }





#=================================================================
#                           TMB modelling
#=================================================================
data$year <- as.factor(data$year)
data$trim <- as.factor(data$trim)
data$destine <- as.factor(data$destine)


#======================================
#                TMB data
#======================================
tmb_data = list(likelihood = 3,
                #cpue     = as.vector(data$cpue),          # Response
                cpue     = sqrt(as.vector(data$cpue)),          # Response
                year = as.numeric(data$year) - 1,
                trim = as.numeric(data$trim) - 1,
                destine = as.numeric(data$destine) - 1,
                depth = as.numeric(data$depth),
                #X    = model.matrix(~1 + data$depth + as.numeric(data$year) - 1),
                TPS     = Xtps,                                   # Design matrix, without intercept and betas
                S     = as(S_combined, "dgTMatrix"),                          # Combined penalty matrix
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
               logalpha = exp(0.1),
               #logalpha = rep(rep(0.1,length(Sdims))),       # Log spline penalization coefficients
               logsigma = exp(0.1),
               loglambda = exp(0.1))
               




#=====================================
#             SKEW NORMAL MODEL
#=====================================
obj <- MakeADFun(tmb_data, random = c("smoothCoefs"), tmb_par, DLL="semipar_SNTPS_tmbstan")


# Lognormal
startTime <- Sys.time()
opt <- nlminb(obj$par,obj$fn,obj$gr)
rep <- sdreport(obj)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)

# Gamma
startTime <- Sys.time()
opt2 <- nlminb(obj$par,obj$fn,obj$gr)
rep2 <- sdreport(obj)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)


# Skew-normal
startTime <- Sys.time()
opt3 <- nlminb(obj$par,obj$fn,obj$gr)
rep3 <- sdreport(obj)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)

obj$report()$log_lik
#-------------------------------------------

library(TMBhelper)
TMBAIC(opt)
TMBAIC(opt2)
TMBAIC(opt3)


#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================
# Calculate the number of cores
no_cores <- parallelly::availableCores()  - 1

init.fn <- function()
  split(unname(obj$env$last.par.best),names(obj$env$last.par.best))


startTime <- Sys.time()
fit = tmbstan(obj, chains = 3, open_progress = FALSE, 
              init = init.fn, control = list(max_treedepth = 13, adapt_delta = 0.9), 
              iter = 3000, warmup=1000, cores = no_cores, seed=483892929)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)


## 1.9 hours minutes!!!


s <- summary(fit, probs = c(0.25, 0.75))
s$summary  # all chaines merged



setwd("C:/Users/joaquin/Desktop/new_models/output_real_application")

posterior_tps_tmbstan <- as.matrix(fit)

saveRDS(posterior_tps_tmbstan, file='posterior_tps_tmbstan.RDS')
saveRDS(fit, file='fit_tps_tmbstan.RDS')

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
pairs(fit, pars = c("logsigma", "logalpha", "lp__"), las = 1) # below the diagonal
pairs(fit, pars = c("logsigma", "logalpha", "loglambda")) # below the diagonal


pairs(fit, pars = c("logsigma", "logalpha", "lp__"), las = 1) # below the diagonal


params_cp <- as.data.frame(fit)
names(params_cp) <- gsub("chain:1.", "", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("[", ".", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("]", "", names(params_cp), fixed = TRUE)
params_cp$iter <- 1:6000


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





# Simulating from the SIMULATE function
mat_sim = matrix(data=NA, nrow=length(obj$simulate()$cpue_sim), ncol=100)
mat_sim


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = obj$simulate()$cpue_sim
  }
}
mat_sim

# Plot TPS model using Gamma likelihood
hist(data$cpue, col = "gray90", prob = TRUE, main = "100 simulated samples", cex.main = 2, xlab = "", col.lab = 'blue', col.main="blue", cex.lab = 1.4, cex.axis = 1.2, ylim = c(0, 0.008))
for (j in 1: ncol(mat_sim)){
  lines(x = density(x = mat_sim[, j]),  lty="dotted", col="azure4", lwd=1)
}
#legend("topright", "A", bty = "n", text.col="blue", cex = 1.4)
legend(6, 0.3, c("response", "simulations"), lwd=4, col=c("gray90", "azure4"), cex = 1.5)
