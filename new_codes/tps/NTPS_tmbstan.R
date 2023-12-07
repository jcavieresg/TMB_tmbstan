rm(list = ls())
setwd("C:/Users/Usuario/Desktop/tps_vs_spde/gaussian_model/tps")

library(pacman)
pacman::p_load(geoR, fields, prodlim, TMB, TMBhelper, mgcv, dplyr, tmbstan, parallel, MASS, Matrix,
               raster, ggplot2, gridExtra, bayesplot)



options(scipen=999)
# Calculate the number of cores
no_cores <- detectCores() - 1

#=================================================
# Compilamos el modelo y lo cargamos en R
TMB::compile('NTPS_tmbstan.cpp')
dyn.load(dynlib("NTPS_tmbstan"))
# #=================================================



#==================================================
#                 SIMULATED DATA
#==================================================

# MAIN FUNCTION to run models for distincs "n" of the grid
set.seed(1234)
run_tmb <- function(nlength){
coords = matrix(runif(2*nlength), nlength)*10 # coordinates

# Creating a response variable only for simulation
z = 0

# Simulated continuous covariate
x0 = runif(length(coords[,1]), 0, 1)


# Build this data frame ONLY for mgcv setup
df <- data.frame(loc.data, z, x0)
names(df) <- c("s1", "s2", "z", "x0")
xt <- unique(df[c("s1", "s2")]) 
nrow(xt)


#====================================================
#      Getting matrices from mgcv
#====================================================
# mgcv setup
m0 <- gam(z ~ x0 + s(s1, s2, bs = "tp", k = 50), data = df, fit = FALSE) 


Xtps <- m0$X[, c(-1, -2)]            # Matricial form without intercept and the parameters asociated with the covariates
Stps <- m0$smooth[[1]]$S[[1]]        # Extrtact penelization matrices
dim(Stps)

S_list = list(Stps)
S_combined = .bdiag(S_list)          # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow))  # Find dimension of each S (in this case only one)
#----------------------------------------



#===========================================================================
# For report, used for constructing plots
s1 = seq(min(df$s1), max(df$s1), by = 0.5)
s2 = seq(min(df$s2), max(df$s2), by = 0.5)

tpsReport = PredictMat(m0$smooth[[1]], data = data.frame(s1, s2))
tpsReport = list(tpsReport)
#-------------------------------------------


#===========================================================================
#                          Simulate data
#===========================================================================
n <- length(df$z)
beta0_ini  <- 1
beta1_ini  <- 2
sigma_ini   <- 0.3                    
alpha_ini  <- 1
E = ginv(Stps)


smoothCoefs <- rmvn(n = 1, mu = rep(0, nrow(Stps)),  V = alpha_ini*E) # From Wood 2011 (apendix)

# Simulated covariate
x1 = runif(n, 0, 1)                     

# mu simulated
mu_sim <- beta0_ini + beta1_ini*x1 + Xtps %*% smoothCoefs # mean value in each location

# y (response variable) simulated
y_sim <- rnorm(n, mean = mu_sim, sd = sigma_ini)


# Add the new variables to the data.frame
df$y_sim <- y_sim
df$x1 <- x1
head(df)






#=================================================================
#                           TMB modelling
#=================================================================

#======================================
#                TMB data
#======================================
tmb_data = list(#space = 2,                              # 0 = No space, 1 = Space
                likelihood = 1,                          # Only Gaussian for now
                y     = as.vector(df$y_sim),             # Response
                x1    = as.vector(df$x1),                # Covariate
                TPS     = Xtps,                          # Design matrix, without intercept and betas
                S     = as(S_combined, "sparseMatrix"),  # Combined penalty matrix
                Sdims = Sdims,
                tpsReport = .bdiag(tpsReport))





#=====================================
#            TMB parameters
#=====================================
tmb_par = list(beta0 = 0.1,                                  # Intercept
               beta1 = 0.1,                                  # beta associated with x1
               logsigma = 0.1, 
               loglambda = 0.1,
               smoothCoefs = rep(0, ncol(Xtps)),             # Spline coefficients
               logalpha = rep(rep(0.1,length(Sdims))))       # Log spline penalization coefficients (Normal prior)




#=====================================
#             Run the TMB model
#=====================================
obj <- MakeADFun(tmb_data, random = c("smoothCoefs"), tmb_par, hessian = TRUE, DLL="semipar_NTPS_tmbstan")
  return(obj)
}

obj1 <- run_tmb(100)
obj2 <- run_tmb(200)
obj3 <- run_tmb(300)
obj4 <- run_tmb(400)





#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================


M = list()
M[[1]] = list()
M[[1]]$model = "spatial_n10"
M[[2]] = list()
M[[2]]$model = "spatial_n20"
M[[3]] = list()
M[[3]]$model = "spatial_n30"
M[[4]] = list()
M[[4]]$model = "spatial_n40"



M[[1]]$formula = obj1
M[[2]]$formula = obj2
M[[3]]$formula = obj3
M[[4]]$formula = obj4


#===========================================
#           Run the models
#===========================================
  for (i in 1:length(M)){
    startTime <- Sys.time()
    print(paste("Running:  ", M[[i]]$model))
    fit <- tmbstan(M[[i]]$formula,
                   chains= 3, open_progress = FALSE,
                   control = list(max_treedepth= 12,  adapt_delta = 0.95),
                   iter = 4000, warmup= 700, cores=no_cores,
                   init = 'last.par.best', seed = 12345)
    endTime <- Sys.time()
    timeUsed = endTime - startTime
    print(timeUsed)
    saveRDS(fit, file=paste0('tps_fit_', i,'.RDS'))
}



c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")


M[[1]]$res <- readRDS('tps_fit_1.RDS')
M[[2]]$res <- readRDS('tps_fit_2.RDS')
M[[3]]$res <- readRDS('tps_fit_3.RDS')
M[[4]]$res <- readRDS('tps_fit_4.RDS')

summary(M[[1]]$res)


posterior_1 = as.matrix(M[[1]]$res)
posterior_2 = as.matrix(M[[2]]$res)
posterior_3 = as.matrix(M[[3]]$res)
posterior_4 = as.matrix(M[[4]]$res)


traceplot(M[[1]]$res, pars=names(obj1$par), inc_warmup=TRUE)
traceplot(M[[2]]$res, pars=names(obj2$par), inc_warmup=TRUE)



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
params_cp$iter <- 1:6900


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

