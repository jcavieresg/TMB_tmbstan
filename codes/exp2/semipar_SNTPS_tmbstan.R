rm(list = ls())
setwd("C:/Users/joaquin/Desktop/TPS")
library(TMB)
library(mgcv)
library(sn)
library(MASS)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(Matrix)
library(gridExtra)
library(grid)
library(lattice)
library(tictoc)
library(tmbstan)
library(rstan)
library(tmbstan)
library(rstan)
library(bayesplot)
library(ggplot2)
library(geoR)
library(raster)


options(scipen=999)

#=================================================
# Compilamos el modelo y lo cargamos en R
TMB::compile('semipar_SNTPS_tmbstan.cpp')
dyn.load(dynlib("semipar_SNTPS_tmbstan"))
# #=================================================


set.seed(456)
to_plot = TRUE   # generate plots

# Build a grid of locations
x <- seq(1, 10, length = 20)
grid.pts <- expand.grid(x, x)
dim(grid.pts)

# Creamos una variable respuesta 0 inicial
z = 0

# Simulamos una variable continua x1
x0 = runif(length(x), 0, 1)


# Build this data frame ONLY for mgcv setup
df <- data.frame(grid.pts, z, x0)
names(df) <- c("s1", "s2", "z", "x0")
head(df)
dim(df)

xt <- unique(df[c("s1", "s2")]) 
nrow(xt)

#====================================================
#      Obtenemos matriz penalizada desde mgcv
#====================================================
# mgcv setup
m0 <- gam(z ~ x0 + s(s1, s2, bs = "tp", k = 50), data = df, fit = FALSE) 

Xtps <- m0$X[, c(-1, -2)]                    # Matricial form without intercept and the parameters asociated with the covariates
Stps <- m0$smooth[[1]]$S[[1]]            # Extrtact penelization matrices
dim(Stps)

S_list = list(Stps)
S_combined = .bdiag(S_list)                               # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow))                       # Find dimension of each S (in this case only one)
#----------------------------------------

#===========================================================================
# For report, used for constructing plots
s1 = seq(min(df$s1), max(df$s1), by = 0.005)
s2 = seq(min(df$s2), max(df$s2), by = 0.005)


tpsReport = PredictMat(m0$smooth[[1]], data = data.frame(s1, s2))
tpsReport = list(tpsReport)
#-------------------------------------------


#===========================================================================
#                          Simulated data
#===========================================================================
# Simule new data
n <- length(df$z)
beta0_ini  <- 1.5
beta1_ini  <- 2
sigma_ini   <- 1                          # Desvest
lambda_ini <- 3.0                          # slant parameter of the Skew Normal
alpha_ini  <- 1
E = ginv(Stps)


smoothCoefs <- rmvn(n = 1, mu = rep(0, nrow(Stps)),  V = alpha_ini*E) # From Wood 2011 (apendix)

# Simulated covariate
x1 = runif(n, 0, 1)                     

# yhat
mu_sim <- beta0_ini + beta1_ini*x1 + Xtps %*% smoothCoefs # mean value in each location


# y real Skew Normal
y_sim <- rsn(n, xi = mu_sim, omega = sigma_ini^2, alpha = lambda_ini)

# Histogram
hist(y_sim)

# Tps simulated
tps_sim = Xtps %*% smoothCoefs


tps.rast <- raster(nrows=length(x), ncols=length(x),
                   xmn=1, xmx=10, ymn=1, ymx=10,
                   vals=(tps_sim))

tps.rast@data


## Define cluster locations and sample size at each
n.clust <- 400                   # This number must be equal to square grid
clust.mean.ss <- 50              # Mean of observations by cluster (site)
dat <- data.frame(x = runif(n.clust, min = min(df$s1), max = max(df$s1)),
                  y = runif(n.clust, min = min(df$s2), max = max(df$s2)),
                  n = rpois(n.clust, clust.mean.ss))

latent.truth <- raster::extract(x = tps.rast, y = cbind(dat$x, dat$y))
p.truth <- plogis(latent.truth)

dat$latent.truth <- latent.truth
dat$p.truth <- p.truth

dat$y_sim <- y_sim
dat$x1 <- x1


# Head of main dataset simulated
head(dat)


# Summary of the simulated data
ggplot(dat, aes(x=x, y=y)) + geom_point(size=3, shape = 1) +
  theme(legend.title = element_text(size=18),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18, face = "bold")) +
  theme(plot.margin = unit(c(0.5,4,1,3), "cm"))


#par(mfrow = c(3, 1))
plot(tps.rast, maxpixels = length(x) ^ 2,
     xlim = range(dat$x), ylim = range(dat$y), main = 'latent truth')

fields::quilt.plot(dat$x, dat$y, dat$y_sim, main = 'empirical Skew Normal observations')





#=================================================================
#                           TMB modelling
#=================================================================

#======================================
#                TMB data
#======================================
tmb_data = list(likelihood = 1,                              # 1 = Skew-Normal, 2 = Normal, 3 = tweddie
                y     = as.vector(dat$y_sim),                 # Response
                X    = model.matrix(~ 1 + x1),
                TPS     = Xtps,                                # Design matrix, without intercept and betas
                S     = as(S_combined, "dgTMatrix"),                          # Combined penalty matrix
                Sdims = Sdims,
                tpsReport = .bdiag(tpsReport))



#=====================================
#            TMB parameters
#=====================================
tmb_par = list(beta = c(0.1, 0.1),
               smoothCoefs = rep(0, ncol(Xtps)),             # Spline coefficients
               logalpha = rep(rep(0.1,length(Sdims))),      # Log spline penalization coefficients (Normal prior)
               logsigma = exp(0.1),
               loglambda = exp(0.1))



#=====================================
#             Run the model
#=====================================
# Skew normal
obj <- MakeADFun(tmb_data, random = c("smoothCoefs"), tmb_par, DLL="semipar_SNTPS_tmbstan", hessian = TRUE)
opt = nlminb(obj$par,obj$fn,obj$gr)
rep = sdreport(obj)
#-------------------------------------------

opt$convergence==0 # El modelo esta convergiendo

# Report
obj$report()

param = cbind(opt$par[1], opt$par[2], exp(opt$par[3]), exp(opt$par[4]), exp(opt$par[5]))
colnames(param) = c("beta0", "beta1", "alpha", "sigma", "lambda")
param



opt$objective
obj$fn(obj$env$last.par.best[-obj$env$random])


# Likelihood profile for 
prof_beta0 <- tmbprofile(obj, "beta0")
prof_beta1 <- tmbprofile(obj, "beta1")
prof_alpha <- tmbprofile(obj, "logalpha")
prof_sigma2 <- tmbprofile(obj, "logsigma")
prof_lambda <- tmbprofile(obj, "loglambda")

par(mfrow=c(2,2), mar = c(3, 5, 4, 4))
plot(prof_beta0, main = expression(paste("Profile likelihood ", beta[0])), col = "blue",  cex.lab = 1.6,
     cex.axis = 1.4, cex.main = 1.6, cex = 1.6, lwd = 2, xlab = "")
# plot(prof_beta1, main = expression(paste("Profile likelihood ", beta[0])), col = "blue", cex.lab = 1.6,
#      cex.axis = 1.4, cex.main = 1.6, lwd = 2, xlab = "")
plot(prof_alpha, main = expression(paste("Profile likelihood ", alpha)), col = "blue", cex.lab = 1.6,
     cex.axis = 1.4, cex.main = 1.6, lwd = 2, xlab = "", ylab = "")
plot(prof_sigma2, main = expression(paste("Profile likelihood ", sigma^2)), col = "blue", cex.lab = 1.6,
     cex.axis = 1.4, cex.main = 1.6, lwd = 2, xlab = "")
plot(prof_lambda, main = expression(paste("Profile likelihood ", lambda)), col = "blue", cex.lab = 1.6,
     cex.axis = 1.4, cex.main = 1.6, lwd = 2, xlab = "", ylab = "")



print(confint(prof, level = 0.8))
print(confint(prof, level = 0.9))
print(confint(prof, level = 0.95))


#========================================
#              Plot results
#========================================
library(plot3D)

muSpline <- rep$value[names(rep$value) == "splines2D"]
sdSpline <- rep$sd[names(rep$value)=="splines2D"]

# #par(mfrow=c(1,2))

# Par?metros obtenidos mediante SNTPS
n <- length(dat$y_sim)
beta0_est <- as.numeric(opt$par[1])
beta1_est <- as.numeric(opt$par[2])
sigma2_est <- as.numeric(exp(opt$par[4])^2)
lambda_est <- as.numeric(exp(opt$par[5]))

vies = sqrt(2/pi)*sqrt(sigma2_est)*lambda_est/sqrt(1+lambda_est^2)

modelmat <- model.matrix(~ dat$x1)
tps_est = dat$y_sim - (modelmat %*% opt$par[1:2] + vies)

scatter3D(dat$x, dat$y, tps_est, bty = "g", pch = 19, cex = 1, ticktype = "detailed", theta = 40, phi = 10)


library(lattice)
wireframe(tps_est ~ s1*s2, data = df,
          xlab = "X Coordinate (feet)", ylab = "Y Coordinate (feet)",
          main = "Surface Thin Plate Spline",
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -80, x = -60)
)


wireframe(tps_est ~ s1*s2, data = df,
          shade = TRUE, aspect = c(1, 1),
          light.source = c(10,10,10), main = "Study 1",
          scales = list(z.ticks=5,arrows=FALSE, col="black", font=10, tck=0.5),
          screen = list(z = 60, x = -70, y = 0))


v <- ggplot(df, aes(s1, s2, z = tps_est))
v + geom_contour_filled()





# Plot prediction
xhat = summary(rep, "random")[,1]
#par_fixed = summary(rep, "fixed")
# beta0_hat = summary(rep, "fixed")["beta0", 1]
df$y_pred = obj$report()$preds
df$tps_sim <- tps_sim
#df$field_est = beta0_hat + Xtps %*% xhat
df$tps_pred = Xtps %*% xhat


#=======================================================
# ggplots


# y_real vs y_pred

# limits for the response variable and prediction
xy_lims = range(c(df$y_sim, df$y_pred))

# y_true
y_true =  ggplot(df) +
  geom_tile(aes(x = s1, y = s2, fill = y_sim)) +
  coord_equal() +
  scale_fill_viridis_c(limits = xy_lims) +
  theme(legend.title = element_text(size=20),
        axis.title.x=element_blank(),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18, face="bold")) +
  theme(legend.text = element_text(size=14))+
  theme(plot.margin = unit(c(0.2,0,0.2,0), "cm"))

# y_pred
y_pred =  ggplot(df) +
  geom_tile(aes(x = s1, y = s2, fill = y_pred)) +
  coord_equal() +
  scale_fill_viridis_c(limits = xy_lims) + 
  theme(legend.title = element_text(size=20),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18, face="bold")) +
  theme(legend.text = element_text(size=14))+
  theme(plot.margin = unit(c(0.2,0,0.2,0), "cm"))



# limits for the simulated and estimated field
scl_lims = range(c(df$tps_sim, df$tps_pred))

# field_true
field_true =  ggplot(df) +
  geom_tile(aes(x = s2, y = s1, fill = tps_sim)) +
  coord_equal() +
  scale_fill_viridis_c(limits = scl_lims) + 
  theme(legend.title = element_text(size=20),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18, face="bold")) +
  theme(legend.text = element_text(size=14))+
  theme(plot.margin = unit(c(0.2,0,0.2,0), "cm"))

# field_est
field_pred =  ggplot(df) +
  geom_tile(aes(x = s2, y = s1, fill = tps_pred)) +
  coord_equal() +
  scale_fill_viridis_c(limits = scl_lims)+ 
  theme(legend.title = element_text(size=20),
        axis.title.y=element_blank(),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18, face="bold")) +
  theme(legend.text = element_text(size=14))+
  theme(plot.margin = unit(c(0.2,0,0.2,0), "cm"))



grid.arrange(y_true, y_pred, field_true, field_pred, ncol = 2,
             top = textGrob("Skew-Normal model for n = 400" , gp=gpar(fontsize=20, font=1))) 

  
  
# grid.arrange(y_true_sim100, y_true_sim900, y_true_sim1600,ncol = 3) 
#   
# 
# grid.arrange(y_true, y_pred, field_true, field_pred, ncol = 2,
#              top = textGrob("Simulation for n = 100" , gp=gpar(fontsize=20, font=1))) 



# Package with others metrics
library(APMtools)
APMtools::error_matrix(validation = as.numeric(dat$y_sim), prediction = obj$report()$preds)

# Bias

library(SimDesign)
bias(as.numeric(dat$y_sim), obj$report()$preds)


#==================
#   Marginal AIC
#==================
AIC <- TMBAIC(opt = opt, p = 2, n = Inf); AIC

#==================
#       BIC
#==================
# La formula es -2*l(\theta), pero como estamos trabajando con la neg-lohlikelihood, le quitamos el -
BIC <- 2*opt$objective + log(length(dat$y_sim))*(length(opt$par)); BIC






#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================
# Calculate the number of cores
no_cores <- parallelly::availableCores()  - 1

init.fn <- function()
  split(unname(obj$env$last.par.best),names(obj$env$last.par.best))


startTime <- Sys.time()
fit = tmbstan(obj, chains = 3, open_progress = FALSE, 
              init = init.fn, control = list(max_treedepth = 12, adapt_delta = 0.9), 
              iter = 3000, warmup=1000, cores = no_cores, seed=483892929)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)

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






