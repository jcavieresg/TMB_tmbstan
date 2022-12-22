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
library(parallelly)
library(tweedie)

options(scipen=999)

# Calculate the number of cores
no_cores <- parallelly::availableCores() 

#=================================================
# Compilamos el modelo y lo cargamos en R
compile("tweedie_spatial.cpp")
dyn.load(dynlib("tweedie_spatial"))
#=================================================

#=============
# EJEMPLO 1 
#=============
set.seed(456)
to_plot = TRUE   # generate plots

# Build a grid of locations
x <- seq(1, 10, length = 30)
grid.pts <- expand.grid(x, x)
dim(grid.pts)


# Creamos una variable respuesta 0 inicial
z = 0


# Build this data frame ONLY for mgcv setup
df <- data.frame(grid.pts, z)
names(df) <- c("s1", "s2", "z")
head(df)
dim(df)

xt <- unique(df[c("s1", "s2")])
nrow(xt)

#====================================================
#      Obtenemos matriz penalizada desde mgcv
#====================================================
# mgcv setup
#m0 <- gam(z ~ x0 + s(s1, s2, bs = "tp", k = 50), data = df, fit = FALSE) 
m0 <- gam(z ~ s(s1, s2, bs = "tp", k = 50), data = df, fit = FALSE) 

Xtps <- m0$X[, -1]
dim(Xtps)
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
beta0_ini  <- 2
phi_ini<- 2                          # Desvest
p_ini  <- 1.5                          # slant parameter of the Skew Normal
alpha_ini  <- 1
E = ginv(Stps)

smoothCoefs <- rmvn(n = 1, mu = rep(0, nrow(Stps)),  V = alpha_ini*E) # From Wood 2011 (apendix)

# mu_sim
mu_sim <- beta0_ini + Xtps %*% smoothCoefs # mean value in each location

# y real Tweddie
y_sim = rTweedie(mu = mu_sim, p = p_ini, phi = p_ini)
hist(y_sim)


# Tps simulated
tps_sim = Xtps %*% smoothCoefs


tps.rast <- raster(nrows=length(x), ncols=length(x),
                   xmn=1, xmx=10, ymn=1, ymx=10,
                   vals=(tps_sim))

tps.rast@data


## Define cluster locations and sample size at each
n.clust <- nrow(xt)                   # This number must be equal to square grid
clust.mean.ss <- 50              # Mean of observations by cluster (site)
dat <- data.frame(x = runif(n.clust, min = min(df$s1), max = max(df$s1)),
                  y = runif(n.clust, min = min(df$s2), max = max(df$s2)),
                  n = rpois(n.clust, clust.mean.ss))

latent.truth <- raster::extract(x = tps.rast, y = cbind(dat$x, dat$y))
p.truth <- plogis(latent.truth)

dat$latent.truth <- latent.truth
dat$p.truth <- p.truth

dat$y_sim <- y_sim

# Head of main dataset simulated
head(dat)

#par(mfrow = c(3, 1))
plot(tps.rast, maxpixels = length(x) ^ 2,
     xlim = range(dat$x), ylim = range(dat$y), main = 'latent truth')

fields::quilt.plot(dat$x, dat$y, dat$y_sim, main = 'empirical Skew Normal observations')

# Summary of the simulated data
ggplot(dat, aes(x=x, y=y)) + geom_point(size=3, shape = 1) +
  theme(legend.title = element_text(size=18),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18, face = "bold")) +
  theme(plot.margin = unit(c(0.5,4,1,3), "cm"))




#=================================================================
#                           TMB modelling
#=================================================================

#======================================
#                TMB data
#======================================
tmb_data = list(y     = dat$y_sim,                           # Response
                X     = Xtps,                                # Design matrix, without intercept and betas
                S     = as(S_combined, "dgTMatrix"),         # Combined penalty matrix
                Sdims = Sdims,
                tpsReport = .bdiag(tpsReport))
                



#=====================================
#            TMB parameters
#=====================================
tmb_par = list(beta0 = 0.5,                                  # Intercept
               phi = 1,
               p = 1.5,
               smoothCoefs = rep(0, ncol(Xtps)),             # Spline coefficients
               logalpha = rep(rep(0.1,length(Sdims))))       # Log spline penalization coefficients (Normal prior)




#=====================================
#             Run the model
#=====================================
obj = MakeADFun(data = tmb_data, parameters = tmb_par, random = "smoothCoefs", DLL = "tweedie_spatial", hessian = TRUE)

# lwr <- c(-Inf, -Inf, 1, -Inf)
# upr <- c(Inf, Inf, 2, Inf)

tic("TMB Tweddie model fitting")
#opt = nlminb(obj$par,obj$fn,obj$gr, lower=lwr, upper=upr)
opt = nlminb(obj$par,obj$fn,obj$gr)
toc()

rep = sdreport(obj)

param = cbind(opt$par[1], exp(opt$par[2]), exp(opt$par[3]), exp(opt$par[4]))
colnames(param) = c("beta0", "phi", "p", "alpha")
param












#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================
init.fn <- function()
  split(unname(obj$env$last.par.best),names(obj$env$last.par.best))

tic("Time of estimation")
fit = tmbstan(obj, chains = 3, open_progress = FALSE,# lower = lwr, upper = upr,
              init = init.fn, control = list(max_treedepth = 12, adapt_delta = 0.9), 
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

pairs(fit, pars = c("phi", "p")) # below the diagonal


params_cp <- as.data.frame(fit)
names(params_cp) <- gsub("chain:1.", "", names(params_cp), fixed = T)
names(params_cp) <- gsub("[", ".", names(params_cp), fixed = T)
names(params_cp) <- gsub("]", "", names(params_cp), fixed = T)
params_cp$iter <- 1:4000


# phi
par(mfrow=c(3,2),mar=c(4,4,0.5,0.5), oma=c(2,3,1,1))
plot(params_cp$iter, params_cp$phi, col=c_dark, pch=16, cex=0.8, type = "l",
     xlab="Iteration", ylab="phi", cex.lab=1.3, cex.axis=1.3)

running_means_phi = sapply(params_cp$iter, function(n) mean(params_cp$phi[1:n]))
plot(params_cp$iter, running_means_phi, col=c_dark, pch=16, cex=0.8,  cex.lab=1.3, cex.axis=1.3,
     xlab="Iteration", ylab="MCMC mean of phi")
abline(h=mean(running_means_phi), col="grey", lty="dashed", lwd=3)



# p
plot(params_cp$iter, params_cp$p, col=c_dark, pch=16, cex=0.8, type = "l",
     xlab="Iteration", ylab="p",  cex.lab=1.3, cex.axis=1.3)

running_means_p = sapply(params_cp$iter, function(n) mean(params_cp$p[1:n]))
plot(params_cp$iter, running_means_p, col=c_dark, pch=16, cex=0.8,  cex.lab=1.3, cex.axis=1.3,
     xlab="Iteration", ylab="MCMC mean of p")
abline(h=mean(running_means_p), col="grey", lty="dashed", lwd=3)


# logalpha
plot(params_cp$iter, params_cp$logalpha, col=c_dark, pch=16, cex=0.8, type = "l",
     xlab="Iteration", ylab="logalpha",  cex.lab=1.3, cex.axis=1.3)

running_means_alpha = sapply(params_cp$iter, function(n) mean(params_cp$logalpha[1:n]))
plot(params_cp$iter, running_means_alpha, col=c_dark, pch=16, cex=0.8,  cex.lab=1.3, cex.axis=1.3,
     xlab="Iteration", ylab="MCMC mean of logalpha")
abline(h=mean(running_means_alpha), col="grey", lty="dashed", lwd=3)
mtext("Convergence of the parameters phi, p and alpha", outer=l,  cex=1, line=-0.5)




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


# Variancean
library(MCMCvis)
color_scheme_set("viridisE")
mcmc_rank_hist(fit, pars = c("logalpha", "logsigma", "loglambda"), ref_line = l)   # spatial random field


## Extract marginal posteriors
posterior <- as.matrix(fit)

mean(posterior[, "beta0"])
mean(posterior[, "phi"])
mean(posterior[, "p"])
exp(mean(posterior[, "logalpha"]))


dyn.unload(dynlib("tweedie_spatial"))
