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
compile("hurdle_spatial_tps.cpp")
dyn.load(dynlib("hurdle_spatial_tps"))
#=================================================




## Load data and prepare for model
load('EBS_pollock_data.rda')
pollock <- EBS_pollock_data
dim(pollock)

# Subset of the year 2014
pollock <- subset(pollock, year==2014)
hist(pollock$catch)


#====================================================
#      Penalized matrices from mgcv
#====================================================
# Artificial response variable
z <- rep(0, length(pollock$catch))
pollock$z <- z
head(pollock)


#===================
#     mgcv setup
#===================
m0 <- gam(z ~ s(lat, long, bs = "tp", k = 75), data = pollock, fit = FALSE) 


Xtps <- m0$X[, c(-1)]                    # Matricial form without intercept and the parameters asociated with the covariates
Stps <- m0$smooth[[1]]$S[[1]]                # Extrtact penelization matrices
dim(Stps)

S_list = list(Stps)
S_combined = .bdiag(S_list)                               # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow))                       # Find dimension of each S (in this case only one)
#----------------------------------------


#===========================================================================
# For report, used for constructing plots
lat = seq(min(pollock$lat), max(pollock$lat), by = 0.38)
long = seq(min(pollock$long), max(pollock$long), by = 0.5)

# by cambia para hacer coincidir los largos de los vectores
tpsReport = PredictMat(m0$smooth[[1]], data = data.frame(lat, long))
tpsReport = list(tpsReport)

#=================================================================
#                           TMB modelling
#=================================================================

#======================================
#                TMB data
#======================================
tps_data = list(likelihood=2,
                y     = pollock$catch,                           # Response
                X     = Xtps,                                # Design matrix, without intercept and betas
                S     = as(S_combined, "dgTMatrix"),         # Combined penalty matrix
                Sdims = Sdims,
                tpsReport = .bdiag(tpsReport))
                



#=====================================
#            TMB parameters
#=====================================
tps_par = list(intercept = 0.5,                                  # Intercept
               theta = -3,
               logsigma = 0.1,
               smoothCoefs = rep(0, ncol(Xtps)),             # Spline coefficients
               logalpha = rep(rep(0.1,length(Sdims))))       # Log spline penalization coefficients (Normal prior)




#=====================================
#             Run the model
#=====================================
obj = MakeADFun(data = tps_data, parameters = tps_par, random = "smoothCoefs", DLL = "hurdle_spatial_tps", hessian = TRUE)

tic("TPS hurdle model")
opt <- nlminb(obj$par,obj$fn,obj$gr)
toc()

# opt <- with(obj, nlminb(par, fn, gr))
# opt <- with(obj, nlminb(opt$par, fn, gr)) # restart
rep <- obj$report()

#rep = sdreport(obj)

param_tps = cbind(opt$par[1], opt$par[2],  exp(opt$par[3]), exp(opt$par[4]))
colnames(param_tps) = c("intercept", "theta", "sigma", "alpha")
param_tps






#=================================================
# Compilamos el modelo y lo cargamos en R
compile("cpue_spatial_spde.cpp")
dyn.load(dynlib("cpue_spatial_spde"))
#=================================================

library(INLA)

#Define mesh and components representing the  precision matrix----
#loc = cbind(pollock$lat, pollock$long)
# boundary = INLA::inla.nonconvex.hull(loc)
# boundary2 = INLA::inla.nonconvex.hull(loc,convex = -0.35)
#mesh = inla.mesh.2d(loc = loc,   boundary = list(boundary,boundary2), max.edge=c(9, 10), cutoff=0.05)
# mesh = inla.mesh.2d(loc = loc,  max.edge=c(5, 10))
# mesh$n
# plot(mesh)
# points(pollock$lat, pollock$long, col = "blue")
# 
# A = inla.spde.make.A(mesh,loc)
# 
# spde = inla.spde2.matern(mesh, alpha=2)


xt <- unique(pollock[c("lat", "long")]) 
nrow(xt)

## The SPDE approach to simplify the model. Using RINLA tools to create
## inputs to TMB. Based of Jim Thorson's function.
create.spde <- function(dat, n_knots, make.plot=FALSE, jitter=.3){
  loc_xy <- data.frame(x=dat$long, y=dat$lat)
  knots <- kmeans( x=loc_xy, centers=n_knots )
  loc_centers <- knots$centers
  ## loc_xy <- cbind(loc_xy, cluster=knots$cluster)
  ## ggplot(loc_xy, aes(x,y, col=factor(cluster))) + geom_point(size=.5)
  #mesh <- inla.mesh.create( loc_centers, refine=TRUE)
  mesh = inla.mesh.2d(loc = loc_centers,  max.edge=c(5, 10))
  plot(mesh)
  spde <- inla.spde2.matern( mesh )
  if(make.plot){
    png(paste0('mesh_', n_knots, '.png'), width=7, height=4, units='in', res=500)
    par(mar=.5*c(1,1,1,1))
    plot(mesh, main=NA, edge.color=gray(.7))
    points( jitter(dat$longitude, amount=jitter), jitter(dat$latitude, amount=jitter), cex=1, pch='.', col=rgb(0,0,0,.3))
    points( loc_centers, cex=.5, pch=20, col="red")
    dev.off()
  }
  return(list(mesh=mesh, spde=spde, cluster=knots$cluster, loc_centers=loc_centers))
}

spde <- create.spde(dat=pollock, n_knots=75, make.plot=T)
spdeMatrices = spde$spde$param.inla[c("M0","M1","M2")]



# mesh = inla.mesh.2d(loc = cbind(pollock$lat, pollock$long),  max.edge=c(5, 10))
# spde = inla.spde2.matern(mesh, alpha=2)
# spdeMatrices = spde$param.inla[c("M0","M1","M2")]


  
# spde$mesh$n
# mesh$n
# 
# plot(spde$mesh)
# plot(mesh)
# 
# 
# loc_xy <- data.frame(x=pollock$long, y=pollock$lat)
# knots <- kmeans( x=loc_xy, centers= 75)


## Make spde inputs for TMB
spde_data <- list(likelihood=2, 
             y=pollock$catch, 
             #site = mesh$idx$loc - 1,
             site = spde$cluster,
             spdeMatrices = spdeMatrices)


spde_par <- list(intercept=3, 
             theta=-3, 
             logsigma=0,
             logtau=1, 
             logkappa=-3,
             u=rnorm(spde$mesh$n,0,1))




## Rerun with spatial SPDE approach
obj2 = MakeADFun(data = spde_data, parameters = spde_par, random = "u", DLL = "cpue_spatial_spde", hessian = TRUE)

tic("SPDE hurdle model")
opt2 <- nlminb(obj2$par,obj2$fn,obj2$gr)
toc()

# opt2 <- with(obj2, nlminb(par, fn, gr))
# opt2 <- with(obj2, nlminb(opt2$par, fn, gr)) # restart
# opt2 <- with(obj2, nlminb(opt2$par, fn, gr)) # restart
rep2 <- obj2$report()


param_spde = cbind(opt2$par[1], opt2$par[2],  exp(opt2$par[3]), exp(opt2$par[4]))
colnames(param_spde) = c("intercept", "theta", "sigma", "kappa")
param_spde


param_tps
param_spde


library(TMBhelper)
AIC_tps <- TMBAIC(opt = opt)
AIC_spde <- TMBAIC(opt = opt2)
        
AIC_tps
AIC_spde


library(reshape2)
## Get residuals for all three models
pollock$resids.tps <- (rep$pred-log(tps_data$y))/exp(rep$logsigma)
pollock$resids.spde <- (rep2$pred-log(spde_data$y))/exp(rep2$logsigma)
resids <- melt(pollock, id.vars = c("lat", "long"), measure.vars = c("resids.tps", "resids.spde"))
ggplot(resids, aes(long, lat, size=abs(value), color=value)) +
  geom_point(alpha=.5)+ scale_color_gradient(low='red', high='black') +
  facet_wrap('variable') + scale_size_area(max_size=4)



plot(log(pollock$catch))
lines(obj$report()$pred, col = "red")
lines(obj2$report()$pred, col = "blue")
















#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================

# TPS MODEL
init.fn <- function()
  split(unname(obj$env$last.par.best),names(obj$env$last.par.best))

tic("Time of estimation TPS model")
fit_tps = tmbstan(obj, chains = 3, open_progress = FALSE,# lower = lwr, upper = upr,
              init = init.fn, control = list(max_treedepth = 12, adapt_delta = 0.95), 
              iter = 5000, warmup=1000, cores=no_cores, seed=483892929)
toc()





# SPDE MODEL
init.fn2 <- function()
  split(unname(obj2$env$last.par.best),names(obj2$env$last.par.best))


tic("Time of estimation SPDE model")
fit_spde = tmbstan(obj2, chains = 3, open_progress = FALSE,# lower = lwr, upper = upr,
              init = init.fn2, control = list(max_treedepth = 12, adapt_delta = 0.95), 
              iter = 5000, warmup=1000, cores=no_cores, seed=483892929)
toc()



# saveRDS(fit_tps, 'fit_tps.RDS')
# saveRDS(fit_spde, 'fit_spde.RDS')


c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")


traceplot(fit_tps, pars=names(obj$par), inc_warmup=TRUE)
traceplot(fit_spde, pars=names(obj2$par), inc_warmup=TRUE)

#load package
require(MCMCvis)
MCMCsummary(fit_tps, round = 2)
MCMCsummary(fit_spde, round = 2)

pairs(fit, pars = c("intercept", "theta", "logsigma")) # below the diagonal


params_cp <- as.data.frame(fit_tps)
names(params_cp) <- gsub("chain:1.", "", names(params_cp), fixed = T)
names(params_cp) <- gsub("[", ".", names(params_cp), fixed = T)
names(params_cp) <- gsub("]", "", names(params_cp), fixed = T)
params_cp$iter <- 1:12000


# Intercept
par(mfrow=c(3,2),mar=c(4,4,0.5,0.5), oma=c(2,3,1,1))
plot(params_cp$iter, params_cp$intercept, col=c_dark, pch=16, cex=0.8, type = "l",
     xlab="Iteration", ylab="Intercept", cex.lab=1.3, cex.axis=1.3)

running_means_intercept = sapply(params_cp$iter, function(n) mean(params_cp$intercept[1:n]))
plot(params_cp$iter, running_means_intercept, col=c_dark, pch=16, cex=0.8,  cex.lab=1.3, cex.axis=1.3,
     xlab="Iteration", ylab="MCMC mean of Intercept")
abline(h=mean(running_means_intercept), col="grey", lty="dashed", lwd=3)



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
divergent = get_sampler_params(fit_tps, inc_warmup=FALSE)[[1]][,'divergent__']
sum(divergent)

divergent = get_sampler_params(fit_spde, inc_warmup=FALSE)[[1]][,'divergent__']
sum(divergent)

## Methods provided by 'rstan'
class(fit)
methods(class ="stanfit")

get_cppo_mode(fit_tps)
get_stancode(fit)
get_stanmodel(fit)
log_prob(fit)
loo(fit)
traceplot(fit_tps)

#launch_shinystan(fit)

## ESS and Rhat from rstan::monitor
mon = monitor(fit_tps)
max(mon$Rhat)
min(mon$Tail_ESS)

# evalaute problem of convergence
sum(mon$Rhat > 1.01)
sum(mon$Tail_ESS < 400)

source('monitornew.R')
source('monitorplot.R')
source('stan_utility.R')

which_min_ess = which.min(mon[1:200, 'Tail_ESS'])
plot_local_ess(fit = fit_tps, par = which_min_ess, nalpha = 10)

plot_quantile_ess(fit = fit_tps, par = which_min_ess, nalpha = 50)

plot_change_ess(fit = fit_tps, par = which_min_ess)

check_rhat(fit_spde)
check_treedepth(fit_spde, 10)
check_energy(fit)  #
check_div(fit)


# Variancean
library(MCMCvis)
library(bayesplot)
color_scheme_set("viridisE")
mcmc_rank_hist(fit_tps, pars = c("intercept", "theta", "logsigma"), ref_line = T)   # spatial random field


## Extract marginal posteriors
posterior_1 <- as.matrix(fit_tps)
posterior_2 <- as.matrix(fit_spde)

mean(posterior[, "intercept"])
mean(posterior[, "theta"])
exp(mean(posterior[, "logsigma"]))
exp(mean(posterior[, "logalpha"]))


# dyn.unload(dynlib("hurdle_spatial_tps"))
# dyn.unload(dynlib("cpue_spatial_spde"))



fit_spde <- readRDS("fit_spde.rds")
posterior_2 <- as.matrix(fit_spde)


# loo for lognormal model without spatial random effect
set.seed(1000)
n_1 = length(pollock$catch)
log_lik_2 <- matrix(NA, nrow=nrow(posterior_2), ncol=n_1)

for(i in 1:nrow(posterior_2)){
  r2 <- obj$report(posterior_2[i,-ncol(posterior_2)])
  log_lik_2[i,] <- r$log_lik
}

loo_1 = loo(log_lik_1, r_eff = NA)

obj$report(posterior_1[1,-ncol(posterior_1)])         # sd0 is only element
log_like1 <- rep(NA, len=nrow(posterior_1))
for(i in 1:nrow(posterior_1)){
  r <- obj$report(posterior_1[i,-ncol(posterior_1)])
  log_like1[i] <- r$sd0
}
hist(sd0)


# loo for lognormal model with spatial random effect
set.seed(1000)
n_2 = length(data_full$cpue)
log_lik_2 <- matrix(NA, nrow=nrow(posterior_2), ncol=n_2)

for(i in 1:nrow(posterior_2)){
  r2 <- obj2$Obj$report(posterior_2[i,-ncol(posterior_2)])
  log_lik_2[i,] <- r2$log_lik
}

loo_2 = loo(log_lik_2, r_eff = NA)
print(loo_2)


# loo for Gamma model without spatial random effect
set.seed(1000)
n_3 = length(data_full$cpue)
log_lik_3 <- matrix(NA, nrow=nrow(posterior_3), ncol=n_3)

for(i in 1:nrow(posterior_3)){
  r3 <- obj3$Obj$report(posterior_3[i,-ncol(posterior_3)])
  log_lik_3[i,] <- r3$log_lik
}

loo_3 = loo(log_lik_3, r_eff = NA)
print(loo_3)


# loo for Gamma model with spatial random effect
set.seed(1000)
n_4 = length(data_full$cpue)
log_lik_4 <- matrix(NA, nrow=nrow(posterior_4), ncol=n_4)

for(i in 1:nrow(posterior_4)){
  r4 <- obj4$Obj$report(posterior_4[i,-ncol(posterior_4)])
  log_lik_4[i,] <- r4$log_lik
}

loo_4 = loo(log_lik_4, r_eff = NA)
print(loo_4)



#================================================================
#               loo comparision for the 6 models
#================================================================
comp <- loo_compare(loo_1, loo_2, loo_3, loo_4)
print(comp, digits = 3)

par(mfrow=c(2,2))
plot(loo_1, main =  expression(bold('Model'["1,1"])), xlab="")
plot(loo_2, main =  expression(bold('Model'["1,2"])))
plot(loo_3, main =  expression(bold('Model'["2,1"])))
plot(loo_4, main =  expression(bold('Model'["2,2"])))


khat <- pareto_k_values(loo_4)
kdata <- data.frame(Data_point = seq_along(khat), khat = khat, bad = khat > 0.7)
ggplot(kdata, aes(x = Data_point, y = khat, label = Data_point)) + 
  geom_hline(yintercept = 0.7, size = 0.25, linetype = 2, color = "blue") +
  geom_point(data = subset(kdata, !bad), size = 0.5) + 
  geom_text(data = subset(kdata, bad), color = "red")

