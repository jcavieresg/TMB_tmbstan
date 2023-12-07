rm(list = ls())
setwd("C:/Users/Usuario/Desktop/tps_vs_spde/gaussian_model/spde")

library(pacman)
pacman::p_load(geoR, fields, prodlim, TMB, INLA, dplyr, tmbstan, rstan, parallel,
               raster, ggplot2)

options(scipen=999)

# Calculate the number of cores
no_cores <- parallelly::availableCores()  - 1




#=================================================
# Compilamos el modelo y lo cargamos en R
compile("NSPDE_tmbstan.cpp")
dyn.load(dynlib("NSPDE_tmbstan"))


#=========================================
#               Simulation
#=========================================
set.seed(1234)

run_tmb <- function(nloc){
inla.seed = sample.int(n=1E6, size=1)
options(width=70, digits=3)

loc_ini = matrix(c(0,0,10,10, 0, 10, 10, 0), nrow = 4, byrow = T)
mesh_sim = inla.mesh.2d(loc = loc_ini, max.edge=c(2, 5))

mesh_sim$n
plot(mesh_sim)

#spde = inla.spde2.pcmatern(mesh.sim, prior.range = c(.5, .5), prior.sigma = c(.5, .5))
spde = inla.spde2.matern(mesh_sim, alpha = 2)
spde_mat = spde$param.inla[c("M0","M1","M2")]

range = 2
kappa = sqrt(8)/ range
tau = 0.5

Qu = inla.spde.precision(spde, theta=c(log(kappa), log(tau)))
u = inla.qsample(n=1, Q=Qu, seed = inla.seed)
u = u[ ,1]

# Number of locations
n = 100
coords = matrix(runif(2*n), n)*10 # coordinates


A = inla.spde.make.A(mesh = mesh_sim, loc = coords)
u = drop(A %*% u)

#===========================================================================
#                          Simulate data
#===========================================================================
sigma_ini = 0.3
x1 = runif(n, 0, 1)
beta0 = 2
beta1 = 0.5

mu_sim = beta0 + beta1*x1 + u

y_sim = mu_sim + rnorm(n, mean = 0, sd = sigma_ini)

# data.frame
df = data.frame(y_sim = y_sim, s1= coords[ ,1], s2 = coords[ ,2], x1 = x1)
summary(df)


#=================================================================
#                           TMB modelling
#=================================================================

#======================================
#                TMB data
#======================================
tmb_data <-  list(model = 2,                  # 1 = (A*u) / tau,    2 = A*u,  3 = u(site)
                  normalization = 1,          # 1 no scaliing the GMRF, 2 = scalling the GMRF (1/ tau)
                  y    = as.vector(df$y_sim),
                  x1 = as.vector(df$x1),
                  site = mesh_sim$idx$loc - 1,
                  spde_mat = spde_mat,
                  A = A,
                  tau0 = exp(spde$param.inla$theta.initial[1]),
                  kappa0 = exp(spde$param.inla$theta.initial[2]))




#=====================================
#            TMB parameters
#=====================================
tmb_par  <- list(#beta = c(0.1, 0.1),
                 beta0 = 0.1, 
                 beta1 = 0.1,
                 logsigma_e = 0.1,
                 logtau   = 0.1,
                 logkappa = 0.1,
                 u = rnorm(spde$mesh$n,0,1))



#dyn.unload(dynlib("spatial"))

## Create the TMB object
obj <- MakeADFun(data = tmb_data, parameters = tmb_par, random="u", DLL="NSPDE_tmbstan", hessian = T)
return(obj)
}


#obj$env$beSilent()
# lwr <- c(-Inf, -Inf, -Inf, 0, 0)
# upr <- c(Inf, Inf, Inf, Inf, Inf)

#opt <- with(obj, nlminb(par, fn, gr), lower=lwr, upper = upr)
# opt <- with(obj, nlminb(par, fn, gr))
# 
# opt$par
# exp(opt$par[3:5])
# rep  <- obj$report()



obj1 <- run_tmb(nloc = 100)
obj2 <- run_tmb(nloc = 200)
obj3 <- run_tmb(nloc = 300)
obj4 <- run_tmb(nloc = 400)





#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================


M = list()
M[[1]] = list()
M[[1]]$model = "spatial_n100"
M[[2]] = list()
M[[2]]$model = "spatial_n200"
M[[3]] = list()
M[[3]]$model = "spatial_n300"
M[[4]] = list()
M[[4]]$model = "spatial_n400"



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
  saveRDS(fit, file=paste0('spde_fit_', i,'.RDS'))
}


s <- summary(fit_tmbstan, probs = c(0.25, 0.75))
s$summary  # all chaines merged

exp(s$summary[3:5, c(1,3)])
s$summary[1:2, c(1,3)]


# Extract posterior draws for later use
posterior_cp <- as.array(fit_tmbstan)


require(MCMCvis)
library(bayesplot)
mcmc_pairs(posterior_cp, pars = c("logsigma_e","logkappa","logtau"),
           off_diag_args = list(size = 0.75))


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
pairs(fit_tmbstan, pars = c("logsigma", "logtau", "logkappa")) # below the diagonal




# Simulating from the SIMULATE function
mat_sim = matrix(data=NA, nrow=length(obj$simulate()$y_sim), ncol=100)
mat_sim


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = obj$simulate()$y_sim
  }
}
mat_sim

# Plot lognormal non spatial simulations
hist(y, col = "gray90", prob = TRUE, main = "100 simulated samples", cex.main = 2, xlab = "", col.lab = 'blue', col.main="blue", cex.lab = 1.4, cex.axis = 1.2, ylim = c(0, 0.5),
     xlim = c(-5, 5))
for (j in 1: ncol(mat_sim)){
  lines(x = density(x = mat_sim[, j]),  lty="dotted", col="azure4", lwd=1)
}
#legend("topright", "A", bty = "n", text.col="blue", cex = 1.4)
legend(6, 0.3, c("response", "simulations"), lwd=4, col=c("gray90", "azure4"), cex = 1.5)




