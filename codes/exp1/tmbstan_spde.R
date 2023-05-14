rm(list = ls())
setwd("C:/Users/joaquin/Desktop/new_models")

library(pacman)
pacman::p_load(RandomFields, geoR, fields, prodlim, TMB, INLA, dplyr, tmbstan, rstan, parallel)

options(scipen=999)

# Calculate the number of cores
no_cores <- parallelly::availableCores()  - 1


## From package 'geoR':
# Matern correlation
matern_fun <- function(h, nu, kappa) {
  ifelse(h > 0, besselK(h * kappa, nu) * (h * kappa)^nu / (gamma(nu) * 2^(nu - 1)), 1)
}

## Test:
set.seed(123)
n <- 50
loc <- matrix(runif(2*n, 0, 10), n)
dim(loc)

# # Spatial parameters
# rho0 = 0.5
# kappa0 = 5
# 
dist_mat <- as.matrix(dist(loc))
# cov_mat <- matern_fun(dist_mat, phi = rho0, kappa = kappa0)
# u <- t(chol(cov_mat)) %*% rnorm(n, 0, 1) 
# 
# beta0 = -0.5
# beta1 <- 0.5
# sigma_e = 1
# x1 <- rnorm(n, 0 ,1)
# 
# 
# y_sim <- rnorm(n, mean = beta0 + beta1*x1 + u, sd = sigma_e)
# hist(y_sim)


sigma_e <- 1.0
sigma_u <- 1.0
kappa <- 5
nu <- 1

mcor <- matern_fun(dist_mat, nu, kappa)

u1 <- rnorm(n,0,sigma_u) %*% chol(mcor) #Spatial random effect
u2 <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = mcor)


plot(loc, cex = u1 - min(u1))
plot(loc, cex = u2 - min(u2))

#sim_grf <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = K)

beta0 <- -0.5
beta1 <-  0.5 
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


# Compile the model and load it
compile("tmbstan_spde.cpp")
dyn.load(dynlib("tmbstan_spde"))


#Define mesh and components representing the  precision matrix----
# boundary = INLA::inla.nonconvex.hull(loc)
# boundary2 = INLA::inla.nonconvex.hull(loc,convex = -0.35)
# mesh = INLA::inla.mesh.2d(loc=loc,
#                           boundary = list(boundary,boundary2),
#                           max.edge=c(0.05, 0.5),
#                           cutoff=0.55)
mesh <- inla.mesh.create(loc = as.matrix(loc), extend = T, refine = T)
mesh$n

plot(mesh)
points(loc, col = "red", pch = 19)
A = inla.spde.make.A(mesh,loc)
spde = inla.spde2.matern(mesh, alpha=2)
#spde = inla.spde2.pcmatern(mesh, alpha = 2, prior.range = c(5, 0.5), prior.sigma = c(5, 0.5))

spde_mat = spde$param.inla[c("M0","M1","M2")]


#=======================================
#           TMB modelling
#=======================================

## TMB data
tmb_data <-  list(model = 1, 
                  y    = y1,
                  X = model.matrix(~ 1 + x1),
                  site = mesh$idx$loc - 1,
                  spde_mat = spde_mat, 
                  A = A)#,
                  #tau0 = exp(spde$param.inla$theta.initial[1]),
                  #kappa0 = exp(spde$param.inla$theta.initial[2]))


## TMB parameters
tmb_par  <- list(beta = c(0.1, 0.1),
                 logsigma = 0.1,
                 logtau   = 0.1,
                 logkappa = 0.1,
                 u = rnorm(spde$mesh$n,0,1))



#dyn.unload(dynlib("spatial"))

## Create the TMB object
obj2 <- MakeADFun(data = tmb_data, parameters = tmb_par, random="u", DLL="tmbstan_spde", hessian = T)



#obj$env$beSilent()
opt2 <- with(obj2, nlminb(par, fn, gr))

opt2$par
exp(opt2$par[3:5])
rep  <- obj2$report()



#=====================================================
#                 Fit with tmbstan
#=====================================================

startTime <- Sys.time()
fit_tmbstan2 = tmbstan(obj2, chains = 3, open_progress = FALSE, 
              init ='last.par.best', 
              control = list(max_treedepth = 12, adapt_delta = 0.95), 
              iter = 4000, warmup=1000, cores = no_cores, seed=483892929)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)

s <- summary(fit_tmbstan2, probs = c(0.25, 0.75))
s$summary  # all chaines merged

exp(s$summary[3:5, c(1,3)])


library(gridExtra)
library(tidyverse)
traceplot(fit_tmbstan2, pars=names(obj2$par), inc_warmup=TRUE) + 
  theme(strip.text = ggplot2::element_text(size = 16, color = "black"),
        axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14))

#load package
require(MCMCvis)
MCMCsummary(fit_tmbstan2, round = 2)
pairs(fit_tmbstan2, pars = c("logsigma", "logtau", "logkappa")) # below the diagonal




# Simulating from the SIMULATE function
mat_sim = matrix(data=NA, nrow=length(obj2$simulate()$y_sim), ncol=100)
mat_sim


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = obj2$simulate()$y_sim
  }
}
mat_sim

# Plot lognormal non spatial simulations
hist(y, col = "gray90", prob = TRUE, main = "100 simulated samples", cex.main = 2, xlab = "", col.lab = 'blue', col.main="blue", cex.lab = 1.4, cex.axis = 1.2, ylim = c(0, 0.3))
for (j in 1: ncol(mat_sim)){
  lines(x = density(x = mat_sim[, j]),  lty="dotted", col="azure4", lwd=1)
}
#legend("topright", "A", bty = "n", text.col="blue", cex = 1.4)
legend(6, 0.3, c("response", "simulations"), lwd=4, col=c("gray90", "azure4"), cex = 1.5)




