rm(list=ls())

setwd("C:/Users/joaquin/Desktop/new_models")

library(pacman)
pacman::p_load(RandomFields, geoR, fields, prodlim, TMB, INLA, dplyr, parallel, rstan)

options(scipen = 999)

# Calculate the number of cores
no_cores <- parallelly::availableCores() - 1  
rstan_options(auto_write = TRUE)


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
nu <- 1

mcor <- matern_fun(dist_mat, nu, kappa)

u1 <- rnorm(n,0, sigma_u) %*% chol(mcor)             #Spatial random effect
u2 <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = mcor)


plot(loc, cex = u1 - min(u1))
plot(loc, cex = u2 - min(u2))

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



# Mesh of INLA
mesh <- inla.mesh.create(loc = as.matrix(loc), extend = T, refine = T)
mesh$n

plot(mesh)
points(loc, col = "red", pch = 19)
A = inla.spde.make.A(mesh,loc)
A <- as(A, "matrix")
spde = inla.spde2.matern(mesh, alpha=2)


# spde_mat <- array(c(as(spde$param.inla$M0, "matrix"),
#              as(spde$param.inla$M1, "matrix"),
#              as(spde$param.inla$M2, "matrix")),
#              dim = c(mesh$n, mesh$n, 3))
# spde_mat <- aperm(spde_mat, c(3, 1, 2))


spdeMatrices = spde$param.inla[c("M0","M1","M2")]


# Stan data
stan_data <- list(n = n,
                  y = as.numeric(y1),
                  k = 2,                           # Number parameters in the design matrix (including the intercept!)
                  X = model.matrix(~ 1 + x1),      # Design matrix
                  n_knots = mesh$n,
                  M0 = as.matrix(spdeMatrices$M0),
                  M1 = as.matrix(spdeMatrices$M1),
                  M2 = as.matrix(spdeMatrices$M2),
                  #spde_mat = spde_mat,
                  A = A)

startTime <- Sys.time()
fit_stan <- stan(file='stan_chol_matern2.stan', data = stan_data, 
                 iter = 2000, chains = 3, warmup = 1000, 
                 control = list(max_treedepth = 12, adapt_delta = 0.9), 
                 cores = no_cores, open_progress=FALSE, seed=483892929)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)


s <- summary(fit_stan, probs = c(0.25, 0.75))
s$summary  # all chaines merged


# Extract posterior draws for later use
posterior_cp <- as.array(fit_stan)


require(MCMCvis)
library(bayesplot)
mcmc_pairs(posterior_cp, pars = c("sigma_e","kappa","tau"),
           off_diag_args = list(size = 0.75))

pairs(fit_stan, pars = c("beta", "sigma_e")) # below the diagonal
pairs(fit_stan, pars = c("tau", "kappa")) # below the diagonal
