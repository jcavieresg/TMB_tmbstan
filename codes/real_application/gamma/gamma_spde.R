rm(list = ls())
setwd("C:/Users/Usuario/Desktop/tps_vs_spde/real_application/codes/gamma")

library(pacman)
pacman::p_load(geoR, fields, prodlim, TMB, INLA, dplyr, tmbstan, rstan, parallel)

options(scipen=999)

# Calculate the number of cores
no_cores <- parallelly::availableCores()  - 1


# Compile the model and load it
compile("gamma_spde.cpp")
dyn.load(dynlib("gamma_spde"))






#=============================================================================
#                              Applied analysis
#=============================================================================
data = read.csv("north2.csv", header = T)
head(data, 3)

loc <- cbind(data$longitude, data$latitude)

#Define mesh and components representing the  precision matrix----
bound1 <- inla.nonconvex.hull(loc)
mesh = inla.mesh.create(loc, plot.delay=NULL, refine=TRUE, boundary = bound1)

plot(mesh)
points(loc, col = "red", pch = 19)
A = inla.spde.make.A(mesh,loc)
spde = inla.spde2.matern(mesh, alpha=2)

spde_mat = spde$param.inla[c("M0","M1","M2")]


#=================================================================
#                           TMB modelling
#=================================================================
data$year <- as.factor(data$year)
data$trim <- as.factor(data$trim)
data$destine <- as.factor(data$destine)


#======================================
#                TMB data
#======================================
tmb_data = list(likelihood = 2,
                cpue     = as.vector(data$cpue),          # Response
                depth = as.numeric(data$depth),
                year = as.numeric(data$year) - 1,
                trim = as.numeric(data$trim) - 1,
                destine = as.numeric(data$destine) - 1,
                spde_mat = spde_mat, 
                A = A,
                tau0 = exp(spde$param.inla$theta.initial[1]),
                kappa0 = exp(spde$param.inla$theta.initial[2]))
                  


## TMB parameters
tmb_par  <- list(beta0 = 0.1,
                 beta_depth = 0.1,
                 beta_year = rep(0, length(levels(data$year))),
                 beta_trim = rep(0, length(levels(data$trim))),
                 beta_destine = rep(0, length(levels(data$destine))),
                 logsigma = 0.1,
                 logtau   = 0.1,
                 logkappa = 0.1,
                 u = rnorm(spde$mesh$n,0,1))



#dyn.unload(dynlib("spatial"))

## Create the TMB object
obj <- MakeADFun(data = tmb_data, parameters = tmb_par, random="u", DLL="gamma_spde", hessian = T)


startTime <- Sys.time()
opt <- nlminb(obj$par,obj$fn,obj$gr)
rep <- sdreport(obj)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)


saveRDS(obj, file='TMB_gamma_spde.RDS')


#=====================================================
#                 Fit with tmbstan
#=====================================================

init.fn <- function()
  split(unname(obj$env$last.par.best),names(obj$env$last.par.best))


startTime <- Sys.time()
stan_gamma_spde = tmbstan(obj, chains = 3, open_progress = FALSE, 
              init = init.fn, 
              control = list(max_treedepth = 13, adapt_delta = 0.95), 
              iter = 3500, warmup=700, cores = no_cores, seed=483892929)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)



posterior_gamma_spde <- as.matrix(stan_gamma_spde)

saveRDS(posterior_gamma_spde, file='posterior_gamma_spde.RDS')
saveRDS(stan_gamma_spde, file='stan_gamma_spde.RDS')

