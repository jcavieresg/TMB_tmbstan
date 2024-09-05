rm(list = ls())
setwd("C:/Users/Usuario/Desktop/tps_vs_spde/answers_reviewers/files_and_codes/sim3")

library(pacman)
pacman::p_load(geoR, fields, prodlim, TMB, TMBhelper, mgcv, dplyr, tmbstan, parallel, MASS, Matrix,
               raster, ggplot2, gridExtra, bayesplot, grid, rSPDE, fmesher, INLA, GeoModels)

options(scipen=999)
# Calculate the number of cores
no_cores <- detectCores() - 1

# #=================================================
# # Compilamos el modelo y lo cargamos en R
TMB::compile('spde_fit_grf_sim3.cpp')
dyn.load(dynlib("spde_fit_grf_sim3"))
# # #=================================================







#==================================================
#                 SIMULATED DATA
#==================================================

run_tmb <- function(nloc, model, grid_model){
set.seed(1234)
  
  # Simulating a GRF
  nloc = nloc
  if(grid_model == 1) {
    sim_grf <- grf(nloc, cov.pars = c(0.5, .15), cov.model = 'matern')
  } else if (grid_model == 2){
    sim_grf <- grf(nloc, grid = "reg", cov.pars = c(0.5, .15), cov.model = 'matern')
    image(sim_grf)
  } else {
    print("Set a correct grid_model")
  }


  # Hyperparameters
  u <- sim_grf$data
  tau = 1/var(u)
  hyper_params <- c(sim_grf$cov.pars, sim_grf$kappa, tau)  # sigmau, rho and kappa respectivelly
  coords <- sim_grf$coords
  

  # Initial values
  beta0 <- 1.0
  beta1 <- 2.0
  sigma_e <- 0.1


if(grid_model == 1) {
  x1 = runif(nloc, 0, 1)
  mu_sim = beta0 + beta1*x1 + u
  y_sim = mu_sim + rnorm(length(u))*sigma_e
} else if (grid_model == 2){
  x1 = runif(length(u), 0, 1)
  mu_sim = beta0 + beta1*x1 + u
  y_sim = mu_sim + rnorm(length(u))*sigma_e
} else {
  print("Set a correct length number for x (must be equal to u)")
}
  

df = data.frame(y_sim = as.numeric(y_sim), s1 = coords[, 1], s2 = coords[, 2], x1 = x1, u = as.numeric(u))

# SPDE quantities
mesh <- fm_mesh_2d(loc = coords, cutoff = 0.06, max.edge = c(nloc*2.4, 5))
mesh$n
plot(mesh, main = "")
points(coords, col = "red")

A <- spde.make.A(mesh = mesh, loc = coords)

spde = inla.spde2.matern(mesh, alpha = 2)
spde_mat = spde$param.inla[c("M0","M1","M2")]


#======================================
#                TMB data
#======================================
tmb_data <-  list(model = model, 
                  y    = as.vector(df$y_sim),
                  x1   = as.vector(df$x1),
                  spde_mat = spde_mat,
                  A = A)

#=====================================
#            TMB parameters
#=====================================
tmb_par  <- list(beta0 = 0.1,
                 beta1 = 0.1,
                 logsigma_e = 0.1,
                 logtau   = 0.1,
                 logkappa = 0.1,
                 u = rnorm(spde$mesh$n,0,1))

## Create the TMB objects
obj <- MakeADFun(data = tmb_data, parameters = tmb_par, random="u", DLL="spde_fit_grf_sim3", hessian = T)
opt = nlminb(obj$par, obj$fn, obj$gr)
sdrep = sdreport(obj, getJointPrecision = TRUE)
res_list = list(obj, opt, sdrep, df, tmb_data, tmb_par, mesh, spde, hyper_params)
return(res_list)
}

#=====================================
#             Run the TMB models
#=====================================
obj1 <- run_tmb(nloc = 100, model = 1, grid_model = 1)
obj2 <- run_tmb(nloc = 200, model = 1, grid_model = 1)
obj3 <- run_tmb(nloc = 300, model = 1, grid_model = 1)
obj4 <- run_tmb(nloc = 400, model = 1, grid_model = 1)


#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================
M = list()
M[[1]] = list()
M[[1]]$model = "spatial_100"
M[[2]] = list()
M[[2]]$model = "spatial_n200"
M[[3]] = list()
M[[3]]$model = "spatial_n300"
M[[4]] = list()
M[[4]]$model = "spatial_n400"


M[[1]]$formula = obj1[[1]]
M[[2]]$formula = obj2[[1]]
M[[3]]$formula = obj3[[1]]
M[[4]]$formula = obj4[[1]]

#===========================================
#           Run the models
#===========================================
for (i in 1:length(M)){
  startTime <- Sys.time()
  print(paste("Running:  ", M[[i]]$model))
  fit <- tmbstan(M[[i]]$formula,
                 chains= 3, open_progress = FALSE,
                 control = list(max_treedepth= 13,  adapt_delta = 0.95),
                 iter = 4000, warmup= 700, cores=no_cores,
                 init = 'last.par.best', seed = 12345)
  endTime <- Sys.time()
  timeUsed = difftime(endTime, startTime, units='mins')
  print(timeUsed)
  saveRDS(fit, file=paste0('spde_fit_grf_sim3_', i,'.RDS'))
}
