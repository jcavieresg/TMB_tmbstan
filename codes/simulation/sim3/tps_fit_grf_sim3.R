rm(list = ls())
setwd("")

library(pacman)
pacman::p_load(geoR, fields, prodlim, TMB, TMBhelper, mgcv, dplyr, tmbstan, parallel, MASS, Matrix,
               raster, ggplot2, gridExtra, bayesplot, grid, fmesher, rSPDE, GeoModels, INLA)

options(scipen=999)

# Calculate the number of cores
no_cores <- parallelly::availableCores()  - 1




#=================================================
compile("tps_fit_grf_sim3.cpp")
dyn.load(dynlib("tps_fit_grf_sim3"))
#=================================================


#=========================================
#               Simulation
#=========================================
set.seed(1234)

run_tmb <- function(nloc, model, grid_model){

set.seed(1234)
  
  # Simulating a GRF
  nloc = nloc
  if(grid_model == 1) {
    sim_grf <- grf(nloc, cov.pars = c(0.5, .15))
  } else if (grid_model == 2){
    sim_grf <- grf(nloc, grid = "reg", cov.pars = c(0.5, .15))
    image(sim_grf)
  } else {
    print("Set a correct grid_model")
  }


  # Hyperparameters
  hyper_params <- c(sim_grf$cov.pars, sim_grf$kappa)  # sigmau, rho and kappa respectivelly
  u <- sim_grf$data
  coords <- sim_grf$coords


  
  #===========================================================================
  #                          Simulate data
  #===========================================================================
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
  
  # Create artificial variables to obtain the penalty smoother
  df$z = 0
  df$x0 = 0
  
  
  #====================================================
  #      Getting matrices from mgcv
  #====================================================
  tp_setup = gam(z ~ x0 + s(s1, s2, bs = "tp", 
                 k = ifelse(round(length(df$s1)*0.3, digits = 0) < 30, 30, 
                            round(length(df$s1)*0.3, digits = 0))),
                 data = df,
                 fit = FALSE)
  
  Stp = tp_setup$smooth[[1]]$S[[1]]
  Xtp = tp_setup$X[, c(-1, -2)] # Igualmente debo remover el intercepto
  E = ginv(Stp)
  log_lam = 0
  
  x = rmvn(n = 1, mu = rep(0, nrow(E)), V = exp(log_lam) * E) 
  knots_app <- tp_setup$smooth[[1]]$bs.dim
  knots_app
  
  #======================================
  #                TMB data
  #======================================
  tmb_data <- list(model = model, 
                   y = df$y_sim,
                   x1 = df$x1, 
                   X = Xtp,
                   S = as(Stp, "sparseMatrix"))
  
  #=====================================
  #            TMB parameters
  #=====================================
  tmb_par <- list(beta0 = 0.5,
                  beta1 = 0.5,
                  logsigma_e = 0.5,
                  x = rep(0, nrow(Stp)),
                  loglambda = 0.5)

  obj <- MakeADFun(tmb_data, random = c("x"), tmb_par, DLL="tps_fit_grf_sim3", hessian = T)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdrep = sdreport(obj, getJointPrecision = TRUE)
  
  res_list <- list(obj, opt, sdrep, df, tmb_data, tmb_par, knots_app)
  return(res_list)
}


obj1 <- run_tmb(nloc = 100, model = 1, grid_model = 1)
obj2 <- run_tmb(nloc = 200, model = 1, grid_model = 1)
obj3 <- run_tmb(nloc = 300, model = 1, grid_model = 1)
obj4 <- run_tmb(nloc = 400, model = 1, grid_model = 1)


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
  saveRDS(fit, file=paste0('tps_fit_grf_sim3_', i,'.RDS'))
}
