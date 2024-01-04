#================================================================
#                        Run models in tmbstan
#================================================================


M = list()
M[[1]] = list()
M[[1]]$model = "lognormal_non-spatial"
M[[2]] = list()
M[[2]]$model = "lognormal_spatial"
M[[3]] = list()
M[[3]]$model = "gamma_non-spatial"
M[[4]] = list()
M[[4]]$model = "gamma_spatial"



M[[1]]$formula = obj1$Obj
M[[2]]$formula = obj2$Obj
M[[3]]$formula = obj3$Obj
M[[4]]$formula = obj4$Obj


## Write each file separately so you can run them in different R
## sessions
if(answer=='YES'){
  for (i in 1:4){
    print(paste("Running:  ", M[[i]]$model))
    fit <- tmbstan(M[[i]]$formula,
                   chains= 4, open_progress = FALSE,
                   control = list(max_treedepth= 13,  adapt_delta = 0.9),
                   iter = 3000, warmup= 700, cores=no_cores,
                   init = 'par')
    saveRDS(fit, file=paste0('tmbstan_fit_', i,'.RDS'))
  }
}

message("Reading in tmbstan fits from file")
M[[1]]$res <- readRDS('tmbstan_fit_1.RDS')
M[[2]]$res <- readRDS('tmbstan_fit_2.RDS')
M[[3]]$res <- readRDS('tmbstan_fit_3.RDS')
M[[4]]$res <- readRDS('tmbstan_fit_4.RDS')

message("Saving posteriors and model fits")
posterior_1 = as.matrix(M[[1]]$res)
posterior_2 = as.matrix(M[[2]]$res)
posterior_3 = as.matrix(M[[3]]$res)
posterior_4 = as.matrix(M[[4]]$res)
posteriors_tmbstan <- list(posterior_1, posterior_2, posterior_3, posterior_4)
saveRDS(posteriors_tmbstan, file='posteriors_tmbstan.RDS')
saveRDS(M, file='tmbstan_models.RDS')
