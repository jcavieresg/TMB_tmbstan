rm(list = ls())
setwd("")

library(pacman)
pacman::p_load(TMB, INLA, dplyr, tmbstan, rstan, parallel, bayesplot, ggplot2, loo, gridExtra, bayestestR)

options(scipen=999)

# Calculate the number of cores
no_cores <- parallelly::availableCores()  - 1


# TMB objects
obj_TMB  = readRDS("C:/Users/Usuario/Desktop/Spatial/fits_TMB.RDS")
M1_tmb = obj_TMB[[4]]
M2_tmb  = readRDS("C:/Users/Usuario/Desktop/tps_vs_spde/real_application/TMB_sn_spde.RDS")

M3_tmb  = readRDS("C:/Users/Usuario/Desktop/tps_vs_spde/real_application/TMB_gamma_tps2.RDS")
M4_tmb  = readRDS("C:/Users/Usuario/Desktop/tps_vs_spde/real_application/TMB_sn_tps.RDS")


# Posteriors (from tmbstan)
posteriors  = readRDS("C:/Users/Usuario/Desktop/Spatial/posteriors_tmbstan.RDS")
M1 = posteriors[[4]]
M2  = readRDS("C:/Users/Usuario/Desktop/tps_vs_spde/real_application/M2/posterior_sn_spde.RDS")

M3  = readRDS("C:/Users/Usuario/Desktop/tps_vs_spde/real_application/M3/posterior_gamma_tps.RDS")
M4  = readRDS("C:/Users/Usuario/Desktop/tps_vs_spde/real_application/M4/posterior_sn_tps.RDS")


#===========================
# Read the data
#===========================
data = read.csv("north2.csv", header = T)
head(data, 3)


#==============================================================
# M1 (here the executable in TMB is called "spatial_tmbstan")
dyn.load(dynlib("C:/Users/Usuario/Desktop/Spatial/spatial_tmbstan")) # only if it is neccesary!
n_1 = length(data$cpue)
log_lik_1 <- matrix(NA, nrow=nrow(M1), ncol=n_1)

for(i in 1:nrow(M1)){
  r1 <- M1_tmb$Obj$report(M1[i,-ncol(M1)])
  log_lik_1[i,] <- r1$log_lik
}

loo_1 = loo(log_lik_1, cores = no_cores)
print(loo_1)


#====================================================================
# M2 (here the executable in TMB is called "semipar_SNSPDE_tmbstan")
dyn.load(dynlib("semipar_SNSPDE_tmbstan")) # only if it is neccesary!
n_2 = length(data$cpue)
log_lik_2 <- matrix(NA, nrow=nrow(M2), ncol=n_2)

for(i in 1:nrow(M2)){
  r2 <- M2_tmb$report(M2[i,-ncol(M2)])
  log_lik_2[i,] <- r2$log_lik
}

loo_2 = loo(log_lik_2, cores = no_cores)
print(loo_2)



#====================================================================
# M3 (here the executable in TMB is called "gamma_tps2")
dyn.load(dynlib("gamma_tps2"))
n_3 = length(data$cpue)
log_lik_3 <- matrix(NA, nrow=nrow(M3), ncol=n_3)

for(i in 1:nrow(M3)){
  r3 <- M3_tmb$report(M3[i,-ncol(M3)])
  log_lik_3[i,] <- r3$log_lik
}

loo_3 = loo(log_lik_3, cores = no_cores)
print(loo_3)



#====================================================================
# M4 (here the executable in TMB is called "semipar_SNTPS_tmbstan")
dyn.load(dynlib("semipar_SNTPS_tmbstan"))
n_4 = length(data$cpue)
log_lik_4 <- matrix(NA, nrow=nrow(M4), ncol=n_4)

for(i in 1:nrow(M4)){
  r4 <- M4_tmb$report(M4[i,-ncol(M4)])
  log_lik_4[i,] <- r4$log_lik
}

loo_4 = loo(log_lik_4, cores = no_cores)
print(loo_4)


#================================================================
#               loo comparision for the 4 models
#================================================================
comp <- loo_compare(loo_1, loo_2, loo_3, loo_4)
print(comp, digits = 2, simplify = FALSE)

par(mfrow=c(2,2))
plot(loo_1, main =  expression(bold('Gamma model')), cex.lab = 1.5)
plot(loo_2, main =  expression(bold('Skew normal model')), xlab="")
plot(loo_3, main =  expression(bold('Gamma model')), cex.lab = 1.5)
plot(loo_4, main =  expression(bold('Skew normal model')), xlab="")


#====================================================================
# JACOBIAN correction since the Skew normal model used sqrt(cpue)
#====================================================================

loo_2_with_jacobian <- loo_2
loo_2_with_jacobian$pointwise[,1] <- loo_2_with_jacobian$pointwise[,1] - log(2*sqrt(data$cpue))
loo_2_with_jacobian$estimates["elpd_loo",] <- loo:::table_of_estimates(loo_2_with_jacobian$pointwise[,"elpd_loo", drop=FALSE])

loo_4_with_jacobian <- loo_4
loo_4_with_jacobian$pointwise[,1] <- loo_4_with_jacobian$pointwise[,1] - log(2*sqrt(data$cpue))
loo_4_with_jacobian$estimates["elpd_loo",] <- loo:::table_of_estimates(loo_4_with_jacobian$pointwise[,"elpd_loo", drop=FALSE])

comp2 <- loo_compare(loo_1, loo_2_with_jacobian, loo_3, loo_4_with_jacobian)
print(comp2, digits = 2, simplify = FALSE)
