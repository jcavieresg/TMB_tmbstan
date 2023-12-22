# if (!require(devtools)) {
#   install.packages("devtools")
# }
# devtools::install_github("stan-dev/loo", build_vignettes = FALSE)
# 

rm(list = ls())
setwd("C:/Users/Usuario/Desktop/semipar_TPS/tmbstan_real_model")

library(loo)
library(dplyr)
library(ggplot2)
library(rstan)
set.seed(12345)

obj_tps  = readRDS("C:/Users/Usuario/Desktop/semipar_TPS/tmbstan_real_model/fit_tps_TMB.RDS")
obj_TMB  = readRDS("C:/Users/Usuario/Desktop/Spatial/fits_TMB.RDS")
obj_spde = obj_TMB[[4]]

posterior_tps  = readRDS("C:/Users/Usuario/Desktop/semipar_TPS/tmbstan_real_model/posterior_tps_tmbstan.RDS")
posteriors  = readRDS("C:/Users/Usuario/Desktop/Spatial/posteriors_tmbstan.RDS")
posterior_spde = posteriors[[4]]

#======================================================================================================
#                                   load the model created in TMB
#                                 dyn.load(dynlib("spatial_tmbstan"))
#======================================================================================================


data = read.csv("north2.csv", header = T)
head(data, 3)

#========================
# loo for the tps model
# dyn.load(dynlib("semipar_SNTPS_tmbstan"))
set.seed(1000)
n_1 = length(data$cpue)
log_lik_1 <- matrix(NA, nrow=nrow(posterior_tps), ncol=n_1)

for(i in 1:nrow(posterior_tps)){
  r1 <- obj_tps$report(posterior_tps[i,-ncol(posterior_tps)])
  log_lik_1[i,] <- r1$log_lik
  }

loo_1 = loo(log_lik_1, r_eff = NA)


#=========================
# loo for the spde model
# dyn.load(dynlib("spatial_tmbstan")) # only if it is neccesary!
n_2 = length(data$cpue)
log_lik_2 <- matrix(NA, nrow=nrow(posterior_spde), ncol=n_2)

for(i in 1:nrow(posterior_spde)){
  r2 <- obj_spde$Obj$report(posterior_spde[i,-ncol(posterior_spde)])
  log_lik_2[i,] <- r2$log_lik
}

loo_2 = loo(log_lik_2, r_eff = NA)
#print(loo_2)






#================================================================
#               loo comparision for the 6 models
#================================================================
comp <- loo_compare(loo_1, loo_2)
print(comp, digits = 3)

par(mfrow=c(1,2))
plot(loo_2, main =  expression(bold('Gamma model')), cex.lab = 1.5)
plot(loo_1, main =  expression(bold('Skew normal model')), xlab="")


khat <- pareto_k_values(loo_1)
kdata <- data.frame(Data_point = seq_along(khat), khat = khat, bad = khat > 0.7)
plot1 = ggplot(kdata, aes(x = Data_point, y = khat, label = Data_point)) + 
  #geom_hline(yintercept = 0.7, size = 1.0, linetype = 2, color = "blue") +
  geom_point(shape = 1) +
  geom_hline(yintercept = 0.7, size = 1.0, linetype = 2, color = "blue") + 
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  ggtitle("Skew normal ST model") + 
  theme(plot.title = element_text(size = 22, hjust = 0.5))


khat <- pareto_k_values(loo_2)
kdata <- data.frame(Data_point = seq_along(khat), khat = khat)
plot2 = ggplot(kdata, aes(x = Data_point, y = khat)) + 
  #geom_point(data = subset(kdata, !bad), size = 1.0, shape = 1) + ylim(-0.25, 0.9) +
  geom_point(shape = 1) +
  geom_hline(yintercept = 0.7, size = 1.0, linetype = 2, color = "blue") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  ggtitle("Gamma ST model") + 
  theme(plot.title = element_text(size = 22, hjust = 0.5))

library(gridExtra)
grid.arrange(plot2, plot1, ncol = 2)




mcmc_areas(
  posterior_tps, 
  pars = c("beta0", "sigma", "alpha", "sigma"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)


library(loo)
library(rstan)
mon1 = monitor(posterior_tps)

data.frame(posterior_tps[, 5])

ggplot(pl, aes(Time)) + 
  geom_line(aes(y=menle), colour="blue") + 
  geom_ribbon(aes(ymin=menlelb, ymax=menleub), alpha=0.2)





#======================================================================
#      Additionaly comparitions for Gamma model (sensivity analysis)
#======================================================================

posteriors  = readRDS("C:/Users/Usuario/Desktop/Spatial/posteriors_tmbstan.RDS")

# posterior_tps = posteriors[[4]]
# posterior_spde = posteriors[[2]]
# posterior_3 = posteriors[[3]]
posterior_4 = posteriors[[4]]



# loo for log(tau) = model base (Spatial Gamma model)----> -3.78
set.seed(1000)
n_1 = length(data$cpue)
log_lik_1 <- matrix(NA, nrow=nrow(posterior_4), ncol=n_1)

for(i in 1:nrow(posterior_4)){
  r1 <- obj4$Obj$report(posterior_4[i,-ncol(posterior_4)])
  log_lik_1[i,] <- r1$log_lik
}

loo_1 = loo(log_lik_1, r_eff = NA)




# loo for log(tau) = model base (Spatial Gamma model)----> -0.975
set.seed(1000)
n_2 = length(data$cpue)
log_lik_2 <- matrix(NA, nrow=nrow(posterior_4), ncol=n_2)

for(i in 1:nrow(posterior_4)){
  r2 <- obj4$Obj$report(posterior_4[i,-ncol(posterior_4)])
  log_lik_2[i,] <- r2$log_lik
}

loo_2 = loo(log_lik_2, r_eff = NA)
print(loo_2)


# loo for log(tau) = model base (Spatial Gamma model)----> -15.12
set.seed(1000)
n_3 = length(data$cpue)
log_lik_3 <- matrix(NA, nrow=nrow(posterior_4), ncol=n_3)

for(i in 1:nrow(posterior_4)){
  r3 <- obj4$Obj$report(posterior_4[i,-ncol(posterior_4)])
  log_lik_3[i,] <- r3$log_lik
}

loo_3 = loo(log_lik_3, r_eff = NA)
print(loo_3)




#================================================================
#               loo comparision for the 6 models
#================================================================
comp <- loo_compare(loo_1, loo_2, loo_3)
print(comp, digits = 3)


