rm(list = ls())
setwd("")

library(pacman)
pacman::p_load(TMB, rSPDE, fmesher, INLA, dplyr, tmbstan, rstan, parallel,
               raster, ggplot2, grid, gridExtra)

options(scipen=999)

# Calculate the number of cores
no_cores <- parallelly::availableCores()  - 1

#=================================================
compile("normal_spde.cpp")
dyn.load(dynlib("normal_spde"))
#=================================================


#=========================================
#               Simulation
#=========================================
set.seed(1234)

run_tmb <- function(nloc){
  set.seed(1234)
  
  nloc <- nloc
  x = runif(nloc, 0, 10)
  y = runif(nloc, 0, 10)
  coords = cbind(x, y)
  mesh <- fm_mesh_2d(loc = coords, cutoff = 0.5, max.edge = c(nloc*2.4, 5))
  mesh$n
  plot(mesh, main = "")
  
  sigma_u0 <- 1
  range0 <- 0.2
  nu <- 0.5
  kappa0 <- sqrt(8 * nu) / range0
  tau0 <- 1 / (sqrt(4*pi)*kappa0*sigma_u0)
  op <- matern.operators(mesh = mesh, alpha = 2, kappa = kappa0, tau = tau0, parameterization = "spde")
  
  u <- simulate(op)
  A <- spde.make.A(mesh = mesh, loc = coords)
  u = A %*% u
  
  
  beta0 = 1.0
  beta1 = 2.0
  sigma_e <- 0.1
  x1 = runif(length(coords[, 1]), 0, 1)
  y <- beta0 + beta1*x1 + u + rnorm(nloc) * sigma_e
  
  df = data.frame(y_sim = as.numeric(y), s1 = coords[, 1], s2 = coords[, 2], x1 = x1, u = as.numeric(u))
  
  spde = inla.spde2.matern(mesh, alpha = 2)
  spde_mat = spde$param.inla[c("M0","M1","M2")]
  
  sp_param = data.frame(sigma_u0, range0, kappa0, tau0)
  
  #======================================
  #                TMB data
  #======================================
  tmb_data <-  list(y    = as.vector(y),
                    x1   = as.vector(x1),
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
  obj <- MakeADFun(data = tmb_data, parameters = tmb_par, random="u", DLL="normal_spde", hessian = T)
  opt = nlminb(obj$par, obj$fn, obj$gr)
  sdrep = sdreport(obj, getJointPrecision = TRUE)
  res_list = list(obj, opt, sdrep, df, tmb_data, tmb_par, mesh, spde, sp_param)
  return(res_list)
}


obj1 <- run_tmb(nloc = 100)
obj2 <- run_tmb(nloc = 200)
obj3 <- run_tmb(nloc = 300)
obj4 <- run_tmb(nloc = 400)
obj5 <- run_tmb(nloc = 500)
obj6 <- run_tmb(nloc = 600)

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
M[[5]] = list()
M[[5]]$model = "spatial_n500"
M[[6]] = list()
M[[6]]$model = "spatial_n600"


M[[1]]$formula = obj1[[1]]
M[[2]]$formula = obj2[[1]]
M[[3]]$formula = obj3[[1]]
M[[4]]$formula = obj4[[1]]
M[[5]]$formula = obj5[[1]]
M[[6]]$formula = obj6[[1]]




#===========================================
#           Run the models
#===========================================

for (i in 1:length(M)){
  startTime <- Sys.time()
  print(paste("Running:  ", M[[i]]$model))
  fit <- tmbstan(M[[i]]$formula,
                 chains= 3, open_progress = FALSE,
                 control = list(max_treedepth= 13,  adapt_delta = 0.95),
                 iter = 4500, warmup= 700, cores=no_cores,
                 init = 'last.par.best', seed = 12345)
  endTime <- Sys.time()
  timeUsed = difftime(endTime, startTime, units='mins')
  print(timeUsed)
  saveRDS(fit, file=paste0('stan_spde_', i,'.RDS'))
}

warnings()



c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")


m1 <- readRDS('stan_spde_1.RDS')
m2 <- readRDS('stan_spde_2.RDS')
m3 <- readRDS('stan_spde_3.RDS')
m4 <- readRDS('stan_spde_4.RDS')
m5 <- readRDS('stan_spde_5.RDS')
m6 <- readRDS('stan_spde_6.RDS')

mon1 = monitor(m1)
mon2 = monitor(m2)
mon3 = monitor(m3)
mon4 = monitor(m4)
mon5 = monitor(m5)
mon6 = monitor(m6)


sum1 = data.frame(mon1$mean, mon1$`25%`, mon1$`75%`)
colnames(sum1) = c("mean", "25%", "75%")
head(sum1)

sum2 = data.frame(mon2$mean, mon2$`25%`, mon2$`75%`)
colnames(sum2) = c("mean", "25%", "75%")
head(sum2)

sum3 = data.frame(mon3$mean, mon3$`25%`, mon3$`75%`)
colnames(sum3) = c("mean", "25%", "75%")
head(sum3)

sum4 = data.frame(mon4$mean, mon4$`25%`, mon4$`75%`)
colnames(sum4) = c("mean", "25%", "75%")
head(sum4)

sum5 = data.frame(mon5$mean, mon5$`25%`, mon5$`75%`)
colnames(sum5) = c("mean", "25%", "75%")
head(sum5)

sum6 = data.frame(mon6$mean, mon6$`25%`, mon6$`75%`)
colnames(sum6) = c("mean", "25%", "75%")
head(sum6)

posterior_1 = as.matrix(M[[1]]$res)
posterior_2 = as.matrix(M[[2]]$res)
posterior_3 = as.matrix(M[[3]]$res)
posterior_4 = as.matrix(M[[4]]$res)
posterior_5 = as.matrix(M[[4]]$res)
posterior_6 = as.matrix(M[[6]]$res)

# Simulating from the SIMULATE function
mat_sim = matrix(data=NA, nrow=length(obj6[[1]]$simulate()$y_sim), ncol=100)
mat_sim

for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = obj6[[1]]$simulate()$y_sim
  }
}
mat_sim


df1 <- data.frame(obj1[[4]]$y_sim, mat_sim)
names(df1)[names(df1) == 'obj1..4...y_sim'] <- 'y_sim'

df2 <- data.frame(obj2[[4]]$y_sim, mat_sim)
names(df2)[names(df2) == 'obj2..4...y_sim'] <- 'y_sim'

df3 <- data.frame(obj3[[4]]$y_sim, mat_sim)
names(df3)[names(df3) == 'obj3..4...y_sim'] <- 'y_sim'

df4 <- data.frame(obj4[[4]]$y_sim, mat_sim)
names(df4)[names(df4) == 'obj4..4...y_sim'] <- 'y_sim'

df5 <- data.frame(obj5[[4]]$y_sim, mat_sim)
names(df5)[names(df5) == 'obj5..4...y_sim'] <- 'y_sim'

df6 <- data.frame(obj6[[4]]$y_sim, mat_sim)
names(df6)[names(df6) == 'obj6..4...y_sim'] <- 'y_sim'



# Histogram with kernel density
p1 <- ggplot(df1, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, label = "SP1", colour = "blue", size = 10)

for (i in 2:ncol(df1)) {
  p1 <- p1 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df1[, i]), 
                                       sd = sd(df1[, i])), lwd = 1, col = 'orange')}




# Histogram with kernel density
p2 <- ggplot(df2, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("SP2") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, label = "SP2", colour = "blue", size = 10)

for (i in 2:ncol(df2)) {
  p2 <- p2 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df2[, i]), 
                                       sd = sd(df2[, i])), lwd = 1, col = 'orange')}



# Histogram with kernel density
p3 <- ggplot(df3, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        #axis.title.y = element_blank(),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("SP3") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, label = "SP3", colour = "blue", size = 10)

for (i in 2:ncol(df3)) {
  p3 <- p3 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df3[, i]), 
                                       sd = sd(df3[, i])), lwd = 1, col = 'orange')}



# Histogram with kernel density
p4 <- ggplot(df4, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("SP4") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, label = "SP4", colour = "blue", size = 10)

for (i in 2:ncol(df4)) {
  p4 <- p4 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df4[, i]), 
                                       sd = sd(df4[, i])), lwd = 1, col = 'orange')}


# Histogram with kernel density
p5 <- ggplot(df5, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("SP5") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, label = "SP5", colour = "blue", size = 10)

for (i in 2:ncol(df5)) {
  p5 <- p5 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df5[, i]), 
                                       sd = sd(df5[, i])), lwd = 1, col = 'orange')}


# Histogram with kernel density
p6 <- ggplot(df6, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("SP6") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, label = "SP6", colour = "blue", size = 10)

for (i in 2:ncol(df6)) {
  p6 <- p6 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df6[, i]), 
                                       sd = sd(df6[, i])), lwd = 1, col = 'orange')}

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3,
             top = textGrob("M-GRF", gp=gpar(fontsize=28,font=1)))

