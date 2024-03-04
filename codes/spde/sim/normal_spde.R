rm(list = ls())
setwd("")

library(pacman)
pacman::p_load(geoR, fields, prodlim, TMB, INLA, dplyr, tmbstan, rstan, parallel,
               raster, ggplot2)

options(scipen=999)

# Calculate the number of cores
no_cores <- parallelly::availableCores()  - 1




#=================================================
# Compilamos el modelo y lo cargamos en R
compile("normal_spde.cpp")
dyn.load(dynlib("normal_spde"))


#=========================================
#               Simulation
#=========================================
set.seed(1234)

run_tmb <- function(nloc){
  
  inla.seed = sample.int(n=1E6, size=1)
  options(width=70, digits=3)
  
  loc_ini = matrix(c(0,0,10,10, 0, 10, 10, 0), nrow = 4, byrow = T)
  mesh_sim = inla.mesh.2d(loc = loc_ini, max.edge=c(1.9, 2.9))
  spde = inla.spde2.matern(mesh_sim, alpha = 2)

  range = 2
  kappa = sqrt(8)/ range
  tau = 0.5
  
  Qgrf = inla.spde.precision(spde, theta=c(log(kappa), log(tau)))
  grf = inla.qsample(n=1, Q=Qgrf)
  grf = grf[ ,1]
  
#=============================================
#            Number of locations
#=============================================
x <- runif(nloc, 0, 10)
y <- runif(nloc, 0, 10)
  
coords = cbind(x, y)
  
A = inla.spde.make.A(mesh = mesh_sim, loc = coords)
grf = drop(A %*% grf)
  
#===========================================================================
#                          Simulate data
#===========================================================================
x1 = runif(length(coords[, 1]), 0, 1)
beta0 = 1.0
beta1 = 2.0
sigma_ini = 0.3
  
mu_sim = beta0 + beta1*x1 + grf
y_sim = mu_sim + rnorm(length(coords[, 1]), mean = 0, sd = sigma_ini)
  
# data.frame
df = data.frame(y_sim = y_sim, s1= coords[ ,1], s2 = coords[ ,2], x1 = x1, grf = grf)
summary(df)
  
mesh = inla.mesh.2d(loc = loc_ini, max.edge=c(1.9, 2.9))
A = inla.spde.make.A(mesh=mesh, loc=data.matrix(df[ , c('s1', 's2')]))
  
spde = inla.spde2.matern(mesh, alpha = 2)
spde_mat = spde$param.inla[c("M0","M1","M2")]
  
#=================================================================
#                           TMB modelling
#=================================================================
  
#======================================
#                TMB data
#======================================
tmb_data <-  list(y    = as.vector(df$y_sim),
                  x1 = as.vector(df$x1),
                  site = mesh_sim$idx$loc - 1,
                  spde_mat = spde_mat,
                  A = A)#,
 
#=====================================
#            TMB parameters
#=====================================
tmb_par  <- list(beta0 = 0.1, 
                 beta1 = 0.1,
                 logsigma_e = 0.1,
                 logtau   = 0.1,
                 logkappa = 0.1,
                 u = rnorm(spde$mesh$n,0,1))
  
  
  
#dyn.unload(dynlib("spatial"))
  
## Create the TMB object
obj <- MakeADFun(data = tmb_data, parameters = tmb_par, random="u", DLL="normal_spde3", hessian = T)
opt = nlminb(obj$par, obj$fn, obj$gr)
sdrep = sdreport(obj, getJointPrecision = TRUE)
res_list = list(obj, opt, sdrep, df, tmb_data, tmb_par, mesh_sim)
return(res_list)
}


obj1 <- run_tmb(nloc = 50)
obj2 <- run_tmb(nloc = 100)
obj3 <- run_tmb(nloc = 150)
obj4 <- run_tmb(nloc = 200)
obj5 <- run_tmb(nloc = 250)
obj6 <- run_tmb(nloc = 300)

#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================


M = list()
M[[1]] = list()
M[[1]]$model = "spatial_n50"
M[[2]] = list()
M[[2]]$model = "spatial_n75"
M[[3]] = list()
M[[3]]$model = "spatial_n100"
M[[4]] = list()
M[[4]]$model = "spatial_n125"
M[[5]] = list()
M[[5]]$model = "spatial_n150"
M[[6]] = list()
M[[6]]$model = "spatial_n175"



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
                 iter = 4000, warmup= 700, cores=no_cores,
                 #lower = lwr, upper = upr, 
                 init = 'last.par.best', seed = 12345)
  endTime <- Sys.time()
  timeUsed = difftime(endTime, startTime, units='mins')
  print(timeUsed)
  saveRDS(fit, file=paste0('stan_spde_', i,'.RDS'))
}

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
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP 1", colour = "blue", size = 10)

for (i in 2:ncol(df1)) {
  p1 <- p1 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df1[, i]), 
                                       sd = sd(df1[, i])), lwd = 1, col = 'orange')}




# Histogram with kernel density
p2 <- ggplot(df2, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 2") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP 2", colour = "blue", size = 10)

for (i in 2:ncol(df2)) {
  p2 <- p2 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df2[, i]), 
                                       sd = sd(df2[, i])), lwd = 1, col = 'orange')}



# Histogram with kernel density
p3 <- ggplot(df3, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        #axis.title.y = element_blank(),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 3") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP 3", colour = "blue", size = 10)

for (i in 2:ncol(df3)) {
  p3 <- p3 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df3[, i]), 
                                       sd = sd(df3[, i])), lwd = 1, col = 'orange')}



# Histogram with kernel density
p4 <- ggplot(df4, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 4") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP 4", colour = "blue", size = 10)

for (i in 2:ncol(df4)) {
  p4 <- p4 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df4[, i]), 
                                       sd = sd(df4[, i])), lwd = 1, col = 'orange')}


# Histogram with kernel density
p5 <- ggplot(df5, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 5") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP 5", colour = "blue", size = 10)

for (i in 2:ncol(df5)) {
  p5 <- p5 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df5[, i]), 
                                       sd = sd(df5[, i])), lwd = 1, col = 'orange')}


# Histogram with kernel density
p6 <- ggplot(df6, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.35, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 6") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Grid 6", colour = "blue", size = 10)

for (i in 2:ncol(df6)) {
  p6 <- p6 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df6[, i]), 
                                       sd = sd(df6[, i])), lwd = 1, col = 'orange')}

library(grid)
library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3,
             top = textGrob("M-GRF", gp=gpar(fontsize=28,font=1)))

# grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3,
#              top = textGrob(expression(M%~~%GRF), gp=gpar(fontsize=28,font=1)))




mean(obj1[[4]]$y_sim)
mean(obj2[[4]]$y_sim)
mean(obj3[[4]]$y_sim)
mean(obj4[[4]]$y_sim)
mean(obj5[[4]]$y_sim)
mean(obj6[[4]]$y_sim)



mean(colMeans(df1[, c(-1)]))
mean(colMeans(df2[, c(-1)]))
mean(colMeans(df3[, c(-1)]))
mean(colMeans(df4[, c(-1)]))
mean(colMeans(df5[, c(-1)]))
mean(colMeans(df6[, c(-1)]))

mean(sapply(df1[, c(-1)], sd))
mean(sapply(df2[, c(-1)], sd))
mean(sapply(df3[, c(-1)], sd))
mean(sapply(df4[, c(-1)], sd))
mean(sapply(df5[, c(-1)], sd))
mean(sapply(df6[, c(-1)], sd))

