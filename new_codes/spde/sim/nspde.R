rm(list = ls())
setwd("C:/Users/Usuario/Desktop/tps_vs_spde/example")

library(pacman)
pacman::p_load(geoR, fields, prodlim, TMB, INLA, dplyr, tmbstan, rstan, parallel,
               raster, ggplot2)

options(scipen=999)

# Calculate the number of cores
no_cores <- parallelly::availableCores()  - 1




#=================================================
# Compilamos el modelo y lo cargamos en R
compile("NSPDE_tmbstan.cpp")
dyn.load(dynlib("NSPDE_tmbstan"))


#=========================================
#               Simulation
#=========================================
set.seed(1234)

run_tmb <- function(nloc){
  
local.plot.field = function(field, mesh, xlim=c(0,10), ylim=c(0,10), ...){
    stopifnot(length(field) == mesh$n)
    proj = inla.mesh.projector(mesh, xlim = xlim, 
                               ylim = ylim, dims=c(300, 300))
    field.proj = inla.mesh.project(proj, field)
    n.col = 20
    image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
               xlim = xlim, ylim = ylim, col = plasma(n.col), nlevel=n.col+1, ...)
}  
  
#inla.seed = sample.int(n=1E6, size=1)
options(width=70, digits=3)

loc_ini = matrix(c(0,0,10,10, 0, 10, 10, 0), nrow = 4, byrow = T)
#loc_ini = cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
mesh_sim = inla.mesh.2d(loc = loc_ini, max.edge=c(1.7, 2))

mesh_sim$n


spde = inla.spde2.matern(mesh_sim, alpha = 2)
spde_mat = spde$param.inla[c("M0","M1","M2")]

range = 2
kappa = sqrt(8)/ range
tau = 0.5

Qgrf = inla.spde.precision(spde, theta=c(log(kappa), log(tau)))
grf = inla.qsample(n=1, Q=Qgrf)
grf = grf[ ,1]


local.plot.field(grf, mesh_sim)


#=============================================
#            Number of locations
#=============================================
# n = nloc
# coords = matrix(runif(2*n), n)*10 # coordinates

nx = nloc
ny = nloc
df = expand.grid(x = seq(0, 10, length.out = nloc),
                 y = seq(0, 10, length.out = nloc))

coords = as.matrix(df)

A = inla.spde.make.A(mesh = mesh_sim, loc = coords)
grf = drop(A %*% grf)

quilt.plot(x=coords[, 1],y=coords[, 2],z=grf,nx=80,ny=80, 
           col = plasma(101), main="Field projected to data locations", 
           zlim = range(grf))


#===========================================================================
#                          Simulate data
#===========================================================================
x1 = runif(df$x, 0, 1)
beta0 = 2.0
beta1 = 1.0
sigma_ini = 0.3

mu_sim = beta0 + beta1*x1 + grf

#y_sim = rnorm(n, mean = mu_sim, sd = sigma_ini)
y_sim = mu_sim + rnorm(df$x, mean = 0, sd = sigma_ini)

# data.frame
df = data.frame(y_sim = y_sim, s1= coords[ ,1], s2 = coords[ ,2], x1 = x1, grf = grf)
summary(df)


#=================================================================
#                           TMB modelling
#=================================================================

#======================================
#                TMB data
#======================================
tmb_data <-  list(#model = 2,                  # 1 = (A*u) / tau,    2 = A*u,  3 = u(site)
                  #normalization = 1,          # 1 no scaliing the GMRF, 2 = scalling the GMRF (1/ tau)
                  y    = as.vector(df$y_sim),
                  x1 = as.vector(df$x1),
                  site = mesh_sim$idx$loc - 1,
                  spde_mat = spde_mat,
                  A = A,
                  tau0 = exp(spde$param.inla$theta.initial[1]),
                  kappa0 = exp(spde$param.inla$theta.initial[2]))




#=====================================
#            TMB parameters
#=====================================
tmb_par  <- list(#beta = c(0.1, 0.1),
                 beta0 = 0.1, 
                 beta1 = 0.1,
                 logsigma_e = 0.1,
                 logtau   = 0.1,
                 logkappa = 0.1,
                 u = rnorm(spde$mesh$n,0,1))



#dyn.unload(dynlib("spatial"))

## Create the TMB object
obj <- MakeADFun(data = tmb_data, parameters = tmb_par, random="u", DLL="NSPDE_tmbstan", hessian = T)
opt = nlminb(obj$par, obj$fn, obj$gr)
sdrep = sdreport(obj, getJointPrecision = TRUE)
res_list = list(obj, opt, sdrep, df, tmb_data, tmb_par, mesh_sim)
return(res_list)
}


obj1 <- run_tmb(nloc = 10)
obj2 <- run_tmb(nloc = 15)
obj3 <- run_tmb(nloc = 20)
obj4 <- run_tmb(nloc = 25)





#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================


M = list()
M[[1]] = list()
M[[1]]$model = "spatial_n10"
M[[2]] = list()
M[[2]]$model = "spatial_n15"
M[[3]] = list()
M[[3]]$model = "spatial_n20"
M[[4]] = list()
M[[4]]$model = "spatial_n25"



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
  timeUsed = endTime - startTime
  print(timeUsed)
  saveRDS(fit, file=paste0('spde_fit_', i,'.RDS'))
}

warnings()



c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")


# M[[1]]$res <- readRDS('spde_fit_1.RDS')
# M[[2]]$res <- readRDS('spde_fit_2.RDS')
# M[[3]]$res <- readRDS('spde_fit_3.RDS')
# M[[4]]$res <- readRDS('spde_fit_4.RDS')


m1 <- readRDS('spde_fit_1.RDS')
m2 <- readRDS('spde_fit_2.RDS')
m3 <- readRDS('spde_fit_3.RDS')
m4 <- readRDS('spde_fit_4.RDS')

mon1 = monitor(m1)
mon2 = monitor(m2)
mon3 = monitor(m3)
mon4 = monitor(m4)


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

summary(M[[1]]$res)
summary(M[[2]]$res)


posterior_1 = as.matrix(M[[1]]$res)
posterior_2 = as.matrix(M[[2]]$res)
posterior_3 = as.matrix(M[[3]]$res)
posterior_4 = as.matrix(M[[4]]$res)


traceplot(M[[1]]$res, pars=names(obj1$par), inc_warmup=TRUE)
traceplot(M[[2]]$res, pars=names(obj2$par), inc_warmup=TRUE)





# Simulating from the SIMULATE function
mat_sim = matrix(data=NA, nrow=length(obj4[[1]]$simulate()$y_sim), ncol=150)
mat_sim


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = obj4[[1]]$simulate()$y_sim
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


# Histogram with kernel density
p1 <- ggplot(df1, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        #axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  ggtitle("Grid 1") + 
  theme(plot.title = element_text(size = 22, hjust = 0.5))

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
  ggtitle("Grid 2") + 
  theme(plot.title = element_text(size = 22, hjust = 0.5))

for (i in 2:ncol(df2)) {
  p2 <- p2 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df2[, i]), 
                                       sd = sd(df2[, i])), lwd = 1, col = 'orange')}



# Histogram with kernel density
p3 <- ggplot(df3, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  ggtitle("Grid 3") + 
  theme(plot.title = element_text(size = 22, hjust = 0.5))

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
  ggtitle("Grid 4") + 
  theme(plot.title = element_text(size = 22, hjust = 0.5))

for (i in 2:ncol(df4)) {
  p4 <- p4 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df4[, i]), 
                                       sd = sd(df4[, i])), lwd = 1, col = 'orange')}

library(grid)
library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol = 2,
             top = textGrob("M1",gp=gpar(fontsize=28,font=1)))



















#load package
require(MCMCvis)
MCMCsummary(fit, round = 2)

#pairs(fit, pars=names(fit)[-grep('beta0',names(fit))][1:10])
pairs(fit, pars = c("logsigma_e", "logtau", "lp__"), las = 1) # below the diagonal
pairs(fit, pars = c("logsigma_e", "logtau", "logkappa")) # below the diagonal


pairs(fit, pars = c("logsigma", "logalpha", "lp__"), las = 1) # below the diagonal













# s <- summary(fit_tmbstan, probs = c(0.25, 0.75))
# s$summary  # all chaines merged
# 
# exp(s$summary[3:5, c(1,3)])
# s$summary[1:2, c(1,3)]
# 
# 
# # Extract posterior draws for later use
# posterior_cp <- as.array(fit_tmbstan)


require(MCMCvis)
library(bayesplot)
mcmc_pairs(posterior_1, pars = c("logsigma_e","logkappa","logtau"),
           off_diag_args = list(size = 0.75))


library(gridExtra)
library(tidyverse)
traceplot(M[[1]]$res, pars=names(obj1$par), inc_warmup=TRUE) + 
  theme(strip.text = ggplot2::element_text(size = 16, color = "black"),
        axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14))

#load package
require(MCMCvis)
MCMCsummary(fit_tmbstan, round = 2)
pairs(fit_tmbstan, pars = c("logsigma", "logtau", "logkappa")) # below the diagonal





grfhat = summary(obj4[[3]], "random")[,1]
par_fixed = summary(obj4[[3]], "fixed")
betahat = summary(obj4[[3]], "fixed")["beta0", 1]
obj4[[4]]$grf_est = obj4[[5]]$A %*% grfhat


scl_lims = range(c(obj4[[4]]$grf, as.numeric(obj4[[4]]$grf_est)))
# truth
g_true =  ggplot(obj4[[4]]) +
  geom_tile(aes(x = s1, y = s2, fill = grf)) +
  coord_equal() +
  scale_fill_viridis_c(limits = scl_lims) +
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# estimate
g_est =  ggplot(obj4[[4]]) +
  geom_tile(aes(x = s1, y = s2, fill = as.matrix(grf_est))) +
  coord_equal() +
  scale_fill_viridis_c(limits = scl_lims) +
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

grid.arrange(g_true, g_est, ncol = 1)


# Simulating from the SIMULATE function
mat_sim = matrix(data=NA, nrow=length(obj1[[1]]$simulate()$y_sim), ncol=100)
mat_sim


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = obj1[[1]]$simulate()$y_sim
  }
}
mat_sim


# Plot lognormal non spatial simulations
hist(obj1[[4]]$y_sim, col = "gray90", prob = TRUE, main = "100 simulated samples", cex.main = 2, xlab = "", col.lab = 'blue', col.main="blue", cex.lab = 1.4, cex.axis = 1.2, ylim = c(0, 0.5),
     xlim = c(0, 5))
for (j in 1: ncol(mat_sim)){
  lines(x = density(x = mat_sim[, j]),  lty="dotted", col="azure4", lwd=1)
}
#legend("topright", "A", bty = "n", text.col="blue", cex = 1.4)
legend(6, 0.3, c("response", "simulations"), lwd=4, col=c("gray90", "azure4"), cex = 1.5)


hist(obj3[[2]]$y_sim)
hist(obj2[[1]]$simulate()$y_sim)





rep  = obj4[[3]]

rangeIndex = which(row.names(summary(rep,"report"))=="rho")
fieldIndex = which(row.names(summary(rep,"report"))=="u")
range = summary(rep,"report")[rangeIndex,]

proj = inla.mesh.projector(obj4[[7]])
latentFieldMAP = rep$par.random[names(rep$par.random)=="u"]/exp(rep$par.fixed[which(names(rep$par.fixed)=="logtau")])
#x11()
par(mfrow=c(3,2), mar = c(5, 4, 2.5, 0.5), oma = c(0.5, 0.5, 0.2, 0.2))
image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldMAP),col =  colorRampPalette(c("white","gold", "blue"))(12),
           xlab = '', ylab = 'Latitude',
           main = "Mean of the spatial random effect",
           cex.lab = 1.6,cex.axis = 1.3, cex.main=1.4,
           cex.sub= 1.1,
           axis.args=list(cex.axis=1.3),
           zlim = range(-9, 9))
bnd <- inla.mesh.boundary(obj4[[7]])
inter <- inla.mesh.interior(obj4[[7]])
plot(obj4[[7]], draw.segments=FALSE, main = '', add = T, col = "black")
#points(cbind(obj1[[4]]$s1,obj1[[4]]$s2) ,type = 'p',lwd = 2, pch = 19, cex = 1.2)
legend(-74.01, -41.956, legend=c("Spatial Gamma model"),
       col=c("NA"), cex=1.4, box.lty=1)
#contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldMAP) ,add = T,labcex  = 1,cex = 1)
text(x = -73.9, y = -41.795, expression(bold("S"[5])), cex = 1.8, font=4, col = "black")
text(x = -73.85, y = -41.725, expression(bold("S"[9])), cex = 1.8)
text(x = -73.79, y = -41.82,  expression(bold("S"[3])), cex = 1.8,  col = "black", font =  2)
text(x = -73.76, y = -41.71, expression(bold("S"[7])), cex = 1.8)
text(x = -73.72, y = -41.835, expression(bold("S"[1])), cex = 1.8, col = "black", font = 2)
text(x = -73.70, y = -41.74, expression(bold("S"[8])), cex = 1.8,  col = "black", font = 2)
text(x = -73.66, y = -41.815, expression(bold("S"[2])), cex = 1.8)
text(x = -73.64, y = -41.76, expression(bold("S"[6])), cex = 1.8, col = "black", font  = 2)

text(x = -73.5, y = -41.86, expression(bold("S"[10])), cex = 1.8)

text(x = -72.98, y = -41.64, expression(bold("S"[13])), cex = 1.8)
text(x = -72.98, y = -41.74, expression(bold("S"[4])), cex = 1.8)
text(x = -73, y = -41.82, expression(bold("S"[12])), cex = 1.8)
text(x = -73.1, y = -41.93, expression(bold("S"[11])), cex = 1.8)
