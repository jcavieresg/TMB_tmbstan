rm(list = ls())
setwd("C:/Users/joaquin/Desktop/stan_tmb")

library(TMB)
library(ggplot2)
library(INLA)
library(reshape2)
library(mgcv)
library(tmbstan)
library(MASS)
library(tidyverse)
library(Matrix)
library(gridExtra)
library(sp)
library(rgdal)
library(raster)
library(leaflet)
library(mapview)
library(tictoc)
library(parallelly)
library(tweedie)


options(scipen=999)

# Calculate the number of cores
no_cores <- parallelly::availableCores() 

#=================================================
# Compilamos el modelo y lo cargamos en R
compile("cpue_spatial_spde.cpp")
dyn.load(dynlib("cpue_spatial_spde"))
#=================================================





## Load data and prepare for model
load('EBS_pollock_data.rda')
pollock <- EBS_pollock_data
dim(pollock)

pollock <- subset(pollock, year==2014)
hist(pollock$catch)

# ## Matrix of pairwise differences
# dd <- as.matrix(dist(pollock[,c('lat', 'long')]))
# n <- length(pollock$catch)
# pollock$lat <- pollock$lat-mean(pollock$lat)
# pollock$long <- pollock$long-mean(pollock$long)

ggplot(pollock, aes(long, lat, size=(catch), color=catch==0)) + geom_jitter(alpha=.5)
# compile("tmb_models/cpue_spatial.cpp")
# dyn.load(dynlib("tmb_models/cpue_spatial"))
# pars = list(intercept=3, theta=-5, beta_lat=0, beta_lon=0,
#             logsigma=0, logsigma_space=1, a=1, u=rep(0, n))
# data = list(likelihood=2, y=pollock$catch, lat=pollock$lat,
#             lon=pollock$long, dd=dd)
# obj = MakeADFun(data=data, parameters=pars, random='u', DLL="cpue_spatial")
# #obj$env$beSilent()
# opt <- with(obj, nlminb(par, fn, gr))
# opt <- with(obj, nlminb(opt$par, fn, gr)) # restart optimizer
# rep1 <- obj$report()

### The SPDE approach to simplify the model. Using RINLA tools to create
### inputs to TMB. Based of Jim Thorson's function.
create.spde <- function(dat, n_knots, make.plot=FALSE, jitter=.3){
  loc_xy <- data.frame(x=dat$long, y=dat$lat)
  knots <- kmeans( x=loc_xy, centers=n_knots )
  loc_centers <- knots$centers
  ## loc_xy <- cbind(loc_xy, cluster=knots$cluster)
  ## ggplot(loc_xy, aes(x,y, col=factor(cluster))) + geom_point(size=.5)
  mesh <- inla.mesh.create( loc_centers, refine=TRUE)
  spde <- inla.spde2.matern( mesh )
  if(make.plot){
    png(paste0('mesh_', n_knots, '.png'), width=7, height=4, units='in', res=500)
    par(mar=.5*c(1,1,1,1))
    plot(mesh, main=NA, edge.color=gray(.7))
    points( jitter(dat$longitude, amount=jitter), jitter(dat$latitude, amount=jitter), cex=1, pch='.', col=rgb(0,0,0,.3))
    points( loc_centers, cex=.5, pch=20, col="red")
    dev.off()
  }
  return(list(mesh=mesh, spde=spde, cluster=knots$cluster, loc_centers=loc_centers))
}
spde <- create.spde(dat=pollock, n_knots=100, make.plot=FALSE)


#Define mesh and components representing the  precision matrix----
#loc = cbind(pollock$lat, pollock$long)
# boundary = INLA::inla.nonconvex.hull(loc)
# boundary2 = INLA::inla.nonconvex.hull(loc,convex = -0.35)
#mesh = inla.mesh.2d(loc = loc,   boundary = list(boundary,boundary2), max.edge=c(5, 10), cutoff=0.05)
#mesh = inla.mesh.2d(loc = loc,  max.edge=c(5, 10))
# mesh$n
# plot(mesh)
# points(pollock$lat, pollock$long, col = "blue")
# 
# A = inla.spde.make.A(mesh,loc)

#spde2 = inla.spde2.matern(mesh, alpha=2)
#spdeMatrices = spde2$param.inla[c("M0","M1","M2")]

spdeMatrices = spde$spde$param.inla[c("M0","M1","M2")]

## Make spde inputs for TMB
data <- list(likelihood=2, 
             y=pollock$catch, 
             # lat=pollock$lat,
             # lon=pollock$long, 
             site=spde$cluster,
             #site = mesh$idx$loc - 1,
             spdeMatrices = spdeMatrices)


pars <- list(intercept=3, 
             theta=-5, 
             logsigma=0,
             logtau=1, 
             logkappa=-3,
             u=rnorm(spde$mesh$n,0,1))
             


# compile("tmb_models/cpue_spatial_spde.cpp")
# dyn.load(dynlib("tmb_models/cpue_spatial_spde"))

## We can use map to turn off spatial component
map <- list(logkappa=factor(NA), logsigma_s=factor(NA), u=factor(rep(NA, length(pars$u))))
obj <- MakeADFun(data=data, parameters=pars, map=map, DLL="cpue_spatial_spde")
#obj$env$beSilent()
opt2 <- with(obj, nlminb(par, fn, gr))
opt2 <- with(obj, nlminb(opt2$par, fn, gr)) # restart
rep2 <- obj$report()

## Rerun with spatial SPDE approach
obj3 = MakeADFun(data=data, parameters=pars, random='u', DLL="cpue_spatial_spde", hessian = TRUE)
#obj$env$beSilent()
opt3 <- with(obj3, nlminb(par, fn, gr))
opt3 <- with(obj3, nlminb(opt3$par, fn, gr)) # restart
opt3 <- with(obj3, nlminb(opt3$par, fn, gr)) # restart
rep3 <- obj3$report()
rep3$u

library(reshape2)
## Get residuals for all three models
pollock$resids.nospace <- (rep2$pred-log(data$y))/rep2$sigma
pollock$resids.spde <- (rep3$pred-log(data$y))/rep3$sigma
resids <- melt(pollock, id.vars = c("lat", "long"), measure.vars = c("resids.nospace", "resids.spde"))
ggplot(resids, aes(long, lat, size=abs(value), color=value)) +
  geom_point(alpha=.5)+ scale_color_gradient(low='red', high='black') +
  facet_wrap('variable') + scale_size_area(max_size=4)

param_spde = cbind(opt3$par[1], opt3$par[2],  exp(opt3$par[3]), exp(opt3$par[4]))
colnames(param_spde) = c("intercept", "theta", "sigma", "kappa")
param_spde
