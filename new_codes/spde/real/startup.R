## Libraries used
library(TMB)
library(INLA)
library(tmbstan)
library(rstan)
library(plyr)
library(reshape2)
library(maps)
library(tidyverse)
library(dplyr)
library(shinystan)
library(rstan)
library(bayesplot)
library(tidyverse)
library(utils)


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

## This function will clean up the compiled files from a TMB object. Useful
## for development when need to ensure all change are propogated.
clean.TMB.files <- function(model.path){
  o <- paste0(model.path,'.o')
  dll <- paste0(model.path, '.dll')
  tryCatch(expr=dyn.unload(dynlib(model.path)), error=function(cond){x=3})
  if(file.exists(dll)) trash <- file.remove(dll)
  if(file.exists(o)) trash <- file.remove(o)
}

## Forces elements of x to be unique by adding number to any duplicated
## entries. Helper function for plotting functions below
add.unique.names <- function(x){
  as.vector(unlist(llply(unique(x), function(y){
    z <- x[x %in% y]
    if(length(z)>1) z <- paste0(z,"_", seq_along(z))
    z})))
}

## This function generates an input list from data and other arguments,
## which can then be used in TMB::MakeADFun to build the model. The list
## includes Data, Params (inits), Map, and Random elements.
make.inputs <- function(dat, model, likelihood=1, ...){
  ## Make the mesh grid using INLA
    coords = cbind(dat$longitude, dat$latitude)
    bound1 <- inla.nonconvex.hull(coords)
    mesh = inla.mesh.create(coords, plot.delay=NULL, refine=TRUE, boundary = bound1)
  ##mesh  = readRDS("mesh.RDS")
    spde  =  inla.spde2.matern(mesh, alpha=2)
    loc = spde$mesh$loc[,1:2]
    n_s = nrow(loc)

  ## Make inputs depending on the model and spacing form
  model <- match.arg(model, choices=c('NS', "S"))
  if(model=='NS') space <- 0 # no space
  if(model=='S') space <- 1 # spatial only

    Data <- list(likelihood=likelihood, space=space,
               n_t=length(unique(dat$year)),
               n_s = n_s,
               #length(unique(dat$site)),
               cpue=dat$cpue,
               year=as.numeric(dat$year) - 1,
               depth=dat$depth,
               trim=as.numeric(dat$trim) -1 ,
               destine=as.numeric(dat$destine) -1 ,
               logtaumeanO = -1.7,
               logkappamean = 2,
               site = mesh$idx$loc - 1,
               M0=spde$param.inla$M0, M1=spde$param.inla$M1,
               M2=spde$param.inla$M2)

  Params <- list(intercept = 2,
                 beta_depth = 0,
                 beta_year=rep(0, length(levels(dat$year))),
                 beta_trim=rep(0, length(levels(dat$trim))),
                 beta_destine=rep(0, length(levels(dat$destine))),
                 logsigma = -0.9,
                 #logtauO = spde$param.inla$theta.initial[1],
                 #logtauO = -0.945,
                 logtauO = -15.12,
                 logkappa = spde$param.inla$theta.initial[2],
                 omega_s=rep(0,spde$mesh$n))


  ## Need to fix first level of each factor at 0 so they are
  ## identifiable. Get merged into the intercept. I.e., contrasts in R.
  list.factors <- list(beta_year =    factor(c(0, 1:(length(levels(dat$year)) -1))),
                       beta_trim =     factor(c(0, 1:(length(levels(dat$trim)) -1))),
                       beta_destine = factor(c(0, 1:(length(levels(dat$destine)) -1))))

  Map = list(NS=c(list.factors,
                  list(logtauO = factor(NA), logkappa=factor(NA), omega_s=factor(NA*Params$omega_s))),

    ## Turn off just spatial effects
           S = c(list.factors, list(logtauO=factor(NA))))

    return(list(Data=Data, Params=Params, Random=NULL[[model]], Map=Map[[model]]))
}

model.name <- function(model) switch(model, NS="No Space", S='Space')

## This function runs the spatiotemporal model for the logbook data.
run_tmb <- function(dat, model, likelihood=1, trace=10){

  message(paste0('Starting model ', model))
  model.name <- model.name(model)
  start <- Sys.time()
  Inputs <- make.inputs(dat=dat, model=model, likelihood=likelihood)
  Obj <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map)
  Obj$env$beSilent()
  report.temp <- Obj$report();

runtime = as.numeric(difftime(Sys.time(),start, units='mins'))
x = list(model=model, likelihood=likelihood,
         model.name=model.name,
         runtime=runtime, report=report.temp, Obj=Obj,
         Inits=Inputs$Params, Map=Inputs$Map, Random=Inputs$Random)
return(x)
}
