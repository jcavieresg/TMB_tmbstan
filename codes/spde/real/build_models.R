### ------------------------------------------------------------

## setwd("C:/Users/Usuario/Desktop/Prueba")

library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

#options(contrasts=c('contr.sum','contr.poly'))


## Step 1: read in and prep the data, and compile and link model
data_full = read.csv("north2.csv", header=T)

data_full$year = as.factor(data_full$year)
data_full$trim = as.factor(data_full$trim)
data_full$destine = as.factor(data_full$destine)
data_full$site = as.factor(data_full$site)


# contrasts(data_full$year) = contr.sum
# contrasts(data_full$trim) = contr.sum
# contrasts(data_full$destine) = contr.sum
# contrasts(data_full$site) = contr.sum

#
#
#d <- subset(data.full)
# glimpse(d) para comprobar


## ## Otherwise can simulate a set like this. See function for more control.
## Compile and link TMB model
m <- "spatial_tmbstan"
#clean.TMB.files(m)
TMB::compile( paste0(m,".cpp"))
dyn.load( dynlib(m))


### ------------------------------------------------------------
## Step 2. Models= no spatial effect (NS), spatial model (S)
#================================================================
#                        Run models in TMB
#================================================================

obj1 = run_tmb(data_full, model='NS', likelihood = 1)
obj2 = run_tmb(data_full, model='S', likelihood = 1)
obj3 = run_tmb(data_full, model='NS', likelihood = 2)
obj4 = run_tmb(data_full, model='S', likelihood = 2)

fits_TMB <- list(obj1, obj2, obj3, obj4)

saveRDS(fits_TMB, file='fits_TMB.RDS')





##==================
##     Cleanup
##==================
dyn.unload(dynlib(m))

