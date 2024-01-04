### Directory
setwd("C:/Users/Usuario/Desktop/Spatial")


### ------------------------------------------------------------
## Step 1: prepare workspace and load data
source('startup.R')

### ------------------------------------------------------------
## Step 2: Build the TMB models
source('build_models.R')

### ------------------------------------------------------------

## Step 3: Run Bayesian integration on all models
## answer = winDialog("yesno", "Run models with 'tmbstan'?")
answer <- "NO"
source("run_tmbstan.R")

## Step 4: Compare models with loo
source("build_models.R") # rebuild TMB objects
source('loo.R')

## Step 5: Diagnostics for the best model
source('diagnostics.R')


# #### ------------------------------------------------------------
# # ## Step 5: Create plots, tables, and figures
# source('make_figures.R')



