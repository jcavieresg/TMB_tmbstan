rm(list = ls())
setwd("")

library(pacman)
pacman::p_load(TMB, rSPDE, fmesher, INLA, dplyr, tmbstan, rstan, parallel,
               raster, ggplot2, grid, gridExtra)

options(scipen=999)
no_cores <- detectCores() - 1

# #=================================================
# # Compilamos el modelo y lo cargamos en R
TMB::compile('normal_tps.cpp')
dyn.load(dynlib("normal_tps"))
# # #=================================================

# Fix the ssed
set.seed(1234)

#==================================================
#                 SIMULATED DATA
#==================================================
run_tmb <- function(nloc){
set.seed(1234)

nloc = nloc
x <- runif(nloc, 0, 10)
y <- runif(nloc, 0, 10)
  
coords = cbind(x, y)
df <- data.frame(coords)

df$z = 0
df$x0 = 0

colnames(df) <- c("s1", "s2", "z", "x0")

#====================================================
#      Getting matrices from mgcv
#====================================================
tp_setup = gam(z ~ x0 + s(s1, s2, bs = "tp", 
                          k = ifelse(round(length(df$s1)*0.1, digits = 0) < 30, 30, 
                              round(length(df$s1)*0.1, digits = 0))),
                data = df,
                fit = FALSE)

Stp = tp_setup$smooth[[1]]$S[[1]]
Xtp = tp_setup$X[, c(-1, -2)]
E = ginv(Stp)
log_lam = 0

x = rmvn(n = 1, mu = rep(0, nrow(E)), V = exp(log_lam) * E) 

knots_app <- tp_setup$smooth[[1]]$bs.dim

#===========================================================================
#                          Simulate data
#===========================================================================
beta0 = 1.0
beta1 = 2.0
sigma_ini = 0.1
x1 = runif(df$x, 0, 1)

mu_sim = beta0 + beta1*x1 + Xtp %*% x
y_sim = mu_sim + rnorm(nrow(df), mean = 0, sd = sigma_ini)

df$tps = Xtp %*% x
df$y_sim = y_sim
df$x1 <- x1

#======================================
#                TMB data
#======================================
tmb_data <- list(y = df$y_sim,
                 x1 = df$x1, 
                 X = Xtp,
                 S = as(Stp, "sparseMatrix"))

#=====================================
#            TMB parameters
#=====================================
tmb_par <- list(beta0 = 0.5,
                beta1 = 0.5,
                logsigma = 0.5,
                x = rep(0, nrow(Stp)),
                loglambda = 0.5)

obj <- MakeADFun(tmb_data, random = c("x"), tmb_par, DLL="normal_tps", hessian = T)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdrep = sdreport(obj, getJointPrecision = TRUE)

res_list <- list(obj, opt, sdrep, df, tmb_data, tmb_par, knots_app)
return(res_list)
}


#=====================================
#             Run the TMB models
#=====================================
obj1 <- run_tmb(100)
obj2 <- run_tmb(200)
obj3 <- run_tmb(300)
obj4 <- run_tmb(400)
obj5 <- run_tmb(500)
obj6 <- run_tmb(600)

#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================
M = list()
M[[1]] = list()
M[[1]]$model = "spatial_100"
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
  saveRDS(fit, file=paste0('stan_tps_', i,'.RDS'))
}

warnings()

#posteriors  = readRDS("C:/Users/Usuario/Desktop/tps_vs_spde/example/posteriors_tmbstan.RDS")

m1 <- readRDS('stan_tps_1.RDS')
m2 <- readRDS('stan_tps_2.RDS')
m3 <- readRDS('stan_tps_3.RDS')
m4 <- readRDS('stan_tps_4.RDS')
m5 <- readRDS('stan_tps_5.RDS')
m6 <- readRDS('stan_tps_6.RDS')

mon1 = monitor(m1)
mon2 = monitor(m2)
mon3 = monitor(m3)
mon4 = monitor(m4)
mon5 = monitor(m5)
mon6 = monitor(m6)

sum1 <- data.frame(mon1$mean, mon1$`25%`, mon1$`75%`)
colnames(sum1) <- c("mean", "25%", "75%")
head(sum1)

sum2 <- data.frame(mon2$mean, mon2$`25%`, mon2$`75%`)
colnames(sum2) <- c("mean", "25%", "75%")
head(sum2)

sum3 <- data.frame(mon3$mean, mon3$`25%`, mon3$`75%`)
colnames(sum3) <- c("mean", "25%", "75%")
head(sum3)

sum4 <- data.frame(mon4$mean, mon4$`25%`, mon4$`75%`)
colnames(sum4) <- c("mean", "25%", "75%")
head(sum4)

sum5 <- data.frame(mon5$mean, mon5$`25%`, mon5$`75%`)
colnames(sum5) <- c("mean", "25%", "75%")
head(sum5)

sum6 <- data.frame(mon6$mean, mon6$`25%`, mon6$`75%`)
colnames(sum6) <- c("mean", "25%", "75%")
head(sum6)

# Simulating from the SIMULATE function
mat_sim = matrix(data=NA, nrow=length(obj6[[1]]$simulate()$y_sim), ncol=100)
mat_sim


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = obj6[[1]]$simulate()$y_sim
    }
}

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



# Histogram with kernel density
p1 <- ggplot(df1, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP1", colour = "blue", size = 10)

for (i in 2:ncol(df1)) {
  p1 <- p1 + stat_function(fun = dnorm, 
                         args = list(mean = mean(df1[, i]), 
                                     sd = sd(df1[, i])), lwd = 1, col = 'orange')}

 


# Histogram with kernel density
p2 <- ggplot(df2, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP2", colour = "blue", size = 10)

for (i in 2:ncol(df2)) {
  p2 <- p2 + stat_function(fun = dnorm, 
                         args = list(mean = mean(df2[, i]), 
                                     sd = sd(df2[, i])), lwd = 1, col = 'orange')}


p3 <- ggplot(df3, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP3", colour = "blue", size = 10)


for (i in 2:ncol(df3)) {
  p3 <- p3 + stat_function(fun = dnorm, 
                         args = list(mean = mean(df3[, i]), 
                                     sd = sd(df3[, i])), lwd = 1, col = 'orange')}


p4 <- ggplot(df4, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP4", colour = "blue", size = 10)

for (i in 2:ncol(df4)) {
  p4 <- p4 + stat_function(fun = dnorm, 
                         args = list(mean = mean(df4[, i]), 
                                     sd = sd(df4[, i])), lwd = 1, col = 'orange')}


p5 <- ggplot(df5, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP5", colour = "blue", size = 10)

for (i in 2:ncol(df5)) {
  p5 <- p5 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df5[, i]), 
                                       sd = sd(df5[, i])), lwd = 1, col = 'orange')}


p6 <- ggplot(df6, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP6", colour = "blue", size = 10)

  
for (i in 2:ncol(df6)) {
  p6 <- p6 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df6[, i]), 
                                       sd = sd(df6[, i])), lwd = 1, col = 'orange')}

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3,
             top = textGrob("M-TPS", gp=gpar(fontsize=28,font=1)))


