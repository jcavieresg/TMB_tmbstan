rm(list = ls())
setwd("C:/Users/Usuario/Desktop/tps_vs_spde/tps_rspde")

library(pacman)
pacman::p_load(geoR, fields, prodlim, TMB, TMBhelper, mgcv, dplyr, tmbstan, parallel, MASS, Matrix,
               raster, ggplot2, gridExtra, bayesplot, grid, posterior)

options(scipen=999)
# Calculate the number of cores
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

k = round(nloc*0.5)
tp_setup = gam(z ~ x0 + s(s1, s2, bs = "tp", k = k),
                data = df,
                fit = FALSE)

Stp = tp_setup$smooth[[1]]$S[[1]]
Xtp = tp_setup$X[, c(-1, -2)]
E = ginv(Stp)
log_lam = 0

x = rmvn(n = 1, mu = rep(0, nrow(E)), V = exp(log_lam) * E) 

knots_app <- tp_setup$smooth[[1]]$bs.dim
knots_app

#===========================================================================
#                          Simulate data
#===========================================================================
beta0 = 1.0
beta1 = 2.0
sigma_ini = 0.3
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
obj1 <- run_tmb(50)
obj2 <- run_tmb(100)
obj3 <- run_tmb(150)
obj4 <- run_tmb(200)
obj5 <- run_tmb(250)
obj6 <- run_tmb(300)


#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================


M = list()
M[[1]] = list()
M[[1]]$model = "spatial_n50"
M[[2]] = list()
M[[2]]$model = "spatial_n100"
M[[3]] = list()
M[[3]]$model = "spatial_n150"
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
                 iter = 4000, warmup= 700, cores=no_cores,
                 init = 'last.par.best', seed = 12345)
  endTime <- Sys.time()
  timeUsed = difftime(endTime, startTime, units='mins')
  print(timeUsed)
  saveRDS(fit, file=paste0('stan_tps_', i,'.RDS'))
}

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


max(mon6$Rhat)
min(mon6$Bulk_ESS)
min(mon6$Tail_ESS)
tps6_df <- as_draws_df(m6)
print(summarise_draws(tps6_df, "mean", "mcse_mean"), n = 100)

params = rstan::extract(m6, permuted=FALSE, inc_warmup=TRUE)

#plot(c(1,3), c(1,3), ty='n', xlab='beta0', ylab='beta0')
#lines(params[,'chain:1','beta0'], params[,'chain:1','beta0'], col='black',ty='o', pch=20)
#lines(params[,'chain:2','beta0'], params[,'chain:2','beta0'], col='orange',ty='o', pch=20)
#lines(params[,'chain:3','beta0'], params[,'chain:3','beta0'], col='red',ty='o', pch=20)
#lines(params[,'chain:4','beta0'], params[,'chain:4','beta0'], col='gray',ty='o', pch=20)
#legend('topright', legend=c('chain1', 'chain2', 'chain3'), col=c('black', 'orange', 'red'), lty='solid', 
#       bty='n')


#plot(c(1.5, 2.5), c(0,2), ty='n', xlab='beta0', ylab='beta1')
#lines(params[,'chain:1','beta0'], params[,'chain:1','beta1'], col='black',ty='o', pch=20)
#lines(params[,'chain:2','beta0'], params[,'chain:2','beta1'], col='orange',ty='o', pch=20)
#lines(params[,'chain:3','beta0'], params[,'chain:3','beta1'], col='grey',ty='o', pch=20)
#legend('topright', legend=c('chain1', 'chain2', 'chain3'), col=c('black', 'orange', 'grey'), lty='solid', bty='n')


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
  # ggtitle("Grid 1") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  #annotate("text", x = -Inf, y = Inf, hjust = 0.05, vjust = 0.8, label = "Grid 1", colour = "black", size = 10)
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP 1", colour = "blue", size = 10)

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
  # ggtitle("Grid 2") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  # annotate("label", x = -Inf, y = Inf, hjust = 0.05, vjust = 0.8, label = "Grid 2", colour = "black", size = 7)
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP 2", colour = "blue", size = 10)

  
  # geom_label(aes(x = -Inf, y = Inf, hjust = 0.05, vjust = 0.8, label = "Grid 2", 
  #                colour = "white", fontface = "bold"), fill = "black") 
for (i in 2:ncol(df2)) {
  p2 <- p2 + stat_function(fun = dnorm, 
                         args = list(mean = mean(df2[, i]), 
                                     sd = sd(df2[, i])), lwd = 1, col = 'orange')}



# Histogram with kernel density
p3 <- ggplot(df3, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  # ggtitle("Grid 3") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  # annotate("label", x = -Inf, y = Inf, hjust = 0.05, vjust = 0.8, label = "Grid 3", colour = "black", size = 7)
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP 3", colour = "blue", size = 10)


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
  # ggtitle("Grid 4") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  #annotate("label", x = -Inf, y = Inf, hjust = 0.05, vjust = 0.8, label = "Grid 4", colour = "black", size = 7)
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP 4", colour = "blue", size = 10)

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
  # ggtitle("Grid 5") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
   #annotate("label", x = -Inf, y = Inf, hjust = 0.05, vjust = 0.8, label = "Grid 5", colour = "black", size = 7)
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "SP 5", colour = "blue", size = 10)

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
  # ggtitle("Grid 6") + 
  # theme(plot.title = element_text(size = 22, hjust = 0.5))
  #annotate("label", x = -Inf, y = Inf, hjust = 0.05, vjust = 0.8, label = "Grid 6", colour = "black", size = 7)
  #annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, label = "Grid 1", colour = "blue", size = 10)
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Grid 6", colour = "blue", size = 10)

  
for (i in 2:ncol(df6)) {
  p6 <- p6 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df6[, i]), 
                                       sd = sd(df6[, i])), lwd = 1, col = 'orange')}


# title <- textGrob("M-TPS", gp=gpar(fontsize=28), vjust = -0.5)
# 
# title <- arrangeGrob(zeroGrob(), title, 
#                      widths = unit(1, 'npc'), 
#                      heights = unit(c(1.5, 1), c('cm', 'npc')),
#                      as.table = FALSE)
# 
# grid.arrange(p1, p2, p3, p4, ncol = 2, top = title)

library(grid)
library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3,
             top = textGrob("M-TPS", gp=gpar(fontsize=28,font=1)))


# Means of the simulated data
mean(colMeans(df1[, c(-1)]))
mean(colMeans(df2[, c(-1)]))
mean(colMeans(df3[, c(-1)]))
mean(colMeans(df4[, c(-1)]))
mean(colMeans(df5[, c(-1)]))
mean(colMeans(df6[, c(-1)]))

# Posterior means
mean(obj1[[4]]$y_sim)
mean(obj2[[4]]$y_sim)
mean(obj3[[4]]$y_sim)
mean(obj4[[4]]$y_sim)
mean(obj5[[4]]$y_sim)
mean(obj6[[4]]$y_sim)

# Posterior sd
mean(sapply(df1[, c(-1)], sd))
mean(sapply(df2[, c(-1)], sd))
mean(sapply(df3[, c(-1)], sd))
mean(sapply(df4[, c(-1)], sd))
mean(sapply(df5[, c(-1)], sd))
mean(sapply(df6[, c(-1)], sd))


# grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3,
#              top = textGrob("M-TPS",gp=gpar(fontsize=28,font=1), vjust = 0.5))

#===============================================
#          With face_wrap() function!
#===============================================


df1$Grid <- ifelse(df1$y_sim<6.0, "Grid1", "NA")
df2$Grid <- ifelse(df2$y_sim<6.0, "Grid2", "NA")
df3$Grid <- ifelse(df3$y_sim<6.0, "Grid3", "NA")
df4$Grid <- ifelse(df4$y_sim<6.0, "Grid4", "NA")
df5$Grid <- ifelse(df5$y_sim<6.0, "Grid5", "NA")
df6$Grid <- ifelse(df6$y_sim<6.0, "Grid6", "NA")

df_plot1 <- bind_rows(df1, df2, df3, df4, df5, df6)

ptot1 <- ggplot(df_plot1, aes(y_sim)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  facet_wrap(~Grid) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 18)) 

for (i in 2:ncol(df_plot1)-1) {
  ptot1 <- ptot1 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df_plot1[, i]), 
                                       sd = sd(df_plot1[, i])), lwd = 1, col = 'orange')}

# ptot1 + ggtitle("M-TPS") + 
#   theme(plot.title = element_text(color="black", size=22, hjust=0.5)) + facet_wrap(~Grid)


