rm(list = ls())
setwd("")

library(pacman)
pacman::p_load(geoR, fields, prodlim, TMB, TMBhelper, mgcv, dplyr, tmbstan, parallel, MASS, Matrix,
               raster, ggplot2, gridExtra, bayesplot, grid)

options(scipen=999)
# Calculate the number of cores
no_cores <- detectCores() - 1

# #=================================================
TMB::compile('tps_sim1.cpp')
dyn.load(dynlib("tps_sim1"))
# # #=================================================

# Fix the ssed
set.seed(1234)

#==================================================
#                 SIMULATED DATA
#==================================================

run_tmb <- function(nloc){

set.seed(1234)

nloc = nloc
s1 <- runif(nloc, 0, 10)
s2 <- runif(nloc, 0, 10)
  
coords = cbind(s1, s2)
df <- data.frame(coords)

df$z = 0
df$x0 = 0

colnames(df) <- c("s1", "s2", "z", "x0")

#====================================================
#      Getting matrices from mgcv
#====================================================
#k = round(nloc*0.5)
tp_setup = gam(z ~ x0 + s(s1, s2, bs = "tp"), 
                          #k = ifelse(round(length(df$s1)*0.1, digits = 0) < 30, 30, 
                          #           round(length(df$s1)*0.1, digits = 0))),
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
sigma_e = 0.1
x1 = runif(df$s1, 0, 1)

y_sim = beta0 + beta1*x1 + Xtp %*% x + rnorm(nrow(df))*sigma_e


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

obj <- MakeADFun(tmb_data, random = c("x"), tmb_par, DLL="tps_sim1", hessian = T)
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
  saveRDS(fit, file=paste0('tps_sim1_', i,'.RDS'))
}


m1 <- readRDS('tps_sim1_1.RDS')
m2 <- readRDS('tps_sim1_2.RDS')
m3 <- readRDS('tps_sim1_3.RDS')
m4 <- readRDS('tps_sim1_4.RDS')
m5 <- readRDS('tps_sim1_5.RDS')
m6 <- readRDS('tps_sim1_6.RDS')

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


for(j in 1:ncol(mat_sim)){
  for(i in 1:nrow(mat_sim)){
    mat_sim[, j] = obj6[[1]]$simulate()$y_sim
    }
}


df1_tps <- data.frame(obj1[[4]]$y_sim, mat_sim)
names(df1_tps)[names(df1_tps) == 'obj1..4...y_sim'] <- 'data'

df2_tps <- data.frame(obj2[[4]]$y_sim, mat_sim)
names(df2_tps)[names(df2_tps) == 'obj2..4...y_sim'] <- 'data'

df3_tps <- data.frame(obj3[[4]]$y_sim, mat_sim)
names(df3_tps)[names(df3_tps) == 'obj3..4...y_sim'] <- 'data'

df4_tps <- data.frame(obj4[[4]]$y_sim, mat_sim)
names(df4_tps)[names(df4_tps) == 'obj4..4...y_sim'] <- 'data'

df5_tps <- data.frame(obj5[[4]]$y_sim, mat_sim)
names(df5_tps)[names(df5_tps) == 'obj5..4...y_sim'] <- 'data'

df6_tps <- data.frame(obj6[[4]]$y_sim, mat_sim)
names(df6_tps)[names(df6_tps) == 'obj6..4...y_sim'] <- 'data'



mean(obj1[[4]]$y_sim)
mean(obj2[[4]]$y_sim)
mean(obj3[[4]]$y_sim)
mean(obj4[[4]]$y_sim)
mean(obj5[[4]]$y_sim)
mean(obj6[[4]]$y_sim)

sd(obj1[[4]]$y_sim)
sd(obj2[[4]]$y_sim)
sd(obj3[[4]]$y_sim)
sd(obj4[[4]]$y_sim)
sd(obj5[[4]]$y_sim)
sd(obj6[[4]]$y_sim)

mean(colMeans(df1_tps[, c(-1)]))
mean(colMeans(df2_tps[, c(-1)]))
mean(colMeans(df3_tps[, c(-1)]))
mean(colMeans(df4_tps[, c(-1)]))
mean(colMeans(df5_tps[, c(-1)]))
mean(colMeans(df6_tps[, c(-1)]))

mean(sapply(df1_tps[, c(-1)], sd))
mean(sapply(df2_tps[, c(-1)], sd))
mean(sapply(df3_tps[, c(-1)], sd))
mean(sapply(df4_tps[, c(-1)], sd))
mean(sapply(df5_tps[, c(-1)], sd))
mean(sapply(df6_tps[, c(-1)], sd))



# Histogram with kernel density
p1 <- ggplot(df1_tps, aes(data)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "a)", colour = "black", size = 7)

for (i in 2:ncol(df1_tps)) {
  p1 <- p1 + stat_function(fun = dnorm, 
                         args = list(mean = mean(df1_tps[, i]), 
                                     sd = sd(df1_tps[, i])), lwd = 1, col = 'orange')}

 


# Histogram with kernel density
p2 <- ggplot(df2_tps, aes(data)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "b)", colour = "black", size = 7)


for (i in 2:ncol(df2_tps)) {
  p2 <- p2 + stat_function(fun = dnorm, 
                         args = list(mean = mean(df2_tps[, i]), 
                                     sd = sd(df2_tps[, i])), lwd = 1, col = 'orange')}



# Histogram with kernel density
p3 <- ggplot(df3_tps, aes(data)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "c)", colour = "black", size = 7)


for (i in 2:ncol(df3_tps)) {
  p3 <- p3 + stat_function(fun = dnorm, 
                         args = list(mean = mean(df3_tps[, i]), 
                                     sd = sd(df3_tps[, i])), lwd = 1, col = 'orange')}



# Histogram with kernel density
p4 <- ggplot(df4_tps, aes(data)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "d)", colour = "black", size = 7)

for (i in 2:ncol(df4_tps)) {
  p4 <- p4 + stat_function(fun = dnorm, 
                         args = list(mean = mean(df4_tps[, i]), 
                                     sd = sd(df4_tps[, i])), lwd = 1, col = 'orange')}


# Histogram with kernel density
p5 <- ggplot(df5_tps, aes(data)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "e)", colour = "black", size = 7)

for (i in 2:ncol(df5_tps)) {
  p5 <- p5 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df5_tps[, i]), 
                                       sd = sd(df5_tps[, i])), lwd = 1, col = 'orange')}



# Histogram with kernel density
p6 <- ggplot(df6_tps, aes(data)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, color="black", fill="grey") + theme_bw() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "f)", colour = "black", size = 7)

  
for (i in 2:ncol(df6_tps)) {
  p6 <- p6 + stat_function(fun = dnorm, 
                           args = list(mean = mean(df6_tps[, i]), 
                                       sd = sd(df6_tps[, i])), lwd = 1, col = 'orange')}


library(grid)
library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3,
             top = textGrob(expression(MTPS[app]), gp=gpar(fontsize=24,font=1)))




df1_tps$Grid <- ifelse(df1_tps$y_sim<6.0, "Grid1", "NA")
df2_tps$Grid <- ifelse(df2_tps$y_sim<6.0, "Grid2", "NA")
df3_tps$Grid <- ifelse(df3_tps$y_sim<6.0, "Grid3", "NA")
df4_tps$Grid <- ifelse(df4_tps$y_sim<6.0, "Grid4", "NA")
df5_tps$Grid <- ifelse(df5_tps$y_sim<6.0, "Grid5", "NA")
df6_tps$Grid <- ifelse(df6_tps$y_sim<6.0, "Grid6", "NA")

df_plot1 <- bind_rows(df1_tps, df2_tps, df3_tps, df4_tps, df5_tps, df6_tps)

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






#==========================================================
#                    L9 versus ~L9
#==========================================================

# L9

# To get the posterior samples of the spatial field
SD0 <- TMB::sdreport(obj9[[1]], getJointPrecision=TRUE,
                     bias.correct = TRUE,
                     bias.correct.control = list(sd = TRUE))
## take samples from fitted model
mu <- c(SD0$par.fixed,SD0$par.random)

## simulate draws
rmvnorm_prec <- function(mu, chol_prec, n.sims) {
  z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L <- chol_prec #Cholesky(prec, super=TRUE)
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  mu + z
}


L <- Cholesky(SD0[['jointPrecision']], super = T)
t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = 100)

## summarize the draws
parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
smoothCoefs_tmb_draws  <- t.draws[parnames == "x",]


# project from mesh to raster, add intercept
pred_tmb <- as.matrix(obj9[[5]]$X %*% smoothCoefs_tmb_draws)

## find the median and sd across draws, as well as 90% intervals
summ_tmb <- cbind(median = (apply(pred_tmb, 1, median)),
                  sd     = (apply(pred_tmb, 1, sd)),
                  lower = (apply(pred_tmb, 1, quantile, .05)),
                  upper = (apply(pred_tmb, 1, quantile, .95)))

head(summ_tmb)


# Tps simulated
tps_sim = obj9[[4]]$tps

tps.rast <- raster(nrows=30, ncols=30,
                   #xmn=0, xmx=10, ymn=0, ymx=10,
                   vals= as.numeric(tps_sim))


## make summary rasters
ras_med_tmb <- ras_sdv_tmb <- ras_lower_tmb <- ras_upper_tmb <- ras_inInt_tmb <- tps.rast
values(ras_med_tmb)   <- summ_tmb[, 1]
values(ras_sdv_tmb)   <- summ_tmb[, 2]
values(ras_lower_tmb) <- summ_tmb[, 3]
values(ras_upper_tmb) <- summ_tmb[, 4]
values(ras_inInt_tmb) <- 0
ras_inInt_tmb[tps.rast < ras_lower_tmb | ras_upper_tmb < tps.rast] <- 1



# Plot prediction
xhat = summary(obj9[[3]], "random")[,1]
tps_sim = obj9[[4]]$tps

y_pred = obj9[[2]]$par[1] + obj9[[2]]$par[2]*obj9[[4]]$x1 + obj9[[5]]$X %*% as.numeric(xhat)

obj9[[4]]$y_pred = y_pred
obj9[[4]]$tps_pred = obj9[[5]]$X %*% xhat
obj9[[4]]$tps_sim <- tps_sim 

# From TMB draws
obj9[[4]]$tps.rast <- tps.rast@data@values
obj9[[4]]$post_mean_tmb <- ras_med_tmb@data@values
obj9[[4]]$post_sd_tmb <- ras_sdv_tmb@data@values


w= 0.02 * (max(obj9[[4]]$s1)-min(obj9[[4]]$s1))
h= 0.02 * (max(obj9[[4]]$s2)-min(obj9[[4]]$s2))


# y_sim
y_sim =  ggplot(obj9[[4]], aes(s1,s2,fill=y_sim)) +
  geom_tile(width=w, height=h) + scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + 
  ggtitle(label = expression(paste("a)"))) + 
  theme(plot.title = element_text(color = "black", size = 22, face = "bold", hjust = 0.5), 
        legend.title = element_text(size=18),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=18), 
        axis.title.y=element_text(size=18, face = "bold")) +
  theme(legend.text = element_text(size=12))+
  theme(plot.margin = unit(c(0.2,0,0.2,0), "cm")) +
  labs(fill = "data")

# y_pred
y_pred =  ggplot(obj9[[4]], aes(s1,s2,fill=y_pred)) +
  geom_tile(width=w, height=h) + scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + 
  ggtitle(label = expression(paste("b)"))) + 
  theme(plot.title = element_text(color = "black", size = 22, face = "bold", hjust = 0.5),
        legend.title = element_text(size=18),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=18), 
        axis.title.y=element_blank()) +
  theme(legend.text = element_text(size=12))+
  theme(plot.margin = unit(c(0.2,0,0.2,0), "cm"))+
  labs(fill = "pred")










# To get the posterior samples of the spatial field
SD02 <- TMB::sdreport(obj9new[[1]], getJointPrecision=TRUE,
                      bias.correct = TRUE,
                      bias.correct.control = list(sd = TRUE))
## take samples from fitted model
mu2 <- c(SD02$par.fixed,SD02$par.random)


L2 <- Cholesky(SD02[['jointPrecision']], super = T)
t.draws2 <- rmvnorm_prec(mu = mu2 , chol_prec = L2, n.sims = 100)

## summarize the draws
parnames2 <- c(names(SD02[['par.fixed']]), names(SD02[['par.random']]))
smoothCoefs_tmb_draws2  <- t.draws2[parnames2 == "x",]


# project from mesh to raster, add intercept
pred_tmb2 <- as.matrix(obj9new[[5]]$X %*% smoothCoefs_tmb_draws2)

## find the median and sd across draws, as well as 90% intervals
summ_tmb2 <- cbind(median = (apply(pred_tmb2, 1, median)),
                   sd     = (apply(pred_tmb2, 1, sd)),
                   lower = (apply(pred_tmb2, 1, quantile, .05)),
                   upper = (apply(pred_tmb2, 1, quantile, .95)))

head(summ_tmb2)


# Tps simulated
tps_sim2 = obj9new[[4]]$tps

tps.rast2 <- raster(nrows=30, ncols=30, vals= as.numeric(tps_sim2))


## make summary rasters
ras_med_tmb2 <- ras_sdv_tmb2 <- ras_lower_tmb2 <- ras_upper_tmb2 <- ras_inInt_tmb2 <- tps.rast2
values(ras_med_tmb2)   <- summ_tmb2[, 1]
values(ras_sdv_tmb2)   <- summ_tmb2[, 2]
values(ras_lower_tmb2) <- summ_tmb2[, 3]
values(ras_upper_tmb2) <- summ_tmb2[, 4]
values(ras_inInt_tmb2) <- 0
ras_inInt_tmb2[tps.rast2 < ras_lower_tmb2 | ras_upper_tmb2 < tps.rast2] <- 1

# Plot prediction
xhat2 = summary(obj9new[[3]], "random")[,1]
y_pred2 = obj9new[[2]]$par[1] + obj9new[[2]]$par[2]*obj9new[[4]]$x1 + obj9new[[5]]$X %*% as.numeric(xhat2)
obj9new[[4]]$y_pred = y_pred2
obj9new[[4]]$tps_pred = obj9new[[5]]$X %*% xhat2

# From TMB draws
obj9new[[4]]$tps.rast2 <- tps.rast2@data@values
obj9new[[4]]$post_mean_tmb2 <- ras_med_tmb2@data@values
obj9new[[4]]$post_sd_tmb2 <- ras_sdv_tmb2@data@values


w2= 0.02 * (max(obj9new[[4]]$s1)-min(obj9new[[4]]$s1))
h2= 0.02 * (max(obj9new[[4]]$s2)-min(obj9new[[4]]$s2))


# y_sim
y_sim2 =  ggplot(obj9new[[4]], aes(s1,s2,fill=y_sim)) +
  geom_tile(width=w2, height=h2) + scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + 
  #ggtitle(label = expression(paste(""%~~%SL9))) +
  ggtitle(label = expression(paste("c)"))) +
  theme(plot.title = element_text(color = "black", size = 22, face = "bold", hjust = 0.5), 
        legend.title = element_text(size=18),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=18), 
        axis.title.y=element_text(size=18, face = "bold")) +
  theme(legend.text = element_text(size=12))+
  theme(plot.margin = unit(c(0.2,0,0.2,0), "cm")) +
  labs(fill = "data")


# y_pred
y_pred2 =  ggplot(obj9new[[4]], aes(s1,s2,fill=y_pred)) +
  geom_tile(width=w2, height=h2) + scale_fill_gradient(low = "white", high = "red") +
  coord_fixed() + 
  #ggtitle(label = expression(paste(""%~~%SL9))) + 
  ggtitle(label = expression(paste("d)"))) +
  theme(plot.title = element_text(color = "black", size = 22, face = "bold", hjust = 0.5),
        legend.title = element_text(size=18),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=18), 
        axis.title.y=element_blank()) +
  theme(legend.text = element_text(size=12))+
  theme(plot.margin = unit(c(0.2,0,0.2,0), "cm"))+
  labs(fill = "pred")



grid.arrange(y_sim, y_pred, y_sim2, y_pred2, ncol = 2)#,

# setwd("C:/Users/Usuario/Desktop/eps_figures")
# setEPS()
# postscript("Fig11.eps")
# grid.arrange(y_sim, y_pred, y_sim2, y_pred2, ncol = 2)#,
# dev.off()



setwd("C:/Users/Usuario/Desktop/eps_figures")
setEPS()
postscript("Fig11.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 14, height = 8)
grid.arrange(y_sim, y_pred, y_sim2, y_pred2, ncol = 2)#,
dev.off()






# RSME




scientific(sqrt(mean((obj9[[4]]$y_sim - obj9[[4]]$y_pred))), digits = 3)
scientific(sqrt(mean((obj9new[[4]]$y_sim - obj9new[[4]]$y_pred))), digits = 3)
