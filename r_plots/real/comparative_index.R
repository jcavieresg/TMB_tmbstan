setwd("")

library(pacman)
pacman::p_load(TMB, INLA, dplyr, tmbstan, rstan, parallel, bayesplot, ggplot2, loo, gridExtra,
               fitdistrplus)

options(scipen=999)

M1    = readRDS("M1/stan_gamma_spde.RDS")
M2    = readRDS("M2/stan_sn_spde.RDS")
M3    = readRDS("M3/stan_gamma_tps.RDS")
M4    = readRDS("M4/stan_sn_tps.RDS")

mon1 = monitor(M1)
mon2 = monitor(M2)
mon3 = monitor(M3)
mon4 = monitor(M4)

options(max.print = 1000000)




#========================================
#       Gamma models (spde and tps)
#========================================

summary.df <- as.data.frame(summary(M1, 
                                    pars=c("beta_year[1]", "beta_year[2]", 
                                           "beta_year[3]", "beta_year[4]", "beta_year[5]","beta_year[6]", "beta_year[7]", 
                                           "beta_year[8]", "beta_year[9]", "beta_year[10]", "beta_year[11]", "beta_year[12]", 
                                           "beta_year[13]", "beta_year[14]", "beta_year[15]", "beta_year[16]", "beta_year[17]", 
                                           "beta_year[18]", "beta_year[19]", "beta_year[20]", "beta_year[21]")))[, c(1, 5, 6, 7)]

years <- 1996:2016
beta0_1 <- summary(M1)$summary[1, 1] 
index <- summary(M1)$summary[2:22, 1] 

all.df <- cbind(summary.df, index, years) 
colnames(all.df) <- make.names(colnames(all.df)) # to remove % in the col names

head(all.df)

p1 <- ggplot(all.df, mapping = aes(x = years)) + theme_bw() +
  geom_ribbon(aes(ymin = summary.25., ymax = summary.75.), fill = "azure3", alpha = 0.6) +
  geom_line(mapping = aes(y = summary.50.)) +
  geom_point(mapping = aes(y = summary.50.)) +
  scale_x_continuous("", labels = as.character(years), breaks = years) + 
  labs(x = "", y = "Index") + # theme_bw() + 
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "M1", colour = "black", size = 7, fontface = "bold") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(plot.title = element_text(size = 22, hjust = 0.5))

p1 


summary.df2 <- as.data.frame(summary(M2,pars=c("beta_year[1]", "beta_year[2]", 
                                                              "beta_year[3]", "beta_year[4]", "beta_year[5]","beta_year[6]", "beta_year[7]", 
                                                              "beta_year[8]", "beta_year[9]", "beta_year[10]", "beta_year[11]", "beta_year[12]", 
                                                              "beta_year[13]", "beta_year[14]", "beta_year[15]", "beta_year[16]", "beta_year[17]", 
                                                              "beta_year[18]", "beta_year[19]", "beta_year[20]", "beta_year[21]")))[, c(1, 5, 6, 7)]

years2 <- 1996:2016
beta0_2 <- summary(M2)$summary[1, 1] 
index2 <- summary(M2)$summary[2:22, 1]

all.df2 <- cbind(summary.df2, index2, years2) 
colnames(all.df2) <- make.names(colnames(all.df2)) # to remove % in the col names

head(all.df2)

p2 <- ggplot(all.df2, mapping = aes(x = years2)) + theme_bw() +
  geom_ribbon(aes(ymin = summary.25., ymax = summary.75.), fill = "azure3", alpha = 0.6) +
  geom_line(mapping = aes(y = summary.50.)) +
  geom_point(mapping = aes(y = summary.50.)) +
  scale_x_continuous("", labels = as.character(years2), breaks = years2) + 
  labs(x = "", y = "") +# theme_bw() + 
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "M2", colour = "black", size = 7, fontface = "bold") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(plot.title = element_text(size = 22, hjust = 0.5))
p2


summary.df3 <- as.data.frame(summary(M3,
                                     pars=c("beta_year[1]", "beta_year[2]", 
                                     "beta_year[3]", "beta_year[4]", "beta_year[5]","beta_year[6]", "beta_year[7]", 
                                     "beta_year[8]", "beta_year[9]", "beta_year[10]", "beta_year[11]", "beta_year[12]", 
                                     "beta_year[13]", "beta_year[14]", "beta_year[15]", "beta_year[16]", "beta_year[17]", 
                                     "beta_year[18]", "beta_year[19]", "beta_year[20]", "beta_year[21]")))[, c(1, 5, 6, 7)]



years3 <- 1996:2016
index3 <- summary(M3)$summary[2:22, 1]

all.df3 <- cbind(summary.df3, index3, years3) 
colnames(all.df3) <- make.names(colnames(all.df3)) # to remove % in the col names

head(all.df3)

p3 <- ggplot(all.df3, mapping = aes(x = years3)) + theme_bw() +
  geom_ribbon(aes(ymin = summary.25., ymax = summary.75.), fill = "azure3", alpha = 0.6) +
  geom_line(mapping = aes(y = summary.50.)) +
  geom_point(mapping = aes(y = index3)) +
  #scale_y_continuous(limits = c(-2, 2), labels = scales::label_number(accuracy = 0)) +
  scale_x_continuous("Years", labels = as.character(years3), breaks = years3) + 
  labs(x = "Years", y = "Index") +# theme_bw() + 
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "M3", colour = "black", size = 7, fontface = "bold") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #ggtitle("M3") + 
  theme(plot.title = element_text(size = 22, hjust = 0.5)) + 
  scale_y_continuous(limits = c(-0.5, 0.8), labels = label_number(accuracy = 0.1))
p3




summary.df4 <- as.data.frame(summary(M4,
                                     pars=c("beta_year[1]", "beta_year[2]", 
                                     "beta_year[3]", "beta_year[4]", "beta_year[5]","beta_year[6]", "beta_year[7]", 
                                     "beta_year[8]", "beta_year[9]", "beta_year[10]", "beta_year[11]", "beta_year[12]", 
                                     "beta_year[13]", "beta_year[14]", "beta_year[15]", "beta_year[16]", "beta_year[17]", 
                                     "beta_year[18]", "beta_year[19]", "beta_year[20]", "beta_year[21]")))[, c(1, 5, 6, 7)]



years4 <- 1996:2016
index4 <- summary(M4)$summary[2:22, 1]

all.df4 <- cbind(summary.df4, index4, years4) 
colnames(all.df4) <- make.names(colnames(all.df4)) # to remove % in the col names

head(all.df4)

p4 <- ggplot(all.df4, mapping = aes(x = years4)) + theme_bw() +
  geom_ribbon(aes(ymin = summary.25., ymax = summary.75.), fill = "azure3", alpha = 0.6) +
  geom_line(mapping = aes(y = summary.50.)) +
  geom_point(mapping = aes(y = index4)) +
  scale_x_continuous("Years", labels = as.character(years4), breaks = years4) + 
  labs(x = "Years", y = "") +# theme_bw() + 
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "M4", colour = "black", size = 7, fontface = "bold") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(plot.title = element_text(size = 22, hjust = 0.5))
p4



cairo_ps(filename='Fig14_p2.eps', width = 10, height = 7)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()
