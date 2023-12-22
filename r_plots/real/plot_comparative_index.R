rm(list = ls())
setwd("C:/Users/Usuario/Desktop/semipar_TPS/tmbstan_real_model")

library(loo)
library(dplyr)
library(ggplot2)
library(rstan)
set.seed(12345)


rs_fit_spde  = readRDS("fit_spde_tmbstan.RDS")
summary(rs_fit_spde)

rs_fit_tps  = readRDS("fit_tps_tmbstan.RDS")
rs_fit_tps


#rs_fit_1_draws <- posterior_tps
# rs_fit_1_draws <- rstan::extract(rs_fit_1)

# summary.df <- as.data.frame(summary(rs_fit_1,pars=c("beta_year[1]", "beta_year[2]", 
#               "beta_year[3]", "beta_year[4]", "beta_year[5]","beta_year[6]", "beta_year[7]", 
#               "beta_year[8]", "beta_year[9]", "beta_year[10]", "beta_year[11]", "beta_year[12]", 
#               "beta_year[13]", "beta_year[14]", "beta_year[15]", "beta_year[16]", "beta_year[17]", 
#               "beta_year[18]", "beta_year[19]", "beta_year[20]", "beta_year[21]"), probs=c(.05,.5,.95))$summary)

summary.df <- as.data.frame(summary(rs_fit_spde,pars=c("beta_year[1]", "beta_year[2]", 
                            "beta_year[3]", "beta_year[4]", "beta_year[5]","beta_year[6]", "beta_year[7]", 
                            "beta_year[8]", "beta_year[9]", "beta_year[10]", "beta_year[11]", "beta_year[12]", 
                            "beta_year[13]", "beta_year[14]", "beta_year[15]", "beta_year[16]", "beta_year[17]", 
                            "beta_year[18]", "beta_year[19]", "beta_year[20]", "beta_year[21]")))[, c(1, 5, 6, 7)]



years <- 1996:2016
index <- summary(rs_fit_spde)$summary[2:22, 1]

all.df <- cbind(summary.df, index, years) 
colnames(all.df) <- make.names(colnames(all.df)) # to remove % in the col names

head(all.df)

p1 <- ggplot(all.df, mapping = aes(x = years)) +
  geom_ribbon(aes(ymin = summary.25., ymax = summary.75.), fill = "azure3", alpha = 0.6) +
  geom_line(mapping = aes(y = summary.50.)) +
  geom_point(mapping = aes(y = index)) +
  scale_x_continuous("Years", labels = as.character(years), breaks = years) + 
  # scale_x_discrete(labels = c("1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", 
  #                             "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) + 
  labs(x = "Years", y = "Index") + # theme_bw() + 
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("Gamma ST model") + 
  theme(plot.title = element_text(size = 22, hjust = 0.5))
p1 


#rs_fit_1_draws <- posterior_tps
# rs_fit_1_draws <- rstan::extract(rs_fit_1)

# summary.df <- as.data.frame(summary(rs_fit_1,pars=c("beta_year[1]", "beta_year[2]", 
#               "beta_year[3]", "beta_year[4]", "beta_year[5]","beta_year[6]", "beta_year[7]", 
#               "beta_year[8]", "beta_year[9]", "beta_year[10]", "beta_year[11]", "beta_year[12]", 
#               "beta_year[13]", "beta_year[14]", "beta_year[15]", "beta_year[16]", "beta_year[17]", 
#               "beta_year[18]", "beta_year[19]", "beta_year[20]", "beta_year[21]"), probs=c(.05,.5,.95))$summary)

summary.df2 <- as.data.frame(summary(rs_fit_tps,pars=c("beta_year[1]", "beta_year[2]", 
              "beta_year[3]", "beta_year[4]", "beta_year[5]","beta_year[6]", "beta_year[7]", 
              "beta_year[8]", "beta_year[9]", "beta_year[10]", "beta_year[11]", "beta_year[12]", 
              "beta_year[13]", "beta_year[14]", "beta_year[15]", "beta_year[16]", "beta_year[17]", 
              "beta_year[18]", "beta_year[19]", "beta_year[20]", "beta_year[21]")))[, c(1, 5, 6, 7)]



years2 <- 1996:2016
index2 <- summary(rs_fit_tps)$summary[2:22, 1]

all.df2 <- cbind(summary.df2, index2, years2) 
colnames(all.df2) <- make.names(colnames(all.df2)) # to remove % in the col names

head(all.df2)

p2 <- ggplot(all.df2, mapping = aes(x = years2)) +
  geom_ribbon(aes(ymin = summary.25., ymax = summary.75.), fill = "azure3", alpha = 0.6) +
  geom_line(mapping = aes(y = summary.50.)) +
  geom_point(mapping = aes(y = index2)) +
  scale_x_continuous("Years", labels = as.character(years2), breaks = years2) + 
  # scale_x_discrete(labels = c("1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", 
  #                             "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")) + 
  labs(x = "Years", y = "") +# theme_bw() + 
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("Skew normal ST model") + 
  theme(plot.title = element_text(size = 22, hjust = 0.5))
p2


library(gridExtra)
grid.arrange(p1, p2, ncol = 2)
