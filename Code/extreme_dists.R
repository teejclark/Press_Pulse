# set directory
rm(list = ls())
setwd("~/Post-Doc/Data/Punta Tombo Data")

library(tidyverse)
library(fitdistrplus)
library(ismev)
library(evd)
library(sensitivity)
library(lubridate)

# FIT distributions to weather data...just to explore
# NOTE: might make more sense to use raw data instead of the yearly means...

##################################################################################################################

# load weather/ocean data as a test

weather <- read.csv("Site Weather.csv")

# replace NAs in precipitation
weather1 <- weather %>% 
  mutate(Precipitation = replace(Precipitation, which(is.na(Precipitation)), 0.0))

# convert character dates to actual dates
weather1$WeathDate <- as.Date(weather1$WeathDate)

# calculate total precipitation between Oct 15 and Dec 15
rain_60 <-  weather1 %>%
  filter(WeathDate >= as.Date(paste(year(WeathDate), 10, 15, sep = "-")),
         WeathDate <= as.Date(paste(year(WeathDate), 12, 15, sep = "-"))) %>%
  group_by(BookYear) %>%
  filter(!is.na(Precipitation)) %>%
  summarize(amt_rain = sum(Precipitation))

# TEMP

# number of days with max temp > 25C per breeding season
num_days_25 <- weather1 %>% group_by(BookYear) %>%
  filter(!is.na(MaxTemp)) %>%
  summarize(n_days = n(),
            n_gt25 = sum(MaxTemp > 25),
            p_gt25 = n_gt25/n_days)

# OCEAN - MIGRATION - SSTA
ssta_migration <- read.csv("ssta.PC1.migration.11.12.21.csv")

# OCEAN - BREEDING - SSTA
ssta_breeding <- read.csv("ssta_summarized.12.2.21.csv")
#ssta_breeding <- read.csv("ssta.PC1.breeding.11.15.21.csv")

##################################################################################################################

#### Fit Normal Curves to Data ####

# temp
hist(num_days_25$p_gt25)
fit <- fitdist(num_days_25$p_gt25, "norm")
plot(fit) # mean is 0.383, sd = 0.111
gofstat(fit)

# ppt = clearly non-normal, for now working with a gamma distribution
hist(rain_60$amt_rain)
fit <- fitdist(rain_60$amt_rain, "gamma")
plot(fit)
gofstat(fit)

# plume
hist(ssta_migration$mean.2)
fit <- fitdist(ssta_migration$mean.2, "norm")
plot(fit)
gofstat(fit)

# ssta-breeding
hist(ssta_breeding$av)
fit <- fitdist(ssta_breeding$av, "norm")
plot(fit)
gofstat(fit)

##################################################################################################################

# Calculate Autocorrelation?

acf(num_days_25$p_gt25) # not statistically significant
acf(rain_60$amt_rain) # not statistically significant
acf(ssta_migration$mean.2) # not statistically significant
acf(ssta_breeding$av) # not statistically significant

##################################################################################################################

#### Fit GEV distribution to Data ####

# ppt
fit <- gev.fit(rain_60$amt_rain)
gev.diag(fit)

# plume
fit <- gev.fit(num_days_25$p_gt25)
gev.diag(fit)

##################################################################################################################

# Run Normal Dist. Sensitivity Analyses, see how Extreme Values are Affected by distributions
# frequency of ECE = >=95% percentile; magnitude of ECE = most extreme temp

# TEMPERATURE
# fit normal distribution to get arguments for sensitivity analysis = just testing here
fit <- fitdist(num_days_25$p_gt25, "norm"); fit # mean = 0.383, sd = 0.111
vals <- rnorm(10000, fit$estimate[1], fit$estimate[2])
ECE <- quantile(vals, 0.95) # 0.566
freqECE <- sum(vals >= ECE)/10000 # baseline frequency of warm ECEs = 5%
magECE <- max(vals) # baseline magnitude of warm ECEs = 0.816

# calculate original values for extremes
ECE <- rep(NA,10000); freqECEval <- rep(NA, 10000); magECEval <- rep(NA,10000)
for (i in 1:length(ECE)){
  vals <- rnorm(10000, fit$estimate[1], fit$estimate[2])
  ECE[i] <- quantile(vals, 0.95)
  freqECEval[i] <- sum(vals >= ECE[i])/10000
  magECEval[i] <- max(vals)
}
ECE <- mean(ECE); ECE #0.5676679
freqECEval <- mean(freqECEval); freqECEval #5%
magECEval <- mean(magECEval); magECEval #0.8148279

# local sensitivity analysis
nsims <- 100
#deviation <- seq(-0.001,0.001,0.001) # for now, very small deviations
#deviation <- seq(-0.5,0.5,0.01)
#mean <- fit$estimate[1]*(1+deviation)
#sd <- fit$estimate[2]*(1+deviation)
mean <- rnorm(100, fit$estimate[1], fit$sd[1])
sd <- rnorm(100, fit$estimate[2], fit$sd[2])
freqECE <- matrix(NA, length(mean)^2, nsims+2)
magECE <- matrix(NA, length(mean)^2, nsims+2)
freqECE[,1] <- expand.grid(mean,sd)[,1]; freqECE[,2] <- expand.grid(mean,sd)[,2]
magECE[,1] <- expand.grid(mean,sd)[,1]; magECE[,2] <- expand.grid(mean,sd)[,2]
colnames(freqECE) <- c("mean", "sd", rep("val", nsims))
colnames(magECE) <- c("mean", "sd", rep("val", nsims))

for (i in 1:nrow(freqECE)){
  for (j in 1:nsims){
    vals <- rnorm(10000, freqECE[i,1], freqECE[i,2])
    freqECE[i,j+2] <- sum(vals >= ECE)/10000 # ECE is 0.5676679
    magECE[i,j+2] <- max(vals)
  }}

# run PCC
freqECE_new <- data.frame(mean = rep(freqECE[,1],nsims),
                          sd = rep(freqECE[,2],nsims),
                          val = data.table::melt(freqECE[,3:ncol(freqECE)])[,3])
magECE_new <- data.frame(mean = rep(magECE[,1],nsims),
                         sd = rep(magECE[,2],nsims),
                         val = data.table::melt(magECE[,3:ncol(magECE)])[,3])

freqECE_pcc <- pcc(X = freqECE_new[,1:2], y = freqECE_new[,3], conf = 0.95, nboot = 100)

magECE_pcc <- pcc(X = magECE_new[,1:2], y = magECE_new[,3], conf = 0.95, nboot = 100)

ECE_frame_a <- data.frame(rbind(freqECE_pcc$PCC, magECE_pcc$PCC),
                        type = c("Freq", "Freq", "Mag", "Mag"),
                        meas = c("Mean", "SD", "Mean", "SD"))
ECE_frame_a <- rownames_to_column(ECE_frame)

ggplot(ECE_frame_a, aes(x=meas, y = original, fill = type))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin = min..c.i., ymax = max..c.i.), width = 0.2, position = position_dodge(0.9))+
  xlab("Parameters")+
  ylab("Partial Correlation Coefficients (PCC)")+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

##############

# PLUME
# NOTE: probably need to recalculate plume...
# fit normal distribution to get arguments for sensitivity analysis = just testing here
fit <- fitdist(ssta_migration$mean.2, "norm"); fit # mean = 0.383, sd = 0.111
vals <- rnorm(10000, fit$estimate[1], fit$estimate[2])
ECE <- quantile(vals, 0.95) # 1.29
freqECE <- sum(vals >= ECE)/10000 # baseline frequency of warm ECEs = 5%
magECE <- max(vals) # baseline magnitude of warm ECEs = 3.2

# calculate original values for extremes
ECE <- rep(NA,10000); freqECEval <- rep(NA, 10000); magECEval <- rep(NA,10000)
for (i in 1:length(ECE)){
  vals <- rnorm(10000, fit$estimate[1], fit$estimate[2])
  ECE[i] <- quantile(vals, 0.95)
  freqECEval[i] <- sum(vals >= ECE[i])/10000
  magECEval[i] <- max(vals)
}
ECE <- mean(ECE); ECE #1.29
freqECEval <- mean(freqECEval); freqECEval #5%
magECEval <- mean(magECEval); magECEval #3.02

# local sensitivity analysis
nsims <- 100
#deviation <- seq(-0.001,0.001,0.001) # for now, very small deviations
#deviation <- seq(-0.5,0.5,0.01)
#mean <- fit$estimate[1]*(1+deviation)
#sd <- fit$estimate[2]*(1+deviation)
mean <- rnorm(100, fit$estimate[1], fit$sd[1])
sd <- rnorm(100, fit$estimate[2], fit$sd[2])
freqECE <- matrix(NA, length(mean)^2, nsims+2)
magECE <- matrix(NA, length(mean)^2, nsims+2)
freqECE[,1] <- expand.grid(mean,sd)[,1]; freqECE[,2] <- expand.grid(mean,sd)[,2]
magECE[,1] <- expand.grid(mean,sd)[,1]; magECE[,2] <- expand.grid(mean,sd)[,2]
colnames(freqECE) <- c("mean", "sd", rep("val", nsims))
colnames(magECE) <- c("mean", "sd", rep("val", nsims))

for (i in 1:nrow(freqECE)){
  for (j in 1:nsims){
    vals <- rnorm(10000, freqECE[i,1], freqECE[i,2])
    freqECE[i,j+2] <- sum(vals >= ECE)/10000 # ECE is 0.5676679
    magECE[i,j+2] <- max(vals)
  }}

# run PCC
freqECE_new <- data.frame(mean = rep(freqECE[,1],nsims),
                          sd = rep(freqECE[,2],nsims),
                          val = data.table::melt(freqECE[,3:ncol(freqECE)])[,3])
magECE_new <- data.frame(mean = rep(magECE[,1],nsims),
                         sd = rep(magECE[,2],nsims),
                         val = data.table::melt(magECE[,3:ncol(magECE)])[,3])

freqECE_pcc <- pcc(X = freqECE_new[,1:2], y = freqECE_new[,3], conf = 0.95, nboot = 100)

magECE_pcc <- pcc(X = magECE_new[,1:2], y = magECE_new[,3], conf = 0.95, nboot = 100)

ECE_frame_b <- data.frame(rbind(freqECE_pcc$PCC, magECE_pcc$PCC),
                          type = c("Freq", "Freq", "Mag", "Mag"),
                          meas = c("Mean", "SD", "Mean", "SD"))
ECE_frame_b <- rownames_to_column(ECE_frame)

ggplot(ECE_frame_b, aes(x=meas, y = original, fill = type))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin = min..c.i., ymax = max..c.i.), width = 0.2, position = position_dodge(0.9))+
  xlab("Parameters")+
  ylab("Partial Correlation Coefficients (PCC)")+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

##############

# SSTA - BREEDING
# fit normal distribution to get arguments for sensitivity analysis = just testing here
fit <- fitdist(ssta_breeding$av, "norm"); fit # mean = 0.383, sd = 0.111
vals <- rnorm(10000, fit$estimate[1], fit$estimate[2])
ECE <- quantile(vals, 0.95) # 0.777
freqECE <- sum(vals >= ECE)/10000 # baseline frequency of warm ECEs = 5%
magECE <- max(vals) # baseline magnitude of warm ECEs = 1.84

# calculate original values for extremes
ECE <- rep(NA,10000); freqECEval <- rep(NA, 10000); magECEval <- rep(NA,10000)
for (i in 1:length(ECE)){
  vals <- rnorm(10000, fit$estimate[1], fit$estimate[2])
  ECE[i] <- quantile(vals, 0.95)
  freqECEval[i] <- sum(vals >= ECE[i])/10000
  magECEval[i] <- max(vals)
}
ECE <- mean(ECE); ECE #0.77
freqECEval <- mean(freqECEval); freqECEval #5%
magECEval <- mean(magECEval); magECEval #1.829

# local sensitivity analysis
nsims <- 100
#deviation <- seq(-0.001,0.001,0.001) # for now, very small deviations
#deviation <- seq(-0.5,0.5,0.01)
#mean <- fit$estimate[1]*(1+deviation)
#sd <- fit$estimate[2]*(1+deviation)
mean <- rnorm(100, fit$estimate[1], fit$sd[1])
sd <- rnorm(100, fit$estimate[2], fit$sd[2])
freqECE <- matrix(NA, length(mean)^2, nsims+2)
magECE <- matrix(NA, length(mean)^2, nsims+2)
freqECE[,1] <- expand.grid(mean,sd)[,1]; freqECE[,2] <- expand.grid(mean,sd)[,2]
magECE[,1] <- expand.grid(mean,sd)[,1]; magECE[,2] <- expand.grid(mean,sd)[,2]
colnames(freqECE) <- c("mean", "sd", rep("val", nsims))
colnames(magECE) <- c("mean", "sd", rep("val", nsims))

for (i in 1:nrow(freqECE)){
  for (j in 1:nsims){
    vals <- rnorm(10000, freqECE[i,1], freqECE[i,2])
    freqECE[i,j+2] <- sum(vals >= ECE)/10000 
    magECE[i,j+2] <- max(vals)
  }}

# run PCC
freqECE_new <- data.frame(mean = rep(freqECE[,1],nsims),
                          sd = rep(freqECE[,2],nsims),
                          val = data.table::melt(freqECE[,3:ncol(freqECE)])[,3])
magECE_new <- data.frame(mean = rep(magECE[,1],nsims),
                         sd = rep(magECE[,2],nsims),
                         val = data.table::melt(magECE[,3:ncol(magECE)])[,3])

freqECE_pcc <- pcc(X = freqECE_new[,1:2], y = freqECE_new[,3], conf = 0.95, nboot = 100)

magECE_pcc <- pcc(X = magECE_new[,1:2], y = magECE_new[,3], conf = 0.95, nboot = 100)

ECE_frame_c <- data.frame(rbind(freqECE_pcc$PCC, magECE_pcc$PCC),
                          type = c("Freq", "Freq", "Mag", "Mag"),
                          meas = c("Mean", "SD", "Mean", "SD"))
ECE_frame_c <- rownames_to_column(ECE_frame)

a <- ggplot(ECE_frame_c, aes(x=meas, y = original, fill = type))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin = min..c.i., ymax = max..c.i.), width = 0.2, position = position_dodge(0.9))+
  xlab("Parameters")+
  ylab("Partial Correlation Coefficients (PCC)")+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.position = "none")

##############

# RAIN
# fit gamma distribution to get arguments for sensitivity analysis = just testing here
fit <- fitdist(rain_60$amt_rain, "gamma"); fit 
vals <- rgamma(10000, fit$estimate[1], fit$estimate[2])
ECE <- quantile(vals, 0.95) # 94.73
freqECE <- sum(vals >= ECE)/10000 # baseline frequency of warm ECEs = 5%
magECE <- max(vals) # baseline magnitude of warm ECEs = 417.75

# calculate original values for extremes
ECE <- rep(NA,10000); freqECEval <- rep(NA, 10000); magECEval <- rep(NA,10000)
for (i in 1:length(ECE)){
  vals <- rgamma(10000, fit$estimate[1], fit$estimate[2])
  ECE[i] <- quantile(vals, 0.95)
  freqECEval[i] <- sum(vals >= ECE[i])/10000
  magECEval[i] <- max(vals)
}
ECE <- mean(ECE); ECE #93.8
freqECEval <- mean(freqECEval); freqECEval #5%
magECEval <- mean(magECEval); magECEval #263.9

# local sensitivity analysis
nsims <- 100
shape <- rnorm(100, fit$estimate[1], fit$sd[1])
rate <- rnorm(100, fit$estimate[2], fit$sd[2])
freqECE <- matrix(NA, length(shape)^2, nsims+2)
magECE <- matrix(NA, length(shape)^2, nsims+2)
freqECE[,1] <- expand.grid(shape,rate)[,1]; freqECE[,2] <- expand.grid(shape,rate)[,2]
magECE[,1] <- expand.grid(shape,rate)[,1]; magECE[,2] <- expand.grid(shape,rate)[,2]
colnames(freqECE) <- c("shape", "rate", rep("val", nsims))
colnames(magECE) <- c("shape", "rate", rep("val", nsims))

for (i in 1:nrow(freqECE)){
  for (j in 1:nsims){
    vals <- rgamma(10000, freqECE[i,1], freqECE[i,2])
    freqECE[i,j+2] <- sum(vals >= ECE)/10000 
    magECE[i,j+2] <- max(vals)
  }}

# run PCC
freqECE_new <- data.frame(shape = rep(freqECE[,1],nsims),
                          scale = 1/rep(freqECE[,2],nsims),
                          val = data.table::melt(freqECE[,3:ncol(freqECE)])[,3])
magECE_new <- data.frame(shape = rep(magECE[,1],nsims),
                         scale = 1/rep(magECE[,2],nsims),
                         val = data.table::melt(magECE[,3:ncol(magECE)])[,3])

freqECE_pcc <- pcc(X = freqECE_new[,1:2], y = freqECE_new[,3], conf = 0.95, nboot = 100)

magECE_pcc <- pcc(X = magECE_new[,1:2], y = magECE_new[,3], conf = 0.95, nboot = 100)

ECE_frame_d <- data.frame(rbind(freqECE_pcc$PCC, magECE_pcc$PCC),
                          type = c("Freq", "Freq", "Mag", "Mag"),
                          meas = c("Shape", "Scale", "Shape", "Scale"))
ECE_frame_d <- rownames_to_column(ECE_frame)

b <- ggplot(ECE_frame_d, aes(x=meas, y = original, fill = type))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin = min..c.i., ymax = max..c.i.), width = 0.2, position = position_dodge(0.9))+
  xlab("Parameters")+
  ylab("Partial Correlation Coefficients (PCC)")+
  ylim(0,1)+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

a+b

ggsave("extreme_sensitivity_12.6.21.png", dpi = 600, width = 13)
