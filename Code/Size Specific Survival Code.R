# set directory
rm(list = ls())
setwd("~/Post-Doc/Data/Punta Tombo Data")

library(tidyverse)
library(lubridate)
library(patchwork)
library(ggplot2)

# GOAL: explore size-specific differences for survival in penguins!

##################################################################################################################

# load data
#CRdat <- read.csv("PresenceAbsenceTOMwithagewithsizeExport_10.20.21.csv")
CRdat <- read.csv("PresenceAbsenceTOMwithagewithsizeExport_weight_10.20.21.csv") # weight included
colnames(CRdat) <- c("sex", "bookyear", "pengid", "age","bl", "bd", "fl", "f", "w")

# keep only unique rows? lots of blank ones...
CRdat <- CRdat %>% distinct(bookyear, pengid, .keep_all = T)

# extract chicks?
CRdat_chick <- CRdat %>% filter(age == "0") #48,538 chicks banded and followed...

# and now merge later data...
CRdat_chick2 <- inner_join(CRdat_chick, CRdat, by = "pengid") %>%
  dplyr::select(-bookyear.x, -age.x, -sex.x, -age.y, -sex.y,
         -bl.x, -bl.y, -bd.x, -bd.y, -fl.x, -fl.y, -f.x, -f.y, -w.x, -w.y)

CRdat2 <- CRdat_chick2 %>%
  mutate(blah = 1) %>% # adds observation for every year
  pivot_wider(names_from = bookyear.y, values_from = blah, names_sort = T) 

# change NA values to 0s
CRdat2[is.na(CRdat2)] <- 0

# now re-merge with sex data
CRdat3 <- inner_join(CRdat2, CRdat_chick, by = "pengid")

#### ADULT DATA
CRdat_adult <- CRdat %>% filter(age != "0")

# get rid of juveniles. turn into adults!
CRdat_adult$age[CRdat_adult$age == 1] <- 2

# REMOVE DUPLICATES; ID = 2281, 9070 
# keep only distinct bookyears for each pengid; removes 94 rows...
CRdat_adult <- CRdat_adult %>% distinct(bookyear, pengid, .keep_all = T)

# remove all entries that are included in CR_chick? = 19,094 penguins not banded as chicks!
CRdat_adult2 <- anti_join(CRdat_adult, CRdat_chick, by = "pengid") %>%
  mutate(blah = 1) %>%
  pivot_wider(names_from = bookyear, values_from = blah, names_sort = T)

# change NA values to 0s
CRdat_adult2[is.na(CRdat_adult2)] <- 0

# lastly, combine datasets (for now)
CRdat3_chick <- CRdat3 %>% dplyr::select(- bookyear) %>% mutate('1982' = 0)#, '2011' = 0)#; use for non-fledged
combo <- rbind(CRdat3_chick, CRdat_adult2) %>%
  dplyr::select(pengid, sex, age, bd, bl, fl, f, w, '1982','1983','1984','1985','1986','1987','1988','1989','1990',
         '1991','1992','1993','1994','1995','1996','1997','1998','1999','2000','2001','2002',
         '2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014',
         '2015','2016','2017','2018','2019')

# remove some errors
combo <- combo %>%
  filter(age != "3") %>%
  filter( age != "NULL")

##################################################################################################################

# let's do some stats!!! n = 15,871 individuals
combo_size <- combo %>%
  filter(bl != "NULL") %>%
  filter(sex != "NULL") # for now, remove the 7 unsexed individuals

combo_size$bd <- as.numeric(combo_size$bd)
combo_size$bl <- as.numeric(combo_size$bl)
combo_size$f <- as.numeric(combo_size$f)
combo_size$fl <- as.numeric(combo_size$fl)
combo_size$w <- as.numeric(combo_size$w)

# cut out outliers
combo_size <- combo_size %>%
  filter(bd < 6) %>%
  filter(bl > 3) %>%
  filter(w > 2) %>%
  filter(f > 7) # no outliers in flipper length

# graph measurements by sex!

# bill length
cols <- c("#F76D5E","#72D8FF")
ggplot(combo_size, aes(x=bl, fill = sex))+
  geom_density(alpha = 0.8, color = NA)+
  scale_fill_manual(values = cols)+
  xlim(4.5,7)+
  xlab("Bill Length (cm)")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# bill depth
ggplot(combo_size, aes(x=bd, fill = sex))+
  geom_density(alpha = 0.8, color = NA)+
  scale_fill_manual(values = cols)+
  xlim(1.5,3)+
  xlab("Bill Depth (cm)")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# foot length; outlier at foot length = 6.75
ggplot(combo_size, aes(x=f, fill = sex))+
  geom_density(alpha = 0.8, color = NA)+
  scale_fill_manual(values = cols)+
  xlim(9.25,15)+
  xlab("Foot Length (cm)")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# flipper length
ggplot(combo_size, aes(x=fl, fill = sex))+
  geom_density(alpha = 0.8, color = NA)+
  scale_fill_manual(values = cols)+
  #xlim(1.5,3)+
  xlab("Flipper Length (cm)")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# weight
ggplot(combo_size, aes(x=w, fill = sex))+
  geom_density(alpha = 0.8, color = NA)+
  scale_fill_manual(values = cols)+
  #xlim(1.5,3)+
  xlab("Weight (kg)")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# look at correlations b/t traits
cor(combo_size[,4:7])
# bill length and depth are most correlated at 0.718
# foot and flipper length aren't all that correlated at 0.602

# let's do an overall pca!!!!!
pca_overall <- prcomp(combo_size[,4:7], center = T, scale. = T) # PC1 = 71.6%, PC2 = 11.4
pca_overall$rotation # when PC1 increases, all traits are decreased!
combo_size$pc1overall <- pca_overall$x[,'PC1']
combo_size$pc2overall <- pca_overall$x[,'PC2']

# graph by sex; NOTE: flipped sign of PC1 to be consistent
ggplot(combo_size, aes(x=(pc1overall)*-1, fill = sex))+
  geom_density(alpha = 0.8, color = NA)+
  scale_fill_manual(values = cols)+
  xlim(-5,5.15)+
  xlab("PC1 - Overall")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

#PC2 is bill length/depth vs. flipper length/foot length
# NOTE: not really any difference between sexes!!!
ggplot(combo_size, aes(x=(pc2overall), fill = sex))+
  geom_density(alpha = 0.8, color = NA)+
  scale_fill_manual(values = cols)+
  xlim(-5,5.15)+
  xlab("PC2")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# let's do an locomotion pca (flipper and foot length)!!!!!
pca_loco <- prcomp(combo_size[,6:7], center = T, scale. = T) # PC1 = 80.1%, PC2 = 19.8
pca_loco$rotation # when PC1 increases, all traits are decreased!
combo_size$pc1loco <- pca_loco$x[,'PC1']

ggplot(combo_size, aes(x=(pc1loco)*-1, fill = sex))+
  geom_density(alpha = 0.8, color = NA)+
  scale_fill_manual(values = cols)+
  xlim(-5,5.15)+
  xlab("PC1 - Locomotion")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# let's do an eating pca (bill depth and bill length)!!!!!
pca_eat <- prcomp(combo_size[,4:5], center = T, scale. = T) # PC1 = 85.9%, PC2 = 14.1
pca_eat$rotation # when PC1 increases, all traits are decreased!
combo_size$pc1eat <- pca_eat$x[,'PC1']

ggplot(combo_size, aes(x=(pc1eat)*-1, fill = sex))+
  geom_density(alpha = 0.8, color = NA)+
  scale_fill_manual(values = cols)+
  xlim(-5,5.15)+
  xlab("PC1 - Eating")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# let's exclude bill depth
pca_nobd <- prcomp(combo_size[,5:7], center = T, scale. = T) # PC1 = 73.1%, PC2 = 13.7
pca_nobd$rotation # when PC1 increases, all traits are decreased!
combo_size$pc1nobd <- pca_nobd$x[,'PC1']

ggplot(combo_size, aes(x=(pc1nobd)*-1, fill = sex))+
  geom_density(alpha = 0.8, color = NA)+
  scale_fill_manual(values = cols)+
  xlim(-5,5.15)+
  xlab("PC1 - Eating")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# explore survival vs. traits!!!!
combo_size$total <- rowSums(combo_size[,8:45])
hist(combo_size$total)

#pc1 - overall
ggplot(combo_size, aes(total, pc1overall*-1))+
  geom_point()+
  facet_grid(~sex)+
  geom_smooth(method = "lm")+
  ylim(-5.5,5.5)+
  xlab("# of Years Detected")+
  ylab("PC1 - Overall")

#pc2 - overall
ggplot(combo_size, aes(total, pc2overall))+
  geom_point()+
  facet_grid(~sex)+
  geom_smooth(method = "lm")+
  ylim(-5.5,5.5)+
  xlab("# of Years Detected")+
  ylab("PC2 - Overall")

#pc1 - locomotion
ggplot(combo_size, aes(total, pc1loco*-1))+
  geom_point()+
  facet_grid(~sex)+
  geom_smooth(method = "lm")+
  ylim(-5.5,5.5)+
  xlab("# of Years Detected")+
  ylab("PC1 - Locomotion")

#pc1 - eating
ggplot(combo_size, aes(total, pc1eat*-1))+
  geom_point()+
  facet_grid(~sex)+
  geom_smooth(method = "lm")+
  ylim(-5.5,5.5)+
  xlab("# of Years Detected")+
  ylab("PC1 - Eating")

#pc1 - nobd
ggplot(combo_size, aes(total, pc1nobd*-1))+
  geom_point()+
  facet_grid(~sex)+
  geom_smooth(method = "lm")+
  ylim(-5.5,5.5)+
  xlab("# of Years Detected")+
  ylab("PC1 - Nobd")

#individual traits...
ggplot(combo_size, aes(total, bd))+
  geom_point()+
  facet_grid(~sex)+
  geom_smooth(method = "lm")+
  xlab("# of Years Detected")+
  ylab("Bill Depth")
ggplot(combo_size, aes(total, bl))+
  geom_point()+
  facet_grid(~sex)+
  geom_smooth(method = "lm")+
  xlab("# of Years Detected")+
  ylab("Bill Length")
ggplot(combo_size, aes(total, fl))+
  geom_point()+
  facet_grid(~sex)+
  geom_smooth(method = "lm")+
  xlab("# of Years Detected")+
  ylab("Flipper Length")
ggplot(combo_size, aes(total, f))+
  geom_point()+
  facet_grid(~sex)+
  geom_smooth(method = "lm")+
  xlab("# of Years Detected")+
  ylab("Foot Length")
ggplot(combo_size, aes(total, w))+
  geom_point()+
  facet_grid(~sex)+
  geom_smooth(method = "lm")+
  xlab("# of Years Detected")+
  ylab("Weight")


##################################################################################################################

#### LET's run some Bayesian code!!!
# NOTE: for now, doing really simple stuff to see which variables are most important!

library(MCMCvis)
#library(parallel)
library(nimble)
library(dplyr)

#this_cluster <- makeCluster(12)
#set.seed(12345)

# parallel function
#run_MCMC_penguin <- function(seed, combo_size){
  
  known.state.cjs <- function(ch){
    state <- ch
    for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
    }
    state[state==0] <- NA
    return(state)
  }
  
  # initial latent state fxn
  cjs.init.z <- function(ch,f){
    for (i in 1:dim(ch)[1]){
      if (sum(ch[i,])==1) next
      n2 <- max(which(ch[i,]==1))
      ch[i,f[i]:n2] <- NA
    }
    for (i in 1:dim(ch)[1]){
      ch[i,1:f[i]] <- NA
    }
    return(ch)
  }

# simple model by sex
peng.surv.sex <- nimbleCode({
  
  # constraints
  for (i in 1:nind){
    for (t in f[i]:(nyears-1)){
      logit(phi[i,t]) <- beta[sex[i]] + beta[3]*trait[i]
      p[i,t] <- p.g[sex[i]]
    }
  }
  
  # priors
  for (u in 1:g){
    beta[u] ~ dnorm(mean = 0, sd = 1.5) #dunif(0,1) # priors for sex specific survival
    p.g[u] ~ dunif(0,1) # priors for sex specific survival
  }
  beta[3] ~ dnorm(mean = 0, sd = 1.5)
  
  # likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1 
    for (t in (f[i]+1):nyears){
      # state process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    }
  }
})

# trim observations to measure stuff
combo_size2 <- combo_size
#combo_size2 <- combo_size[1:1000,]

# capture histories
CH <- combo_size2 %>%
  dplyr::select(-pengid, -sex, -age, -bd, -bl, -f, -fl, -w, -pc1overall, -pc2overall, -pc1loco, -pc1eat, -pc1nobd, -total) %>%
  as.matrix(.)

# my dat
my.data <- list(y = CH, z = known.state.cjs(CH))
  
# first obs
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)
                 
# my constants
my.constants <- list(nind = dim(CH)[1],
                     nyears = dim(CH)[2],
                     g = length(unique(combo_size2$sex)),
                     f = f,
                     sex = if_else(combo_size2$sex == "M", 1, 2),
                     trait = as.vector(combo_size2$w))
                     #trait = as.vector(scale(combo_size2$bd)))

# inits
initial.values <- function(){list(
  z = cjs.init.z(CH,f),
  beta = rnorm(3,0,5),
  p.g = runif(length(unique(combo_size2$sex)),0,1))}

# run model
peng.surv.sex.out <- nimbleMCMC(code = peng.surv.sex ,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values,
                             monitors = c("beta", "p.g"),
                             niter = 5000,
                             nburnin = 1000,
                             nchains = 2)
#return(peng.surv.sex.out)
#}

#peng.surv.sex.out <- parLapply(cl = this_cluster, X = 1:2, # X should be the number of chains to run...
#                            fun = run_MCMC_penguin,
 #                           combo_size = combo_size)
#stopCluster(this_cluster)


MCMCsummary(peng.surv.sex.out, round = 2)
MCMCtrace(peng.surv.sex.out, pdf = F)

# extract data and calculate means/95% CIs to make graphs
beta1 <- c(peng.surv.sex.out$chain1[,'beta[1]'], peng.surv.sex.out$chain2[,'beta[1]'])
beta2 <- c(peng.surv.sex.out$chain1[,'beta[2]'], peng.surv.sex.out$chain2[,'beta[2]'])
beta3 <- c(peng.surv.sex.out$chain1[,'beta[3]'], peng.surv.sex.out$chain2[,'beta[3]'])
#beta <- data.frame(beta1=beta1,beta2=beta2,beta3=beta3)
#write.csv(beta, "beta_results_11.1.21.csv")
predicted_survivalM <- matrix(NA, nrow = length(beta1), ncol = length(my.constants$trait))
predicted_survivalF <- matrix(NA, nrow = length(beta1), ncol = length(my.constants$trait))
for (i in 1:length(beta1)){
  for (j in 1:length(my.constants$trait)){
    predicted_survivalM[i,j] <- plogis(beta1[i] + beta3[i] * my.constants$trait[j]) # males
    predicted_survivalF[i,j] <- plogis(beta2[i] + beta3[i] * my.constants$trait[j]) # females
  }
}
mean_survivalM <- apply(predicted_survivalM, 2, mean)
lciM <- apply(predicted_survivalM, 2, quantile, prob = 2.5/100)
uciM <- apply(predicted_survivalM, 2, quantile, prob = 97.5/100)
mean_survivalF <- apply(predicted_survivalF, 2, mean)
lciF <- apply(predicted_survivalF, 2, quantile, prob = 2.5/100)
uciF <- apply(predicted_survivalF, 2, quantile, prob = 97.5/100)
ord <- order(my.constants$trait)
df <- data.frame(trait = c(my.constants$trait[ord], my.constants$trait[ord]),
                 survival = c(mean_survivalM[ord], mean_survivalF[ord]),
                 lci = c(lciM[ord],lciF[ord]),
                 uci = c(uciM[ord],uciF[ord]),
                 sex = c(rep("male", length(mean_survivalM)), rep("female", length(mean_survivalF))))

# graph!
df %>%
  ggplot() + 
  aes(x = trait, y = survival, color = sex) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = sex), alpha = 0.5) + 
  ylim(0.78,1) + 
  xlim(-5,5)+
  labs(x = "PC1 - Overall", y = "estimated survival", color = "", fill = "")
ggsave("pc1overall_11.1.21.png",dpi=600,width=10)











