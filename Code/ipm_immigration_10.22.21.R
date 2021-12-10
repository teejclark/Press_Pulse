# set directory
rm(list = ls())
setwd("~/Post-Doc/Data/Punta Tombo Data")

library(tidyverse)
library(nimble)
library(ggplot2)
library(MCMCvis)
library(lubridate)

# GOAL: BUILD an IPM and add Environmental/Weather Covariates.
# NOTE: also consider the role of immigration in this system...

##################################################################################################################

#### LOAD POPULATION DATA ####

# REPRODUCTION: NUMNESTS = number of nests, numfledged = number fledged!
repro <- read.csv("ReproSuccess_9.21.21.csv", header = F)
repro$V1[1] <- "1984" # fixed!
colnames(repro) <- c("bookyear", "numnests", "numeggs", "clutchsize", "numhatch", "meanhatch", "numfledged", "rs")
repro$bookyear <- as.Date(as.character(repro$bookyear), "%Y")
repro <- repro %>% arrange(bookyear)
# missing 1982, 2011

# number of surveyed nests
R <- as.numeric(c(NA, repro$numnests[1:28], NA, repro$numnests[29:36]))

# number of fledges
J <- as.numeric(c(NA, repro$numfledged[1:28], NA, repro$numfledged[29:36]))

# POPULATION DATA (using Ginger's normal adjustment, not sliding scale...)
# NOTE: I have 2020 pop data, going to exclude for now...
pop <- read.csv("StakeSurveyTOMACT22_9.21.21.csv", header = F)
pop[1,1] <- 1987 # fix weird stuff
colnames(pop) <- c("bookyear","00N01E", "00N02E", "00N03E", "00N04E", "00N05E", "00N06E", "00N07E", "00N08E",
                   "00N09E", "00N12E", "00N14E", "00N15E", "00N16E", "00N17E", "00N18E", "01N14E", "01S14E",
                   "02N14E", "02S14E", "03S14E", "04S14E", "05S14E", "total")
pop$density <- pop$total/22 # number of 100 meter survey areas = 22

# convert by 1.36
pop$blah <- c(pop$density[1:28]/1.36, pop$density[29:33])
pop$popest2 <- pop$blah*35240 # total habitat reported in the paper...(Table 1, 2012 survey)

y <- c(rep(NA, 5), pop$popest2[1:24], NA, pop$popest2[25:32])

# SURVIVAL DATA - see "CRData.R" to load data
# NOTE: only using verified adult female survival, use full mix of juveniles...
# NOTE: for now, using random sex for the unsexed adults!

# drop banded adults never sexed (we assume all juveniles are not sexed)
# drop banded juveniles that were recaptured and never sexed
combo$pengid <- as.character(combo$pengid)
combo$age <- as.character(combo$age)
combo_b <- combo %>%
  mutate(sum = rowSums(across(where(is.numeric)))) # sum across capture histories

combo2 <- combo_b %>%
  #filter(age == 2, sex != "NULL") # get rid of unsexed adults
  filter(age == 2, sex == "NULL") %>%
  mutate(sex = sample(c("M","F"), n(), prob = c(0.5,0.5), replace = T))

combo2.5 <- combo_b %>%
  filter(age == 2, sex != "NULL")

combo3 <- combo_b %>%
  filter(age == 0) %>%
  bind_rows(combo2) %>%
  bind_rows(combo2.5)

# randomly assign sex to recaptured, unsexed individuals (n = 784)
combo4 <- combo3 %>%
  #filter(age == 0, sum > 1, sex != "NULL") # recaptured, sexed individuals (excludes unsexed recaptured juvs)
  filter(age == 0, sum > 1, sex == "NULL") %>%
  mutate(sex = sample(c("M","F"), n(), prob = c(0.5,0.5), replace = T))

# retain unrecaptured, unsexed individuals
combo4.5 <- combo3 %>%
  filter(age == 0, sum <= 1, sex == "NULL")

# retain unsexed individuals
combo5 <- combo3 %>%
  filter(age == 0, sex != "NULL")

combo6 <- combo3 %>%
  filter(age == 2) %>%
  bind_rows(combo4) %>%
  bind_rows(combo4.5) %>%
  bind_rows(combo5)

# individuals marked as juveniles
CH.J <- combo6 %>% filter(age == 0) %>%
  dplyr::select(-pengid, -sex, -age, -sum) %>%
  as.matrix(.)

# individuals marked as adults, sexed as females
CH.A <- combo6 %>% filter(age == 2) %>%
  filter(sex == "F") %>%
  dplyr::select(-pengid, -sex, -age, -sum) %>%
  as.matrix(.)

# reformat into pure numerics
CH.J2 <- matrix(as.numeric(CH.J), nrow(CH.J), ncol(CH.J))
CH.A2 <- matrix(as.numeric(CH.A), nrow(CH.A), ncol(CH.A))

# set-up
cap <- apply(CH.J2, 1, sum) # total # of juvenile captures
ind <- which(cap >=2) # which inds were recaptured; only 3321
CH.J.R <- CH.J2[ind,] # juveniles recaptured at least once
CH.J.N <- CH.J2[-ind,] # juveniles not recaptured

# remove first capture of re-captured individuals (b/c they are adults with re-capture)
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
  first[i] <- min(which(CH.J.R[i,]==1))
}
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R1[i,first[i]] <- 0
}

# add grown-up juveniles to adults and create m-array
CH.A.m <- rbind(CH.A2, CH.J.R1)
CH.A.marray <- marray(CH.A.m)

# create CH matrix for juveniles, ignoring subsequent recaptures (those who were recaptured)
second <- numeric()
for (i in 1:dim(CH.J.R1)[1]){
  second[i] <- min(which(CH.J.R1[i,]==1))
}
CH.J.R2 <- matrix(0, nrow = dim(CH.J.R)[1], ncol = dim(CH.J.R)[2])
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R2[i,first[i]] <- 1
  CH.J.R2[i,second[i]] <- 1
}

# create m-array for these
CH.J.R.marray <- marray(CH.J.R2)

# the last column ought to show # of juveniles not recaptured again,
# and should all be zeros, since all were released as adults
CH.J.R.marray[,dim(CH.J2)[2]] <- 0

# create m-array for juvs never recaptured and add to previous array
CH.J.N.marray <- marray(CH.J.N)
CH.J.marray <- CH.J.R.marray + CH.J.N.marray

##################################################################################################################

#### LOAD ENVIRONMENTAL DATA ####

# PRECIPITATION

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
#ssta_breeding <- read.csv("ssta.PC1.breeding.11.15.21.csv") # for now, not using the PCA
ssta_breeding <- read.csv("ssta_summarized.12.2.21.csv")

# OCEAN - MIGRATION - CHLORA
# OCEAN - BREEDING - CHLORA
# NOTE: for now, omitting the chlorophyll data b/c it's unclear what to do here.

##################################################################################################################

#### Simple Immigration Model ####
# build a pre-breeding, two-age model for the penguin population! with immigration added!

peng.ipm.imm.1 <- nimbleCode({
  
  ####
  # 1. Define priors for the parameters - very vague priors for everything...
  ####
  
  # initial population sizes
  N1[1] ~ T(dnorm(25000, 4e-08),0,) # from 6500
  Nad[1] ~ T(dnorm(240000, 5e-08),0,) # from 200000
  Nimm[1] ~ T(dnorm(20000, 4e-08),0,) # just guessing here
  
  # mean demographic parameters 
  l.mphij ~ dnorm(0, 0.34)
  l.mphia ~ dnorm(0, 0.34)
  l.mfec ~ dnorm(0, 0.34)
  l.mim ~ dnorm(0, 0.34)
  l.pj ~ dnorm(0, 0.34) 
  l.pa ~ dnorm(0, 0.34)
  
  # observation error
  tauy <- pow(sigma.y, -2) # (precision)
  sigma.y ~ dunif(0, 50000) # (sd) changed from 50000
  
  ####
  # 2. Constrain parameters
  ####
  
  for (t in 1:(nyears-1)){
    logit(phi.juv[t]) <- l.mphij
    logit(phi.ad[t]) <-  l.mphia
    log(fec[t]) <-  l.mfec
    log(omega[t]) <-  l.mim
    #logit(pj[t]) <- l.pj
    #logit(pa[t]) <- l.pa
  }
  
  # detection constraints - juveniles
  for (t in 1:29){logit(pj[t]) <- l.pj}
  pj[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pj[t]) <- l.pj}
  
  # detection constraints - adults
  for (t in 1:29){logit(pa[t]) <- l.pa}
  pa[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pa[t]) <- l.pa}
  
  ####
  # 3. Derived parameters
  ####
  
  mphij <- exp(l.mphij) / (1+exp(l.mphij))
  mphia <- exp(l.mphia) / (1+exp(l.mphia))
  mfec <- exp(l.mfec)
  mim <- exp(l.mim)
  mpj <- exp(l.pj) / (1+exp(l.pj))
  mpa <- exp(l.pa) / (1+exp(l.pa))
  
  # population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1]/Ntot[t]
    l.lambda[t] <- log(lambda[t])
  }
  
  # geometric mean
  geomean.lambda <- exp((1/(nyears-1))*sum(l.lambda[1:(nyears-1)]))
  
  ####
  # 4. Likelihoods
  ####
  
  # 4.1 Likelihoods for population count data
  for (t in 2:nyears){
    mean1[t] <- 0.5 * fec[t-1] * phi.juv[t-1] * Nad[t-1]
    N1[t] ~ dpois(mean1[t])
    Nad[t] ~ dbin(phi.ad[t-1], Ntot[t-1])
    mpo[t] <- Ntot[t-1] * omega[t-1]
    Nimm[t] ~ dpois(mpo[t])
  }
  for (t in 1:nyears){
    Ntot[t] <- N1[t] + Nad[t] + Nimm[t]
  }
  
  # observation process
  for (t in 1:nyears){
    y[t] ~ T(dnorm(Ntot[t], tauy),0,)
  }
  
  # 4.2 Likelihoods for capture-recapture data
  # define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr.j[t,1:nyears] ~ dmulti(pr.j[t,1:nyears], rel.j[t])
    marr.a[t,1:nyears] ~ dmulti(pr.a[t,1:nyears], rel.a[t])
  }
  
  # m-array cell probabilities for juveniles and adults
  for (t in 1:(nyears-1)){
    qa[t] <- 1 - pa[t] # probability of recapture, adults
    qj[t] <- 1 - pj[t] # probability of recapture, juveniles
    
    # main diagonal
    pr.j[t,t] <- phi.juv[t] * pj[t]
    pr.a[t,t] <- phi.ad[t] * pa[t]
    
    # above main diagonal
    for (j in (t+1):(nyears-1)){
      pr.j[t,j] <- phi.juv[t] * prod(phi.ad[(t+1):j]) * qj[t] * prod(qa[t:(j-1)]) * pa[j]/qa[t]
      pr.a[t,j] <- prod(phi.ad[t:j]) * prod(qa[t:(j-1)]) * pa[j]
    }
    
    # below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  
  for (t in 1:(nyears-1)){
    # last column: probability of non-recapture
    pr.j[t,nyears] <- 1 - sum(pr.j[t, 1:(nyears-1)])
    pr.a[t,nyears] <- 1 - sum(pr.a[t, 1:(nyears-1)])
  }
  
  # 4.3 Likelihoods for productivity data
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*fec[t]
  }
})

# my data

# fucking with data
jz <- c(60,J[2],90,J[4:29],50,J[31:38])
#rz <- c(123,R[2],200,R[4:29],170,R[31:38])
rz <- c(140,R[2:29],140,R[31:38]) # edit with av.

my.data <- list(marr.j = CH.J.marray,
                marr.a = CH.A.marray,
                rel.j = rowSums(CH.J.marray),
                rel.a = rowSums(CH.A.marray),
                y = y,
                J = J,#jz,
                R = rz)#R

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2])

# inits
initial.values <- function(){list(l.mphij = rnorm(1, 0.2, 0.5),
                                  l.mphia = rnorm(1, 0.2, 0.5),
                                  l.mfec = rnorm(1, 0.2, 0.5),
                                  l.mim = rnorm(1, 0.2, 0.5),
                                  l.pj = rnorm(1, 0.2, 1),
                                  l.pa = rnorm(1, 0.2, 1),
                                  sigma.y = runif(1,0,10),
                                  N1 = rpois(dim(CH.J.marray)[2], 25000), # from 6500
                                  Nimm = rpois(dim(CH.J.marray)[2], 20000), # from 10000
                                  Nad = rpois(dim(CH.J.marray)[2], 240000))} #changed from 200000

# run code!
peng.ipm1.out <- nimbleMCMC(code = peng.ipm.imm.1,
                            constants = my.constants,
                            data = my.data,
                            inits = initial.values,
                            monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                         "Ntot", "sigma.y", 
                                         "geomean.lambda"),
                            niter = 30000,
                            nburnin = 10000,
                            nchains = 3)
MCMCsummary(peng.ipm1.out, round = 2) 
MCMCtrace(peng.ipm1.out, pdf = F)

pop <- MCMCchains(peng.ipm1.out, params = "Ntot", mcmc.list = F)
pop.med <- as.vector(apply(pop,2,median))
pop.CI <- apply(pop,2,function(x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
pop.data <- rbind(pop.med, pop.CI)
pop.data <- tibble::rownames_to_column(data.frame(t(pop.data)))
pop.data$rowname <- c(seq(1982,2019,1))

ggplot() +
  geom_line(aes(x=pop.data$rowname, y = pop.data$pop.med), color = "red", size = 1.5)+
  #geom_point(aes(x=pop.data$rowname, y = pop.data$pop.med), size = 3.5, color = "red")+
  geom_point(aes(x=pop.data$rowname, y = y), size = 3.5, color = "black")+
  #geom_ribbon(aes(x=pop.data$rowname, ymin=pop.data$X2.5., ymax = pop.data$X97.5.), linetype = 2, alpha = 0.3, fill = "red")+
  geom_ribbon(aes(x=pop.data$rowname, ymin=pop.data$X25., ymax = pop.data$X75.), linetype = 2, alpha = 0.4, fill = "red")+
  xlim(c(1980,2020))+
  xlab("Year")+
  ylab("# of Breeding Pairs")+
  theme_bw() +
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())
ggsave("projected_counts_11.22.21.png", dpi = 600, width = 12, height = 10)

##################################################################################################################

#### Time-varying Immigration Model ####
# build a pre-breeding, two-age model for the penguin population! with immigration added!

peng.ipm.imm.2 <- nimbleCode({
  
  ####
  # 1. Define priors for the parameters - very vague priors for everything...
  ####
  
  # initial population sizes
  N1[1] ~ T(dnorm(25000, 4e-08),0,)
  Nad[1] ~ T(dnorm(240000, 5e-08),0,) # from 190,000
  Nimm[1] ~ T(dnorm(20000, 4e-08),0,) # just guessing here
  
  # mean demographic parameters 
  #l.mphij ~ dnorm(0, 0.34)
  #l.mphia ~ dnorm(0, 0.34)
  #l.mfec ~ dnorm(0, 0.34)
  #l.mim ~ dnorm(0, 0.34)
  l.pj ~ dnorm(0, 0.34) 
  l.pa ~ dnorm(0, 0.34)
  
  # observation error
  tauy <- pow(sigma.y, -2) # (precision)
  sigma.y ~ dunif(0, 50000) # (sd) changed from 50000
  
  ####
  # 2. Constrain parameters
  ####
  
  for (t in 1:(nyears-1)){
    logit(phi.juv[t]) ~ dnorm(0, 0.34)#l.mphij
    logit(phi.ad[t]) ~ dnorm(0, 0.34)#l.mphia
    log(fec[t]) ~ dnorm(0, 0.34)#l.mfec
    log(omega[t]) ~ dnorm(0, 0.34)#l.mim
    #logit(pj[t]) <- l.pj
    #logit(pa[t]) <- l.pa
  }
  
  # detection constraints - juveniles
  for (t in 1:29){logit(pj[t]) <- l.pj}
  pj[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pj[t]) <- l.pj}
  
  # detection constraints - adults
  for (t in 1:29){logit(pa[t]) <- l.pa}
  pa[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pa[t]) <- l.pa}
  
  ####
  # 3. Derived parameters
  ####
  
  #mphij <- exp(l.mphij) / (1+exp(l.mphij))
  #mphia <- exp(l.mphia) / (1+exp(l.mphia))
  #mfec <- exp(l.mfec)
  #mim <- exp(l.mim)
  mpj <- exp(l.pj) / (1+exp(l.pj))
  mpa <- exp(l.pa) / (1+exp(l.pa))
  
  # population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1]/Ntot[t]
    l.lambda[t] <- log(lambda[t])
  }
  
  # geometric mean
  geomean.lambda <- exp((1/(nyears-1))*sum(l.lambda[1:(nyears-1)]))
  
  ####
  # 4. Likelihoods
  ####
  
  # 4.1 Likelihoods for population count data
  for (t in 2:nyears){
    mean1[t] <- 0.5 * fec[t-1] * phi.juv[t-1] * Nad[t-1]
    N1[t] ~ dpois(mean1[t])
    Nad[t] ~ dbin(phi.ad[t-1], Ntot[t-1])
    mpo[t] <- Ntot[t-1] * omega[t-1]
    Nimm[t] ~ dpois(mpo[t])
  }
  for (t in 1:nyears){
    Ntot[t] <- N1[t] + Nad[t] + Nimm[t]
  }
  
  # observation process
  for (t in 1:nyears){
    y[t] ~ T(dnorm(Ntot[t], tauy),0,)
  }
  
  # 4.2 Likelihoods for capture-recapture data
  # define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr.j[t,1:nyears] ~ dmulti(pr.j[t,1:nyears], rel.j[t])
    marr.a[t,1:nyears] ~ dmulti(pr.a[t,1:nyears], rel.a[t])
  }
  
  # m-array cell probabilities for juveniles and adults
  for (t in 1:(nyears-1)){
    qa[t] <- 1 - pa[t] # probability of recapture, adults
    qj[t] <- 1 - pj[t] # probability of recapture, juveniles
    
    # main diagonal
    pr.j[t,t] <- phi.juv[t] * pj[t]
    pr.a[t,t] <- phi.ad[t] * pa[t]
    
    # above main diagonal
    for (j in (t+1):(nyears-1)){
      pr.j[t,j] <- phi.juv[t] * prod(phi.ad[(t+1):j]) * qj[t] * prod(qa[t:(j-1)]) * pa[j]/qa[t]
      pr.a[t,j] <- prod(phi.ad[t:j]) * prod(qa[t:(j-1)]) * pa[j]
    }
    
    # below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  
  for (t in 1:(nyears-1)){
    # last column: probability of non-recapture
    pr.j[t,nyears] <- 1 - sum(pr.j[t, 1:(nyears-1)])
    pr.a[t,nyears] <- 1 - sum(pr.a[t, 1:(nyears-1)])
  }
  
  # 4.3 Likelihoods for productivity data
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*fec[t]
  }
})

# my data

# fucking with data
jz <- c(60,J[2],90,J[4:29],50,J[31:38])
#rz <- c(123,R[2],200,R[4:29],170,R[31:38])
rz <- c(140,R[2:29],140,R[31:38]) # edit with av.

my.data <- list(marr.j = CH.J.marray,
                marr.a = CH.A.marray,
                rel.j = rowSums(CH.J.marray),
                rel.a = rowSums(CH.A.marray),
                y = y,
                J = J,#jz,
                R = rz)#R

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2])

# inits
initial.values <- function(){list(phi.juv = rnorm(my.constants$nyears-1, 0.2, 0.5),
                                  phi.ad = rnorm(my.constants$nyears-1, 0.2, 0.5),
                                  fec = rnorm(my.constants$nyears-1, 0.2, 0.5),
                                  omega = rnorm(my.constants$nyears-1, 0.2, 0.5),
                                  l.pj = rnorm(1, 0.2, 1),
                                  l.pa = rnorm(1, 0.2, 1),
                                  sigma.y = runif(1,0,10),
                                  N1 = rpois(dim(CH.J.marray)[2], 25000),
                                  Nimm = rpois(dim(CH.J.marray)[2], 20000),
                                  Nad = rpois(dim(CH.J.marray)[2], 240000))} #changed from 200000

# run code!
peng.ipm2.out <- nimbleMCMC(code = peng.ipm.imm.2,
                            constants = my.constants,
                            data = my.data,
                            inits = initial.values,
                            monitors = c("mpj", "mpa", #"mphij", "mphia", "mfec", "mim",
                                         "Ntot", "sigma.y", 
                                         "omega", "fec", "phi.ad", "phi.juv",
                                         "geomean.lambda"),
                            niter = 300000,
                            nburnin = 100000,
                            nchains = 3)
MCMCsummary(peng.ipm2.out, round = 2) 
MCMCtrace(peng.ipm2.out, pdf = F)

# graphs!
# population size
pop <- MCMCchains(peng.ipm2.out, params = "Ntot", mcmc.list = F)
pop.med <- as.vector(apply(pop,2,median))
pop.CI <- apply(pop,2,function(x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
pop.data <- rbind(pop.med, pop.CI)
pop.data <- tibble::rownames_to_column(data.frame(t(pop.data)))
pop.data$rowname <- c(seq(1982,2019,1))

a <- ggplot() +
  geom_line(aes(x=pop.data$rowname, y = pop.data$pop.med), color = "red", size = 1.5)+
  #geom_point(aes(x=pop.data$rowname, y = pop.data$pop.med), size = 3.5, color = "red")+
  geom_point(aes(x=pop.data$rowname, y = y), size = 3.5, color = "black")+
  #geom_ribbon(aes(x=pop.data$rowname, ymin=pop.data$X2.5., ymax = pop.data$X97.5.), linetype = 2, alpha = 0.3, fill = "red")+
  geom_ribbon(aes(x=pop.data$rowname, ymin=pop.data$X25., ymax = pop.data$X75.), linetype = 2, alpha = 0.4, fill = "red")+
  xlim(c(1980,2020))+
  xlab("Year")+
  ylab("# of Breeding Pairs")+
  theme_bw() +
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# adult survival
asurv <- MCMCchains(peng.ipm2.out, params = "phi.ad", mcmc.list = F)
asurv.med <- as.vector(apply(asurv,2,median))
asurv.CI <- apply(asurv,2,function(x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
asurv.data <- rbind(asurv.med, asurv.CI)
asurv.data <- tibble::rownames_to_column(data.frame(t(asurv.data)))
asurv.data$rowname <- c(seq(1982,2018,1))

b <- ggplot() +
  geom_line(aes(x=asurv.data$rowname, y = asurv.data$asurv.med), color = "red", size = 1.5)+
  geom_ribbon(aes(x=asurv.data$rowname, ymin=asurv.data$X25., ymax = asurv.data$X75.), linetype = 2, alpha = 0.4, fill = "red")+
  xlim(c(1980,2020))+
  xlab("Year")+
  ylab("Adult Survival")+
  ylim(0,1)+
  theme_bw() +
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# juvenile survival
jsurv <- MCMCchains(peng.ipm2.out, params = "phi.juv", mcmc.list = F)
jsurv.med <- as.vector(apply(jsurv,2,median))
jsurv.CI <- apply(jsurv,2,function(x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
jsurv.data <- rbind(jsurv.med, jsurv.CI)
jsurv.data <- tibble::rownames_to_column(data.frame(t(jsurv.data)))
jsurv.data$rowname <- c(seq(1982,2018,1))
jsurv.data <- jsurv.data[-30,]

c <- ggplot() +
  geom_line(aes(x=jsurv.data$rowname, y = jsurv.data$jsurv.med), color = "red", size = 1.5)+
  geom_ribbon(aes(x=jsurv.data$rowname, ymin=jsurv.data$X25., ymax = jsurv.data$X75.), linetype = 2, alpha = 0.4, fill = "red")+
  xlim(c(1980,2020))+
  xlab("Year")+
  ylab("Juvenile Survival")+
  ylim(0,1)+
  theme_bw() +
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# fecundity
fec <- MCMCchains(peng.ipm2.out, params = "fec", mcmc.list = F)
fec.med <- as.vector(apply(fec,2,median))
fec.CI <- apply(fec,2,function(x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
fec.data <- rbind(fec.med, fec.CI)
fec.data <- tibble::rownames_to_column(data.frame(t(fec.data)))
fec.data$rowname <- c(seq(1982,2018,1))

d <- ggplot() +
  geom_line(aes(x=fec.data$rowname, y = fec.data$fec.med), color = "red", size = 1.5)+
  geom_ribbon(aes(x=fec.data$rowname, ymin=fec.data$X25., ymax = fec.data$X75.), linetype = 2, alpha = 0.4, fill = "red")+
  xlim(c(1980,2020))+
  xlab("Year")+
  ylab("Fecundity")+
  theme_bw() +
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

# immigration
imm <- MCMCchains(peng.ipm2.out, params = "omega", mcmc.list = F)
imm.med <- as.vector(apply(imm,2,median))
imm.CI <- apply(imm,2,function(x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
imm.data <- rbind(imm.med, imm.CI)
imm.data <- tibble::rownames_to_column(data.frame(t(imm.data)))
imm.data$rowname <- c(seq(1982,2018,1))
imm.data <- imm.data[-29,]

e <- ggplot() +
  geom_line(aes(x=imm.data$rowname, y = imm.data$imm.med), color = "red", size = 1.5)+
  geom_ribbon(aes(x=imm.data$rowname, ymin=imm.data$X25., ymax = imm.data$X75.), linetype = 2, alpha = 0.4, fill = "red")+
  xlim(c(1980,2020))+
  xlab("Year")+
  ylab("Immigration")+
  theme_bw() +
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())

library(patchwork)
a/(b|c)/(d|e)

ggsave("demandpopresults_plusimm_11.23.21.png", dpi = 600, width = 12, height = 12)

# plot relationships between variables
# NOTE: re-check this work...but I think I did this correctly?
all <- data.frame(time = seq(1982,2018,1),
                  growth = pop.data$pop.med[2:38]/pop.data$pop.med[1:37],
                  growth2.5 = pop.data$X2.5.[2:38]/pop.data$X2.5.[1:37],
                  growth97.5 = pop.data$X97.5.[2:38]/pop.data$X97.5.[1:37],
                  asurv = asurv.data$asurv.med, asurv2.5 = asurv.data$X2.5., asurv97.5 = asurv.data$X97.5.,
                  jsurv = jsurv.data$jsurv.med, jsurv2.5 = jsurv.data$X2.5., jsurv97.5 = jsurv.data$X97.5.,
                  fec = fec.data$fec.med, fec2.5 = fec.data$X2.5., fec97.5 = fec.data$X97.5.,
                  imm = imm.data$imm.med, imm2.5 = imm.data$X2.5., imm97.5 = imm.data$X97.5.)
# adjust overestimates with averages
all$fec[1] <- 0.534; all$fec2.5[1] <- 0.376; all$fec97.5[1] <- 0.959
all$imm[29] <- 0.105; all$imm2.5[29] <- 0.0287; all$imm97.5[29] <- 0.27
# calculate correlation
cor(all[,-1])

# relationships to population growth
# graph
a <- ggplot(all, aes(jsurv, growth)) + 
  geom_point()+
  #geom_errorbar(aes(xmin = jsurv2.5, xmax = jsurv97.5))+
  #geom_errorbar(aes(ymin = growth2.5, ymax = growth97.5))+
  geom_smooth(method = "lm")+
  geom_text(x = 0.35, y = 1.1, label = "r=0.347", color = "red", size = 8)+
  xlab("Juvenile Survival")+
  ylab("Population Growth")+
  theme_bw() +
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())
b <- ggplot(all, aes(asurv, growth)) + 
  geom_point()+
  geom_smooth(method = "lm")+
  geom_text(x = 0.95, y = 1.1, label = "r=0.561", color = "red", size = 8)+
  xlab("Adult Survival")+
  ylab("Population Growth")+
  theme_bw() +
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())
c <- ggplot(all, aes(fec, growth)) + 
  geom_point()+
  geom_smooth(method = "lm")+
  geom_text(x = 0.90, y = 1.1, label = "r=0.062", color = "red", size = 8)+
  xlab("Fecundity")+
  ylab("Population Growth")+
  theme_bw() +
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())
d <- ggplot(all, aes(imm, growth)) + 
  geom_point()+
  geom_smooth(method = "lm")+
  geom_text(x = 0.17, y = 1.1, label = "r=0.332", color = "red", size = 8)+
  xlab("Immigration")+
  ylab("Population Growth")+
  theme_bw() +
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank())
(a+b)/(c+d)
ggsave("demandpopresults_plusimm_cors_11.23.21.png", dpi = 600, width = 12, height = 12)

##################################################################################################################

#### Simple Immigration Model WITH ENVIRONMENTAL VARS ####
# build a pre-breeding, two-age model for the penguin population! with immigration added! add variables!

peng.ipm.imm.3 <- nimbleCode({
  
  ####
  # 1. Define priors for the parameters - very vague priors for everything...
  ####
  
  # initial population sizes
  N1[1] ~ T(dnorm(25000, 4e-08),0,) # from 6500
  Nad[1] ~ T(dnorm(240000, 5e-08),0,) # from 200000
  Nimm[1] ~ T(dnorm(20000, 4e-08),0,) # just guessing here
  
  # mean demographic parameters 
  l.mphij ~ dnorm(0, 0.34)
  l.mphia ~ dnorm(0, 0.34)
  l.mfec ~ dnorm(0, 0.34)
  l.mim ~ dnorm(0, 0.34)
  l.pj ~ dnorm(0, 0.34) 
  l.pa ~ dnorm(0, 0.34)
  
  # observation error
  tauy <- pow(sigma.y, -2) # (precision)
  sigma.y ~ dunif(0, 50000) # (sd) changed from 50000
  
  # covariates
  for (i in 1:11){alpha[i] ~ dnorm(0, 0.34)}
  
  ####
  # 2. Constrain parameters
  ####
  
  for (t in 1:(nyears-1)){
    log(fec[t]) <-  l.mfec + (alpha[1] * rain[t]) + (alpha[2] * temp[t]) + (alpha[3] * ssta_b[t]) + (alpha[4] * ssta_m[t])
    logit(phi.juv[t]) <- l.mphij + (alpha[5] * ssta_b[t]) + (alpha[6] * ssta_m_l[t])
    logit(phi.ad[t]) <-  l.mphia + (alpha[7] * temp[t]) + (alpha[8] * ssta_b[t]) + (alpha[9] * ssta_m_l[t])
    #log(omega[t]) <-  l.mim + (alpha[10] * ssta_b_l[t]) + (alpha[11] * ssta_m_l[t]) #lagged
    #log(omega[t]) <-  l.mim + (alpha[10] * ssta_b[t]) + (alpha[11] * ssta_m[t]) #non-lagged
    log(omega[t]) <-  l.mim + (alpha[10] * ssta_b[t]) + (alpha[11] * ssta_m_l[t]) #mixed
    #logit(pj[t]) <- l.pj
    #logit(pa[t]) <- l.pa
  }
  
  # detection constraints - juveniles
  for (t in 1:29){logit(pj[t]) <- l.pj}
  pj[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pj[t]) <- l.pj}
  
  # detection constraints - adults
  for (t in 1:29){logit(pa[t]) <- l.pa}
  pa[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pa[t]) <- l.pa}
  
  ####
  # 3. Derived parameters
  ####
  
  mphij <- exp(l.mphij) / (1+exp(l.mphij))
  mphia <- exp(l.mphia) / (1+exp(l.mphia))
  mfec <- exp(l.mfec)
  mim <- exp(l.mim)
  mpj <- exp(l.pj) / (1+exp(l.pj))
  mpa <- exp(l.pa) / (1+exp(l.pa))
  
  # population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1]/Ntot[t]
    l.lambda[t] <- log(lambda[t])
  }
  
  # geometric mean
  geomean.lambda <- exp((1/(nyears-1))*sum(l.lambda[1:(nyears-1)]))
  
  ####
  # 4. Likelihoods
  ####
  
  # 4.1 Likelihoods for population count data
  for (t in 2:nyears){
    mean1[t] <- 0.5 * fec[t-1] * phi.juv[t-1] * Nad[t-1]
    N1[t] ~ dpois(mean1[t])
    Nad[t] ~ dbin(phi.ad[t-1], Ntot[t-1])
    mpo[t] <- Ntot[t-1] * omega[t-1]
    Nimm[t] ~ dpois(mpo[t])
  }
  for (t in 1:nyears){
    Ntot[t] <- N1[t] + Nad[t] + Nimm[t]
  }
  
  # observation process
  for (t in 1:nyears){
    y[t] ~ T(dnorm(Ntot[t], tauy),0,)
  }
  
  # 4.2 Likelihoods for capture-recapture data
  # define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr.j[t,1:nyears] ~ dmulti(pr.j[t,1:nyears], rel.j[t])
    marr.a[t,1:nyears] ~ dmulti(pr.a[t,1:nyears], rel.a[t])
  }
  
  # m-array cell probabilities for juveniles and adults
  for (t in 1:(nyears-1)){
    qa[t] <- 1 - pa[t] # probability of recapture, adults
    qj[t] <- 1 - pj[t] # probability of recapture, juveniles
    
    # main diagonal
    pr.j[t,t] <- phi.juv[t] * pj[t]
    pr.a[t,t] <- phi.ad[t] * pa[t]
    
    # above main diagonal
    for (j in (t+1):(nyears-1)){
      pr.j[t,j] <- phi.juv[t] * prod(phi.ad[(t+1):j]) * qj[t] * prod(qa[t:(j-1)]) * pa[j]/qa[t]
      pr.a[t,j] <- prod(phi.ad[t:j]) * prod(qa[t:(j-1)]) * pa[j]
    }
    
    # below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  
  for (t in 1:(nyears-1)){
    # last column: probability of non-recapture
    pr.j[t,nyears] <- 1 - sum(pr.j[t, 1:(nyears-1)])
    pr.a[t,nyears] <- 1 - sum(pr.a[t, 1:(nyears-1)])
  }
  
  # 4.3 Likelihoods for productivity data
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*fec[t]
  }
})

# my data

# fucking with data
jz <- c(60,J[2],90,J[4:29],50,J[31:38])
#rz <- c(123,R[2],200,R[4:29],170,R[31:38])
rz <- c(140,R[2:29],140,R[31:38]) # edit with av.

my.data <- list(marr.j = CH.J.marray,
                marr.a = CH.A.marray,
                rel.j = rowSums(CH.J.marray),
                rel.a = rowSums(CH.A.marray),
                y = y,
                J = J,#jz,
                R = rz)#R

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = as.vector(scale(c(36.6, rain_60$amt_rain[1:28], 36.6, rain_60$amt_rain[29:36]))),
                     temp = as.vector(scale(c(0.383, num_days_25$p_gt25[1:28], 0.383, num_days_25$p_gt25[29:36]))),
                     ssta_b = as.vector(ssta_breeding$av_s[1:38]),
                     ssta_b_l = as.vector(c(ssta_breeding$av_s[2:38],0)), # lagged (forward) SSTA
                     ssta_m = as.vector(ssta_migration$mean.2[1:38]),
                     ssta_m_l = as.vector(c(ssta_migration$mean.2[2:39]))) # lagged (forward) SSTA!

# inits
initial.values <- function(){list(l.mphij = rnorm(1, 0.2, 0.5),
                                  l.mphia = rnorm(1, 0.2, 0.5),
                                  l.mfec = rnorm(1, 0.2, 0.5),
                                  l.mim = rnorm(1, 0.2, 0.5),
                                  l.pj = rnorm(1, 0.2, 1),
                                  l.pa = rnorm(1, 0.2, 1),
                                  alpha = rnorm(11, 0.2, 1),
                                  sigma.y = runif(1,0,10),
                                  N1 = rpois(dim(CH.J.marray)[2], 25000), # from 6500
                                  Nimm = rpois(dim(CH.J.marray)[2], 20000), # from 10000
                                  Nad = rpois(dim(CH.J.marray)[2], 240000))} #changed from 200000

# run code!
peng.ipm3.out <- nimbleMCMC(code = peng.ipm.imm.3,
                            constants = my.constants,
                            data = my.data,
                            inits = initial.values,
                            monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                         "Ntot", "sigma.y", "alpha",
                                         "geomean.lambda"),
                            niter = 300000,
                            nburnin = 100000,
                            nchains = 3)
MCMCsummary(peng.ipm3.out, round = 2) 
MCMCtrace(peng.ipm3.out, pdf = F) # NOTE: alpha[10] and alpha[11] are not mixing. add more sims and thin

# GRAPHS
# graph parameter values
MCMCplot(peng.ipm3.out, params = 'alpha')

ggsave("ipm_imm_param_results_11.23.21.png", dpi = 600)


