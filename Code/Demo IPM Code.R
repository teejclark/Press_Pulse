
# set directory
rm(list = ls())
setwd("~/Post-Doc/Data/Gownaris et al. 2019/Gownaris&Boersma_DataS1")

library(nimble)
library(tidyverse)
library(MCMCvis)

##################################################################################################################

#### First Goal: Run IPM in Nimble

# load data

# Population counts (from years 1 to 10)
y <- c(45, 48, 44, 59, 62, 62, 55, 51, 46, 42)

# Capture-recapture data (in m-array format, from years 1 to 10)
m <- matrix(c(11,  0,  0,  0,  0,  0,  0,  0,  0,  70,
              0, 12,  0,  1,  0,  0,  0,  0,  0,  52,
              0,  0, 15,  5,  1,  0,  0,  0,  0,  42,
              0,  0,  0,  8,  3,  0,  0,  0,  0,  51,
              0,  0,  0,  0,  4,  3,  0,  0,  0,  61,
              0,  0,  0,  0,  0, 12,  2,  3,  0,  66,
              0,  0,  0,  0,  0,  0, 16,  5,  0,  44,
              0,  0,  0,  0,  0,  0,  0, 12,  0,  46,
              0,  0,  0,  0,  0,  0,  0,  0, 11,  71,
              10,  2,  0,  0,  0,  0,  0,  0,  0,  13,
              0,  7,  0,  1,  0,  0,  0,  0,  0,  27,
              0,  0, 13,  2,  1,  1,  0,  0,  0,  14,
              0,  0,  0, 12,  2,  0,  0,  0,  0,  20,
              0,  0,  0,  0, 10,  2,  0,  0,  0,  21,
              0,  0,  0,  0,  0, 11,  2,  1,  1,  14,
              0,  0,  0,  0,  0,  0, 12,  0,  0,  18,
              0,  0,  0,  0,  0,  0,  0, 11,  1,  21,
              0,  0,  0,  0,  0,  0,  0,  0, 10,  26), ncol = 10, byrow = TRUE)
m.j <- m[1:9,]
m.a <- m[10:18,]

# Productivity data (from years 1 to 9)
J <- c(64, 132,  86, 154, 156, 134, 116, 106, 110)
R <- c(21, 28, 26, 38, 35, 33, 31, 30, 33) 

# build code
simp.IPM <- nimbleCode({
  
  ####
  # 1. Define priors
  ####
  
  # observation error
  tauy <- pow(sigma.y, -2)
  sigma.y ~ dunif(0, 50)
  sigma2.y <- pow(sigma.y,2)
  
  # initial population sizes
  N1[1] ~ T(dnorm(100, 0.0001),0,)
  Nad[1] ~ T(dnorm(100, 0.0001),0,)
  
  # survival and recapture probs and productivity
  for (t in 1:(nyears-1)){
    sjuv[t] <- mean.sjuv
    sad[t] <- mean.sad
    p[t] <- mean.p
    f[t] <- mean.fec
  }
  
  mean.sjuv ~ dunif(0,1)
  mean.sad ~ dunif(0,1)
  mean.p ~ dunif(0,1)
  mean.fec ~ dunif(0,20)
  
  ####
  # 2. Derived parameters
  ####
  
  # pop growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1]/Ntot[t]
  }
  
  ####
  # 3. Likelihoods
  ####
  
  # population count data
  for (t in 2:nyears){
    mean1[t] <- f[t-1]/2 * sjuv[t-1] * Ntot[t-1]
    N1[t] ~ dpois(mean1[t])
    Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
  }
  for (t in 1:nyears){
    Ntot[t] <- Nad[t] + N1[t]
  }
  
  # observation process
  for (t in 1:nyears){
    y[t] ~ dnorm(Ntot[t], tauy)
  }
  
  # define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr.j[t,1:nyears] ~ dmulti(pr.j[t,1:nyears], rel.j[t])
    marr.a[t,1:nyears] ~ dmulti(pr.a[t,1:nyears], rel.a[t])
  }
  
  # m-array cell probabilities for juveniles and adults
  for (t in 1:(nyears-1)){
    q[t] <- 1 - p[t] # probability of recapture
    
    # main diagonal
    pr.j[t,t] <- phi.juv[t] * p[t]
    pr.a[t,t] <- phi.ad[t] * p[t]
    
    # above main diagonal
    for (j in (t+1):(nyears-1)){
      pr.j[t,j] <- phi.juv[t] * prod(phi.ad[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(phi.ad[t:j]) * prod(q[t:(j-1)]) * p[j]
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
  
  # productivity data
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*f[t]
  }
            
})

# my data
my.data <- list(marr.j = m.j,
                marr.a = m.a,
                rel.j = rowSums(m.j),
                rel.a = rowSums(m.a),
                y = y,
                J = J,
                R = R)

# constants
my.constants <- list(nyears = dim(m.j)[2])

# inits
initial.values <- function(){list(mean.sjuv = runif(1,0,1),
                                  mean.sad = runif(1,0,1),
                                  mean.p = runif(1,0,1),
                                  mean.fec = runif(1,0,10),
                                  N1 = rpois(dim(m.j)[2],30),
                                  Nad = rpois(dim(m.j)[2],30),
                                  sigma.y = runif(1,0,10))}

simp.IPM.out <- nimbleMCMC(code = simp.IPM,
                            constants = my.constants,
                            data = my.data,
                            inits = initial.values,
                            monitors = c("mean.sjuv","mean.sad","mean.p","mean.fec",
                                         "N1","Nad","Ntot","lambda"),
                            niter = 30000,
                            nburnin = 10000,
                            nchains = 3)
MCMCsummary(simp.IPM.out, round = 2)

##################################################################################################################

#### Data: 33 years from 1983 to 2015 (last chicks were banded 2010)
# NOTE: lots of data missing, might be useful to add in fake data?

# capture histories from Gownaris et al. 2019
# NOTE: add code from "Multinomial Survival Code.R" sometime later
y <- as.matrix(read.table("Gownaris&Boersma_RecaptureHistories.txt", sep = ",", header = F))
colnames(y) <- NULL
head(y)

# nests monitored from Pozzi et al. 2015 = might not be needed
# real data is from 1987 to 2009
R <- c(NA, NA, NA, NA, 
       183,160,180,161,127,151,151,157,151,150,148,155,140,121,150,131,146,142,140,125,128,134,148,
       NA, NA, NA, NA, NA, NA) # fill in later

# penguin productivity from Boersma 2008 = chicks fledged per nest
# real data is from 1983 to 2006
f <- c(0.51,0.03,0.63,0.44,0.12,0.70,0.51,0.30,0.13,0.59,0.60,0.52,
       0.56,0.95,0.43,0.70,0.13,0.06,0.64,0.32,0.58,0.81,0.85,0.85,
       NA, NA, NA, NA, NA, NA, NA, NA, NA) # fill in later?

# total number of fledglings counted, product of productivity and nests monitored???
J <- pop*f

# mean active nest/100 m2 data is from Rebstock et al. 2016
# real data is 1987 to 2014
# NOTE: this data is weird and should be figured out!
nestdensity <- c(NA, NA, NA, NA,
                 9.67, 8.35, 9.49, 8.34, 5.94, 7.78, 7.32, 7.61, 7.081, 
                 7.20, 6.77, 7.06, 6.01, 4.85, 6.40, 5.22, 6.26, 6.00, 6.14, 
                 5.00, 4.98, 5.54, 6.29, 4.81,NA,7.43,3.89,3.95,
                 NA) # fill in later?
pop <- nestdensity*(3.5*10^4) # this is a guesstimate!

####

peng.ipm.1 <- nimbleCode({
  
  ####
  # 1. Define priors
  ####
  
  # initial population sizes
  N1[1] ~ T(dnorm(12500, 10),0,)
  Nad[1] ~ T(dnorm(200000, 100),0,)
  
  # observation error
  tauy <- pow(sigma.y, -2)
  sigma.y ~ dunif(0, 50)
  taufec <- pow(sigma.fec, -2)
  sigma.fec ~ dunif(0,50)
  
  # constraints
  for (t in 1:(nyears-1)){
    phi.juv[t] <- mean.phijuv
    phi.ad[t] <- mean.phiad
    p[t] <- mean.p
    fec[t] <- mean.fec
  }
  
  # priors
  mean.phijuv ~ dunif(0,1) # prior for mean juv. survival
  mean.phiad ~ dunif(0,1) # prior for mean ad. survival
  mean.p ~ dunif(0,1) # prior for mean recapture
  mean.fec ~ dunif(0,2) # prior for fecundity
  
  ####
  # 2. Derived parameters
  ####
  
  # population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1]/Ntot[t]
    l.lambda[t] <- log(lambda[t])
  }
  
  # geometric mean
  geomean.lambda <- exp((1/(nyears-1))*sum(l.lambda[1:(nyears-1)]))
  
  ####
  # 3. Likelihoods
  ####
  
  # Likelihoods for population count data
  for (t in 2:nyears){
    mean1[t] <- fec[t-1]/2 * phi.juv[t-1] * Nad[t-1]
    N1[t] ~ dpois(mean1[t])
    Nad[t] ~ dbin(phi.ad[t-1], Ntot[t-1])
  }
  
  for (t in 1:nyears){
    Ntot[t] <- N1[t] + Nad[t]
  }
  
  # observation process
  for (t in 1:nyears){
    pop[t] ~ dnorm(Ntot[t], tauy)
  }
  
  # likelihood for survival model
  
  # define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr.j[t,1:nyears] ~ dmulti(pr.j[t,1:nyears], rel.j[t])
    marr.a[t,1:nyears] ~ dmulti(pr.a[t,1:nyears], rel.a[t])
  }
  
  # m-array cell probabilities for juveniles and adults
  for (t in 1:(nyears-1)){
    q[t] <- 1 - p[t] # probability of recapture
    
    # main diagonal
    pr.j[t,t] <- phi.juv[t] * p[t]
    pr.a[t,t] <- phi.ad[t] * p[t]
    
    # above main diagonal
    for (j in (t+1):(nyears-1)){
      pr.j[t,j] <- phi.juv[t] * prod(phi.ad[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(phi.ad[t:j]) * prod(q[t:(j-1)]) * p[j]
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
  
  # Likelihood for productivity data
  for (t in 1:(nyears-1)){
    f[t] ~ dnorm(fec[t], taufec)
  }
})

# my data
my.data <- list(marr.j = CH.J.marray,
                marr.a = CH.A.marray,
                rel.j = rowSums(CH.J.marray),
                rel.a = rowSums(CH.A.marray),
                f = f,
                pop = pop)

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2])

# inits
initial.values <- function(){list(mean.phijuv = runif(1,0,1),
                                  mean.phiad = runif(1,0,1),
                                  mean.p = runif(1,0,1),
                                  mean.fec = runif(1,0,2),
                                  N1 = rpois(dim(CH.J.marray)[2], 12500),
                                  Nad = rpois(dim(CH.J.marray)[2], 200000),
                                  sigma.y = runif(1,0,10),
                                  sigma.fec = runif(1,0,10))}

# run code!
peng.ipm1.out <- nimbleMCMC(code = peng.ipm.1,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values,
                             monitors = c("mean.phijuv", "mean.phiad", "mean.p", 
                                          "mean.fec",
                                          "N1", "Nad", "Ntot", 
                                          "lambda","geomean.lambda"),
                             niter = 300000,
                             nburnin = 100000,
                             nchains = 3)
MCMCsummary(peng.ipm1.out, round = 2) 
MCMCtrace(peng.ipm1.out, pdf = F)




