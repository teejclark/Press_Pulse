# set directory
rm(list = ls())
setwd("~/Post-Doc/Data/Punta Tombo Data")

library(tidyverse)
library(nimble)
library(ggplot2)
library(MCMCvis)

# GOAL: RUN sensitivity analyses to explore the effects of demographic rates, climate, and extremes on lambda

##################################################################################################################

#### 1) Testing Simulating with Nimble
#https://r-nimble.org/nimbleExamples/simulation_from_model.html


# 1) Create NimbleCode describing the model
simCode <- nimbleCode({
  for (i in 1:nsite){
    z[i] ~ dbern(prob = psi) # true occupancy status
    for (j in 1:nvisit){
      y[i,j] ~ dbern(z[i] * p) # observed data
    }
  }
  psi ~ dunif(0,1)
  p ~ dunif(0,1)
})

# 2) Build the model
# provide constants and initial values for p and psi
simModel <- nimbleModel(code = simCode,
                        constants = list(nsite = 10, nvisit = 3),
                        inits = list(psi = 0.7, p = 0.33))

# 3) Identify nodes and simulate
# can be useful to run "model$getDependencies() to gather nodes needed to be simulated
nodesToSim <- simModel$getDependencies(c("psi", "p"),
                                       self = F, # don't want to see psi or p in list
                                       downstream = T)

# 4) Simulate!
simModel$simulate(nodesToSim)
simModel$y # get back values for y!!!!
simModel$z # also get latent state, z!!!!


# OK - let's run our IPM with immigration to prove it works

# run model = simple immigration one
peng.ipm.imm.sim <- nimbleCode({
  
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
    logit(pj[t]) <- l.pj
    logit(pa[t]) <- l.pa
  }
  
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

# constants = change years simulated
my.constants <- list(nyears = 100)

# inits = this is how you mess with the simulations!
# use initial results from parameterized model
initial.values = list(N1 = 12000, 
                      Nimm = 10000,
                      Nad = 100000,
                      sigma.y = 22000,
                      l.mphij = -2.197, # need to logit/log these values
                      l.mphia = 1.658,
                      l.mfec = -0.673,
                      l.mim = -2.120,
                      l.pj = -2.197,
                      l.pa = -0.040)

# build model
peng.simModel <- nimbleModel(code = peng.ipm.imm.sim,
                        constants = my.constants,
                        inits = initial.values)

# identify nodes to simulate
nodesToSim <- peng.simModel$getDependencies(c("N1[1]", "Nimm[1]", "Nad[1]", "sigma.y", "J",
                                         "l.mphij", "l.mphia", "l.mfec", "l.mim", "l.pj", "l.pa"),
                                       self = F, downstream = T)

# simulate?
peng.simModel$simulate(nodesToSim)

peng.simModel$Ntot


# another option: build model by yourself? cause I'm impatient
nyears <- 10
# params - NOTE: can make this stochastic too...
phi.juv <- 0.10
phi.ad <- 0.84
fec <- 0.51
omega <- 0.12
# data.frame
sim <- data.frame(N1 = c(12000, rep(NA,nyears)),
                  Nimm = c(10000, rep(NA, nyears)),
                  Nad = c(100000, rep(NA, nyears)))
sim$Ntot <- sim$N1 + sim$Nimm + sim$Nad

for (i in 2:(nyears+1)){
  mean1 <- 0.5 * fec * phi.juv * sim$Nad[i-1]
  sim$N1[i] <- rpois(1, mean1)
  sim$Nad[i] <- rbinom(1, sim$Ntot[i-1], phi.ad)
  mpo <- sim$Ntot[i-1] * omega
  sim$Nimm[i] <- rpois(1, mpo)
  sim$Ntot[i] <- sim$N1[i] + sim$Nad[i] + sim$Nimm[i]

}

##################################################################################################################

#### 2) Beefing Up the Simulations

# run environmentally varying model

peng.ipm.imm.sim2 <- nimbleCode({
  
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
  alpha1 ~ dnorm(0, 0.34)
  alpha2 ~ dnorm(0, 0.34)
  alpha3 ~ dnorm(0, 0.34)
  alpha4 ~ dnorm(0, 0.34)
  alpha5 ~ dnorm(0, 0.34)
  alpha6 ~ dnorm(0, 0.34)
  alpha7 ~ dnorm(0, 0.34)
  alpha8 ~ dnorm(0, 0.34)
  alpha9 ~ dnorm(0, 0.34)

  
  ####
  # 2. Constrain parameters
  ####
  
  for (t in 1:(nyears-1)){
    log(fec[t]) <- l.mfec + (alpha1 * rain[t]) + (alpha2 * temp[t]) + (alpha3 * ssta_b[t]) + (alpha4 * ssta_m[t])
    logit(phi.juv[t]) <- l.mphij + (alpha5 * ssta_b[t]) + (alpha6 * ssta_m_l[t])
    logit(phi.ad[t]) <-  l.mphia + (alpha7 * temp[t]) + (alpha8 * ssta_b[t]) + (alpha9 * ssta_m_l[t])
    log(omega[t]) <-  l.mim
    logit(pj[t]) <- l.pj
    logit(pa[t]) <- l.pa
  }
  
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

nyears <- 50
# constants = change years simulated
# my.constants <- list(nyears = nyears,
#                      rain = as.vector(scale(rgamma(nyears,1.57396674,0.04298969))),
#                      temp = as.vector(scale(rnorm(nyears,0.3836863,0.1118897))),
#                      ssta_b = as.vector(scale(rnorm(nyears,-0.01683185,0.47938401))),
#                      ssta_m = as.vector(scale(rnorm(nyears,-0.001235519,0.787513205))))
# my.constants$ssta_m_l <- as.vector(c(my.constants$ssta_m[2:nyears],rnorm(1,0,1)))

my.constants <- list(nyears = nyears, # UNSCALED
                     rain = as.vector(rgamma(nyears,1.57396674,0.04298969)),
                     temp = as.vector(rnorm(nyears,0.3836863,0.1118897)),
                     ssta_b = as.vector(rnorm(nyears,-0.01683185,0.47938401)),
                     ssta_m = as.vector(rnorm(nyears,-0.001235519,0.787513205)))
my.constants$ssta_m_l <- as.vector(c(my.constants$ssta_m[2:nyears],rnorm(1,0,1)))

# inits = this is how you mess with the simulations!
# use initial results from parameterized model
initial.values = list(N1 = 12000, 
                      Nimm = 10000,
                      Nad = 100000,
                      sigma.y = 22000,
                      l.mphij = qlogis(0.1), 
                      l.mphia = qlogis(0.84),
                      l.mfec = log(0.51),
                      l.mim = log(0.12),
                      l.pj = qlogis(0.1),
                      l.pa = qlogis(0.49),
                      alpha1 = -0.01, alpha2 = 0.22, alpha3 = 0.14, alpha4 = 0.06, # UNSCALED VALS
                      alpha5 = -0.82, alpha6 = 0.62,
                      alpha7 = -0.74, alpha8 = 0.08, alpha9 = .29)
                      #alpha1 = -0.33, alpha2 = 0.02, alpha3 = 0.07, alpha4 = 0.05, # SCALED VALS
                      #alpha5 = -0.4, alpha6 = 0.61,
                      #alpha7 = -0.10, alpha8 = 0.04, alpha9 = 0.38)

# build model
peng.simModel2 <- nimbleModel(code = peng.ipm.imm.sim2,
                             constants = my.constants,
                             inits = initial.values)

# identify nodes to simulate
nodesToSim <- peng.simModel2$getDependencies(c("N1[1]", "Nimm[1]", "Nad[1]", "sigma.y", "J",
                                              "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", 
                                              "alpha6", "alpha7", "alpha8", "alpha9",
                                              "l.mphij", "l.mphia", "l.mfec", "l.mim", "l.pj", "l.pa"),
                                            self = F, downstream = T)

# simulate?
peng.simModel2$simulate(nodesToSim)

peng.simModel2$Ntot; plot(seq(1,nyears,1),peng.simModel2$Ntot,type="l")
peng.simModel2$geomean.lambda

##################################################################################################################

#### 3) SIMULATE MODELS NOW USING FOR LOOPS

nsims <- 10
nyears <- 50

pop <- data.frame(matrix(NA,nyears,nsims))
lambda <- rep(NA,nsims)

for (i in 1:nsims){

my.constants <- list(nyears = nyears, # UNSCALED
                     rain = as.vector(rgamma(nyears,1.57396674,0.04298969)),
                     temp = as.vector(rnorm(nyears,0.3836863,0.1118897)),
                     ssta_b = as.vector(rnorm(nyears,-0.01683185,0.47938401)),
                     ssta_m = as.vector(rnorm(nyears,-0.001235519,0.787513205))) 
my.constants$ssta_m_l <- as.vector(c(my.constants$ssta_m[2:nyears],rnorm(1,-0.001235519,0.787513205)))

# inits = this is how you mess with the simulations!
# use initial results from parameterized model
initial.values = list(N1 = 12000, 
                      Nimm = 10000,
                      Nad = 120000,
                      sigma.y = 22000,
                      l.mphij = qlogis(0.1), 
                      l.mphia = qlogis(0.84),
                      l.mfec = log(0.51),
                      l.mim = log(0.12),
                      l.pj = qlogis(0.1),
                      l.pa = qlogis(0.49),
                      alpha1 = -0.01, alpha2 = 0.11, alpha3 = 0.16, alpha4 = 0.05, # UNSCALED VALS
                      alpha5 = -0.82, alpha6 = 0.61,
                      alpha7 = -0.78, alpha8 = 0.08, alpha9 = 0.29)

# build model
peng.simModel2 <- nimbleModel(code = peng.ipm.imm.sim2,
                              constants = my.constants,
                              inits = initial.values)

# identify nodes to simulate
nodesToSim <- peng.simModel2$getDependencies(c("N1[1]", "Nimm[1]", "Nad[1]", "sigma.y", "J",
                                               "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", 
                                               "alpha6", "alpha7", "alpha8", "alpha9",
                                               "l.mphij", "l.mphia", "l.mfec", "l.mim", "l.pj", "l.pa"),
                                             self = F, downstream = T)

# simulate
peng.simModel2$simulate(nodesToSim)

pop[,i] <- peng.simModel2$Ntot
lambda[i] <- peng.simModel2$geomean.lambda

}

# plot simulations with real pop
y_adj <- y[6:length(y)]
totals <- data.frame(time = seq(1987,2069,1),
                     orig = c(y_adj,rep(NA,nrow(pop))),
                     new1 = c(rep(NA,length(y_adj)),pop[,1]),
                     new2 = c(rep(NA,length(y_adj)),pop[,2]),
                     new3 = c(rep(NA,length(y_adj)),pop[,3]),
                     new4 = c(rep(NA,length(y_adj)),pop[,4]),
                     new5 = c(rep(NA,length(y_adj)),pop[,5]),
                     new6 = c(rep(NA,length(y_adj)),pop[,6]),
                     new7 = c(rep(NA,length(y_adj)),pop[,7]),
                     new8 = c(rep(NA,length(y_adj)),pop[,8]),
                     new9 = c(rep(NA,length(y_adj)),pop[,9]),
                     new10 = c(rep(NA,length(y_adj)),pop[,10]))

ggplot(totals)+
  geom_point(aes(x = time, y = orig), size = 2)+
  geom_line(aes(x = time, y = new1), color = "red", size = 2)+
  geom_line(aes(x = time, y = new2), color = "blue", size = 2)+
  geom_line(aes(x = time, y = new3), color = "yellow", size = 2)+
  geom_line(aes(x = time, y = new4), color = "green", size = 2)+
  geom_line(aes(x = time, y = new5), color = "purple", size = 2)+
  geom_line(aes(x = time, y = new6), color = "brown", size = 2)+
  geom_line(aes(x = time, y = new7), color = "black", size = 2)+
  geom_line(aes(x = time, y = new8), color = "magenta", size = 2)+
  geom_line(aes(x = time, y = new9), color = "orange", size = 2)+
  geom_line(aes(x = time, y = new10), color = "gray", size = 2)+
  theme_bw()


# plot populations
plot(seq(1,nyears,1),pop[,1], type = "l")
for(i in 2:nsims){
  lines(seq(1,nyears,1),pop[,i], col = i)
}

##################################################################################################################

#### 4) RUN SENSITIVITY ANALYSIS ON DEMOGRAPHIC RATES
# NOTE: it's an open question how to run these. For now, prefer to run over range of natural variations in demography
# NOTE: also does partial rank correlation coefficients make sense? I think so, cause you vary everything all at once.
# NOTE: lastly, 100 parameter combinations, 100 times make sense? maybe up the param combos if possible?

nyears <- 100 # number of years to run stochastic sims
npars <- 4 # number of parameters to analyze
nparsim <- 100 # number of parameter combinations (i.e., unique draws from distributions)
nsims <- 100 # number of simulations per parameter combination

# parameters
dem <- data.frame(matrix(NA,nparsim,npars))
colnames(dem) <- c("asurv", "jsurv", "imm", "fec")

# pop growth
pop <- array(data = NA, dim = c(nyears,nsims,nparsim))
lambda <- data.frame(matrix(NA,nparsim,nsims))

for (i in 1:nparsim){
  for (j in 1:nsims){
  
  my.constants <- list(nyears = nyears, # UNSCALED
                       rain = as.vector(rgamma(nyears,1.57396674,0.04298969)),
                       temp = as.vector(rnorm(nyears,0.3836863,0.1118897)),
                       ssta_b = as.vector(rnorm(nyears,-0.01683185,0.47938401)),
                       ssta_m = as.vector(rnorm(nyears,-0.001235519,0.787513205))) 
  my.constants$ssta_m_l <- as.vector(c(my.constants$ssta_m[2:nyears],rnorm(1,-0.001235519,0.787513205)))
  
  # extract parameter values
  dem$asurv[i] <- runif(1, 0.63, 0.99)
  dem$jsurv[i] <- runif(1, 0.01, 0.42)
  dem$imm[i] <- runif(1, 0.01, 0.2)
  dem$fec[i] <- runif(1, 0.03, 1)
  
  # inits = this is how you mess with the simulations!
  initial.values = list(N1 = 12000, 
                        Nimm = 10000,
                        Nad = 120000,
                        sigma.y = 22000,
                        # random draws from uniform dists
                        l.mphia = qlogis(dem$asurv[i]),
                        l.mphij = qlogis(dem$jsurv[i]),
                        l.mim = log(dem$imm[i]),
                        l.mfec = log(dem$fec[i]),
                        l.pj = qlogis(0.1),
                        l.pa = qlogis(0.49),
                        alpha1 = -0.01, alpha2 = 0.11, alpha3 = 0.16, alpha4 = 0.05, # UNSCALED VALS
                        alpha5 = -0.82, alpha6 = 0.61,
                        alpha7 = -0.78, alpha8 = 0.08, alpha9 = 0.29)
  
  # build model
  peng.simModel2 <- nimbleModel(code = peng.ipm.imm.sim2,
                                constants = my.constants,
                                inits = initial.values)
  
  # identify nodes to simulate
  nodesToSim <- peng.simModel2$getDependencies(c("N1[1]", "Nimm[1]", "Nad[1]", "sigma.y", "J",
                                                 "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", 
                                                 "alpha6", "alpha7", "alpha8", "alpha9",
                                                 "l.mphij", "l.mphia", "l.mfec", "l.mim", "l.pj", "l.pa"),
                                               self = F, downstream = T)
  
  # simulate
  peng.simModel2$simulate(nodesToSim)
  
  pop[,j,i] <- peng.simModel2$Ntot
  
  # need to calculate lambda differently. rename 0s to 1s...
  # and only calculate lambda prior to extinction
  blah <- peng.simModel2$Ntot
  yay <- match(0, blah) # get index of 0 value
  if (is.na(yay)){
    geomean <- peng.simModel2$geomean.lambda
  } else{
    blah[blah == 0] <- 1
    bloop <- rep(NA, yay-1)
    for (t in 1:(yay-1)){bloop[t] <- blah[t+1]/blah[t]}
    l.bloop <- log(bloop)
    geomean <- exp((1/(yay-1))*sum(l.bloop[1:(yay-1)]))
  }
  
  lambda[i,j] <- geomean
  #lambda[i,j] <- peng.simModel2$geomean.lambda
  
}}

##################################################################################################################

#### 4) RUN SENSITIVITY ANALYSIS ON PRESS VS. PULSE
# NOTE: can't be sure about this, but I think I can run these all at once
# and then compare long-term change (mean, sd) to extreme vals (pct, mag)

nyears <- 100 # number of years to run stochastic sims
npars <- 8 # number of parameters to analyze
nparsim <- 5 # number of parameter combinations (i.e., unique draws from distributions)
nsims <- 5 # number of simulations per parameter combination

# parameters
clim <- data.frame(matrix(NA,nparsim,npars))
colnames(clim) <- c("rain_shape", "rain_rate", "temp_mean", "temp_sd",
                    "ssta_b_mean", "ssta_b_sd", "ssta_m_mean", "ssta_m_sd")

ext <- data.frame(matrix(NA,nparsim,npars))
colnames(ext) <- c("rain_freqECE", "rain_magECE", "temp_freqECE", "temp_magECE",
                   "ssta_b_freqECE", "ssta_b_magECE", "ssta_m_freqECE", "ssta_m_magECE")

# pop growth
pop <- array(data = NA, dim = c(nyears,nsims,nparsim))
lambda <- data.frame(matrix(NA,nparsim,nsims))

for (i in 1:nparsim){
  for (j in 1:nsims){
    
    clim$rain_shape[i] <- runif(1,1.3619,1.7859) # VARYING BY +/- Change in real data
    clim$rain_rate[i] <- runif(1,0.0183,0.0675)
    clim$temp_mean[i] <- runif(1,0.3744,0.3928)
    clim$temp_sd[i] <- runif(1,0.11076,0.1130097)
    clim$ssta_b_mean[i] <- runif(1,-0.2768,0.24317)
    clim$ssta_b_sd[i] <- runif(1,0.34338,0.61538)
    clim$ssta_m_mean[i] <- runif(1,-0.605,0.6025)
    clim$ssta_m_sd[i] <- runif(1,0.72651,0.84851)
    
    # calculate extremes
    # NOTE: should these be calculated off of the distribution, or from the local values actually...
    #... experienced by the penguins???
    rain_vals <- rgamma(10000, clim$rain_shape[i], clim$rain_rate[i])
    temp_vals <- rnorm(10000, clim$temp_mean[i], clim$temp_sd[i])
    ssta_b_vals <- rnorm(10000, clim$ssta_b_mean[i], clim$ssta_b_sd[i])
    ssta_m_vals <- rnorm(10000, clim$ssta_m_mean[i], clim$ssta_m_sd[i])
    
    # ECEs calculated from before based on distributions of data...95% quantiles
    ext$rain_freqECE[i] <- sum(rain_vals >= 93.8)/10000
    ext$rain_magECE[i] <- max(rain_vals)
    ext$temp_freqECE[i] <- sum(temp_vals >= 0.5676631)/10000
    ext$temp_magECE[i] <- max(temp_vals)
    ext$ssta_b_freqECE[i] <- sum(ssta_b_vals >= 0.7713634)/10000
    ext$ssta_b_magECE[i] <- max(ssta_b_vals)
    ext$ssta_m_freqECE[i] <- sum(ssta_m_vals >= 1.293878)/10000
    ext$ssta_m_magECE[i] <- max(ssta_m_vals)
  
    my.constants <- list(nyears = nyears, # UNSCALED
                         rain = as.vector(rgamma(nyears,clim$rain_shape[i],clim$rain_rate[i])),
                         temp = as.vector(rnorm(nyears,clim$temp_mean[i],clim$temp_sd[i])),
                         ssta_b = as.vector(rnorm(nyears,clim$ssta_b_mean[i],clim$ssta_b_sd[i])),
                         ssta_m = as.vector(rnorm(nyears,clim$ssta_m_mean[i],clim$ssta_m_sd[i]))) 
    my.constants$ssta_m_l <- as.vector(c(my.constants$ssta_m[2:nyears],rnorm(1,clim$ssta_m_mean[i],clim$ssta_m_sd[i])))
    
    # my.constants <- list(nyears = nyears, # UNSCALED
    #                      rain = as.vector(rgamma(nyears,1.57396674,0.04298969)),
    #                      temp = as.vector(rnorm(nyears,0.3836863,0.1118897)),
    #                      ssta_b = as.vector(rnorm(nyears,-0.01683185,0.47938401)),
    #                      ssta_m = as.vector(rnorm(nyears,-0.001235519,0.787513205))) 
    # my.constants$ssta_m_l <- as.vector(c(my.constants$ssta_m[2:nyears],rnorm(1,-0.001235519,0.787513205)))
    
    # inits = this is how you mess with the simulations!
    initial.values = list(N1 = 12000, 
                          Nimm = 10000,
                          Nad = 120000,
                          sigma.y = 22000,
                          # random draws from uniform dists
                          l.mphij = qlogis(0.1), 
                          l.mphia = qlogis(0.84),
                          l.mfec = log(0.51),
                          l.mim = log(0.12),
                          l.pj = qlogis(0.1),
                          l.pa = qlogis(0.49),
                          alpha1 = -0.01, alpha2 = 0.11, alpha3 = 0.16, alpha4 = 0.05, # UNSCALED VALS
                          alpha5 = -0.82, alpha6 = 0.61,
                          alpha7 = -0.78, alpha8 = 0.08, alpha9 = 0.29)
    
    # build model
    peng.simModel2 <- nimbleModel(code = peng.ipm.imm.sim2,
                                  constants = my.constants,
                                  inits = initial.values)
    
    # identify nodes to simulate
    nodesToSim <- peng.simModel2$getDependencies(c("N1[1]", "Nimm[1]", "Nad[1]", "sigma.y", "J",
                                                   "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", 
                                                   "alpha6", "alpha7", "alpha8", "alpha9",
                                                   "l.mphij", "l.mphia", "l.mfec", "l.mim", "l.pj", "l.pa"),
                                                 self = F, downstream = T)
    
    # simulate
    peng.simModel2$simulate(nodesToSim)
    
    pop[,j,i] <- peng.simModel2$Ntot
    
    # need to calculate lambda differently. rename 0s to 1s...
    # and only calculate lambda prior to extinction
    blah <- peng.simModel2$Ntot
    yay <- match(0, blah) # get index of 0 value
    if (is.na(yay)){
      geomean <- peng.simModel2$geomean.lambda
    } else{
      blah[blah == 0] <- 1
      bloop <- rep(NA, yay-1)
      for (t in 1:(yay-1)){bloop[t] <- blah[t+1]/blah[t]}
      l.bloop <- log(bloop)
      geomean <- exp((1/(yay-1))*sum(l.bloop[1:(yay-1)]))
    }
    
    lambda[i,j] <- geomean
    
  }}

