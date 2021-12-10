
# set directory
rm(list = ls())
setwd("~/Post-Doc/Data/Gownaris et al. 2019/Gownaris&Boersma_DataS1")

library(nimble)
library(nimbleEcology)
library(tidyverse)
library(MCMCvis)

# load data
y <- as.matrix(read.table("Gownaris&Boersma_RecaptureHistories.txt", sep = ",", header = F))
colnames(y) <- NULL
head(y)

##################################################################################################################

#### SIMPLE VERSION: 0s and 1s ####

# first, convert all observations to 1s
y2 <- y
y2[y2 !=0] <- 1
y2 <- matrix(as.numeric(y2), nrow(y2), ncol(y2)) # convert to numbers

#y3 <- y2[1:5000,] # for now, just try first 5000 histories
y3 <- y2 # try all histories

peng.surv1 <- nimbleCode({
  phi ~ dunif(0,1) # survival
  p ~ dunif(0,1) # detection prob

  # transitions
  gamma[1,1] <- phi
  gamma[1,2] <- 1-phi
  gamma[2,1] <- 0
  gamma[2,2] <- 1
  omega[1,1] <- 1-p
  omega[1,2] <- p
  omega[2,1] <- 1
  omega[2,2] <- 0

  # first state propagation - NOTE: constant y
  for (i in 1:N){
    init[i, 1:2] <- gamma[const_y[i, first[i]]-1, 1:2]
  }

  # likelihood
  for (i in 1:N){
    y[i,(first[i]+1):K] ~ dHMMweighted(init = init[i,1:2],
                                       mult = mult[i],
                                       probObs = omega[1:2, 1:2],
                                       probTrans = gamma[1:2, 1:2],
                                       len = K - first[i],
                                       checkRowSums = 0)
  }

})

# pool individual encounter histories by unique histories
y_weighted <- y3 %>%
  as_tibble() %>%
  group_by_all() %>%
  summarize(mult = n()) %>%
  relocate(mult) %>%
  as.matrix()
mult <- y_weighted[,1] # number of individuals w/a particular encounter history
y_weighted2 <- y_weighted[,-1] # pooled data

# get first capture history
first <- apply(y_weighted2, 1, function(x) min(which(x != 0)))

# filter individuals that are captured at last occasion
mask <- which(first!=ncol(y_weighted2))
y_weighted2 <- y_weighted2[mask,]
first <- first[mask]
mult <- mult[mask]

# data
my.data <- list(y = y_weighted2 + 1)

# constants
my.constants <- list(first = first,
                     K = ncol(y_weighted2),
                     N = nrow(y_weighted2),
                     mult = mult,
                     const_y = y_weighted2 + 1)

# initial values
initial.values <- function(){list(phi = runif(1,0,1),
                                  p = runif(1,0,1))}

# run code
peng.surv1.out <- nimbleMCMC(code = peng.surv1,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values(),
                             monitors = c("phi", "p"),
                             niter = 10000,
                             nburnin = 5000,
                             nchains = 2)
MCMCsummary(peng.surv1.out, round = 2) # p = .16, phi = 0.66
MCMCtrace(peng.surv1.out, pdf = F)
MCMCplot(peng.surv1.out)

##################################################################################################################

#### SIMPLE VERSION: 0s and 1s ADD TIME####

# first, convert all observations to 1s
y2 <- y
y2[y2 !=0] <- 1
y2 <- matrix(as.numeric(y2), nrow(y2), ncol(y2)) # convert to numbers

y3 <- y2[1:500,] # for now, just try first 500 histories
y3 <- y2 # try all histories

peng.surv1.t <- nimbleCode({
  for(t in 1:(T-1)){
    phi[t] ~ dunif(0,1)
    gamma[1,1,t] <- phi[t]
    gamma[1,2,t] <- 1 - phi[t]
    gamma[2,1,t] <- 0
    gamma[2,2,t] <- 1
  }
  p ~ dunif(0,1)
  delta[1] <- 1         
  delta[2] <- 0          
  omega[1,1] <- 1 - p   
  omega[1,2] <- p       
  omega[2,1] <- 1       
  omega[2,2] <- 0      
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1]) # same as before, now survival is related to previous time-step
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# get first occasion of capture
first <- apply(y3,1,function(x)min(which(x!=0)))

# constants, data, inits
my.constants <- list(N = nrow(y3),
                     T = ncol(y3),
                     first = first)
my.data <- list(y = y3+1)
zinits = y3+1
zinits[zinits ==2] <- 1
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)

peng.surv1.out <- nimbleMCMC(code = peng.surv1.t,
                        constants = my.constants,
                        data = my.data,
                        inits = initial.values,
                        monitors = c("phi", "p"),
                        niter = 2500,
                        nburnin = 1000,
                        nchains = 2)
MCMCsummary(mcmc.phip, round = 2) # p = 0.9, phi = 0.56



##################################################################################################################

#### MULTI-STATE VERSION: Adults and Juveniles ####

# first, convert all uPs (juveniles) to 1s, everything else to 2s
y2 <- y
y2[y2 == "uP"] <- 1
y2[y2 == "FB"] <- 2
y2[y2 == "FN"] <- 2
y2[y2 == "Fu"] <- 2
y2[y2 == "MB"] <- 2
y2[y2 == "MN"] <- 2
y2[y2 == "Mu"] <- 2
y2[y2 == "uB"] <- 2
y2[y2 == "uN"] <- 2
y2[y2 == "uu"] <- 2

# first, convert all observations to 1s
y2 <- matrix(as.numeric(y2), nrow(y2), ncol(y2)) # convert to numbers

# multi-state with Adults and Juveniles (0 = not seen, 1 = juvs, 2 = adults)
peng.surv2 <- nimbleCode({

  phiJ ~ dunif(0,1)
  phiA ~ dunif(0,1)
  #psiJA ~ dunif (0,1) # transition rate from juvenile to adult (NOT NEEDED)
  pJ ~ dunif(0,1) # I assume this is ~1?
  pA ~ dunif(0,1)

  # transitions
  gamma[1,1] <- 0
  gamma[1,2] <- phiJ #* psiJA
  gamma[1,3] <- 1 - phiJ
  gamma[2,1] <- 0
  gamma[2,2] <- phiA #* (1 - psiJA)
  gamma[2,3] <- 1 - phiA
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- 1

  omega[1,1] <- 1 - pJ
  omega[1,2] <- pJ
  omega[1,3] <- 0
  omega[2,1] <- 1 - pA
  omega[2,2] <- 0
  omega[2,3] <- pA
  omega[3,1] <- 1
  omega[3,2] <- 0
  omega[3,3] <- 0

  # first state propagation - NOTE: constant y
  for (i in 1:N){
    init[i, 1:3] <- gamma[const_y[i, first[i]]-1, 1:3]
  }

  # likelihood
  for (i in 1:N){
    y[i,(first[i]+1):K] ~ dHMMweighted(init = init[i,1:3],
                                       mult = mult[i],
                                       probObs = omega[1:3, 1:3],
                                       probTrans = gamma[1:3, 1:3],
                                       len = K - first[i],
                                       checkRowSums = 0)
  }

})

# pool individual encounter histories by unique histories
y_weighted <- y2 %>%
  as_tibble() %>%
  group_by_all() %>%
  summarize(mult = n()) %>%
  relocate(mult) %>%
  as.matrix()
mult <- y_weighted[,1] # number of individuals w/a particular encounter history
y_weighted2 <- y_weighted[,-1] # pooled data

# get first capture history
first <- apply(y_weighted2, 1, function(x) min(which(x != 0)))

# filter individuals that are captured at last occasion
mask <- which(first!=ncol(y_weighted2))
y_weighted2 <- y_weighted2[mask,]
first <- first[mask]
mult <- mult[mask]

# data
my.data <- list(y = y_weighted2 + 1)

# constants
my.constants <- list(first = first,
                     K = ncol(y_weighted2),
                     N = nrow(y_weighted2),
                     mult = mult,
                     const_y = y_weighted2 + 1)

# initial values
initial.values <- function(){list(phiJ = runif(1,0,1),
                                  phiA = runif(1,0,1),
                                  #psiJA = runif(1,0,1),
                                  pJ = runif(1,0,1),
                                  pA = runif(1,0,1))}

# run code
peng.surv2.out <- nimbleMCMC(code = peng.surv2,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values(),
                             monitors = c("phiJ", "phiA",  "pJ", "pA"),
                             niter = 10000,
                             nburnin = 5000,
                             nchains = 2)
MCMCsummary(peng.surv1.out, round = 2) # phiA = 0.91, phiJ = 0.08 (seems about right!)
MCMCtrace(peng.surv1.out, pdf = F)

##################################################################################################################

#### TIME-VARYING MULTI-STATE VERSION: Adults and Juveniles ####

# first, convert all uPs (juveniles) to 1s, everything else to 2s
y2 <- y
y2 <- y2[1:1000,] # shorten time-series to get things to work
y2[y2 == "uP"] <- 1
y2[y2 == "FB"] <- 2
y2[y2 == "FN"] <- 2
y2[y2 == "Fu"] <- 2
y2[y2 == "MB"] <- 2
y2[y2 == "MN"] <- 2
y2[y2 == "Mu"] <- 2
y2[y2 == "uB"] <- 2
y2[y2 == "uN"] <- 2
y2[y2 == "uu"] <- 2

# first, convert all observations to 1s
y2 <- matrix(as.numeric(y2), nrow(y2), ncol(y2)) # convert to numbers

# multi-state with Adults and Juveniles (0 = not seen, 1 = juvs, 2 = adults)

peng.surv2.t <- nimbleCode({
  
  for (t in 1:(K-1)){ # loop over time
    logit(phiJ[t]) <- beta[1] # boring intercept
    logit(phiA[t]) <- beta[2] # boring intercept
    
    # time-varying transitions
    gamma[1,1,t] <- 0
    gamma[1,2,t] <- phiJ[t]
    gamma[1,3,t] <- 1 - phiJ[t]
    gamma[2,1,t] <- 0
    gamma[2,2,t] <- phiA[t]
    gamma[2,3,t] <- 1 - phiA[t]
    gamma[3,1,t] <- 0
    gamma[3,2,t] <- 0
    gamma[3,3,t] <- 1
  }
  
  # non-time varying observations
  omega[1,1] <- 1 - pJ
  omega[1,2] <- pJ
  omega[1,3] <- 0
  omega[2,1] <- 1 - pA
  omega[2,2] <- 0
  omega[2,3] <- pA
  omega[3,1] <- 1
  omega[3,2] <- 0
  omega[3,3] <- 0
  
  pJ ~ dunif(0,1)
  pA ~ dunif(0,1)
  beta[1] ~ dnorm(0, 1.5) # prior juvenile survival
  beta[2] ~ dnorm(0, 1.5) # prior adult survival
  
  # likelihood
  for (i in 1:N){
    for (j in 1:3)
      init[i,j] <- y_const[i, first[i]] == j
    #init[i, j] <-  gamma[y_const[i, first[i] ] - 1, 1:3]
  }
  
  for (i in 1:N){
    #for (j in (first[i]+1):K){
      y[i,(first[i]):K] ~ dDHMM(init = init[i, 1:3],
                                  probObs = omega[1:3, 1:3],
                                  probTrans = gamma[1:3, 1:3, first[i]:(K-1)],
                                  len = K - first[i] + 1,
                                  checkRowSums = 0)
    #}
    }
})



peng.surv2.t <- nimbleCode({
  delta[1] <- 1 # initial states - always start as juveniles
  delta[2] <- 0
  delta[3] <- 0
  
  for (t in 1:(K-1)){ # loop over time
    logit(phiJ[t]) <- beta[1] # boring intercept
    logit(phiA[t]) <- beta[2] # boring intercept
    
    # time-varying transitions
    gamma[1,1,t] <- 0
    gamma[1,2,t] <- phiJ[t]
    gamma[1,3,t] <- 1 - phiJ[t]
    gamma[2,1,t] <- 0
    gamma[2,2,t] <- phiA[t]
    gamma[2,3,t] <- 1 - phiA[t]
    gamma[3,1,t] <- 0
    gamma[3,2,t] <- 0
    gamma[3,3,t] <- 1
  }
  
  # non-time varying observations
  omega[1,1] <- 1 - pJ
  omega[1,2] <- pJ
  omega[1,3] <- 0
  omega[2,1] <- 1 - pA
  omega[2,2] <- 0
  omega[2,3] <- pA
  omega[3,1] <- 1
  omega[3,2] <- 0
  omega[3,3] <- 0
  
  pJ ~ dunif(0,1)
  pA ~ dunif(0,1)
  beta[1] ~ dnorm(0, 1.5) # prior juvenile survival
  beta[2] ~ dnorm(0, 1.5) # prior adult survival

  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:3])
    for (j in (first[i]+1):K){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:3, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:3])
    }
  }
})


# get first capture history
first <- apply(y2, 1, function(x) min(which(x != 0)))

# constants, data, inits
my.constants <- list(N = nrow(y2),
                     K = ncol(y2),
                     first = first)
                     #y_const = y2+1)
my.data <- list(y = y2+1)
zinits = y2+1
zinits[zinits ==2] <- 1
initial.values <- function() list(beta = rnorm(2,0,1),
                                  pJ = runif(1,0,1),
                                  pA = runif(1,0,1),
                                  z = zinits)

# run code
peng.surv2.t.out <- nimbleMCMC(code = peng.surv2.t,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values(),
                             monitors = c("phiJ", "phiA",  "pJ", "pA", "beta"),
                             niter = 2500,
                             nburnin = 500,
                             nchains = 2)
MCMCsummary(peng.surv2.t.out, round = 2) # 
MCMCtrace(peng.surv1.out, pdf = F)

##################################################################################################################


#### MULTI-STATE VERSION: Unknown Breeding Status, Adults and Pre-Breeders ####
# NOTE: messing with stuff here...

# first, convert
# 0 = not seen,
# 1 = unknown juvenile, 2 = male pre-breeder, 3 = female pre-breeder, 4 = unknown prebreeder,
# 5 = male breeder, # 6 = female breeder, #7 = unknown breeder
y2 <- y

y2[y2 == "uP"] <- 1 # unknown juvenile
y2[y2 == "FB"] <- 6 # female breeder
y2[y2 == "FN"] <- 6 # female breeder
y2[y2 == "MB"] <- 5 # male breeder
y2[y2 == "MN"] <- 5 # male breeder
y2[y2 == "uB"] <- 7 # unknown breeder
y2[y2 == "uN"] <- 7 # unknown breeder

first <- apply(y2, 1, function(x) min(which(x != 0))) # first sighting!

# if uu < 5 = unknown pre-breeder
for (i in 1:nrow(y2)){
    y2[i,(first[i]+which(y[i,first[i]:(first[i]+5)]=="uu")-1)] <- 4
}
# if Mu < 5 = male pre-breeder
for (i in 1:nrow(y2)){
  y2[i,(first[i]+which(y[i,first[i]:(first[i]+5)]=="Mu")-1)] <- 2
}
# if Fu < 5 = female pre-breeder
for (i in 1:nrow(y2)){
  y2[i,(first[i]+which(y[i,first[i]:(first[i]+5)]=="Fu")-1)] <- 3
}

# if uu > 5 = unknown breeder, if Mu > 5 = male breeder, if Fu > 5 = female breeder
y2[y2 == "uu"] <- 7
y2[y2 == "Fu"] <- 6
y2[y2 == "Mu"] <- 5

y2 <- matrix(as.numeric(y2), nrow(y2), ncol(y2)) # convert to numbers

# multistate code?!
peng.surv3 <- nimbleCode({

  # params
  phiMJ ~ dunif(0,1)
  phiFJ ~ dunif(0,1)
  phiMP ~ dunif(0,1)
  phiFP ~ dunif(0,1)
  phiMA ~ dunif(0,1)
  phiFA ~ dunif(0,1)
  psiPA ~ dunif(0,1)
  pP ~ dunif(0,1)
  pA ~ dunif(0,1)
  betaM ~ dunif(0,1)
  betaF ~ dunif(0,1)

  # transitions: 1 = MJ, 2 = FJ, 3 = MP, 4 = FP, 5 = MA, 6 = FA, 7 = Dead (7 x 7)
  gamma[1,1] <- 0
  gamma[1,2] <- 0
  gamma[1,3] <- phiMJ
  gamma[1,4] <- 0
  gamma[1,5] <- 0
  gamma[1,6] <- 0
  gamma[1,7] <- 1 - phiMJ
  gamma[2,1] <- 0
  gamma[2,2] <- 0
  gamma[2,3] <- 0
  gamma[2,4] <- phiFJ
  gamma[2,5] <- 0
  gamma[2,6] <- 0
  gamma[2,7] <- 1 - phiFJ
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- phiMP * (1 - psiPA)
  gamma[3,4] <- 0
  gamma[3,5] <- phiMP * psiPA
  gamma[3,6] <- 0
  gamma[3,7] <- 1 - phiMP
  gamma[4,1] <- 0
  gamma[4,2] <- 0
  gamma[4,3] <- 0
  gamma[4,4] <- phiFP * (1 - psiPA)
  gamma[4,5] <- 0
  gamma[4,6] <- phiFP * psiPA
  gamma[4,7] <- 1 - phiFP
  gamma[5,1] <- 0
  gamma[5,2] <- 0
  gamma[5,3] <- 0
  gamma[5,4] <- 0
  gamma[5,5] <- phiMA
  gamma[5,6] <- 0
  gamma[5,7] <- 1 - phiMA
  gamma[6,1] <- 0
  gamma[6,2] <- 0
  gamma[6,3] <- 0
  gamma[6,4] <- 0
  gamma[6,5] <- 0
  gamma[6,6] <- phiFA
  gamma[6,7] <- 1 - phiFA
  gamma[7,1] <- 0
  gamma[7,2] <- 0
  gamma[7,3] <- 0
  gamma[7,4] <- 0
  gamma[7,5] <- 0
  gamma[7,6] <- 0
  gamma[7,7] <- 1

  # observations  # transitioning from 7 states to 8 observations (7 x 8)
  # states: 1 = MJ, 2 = FJ, 3 = MP, 4 = FP, 5 = MA, 6 = FA, 7 = Dead
  # observations: 1 = NS, 2 = uJ, 3 = MP, 4 = FP, 5 = uP, 6 = MB, 7 = FB, 8 = uB
  omega[1,1] <- 0
  omega[1,2] <- 1
  omega[1,3] <- 0
  omega[1,4] <- 0
  omega[1,5] <- 0
  omega[1,6] <- 0
  omega[1,7] <- 0
  omega[1,8] <- 0
  omega[2,1] <- 1
  omega[2,2] <- 0
  omega[2,3] <- 0
  omega[2,4] <- 0
  omega[2,5] <- 0
  omega[2,6] <- 0
  omega[2,7] <- 0
  omega[2,8] <- 0
  omega[3,1] <- 1 - pP
  omega[3,2] <- 0
  omega[3,3] <- pP * betaM
  omega[3,4] <- 0
  omega[3,5] <- pP * (1 - betaM)
  omega[3,6] <- 0
  omega[3,7] <- 0
  omega[3,8] <- 0
  omega[4,1] <- 1 - pP
  omega[4,2] <- 0
  omega[4,3] <- 0
  omega[4,4] <- pP * betaF
  omega[4,5] <- pP * (1 - betaF)
  omega[4,6] <- 0
  omega[4,7] <- 0
  omega[4,8] <- 0
  omega[5,1] <- 1 - pA
  omega[5,2] <- 0
  omega[5,3] <- 0
  omega[5,4] <- 0
  omega[5,5] <- 0
  omega[5,6] <- pA * betaM
  omega[5,7] <- 0
  omega[5,8] <- pA * (1 - betaM)
  omega[6,1] <- 1 - pA
  omega[6,2] <- 0
  omega[6,3] <- 0
  omega[6,4] <- 0
  omega[6,5] <- 0
  omega[6,6] <- 0
  omega[6,7] <- pA * betaF
  omega[6,8] <- pA * (1 - betaF)
  omega[7,1] <- 1
  omega[7,2] <- 0
  omega[7,3] <- 0
  omega[7,4] <- 0
  omega[7,5] <- 0
  omega[7,6] <- 0
  omega[7,7] <- 0
  omega[7,8] <- 0

  # first state propagation - NOTE: constant y
  for (i in 1:N){
    init[i, 1:7] <- gamma[const_y[i, first[i]]-1, 1:7]
  }

  # likelihood
  for (i in 1:N){
    y[i,(first[i]+1):K] ~ dHMMweighted(init = init[i,1:7],
                                       mult = mult[i],
                                       probObs = omega[1:7, 1:8],
                                       probTrans = gamma[1:7, 1:7],
                                       len = K - first[i],
                                       checkRowSums = 0)
  }

})

# pool individual encounter histories by unique histories
y_weighted <- y2 %>%
  as_tibble() %>%
  group_by_all() %>%
  summarize(mult = n()) %>%
  relocate(mult) %>%
  as.matrix()
mult <- y_weighted[,1] # number of individuals w/a particular encounter history
y_weighted2 <- y_weighted[,-1] # pooled data

# get first capture history
first <- apply(y_weighted2, 1, function(x) min(which(x != 0)))

# filter individuals that are captured at last occasion
mask <- which(first!=ncol(y_weighted2))
y_weighted2 <- y_weighted2[mask,]
first <- first[mask]
mult <- mult[mask]

# data
my.data <- list(y = y_weighted2 + 1)

# constants
my.constants <- list(first = first,
                     K = ncol(y_weighted2),
                     N = nrow(y_weighted2),
                     mult = mult,
                     const_y = y_weighted2 + 1)

# initial values
initial.values <- function(){list(phiMJ = runif(1,0,1),
                                  phiFJ = runif(1,0,1),
                                  phiMP = runif(1,0,1),
                                  phiFP = runif(1,0,1),
                                  phiMA = runif(1,0,1),
                                  phiFA = runif(1,0,1),
                                  psiPA = runif(1,0,1),
                                  pP = runif(1,0,1),
                                  pA = runif(1,0,1),
                                  betaM = runif(1,0,1),
                                  betaF = runif(1,0,1))}

# run code
peng.surv3.out <- nimbleMCMC(code = peng.surv3,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values(),
                             monitors = c("phiMJ", "phiFJ", "phiMP", "phiFP", "phiMA", "phiFA", "psiPA", "pP", "pA", "betaM", "betaF"),
                             niter = 10000,
                             nburnin = 5000,
                             nchains = 2)
MCMCsummary(peng.surv3.out, round = 2)
MCMCtrace(peng.surv3.out, pdf = F)


##################################################################################################################

#### REALLY COMPLICATED VERSION ####
# NOTE: 8/26/21, trying out again! DOES NOT WORK!!!

# convert letters to numbers (0:7)
y[y == "uP"] <- 1
y[y == "MB"] <- 2
y[y == "FB"] <- 3
y[y == "uB"] <- 4
y[y == "MN"] <- 5
y[y == "FN"] <- 6
y[y == "uN"] <- 7
y[y == "Mu"] <- 8
y[y == "Fu"] <- 9
y[y == "uu"] <- 10
y[y == "0"] <- 11

# for now, shorten code to the first 5,000 obs?
y2 <- y
table(y2) # check to make sure all observed states are available (11 obvs)

# code
penguinsurvival <- nimbleCode({

  # -------------------------------------------------
  # Parameters: 18 total parameters!!
  # phiMP: survival probability state MP
  # phiFP: survival probability state FP
  # phiMB: survival probability state MB
  # phiFB: survival probability state FB
  # phiMN: survival probability state MN
  # phiFN: survival probability state FN
  # psiPB: transition probability from P to N
  # psiBN: transition probability from B to N
  # psiNB: transition probability from N to B
  # piMP: prob. of being in initial state MP
  # pP: recapture probability P
  # pB: recapture probability B
  # pN: recapture probability N
  # betaP: prob to ascertain the breeding status of an individual encountered as pre-breeder
  # betaB: prob to ascertain the breeding status of an individual encountered as breeder
  # betaN: prob to ascertain the breeding status of an individual encountered as non-breeder
  # betaM: prob to ascertain the sex of an individual encountered as male
  # betaF: prob to ascertain the sex of an individual encountered as female
  # -------------------------------------------------
  # States (z):
  # 1 male, pre-breeder (MP)
  # 2 female, pre-breeder (FP)
  # 3 male, breeder (MB)
  # 4 female, breeder (FB)
  # 5 male, nonbreeder (MN)
  # 6 female, nonbreeder (FN)
  # 7 dead
  # Observations (y):
  # 1 found, unknown pre-breeding
  # 2 found, male breeding
  # 3 found, female breeding
  # 4 found, unknown breeding
  # 5 found, male nonbreeding
  # 6 found, female nonbreeding
  # 7 found, unknown nonbreeding
  # 8 found, male unknown
  # 9 found, female unknown
  # 10 found, unknown unknown
  # 11 not found
  # -------------------------------------------------

  # priors
  phiMP ~ dunif(0,1)
  phiFP ~ dunif(0,1)
  phiMB ~ dunif(0,1)
  phiFB ~ dunif(0,1)
  phiMN ~ dunif(0,1)
  phiFN ~ dunif(0,1)
  psiPB ~ dunif(0,1)
  psiBN ~ dunif(0,1)
  psiNB ~ dunif(0,1)
  #piMP ~ dunif(0,1)
  pP ~ dunif(0,1)
  pB ~ dunif(0,1)
  pN ~ dunif(0,1)
  betaP ~ dunif(0,1)
  betaB ~ dunif(0,1)
  betaN ~ dunif(0,1)
  betaM ~ dunif(0,1)
  betaF ~ dunif(0,1)

  # vector of initial stats probs - have to start as MP or FP! (1x7)
  # delta[1] <- piMP      # prob. of being in initial state MP
  # delta[2] <- 1 - piMP  # prob. of being in initial state FP
  # delta[3] <- 0         # prob. of being in initial state MB
  # delta[4] <- 0         # prob. of being in initial state FB
  # delta[5] <- 0         # prob. of being in initial state MN
  # delta[6] <- 0         # prob. of being in initial state FN
  # delta[7] <- 0         # prob. of being in initial state dead

  # probabilities of state z(t+1) given z(t) = same as multistate model (7x7)
  gamma[1,1] <- phiMP * (1 - psiPB)
  gamma[1,2] <- 0
  gamma[1,3] <- phiMP * psiPB
  gamma[1,4] <- 0
  gamma[1,5] <- 0
  gamma[1,6] <- 0
  gamma[1,7] <- 1 - phiMP
  gamma[2,1] <- 0
  gamma[2,2] <- phiFP * (1 - psiPB)
  gamma[2,3] <- 0
  gamma[2,4] <- phiFP * psiPB
  gamma[2,5] <- 0
  gamma[2,6] <- 0
  gamma[2,7] <- 1 - phiFP
  gamma[3,1] <- 0
  gamma[3,2] <- 0
  gamma[3,3] <- phiMB * (1 - psiBN)
  gamma[3,4] <- 0
  gamma[3,5] <- phiMB * psiBN
  gamma[3,6] <- 0
  gamma[3,7] <- 1 - phiMB
  gamma[4,1] <- 0
  gamma[4,2] <- 0
  gamma[4,3] <- 0
  gamma[4,4] <- phiFB * (1 - psiBN)
  gamma[4,5] <- 0
  gamma[4,6] <- phiFB * psiBN
  gamma[4,7] <- 1 - phiFB
  gamma[5,1] <- 0
  gamma[5,2] <- 0
  gamma[5,3] <- phiMN * psiNB
  gamma[5,4] <- 0
  gamma[5,5] <- phiMN * (1 - psiNB)
  gamma[5,6] <- 0
  gamma[5,7] <- 1 - phiMN
  gamma[6,1] <- 0
  gamma[6,2] <- 0
  gamma[6,3] <- 0
  gamma[6,4] <- phiFN * psiNB
  gamma[6,5] <- 0
  gamma[6,6] <- phiFN * (1 - psiNB)
  gamma[6,7] <- 1 - phiFN
  gamma[7,1] <- 0
  gamma[7,2] <- 0
  gamma[7,3] <- 0
  gamma[7,4] <- 0
  gamma[7,5] <- 0
  gamma[7,6] <- 0
  gamma[7,7] <- 1

  # probabilities of y(t) given z(t) = new, multi-event model b/c beta (7x11)
  # NOTE: simplify code later (start with 0s matrix)...
  omega[1,1] <- pP * (1 - betaM) * betaP
  omega[1,2] <- 0
  omega[1,3] <- 0
  omega[1,4] <- 0
  omega[1,5] <- 0
  omega[1,6] <- 0
  omega[1,7] <- 0
  omega[1,8] <- pP * betaM * (1 - betaP)
  omega[1,9] <- 0
  omega[1,10] <- pP * (1 - betaM) * (1 - betaP)
  omega[1,11] <- 1 - pP
  omega[2,1] <- pP * (1 - betaF) * betaP
  omega[2,2] <- 0
  omega[2,3] <- 0
  omega[2,4] <- 0
  omega[2,5] <- 0
  omega[2,6] <- 0
  omega[2,7] <- 0
  omega[2,8] <- 0
  omega[2,9] <- pP * betaF * (1 - betaP)
  omega[2,10] <- pP * (1 - betaF) * (1 - betaP)
  omega[2,11] <- 1 - pP
  omega[3,1] <- 0
  omega[3,2] <- pB * betaM * betaB
  omega[3,3] <- 0
  omega[3,4] <- pB * (1 - betaM) * betaB
  omega[3,5] <- 0
  omega[3,6] <- 0
  omega[3,7] <- 0
  omega[3,8] <- pB * betaM * (1 - betaB)
  omega[3,9] <- 0
  omega[3,10] <- pB * (1 - betaM) * (1 - betaB)
  omega[3,11] <- 1 - pB
  omega[4,1] <- 0
  omega[4,2] <- 0
  omega[4,3] <- pB * betaF * betaB
  omega[4,4] <- pB * (1 - betaF) * betaB
  omega[4,5] <- 0
  omega[4,6] <- 0
  omega[4,7] <- 0
  omega[4,8] <- 0
  omega[4,9] <- pB * betaF * (1 - betaB)
  omega[4,10] <- pB * (1 - betaF) * (1 - betaB)
  omega[4,11] <- 1 - pB
  omega[5,1] <- 0
  omega[5,2] <- 0
  omega[5,3] <- 0
  omega[5,4] <- 0
  omega[5,5] <- pN * betaM * betaN
  omega[5,6] <- 0
  omega[5,7] <- pN * (1 - betaM) * betaN
  omega[5,8] <- pN * betaM * (1 - betaN)
  omega[5,9] <- 0
  omega[5,10] <- pN * (1 - betaM) * (1 - betaN)
  omega[5,11] <- 1 - pN
  omega[6,1] <- 0
  omega[6,2] <- 0
  omega[6,3] <- 0
  omega[6,4] <- 0
  omega[6,5] <- 0
  omega[6,6] <- pN * betaF * betaN
  omega[6,7] <- pN * (1 - betaF) * betaN
  omega[6,8] <- 0
  omega[6,9] <- pN * betaF * (1 - betaN)
  omega[6,10] <- pN * (1 - betaF) * (1 - betaN)
  omega[6,11] <- 1 - pN
  omega[7,1] <- 0
  omega[7,2] <- 0
  omega[7,3] <- 0
  omega[7,4] <- 0
  omega[7,5] <- 0
  omega[7,6] <- 0
  omega[7,7] <- 0
  omega[7,8] <- 0
  omega[7,9] <- 0
  omega[7,10] <- 0
  omega[7,11] <- 1

  # omega at t=1, where since we captured them, we HAD to detect them
  # ALL chicks are observed as "uP" -> only MP or FP, collpases most probabilities
  # omega.init[1,1] <- (1 - betaM) * betaP   # Pr(alive MP t = 1 -> detected uP t = 1)
  # omega.init[1,2] <- 0
  # omega.init[1,3] <- 0
  # omega.init[1,4] <- 0
  # omega.init[1,5] <- 0
  # omega.init[1,6] <- 0
  # omega.init[1,7] <- 0
  # omega.init[1,8] <- 0
  # omega.init[1,9] <- 0
  # omega.init[1,10] <- 0
  # omega.init[1,11] <- 0
  # omega.init[2,1] <- (1 - betaF) * betaP   # Pr(alive FP t = 1 -> detected uP t = 1)
  # omega.init[2,2] <- 0
  # omega.init[2,3] <- 0
  # omega.init[2,4] <- 0
  # omega.init[2,5] <- 0
  # omega.init[2,6] <- 0
  # omega.init[2,7] <- 0
  # omega.init[2,8] <- 0
  # omega.init[2,9] <- 0
  # omega.init[2,10] <- 0
  # omega.init[2,11] <- 0
  # omega.init[3,1] <- 0
  # omega.init[3,2] <- 0
  # omega.init[3,3] <- 0
  # omega.init[3,4] <- 0
  # omega.init[3,5] <- 0
  # omega.init[3,6] <- 0
  # omega.init[3,7] <- 0
  # omega.init[3,8] <- 0
  # omega.init[3,9] <- 0
  # omega.init[3,10] <- 0
  # omega.init[3,11] <- 0
  # omega.init[4,1] <- 0
  # omega.init[4,2] <- 0
  # omega.init[4,3] <- 0
  # omega.init[4,4] <- 0
  # omega.init[4,5] <- 0
  # omega.init[4,6] <- 0
  # omega.init[4,7] <- 0
  # omega.init[4,8] <- 0
  # omega.init[4,9] <- 0
  # omega.init[4,10] <- 0
  # omega.init[4,11] <- 0
  # omega.init[5,1] <- 0
  # omega.init[5,2] <- 0
  # omega.init[5,3] <- 0
  # omega.init[5,4] <- 0
  # omega.init[5,5] <- 0
  # omega.init[5,6] <- 0
  # omega.init[5,7] <- 0
  # omega.init[5,8] <- 0
  # omega.init[5,9] <- 0
  # omega.init[5,10] <- 0
  # omega.init[5,11] <- 0
  # omega.init[6,1] <- 0
  # omega.init[6,2] <- 0
  # omega.init[6,3] <- 0
  # omega.init[6,4] <- 0
  # omega.init[6,5] <- 0
  # omega.init[6,6] <- 0
  # omega.init[6,7] <- 0
  # omega.init[6,8] <- 0
  # omega.init[6,9] <- 0
  # omega.init[6,10] <- 0
  # omega.init[6,11] <- 0
  # omega.init[7,1] <- 0
  # omega.init[7,2] <- 0
  # omega.init[7,3] <- 0
  # omega.init[7,4] <- 0
  # omega.init[7,5] <- 0
  # omega.init[7,6] <- 0
  # omega.init[7,7] <- 0
  # omega.init[7,8] <- 0
  # omega.init[7,9] <- 0
  # omega.init[7,10] <- 0
  # omega.init[7,11] <- 0

  # initial state probs
  for (i in 1:N){
    init[i, 1:7] <- gamma[const_y[i, first[i]]-1, 1:7] # first state propagation
  }

  # likelihood
  for (i in 1:N){
    y[i, (first[i]+1):K] ~ dHMMweighted(init = init[i,1:7],
                                mult = mult[i],
                                probObs = omega[1:7,1:11],
                                probTrans = gamma[1:7,1:7],
                                len = K - first[i],
                                checkRowSums = 0)
  }

  # likelihood
  # for (i in 1:N){
  #   # latent state at first capture
  #   z[i,first[i]] ~ dcat(delta[1:7])
  #   y[i,first[i]] ~ dcat(omega.init[z[i,first[i]],1:11])
  #   for (t in (first[i]+1):K){
  #     # z(t) given z(t-1)
  #     z[i,t] ~ dcat(gamma[z[i,t-1],1:7])
  #     # y(t) given z(t)
  #     y[i,t] ~ dcat(omega[z[i,t],1:11])
  #  }
  #}
})

# pool individual encounter histories by unique histories
y_weighted <- y2 %>%
  as_tibble() %>%
  group_by_all() %>%
  summarize(mult = n()) %>%
  relocate(mult) %>%
  as.matrix()
mult <- y_weighted[,1] # number of individuals w/a particular encounter history
y_weighted2 <- y_weighted[,-1] # pooled data

# get first capture history
first <- apply(y_weighted2, 1, function(x) min(which(x != 11)))

# filter individuals that are captured at last occasion
mask <- which(first!=ncol(y_weighted2))
y_weighted2 <- y_weighted2[mask,]
first <- first[mask]
mult <- mult[mask]



# constants
my.constants <- list(first = first,
                     K = ncol(y_weighted2),
                     N = nrow(y_weighted2),
                     mult = mult,
                     const_y = y_weighted2)
my.data <- list(y = y_weighted2)


my.data <- list(y = matrix(as.numeric(y2), nrow(y2), ncol(y2)))


# initial values
initial.values <- function(){list(phiMP = runif(1,0,1),
                                  phiFP = runif(1,0,1),
                                  phiMB = runif(1,0,1),
                                  phiFB = runif(1,0,1),
                                  phiMN = runif(1,0,1),
                                  phiFN = runif(1,0,1),
                                  psiPB = runif(1,0,1),
                                  psiBN = runif(1,0,1),
                                  psiNB = runif(1,0,1),
                                  pP = runif(1,0,1),
                                  pB = runif(1,0,1),
                                  pN = runif(1,0,1),
                                  betaP = runif(1,0,1),
                                  betaB = runif(1,0,1),
                                  betaN = runif(1,0,1),
                                  betaM = runif(1,0,1),
                                  betaF = runif(1,0,1))}

# run code
mcmc.penguin.surv <- nimbleMCMC(code = penguinsurvival,
                                constants = my.constants,
                                data = my.data,
                                inits = initial.values,
                                monitors = c("phiMP", "phiFP", "phiMB", "phiFB", "phiMN", "phiFN",
                                             "psiPB", "psiBN", "psiNB",
                                             "pP", "pB", "pN",
                                             "betaP", "betaB", "betaN", "betaM", "betaF"),
                                niter = 10000,
                                nburnin = 5000,
                                nchains = 2)
MCMCsummary(mcmc.penguin.surv, round = 2)
MCMCtrace(mcmc.penguin.surv, pdf = F)


