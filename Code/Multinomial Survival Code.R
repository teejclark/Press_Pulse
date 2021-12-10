
# set directory
rm(list = ls())
setwd("~/Post-Doc/Data/Gownaris et al. 2019/Gownaris&Boersma_DataS1")

library(nimble)
library(tidyverse)
library(MCMCvis)

##################################################################################################################

#### First Goal: Run an Example, but in Nimble
#### Simulated Data and Capture History: Age-Dependent Models ####

# create two m-arrays for those released as juveniles and as adults
# m-array for juveniles includes transition from j -> a; m-array for adults as before
# simulate data!

# Define parameter values
n.occasions <- 12                    # Number of capture occasions
marked.j <- rep(200, n.occasions-1)  # Annual number of newly marked juveniles
marked.a <- rep(30, n.occasions-1)   # Annual number of newly marked adults
phi.juv <- 0.3                       # Juvenile annual survival
phi.ad <- 0.65                       # Adult annual survival
p <- rep(0.5, n.occasions-1)         # Recapture
phi.j <- c(phi.juv, rep(phi.ad,n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI.J <- matrix(0, ncol = n.occasions-1, nrow = sum(marked.j))
for (i in 1:(length(marked.j)-1)){
  PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:(n.occasions-1)] <- matrix(rep(phi.j[1:(n.occasions-i)],marked.j[i]), ncol = n.occasions-i, byrow = TRUE)
}
P.J <- matrix(rep(p, n.occasions*sum(marked.j)), ncol = n.occasions-1, nrow = sum(marked.j), byrow = TRUE)
PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
    if (mark.occ[i]==n.occasions) next
    for (t in (mark.occ[i]+1):n.occasions){
      # Bernoulli trial: does individual survive occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      if (sur==0) break		# If dead, move to next individual
      # Bernoulli trial: is individual recaptured?
      rp <- rbinom(1, 1, P[i,t-1])
      if (rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
}

# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a)

CH <- rbind(CH.J, CH.A)
age <- c(rep(1,nrow(CH.J)), rep(2,nrow(CH.A))) # designate ages

# function to create basic m-arrays from single-state capture-recapture data!
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}


# function to create age-dependent m-arrays from single-state capture-recapture data!
# Input variables:
# ch = matrix with capture histories = includes inds of all age classes
# age = vector with age class at first capture for each individual
# mAge = maximal number of age classes for m-arrays constructed. Input is optional and only
# required if age matrix has fewer age classes then we want to separate
# (e.g., CH contains only individuals marked as juveniles, but we want 2 age classes)
# Output:
# 3-d array with the m-array. Third dimension is age class.

marray.age <- function(ch, age, mAge = 1){
  
  # 1. Helper functions
  # 1.1. Function to create a m-array based on capture-histories (ch)
  marray <- function(ch){
    ns <- length(table(ch)) - 1
    no <- ncol(ch)
    out <- matrix(0, ncol = ns*(no-1)+1, nrow = ns*(no-1))
    # Remove capture histories of individuals that are marked at last occasion
    get.first <- function(x) min(which(x!=0))
    first <- apply(ch, 1, get.first)
    last <- which(first==no)
    if (length(last) > 0) ch <- ch[-last,]
    # Compute m-array
    for (i in 1:nrow(ch)){
      cap.occ <- which(ch[i,]!=0)
      state <- ch[i,cap.occ]
      if (length(state) == 1) {
        out[state[1]+ns*(cap.occ[1]-1), ns*(no-1)+1] <- out[state[1]+ns*(cap.occ[1]-1), ns*(no-1)+1] + 1
      }
      if (length(state) > 1) {
        for (t in 2:length(cap.occ)){
          out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*ns+state[t]] <- out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*ns+state[t]] + 1
        } # t
        if (max(cap.occ) < no){
          out[(cap.occ[t]-1)*ns+state[t], ns*(no-1)+1] <- out[(cap.occ[t]-1)*ns+state[t], ns*(no-1)+1] + 1
        } # if
      } # if
    } # i
    return(out)
  }
  
  # 1.2. Function to remove histories without any capture from a capture-recapture matrix
  clean.ch <- function(ch){
    incl <- which(rowSums(ch)>=1)
    ch <- ch[incl,]
    return(ch)
  }
  
  # 1.3. Function to remove the first capture in a capture-recapture matrix
  rm.first <- function(ch) {
    get.first <- function(x) min(which(x==1))
    first <- apply(ch, 1, get.first)
    for (i in 1:nrow(ch)){
      ch[i,first[i]] <- 0
    }
    return(ch)
  }
  
  # 1.4. Function to calculate the occasion of first capture
  get.first <- function(x) min(which(x!=0))
  
  
  # 2. Calculations
  if (is.matrix(ch)==FALSE) ch <- matrix(ch, nrow = 1)
  maxAge <- max(c(max(age), mAge))
  nind <- nrow(ch)
  n.occasions <- ncol(ch)
  
  first <- apply(ch, 1, get.first)
  age.matrix <- matrix(0, ncol = n.occasions, nrow = nind)
  for (i in 1:nind){
    age.matrix[i,first[i]:n.occasions] <- 1:(n.occasions-first[i]+1)+(age[i]-1)
  }
  age.matrix[age.matrix > maxAge] <- maxAge
  
  # Recode capture history
  ch.rec <- ch
  for (i in 1:nind){
    h <- which(ch.rec[i,]==1)
    for (j in 1:length(h)){
      ch.rec[i,h[j]] <- j
    } # j
  } # i
  ch.rec[ch.rec > maxAge] <- maxAge
  
  ch.split <- array(0, dim = c(nrow(ch), ncol(ch), maxAge))
  for (a in 1:maxAge){
    for (i in 1:nind){
      j <- which(ch.rec[i,]==a | ch.rec[i,]==(a+1))
      if (length(j)==0) next
      ch.split[i,j[1:2],age.matrix[i,j[1]]] <- 1
      if (length(j)>1){
        ch.split[i,j[2:length(j)],age.matrix[i,j[2]]] <- 1
      }
    } # i
  } # a
  
  marr <- array(0, dim = c(n.occasions-1, n.occasions, maxAge))
  for (a in 1:(maxAge-1)){
    for (i in 1:nind){
      u <- which(ch.split[i,,a]==1)
      if (length(u)==0) next
      if (u[1]==n.occasions) next
      if (length(u)==1) marr[u,n.occasions,a] <- marr[u,n.occasions,a] + 1
      if (length(u)==2) marr[u[1],u[2]-1,a] <- marr[u[1],u[2]-1,a] + 1
    } # i
  } # a
  a <- maxAge
  
  if (is.matrix(ch.split[,,a])==FALSE){
    ch.split1 <- matrix(ch.split[,,a], nrow = 1)
    marr[,,a] <- marray(ch.split1)
  } # if
  else marr[,,a] <- marray(clean.ch(ch.split[,,a]))
  return(marr)
}

# get capture history
marray <- marray.age(CH, age)

##################################################################################################################

#### Running in Nimble: Age Dependent Model ####

sim.surv <- nimbleCode({
  
  # constraints
  for (t in 1:(nyears-1)){
    phi.juv[t] <- mean.phijuv
    phi.ad[t] <- mean.phiad
    p[t] <- mean.p
  }
  
  # priors
  mean.phijuv ~ dunif(0,1) # prior for mean juv. survival
  mean.phiad ~ dunif(0,1) # prior for mean ad. survival
  mean.p ~ dunif(0,1) # prior for mean recapture
  
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
})

# my data
my.data <- list(marr.j = marray[,,1],
                marr.a = marray[,,2],
                rel.j = rowSums(m.j),
                rel.a = rowSums(m.a))

# my constants
my.constants <- list(nyears = dim(marray)[2])

# inits
initial.values <- function(){list(mean.phijuv = runif(1,0,1),
                                  mean.phiad = runif(1,0,1),
                                  mean.p = runif(1,0,1))}

sim.surv.out <- nimbleMCMC(code = sim.surv,
                           constants = my.constants,
                           data = my.data,
                           inits = initial.values,
                           monitors = c("mean.phijuv", "mean.phiad", "mean.p"),
                           niter = 3000,
                           nburnin = 1000,
                           nchains = 2)
MCMCsummary(sim.surv.out, round = 2) # re-creates the data!!!!
MCMCtrace(sim.surv.out, pdf = F)


##################################################################################################################

#### Simple Penguin Survival Model ####
# this will only have juveniles (0-1) and adults (1+)

# load data
y <- as.matrix(read.table("Gownaris&Boersma_RecaptureHistories.txt", sep = ",", header = F))
colnames(y) <- NULL
head(y)

# convert everything to 1s and 0s. Assumption is that everyone turns to adult after 1 year, all start as juvs
y2 <- y
y2[y2 !=0] <- 1

# individuals marked as juveniles (all are)
CH.J2 <- matrix(as.numeric(y2), nrow(y2), ncol(y2))

# set-up
cap <- apply(CH.J2, 1, sum) # total # of captures
ind <- which(cap >=2) # which inds were re-captured
CH.J.R <- CH.J2[ind,] # Juvenile CH recaptured at least once
CH.J.N <- CH.J2[-ind,] # Juvenile CH never recaptured

# remove first capture of re-captured individuals (b/c they are adults with re-capture)
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
  first[i] <- min(which(CH.J.R[i,]==1))
}
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R1[i,first[i]] <- 0
}

# add grown-up juveniles to adults and create m-array (we don't have adults...)
CH.A.marray <- marray(CH.J.R1)

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


#####

# let's try looking at survival!!!!
peng.surv1 <- nimbleCode({
  
  # constraints
  for (t in 1:(nyears-1)){
    phi.juv[t] ~ dunif(0,1)#mean.phijuv
    phi.ad[t] ~ dunif(0,1)#mean.phiad
    #p[t] <- mean.p
  }
  
  # detection constraints
  # remove 2011, cause there was no sampling? (verify later)
  for (t in 1:29){p[t] <- mean.p}
  p[30] <- 0.00001 # OK - setting to 0 is bad, maybe this works?
  for (t in 31:(nyears-1)){p[t] <- mean.p}
  
  # priors
  #mean.phijuv ~ dunif(0,1) # prior for mean juv. survival
  #mean.phiad ~ dunif(0,1) # prior for mean ad. survival
  mean.p ~ dunif(0,1) # prior for mean recapture
  
  # define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr.j[t,1:nyears] ~ dmulti(pr.j[t,1:nyears], rel.j[t])
    marr.a[t,1:nyears] ~ dmulti(pr.a[t,1:nyears], rel.a[t])
  }
  
  # calculate # of birds released each year
  # for (t in 1:(nyears-1)){
  #   rel.j[t] <- sum(marr.j[t,1:nyears])
  #   rel.a[t] <- sum(marr.a[t,1:nyears])
  # }
  
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
})

# my data
my.data <- list(marr.j = CH.J.marray,
                marr.a = CH.A.marray,
                rel.j = rowSums(CH.J.marray),
                rel.a = rowSums(CH.A.marray))

# my constants
my.constants <- list(nyears = dim(CH.J.marray)[2])

# inits
initial.values <- function(){list(#mean.phijuv = runif(1,0,1),
                                  #mean.phiad = runif(1,0,1),
                                  phi.juv = runif(my.constants$nyears-1,0,1),
                                  phi.ad = runif(my.constants$nyears-1,0,1),
                                  mean.p = runif(1,0,1))}

peng.surv1.out <- nimbleMCMC(code = peng.surv1,
                           constants = my.constants,
                           data = my.data,
                           inits = initial.values,
                           monitors = c("phi.juv", "phi.ad", "mean.p"),
                           niter = 30000,
                           nburnin = 10000,
                           nchains = 3)
MCMCsummary(peng.surv1.out, round = 2) # re-creates the data!!!!
MCMCtrace(peng.surv1.out, pdf = F)
MCMCplot(peng.surv1.out)

