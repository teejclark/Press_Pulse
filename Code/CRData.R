# set directory
rm(list = ls())
setwd("~/Post-Doc/Data/Punta Tombo Data")

library(tidyverse)
library(lubridate)
library(patchwork)

# GOAL: manipulate CR data for inclusion into model; also confirm potential bugs in the data

##################################################################################################################

# 1) Capture-Recapture Data = 126,019 rows...
# length of unique PengSeqs = 66,461 penguins between 1982-2019
# Tasha said 44,374 birds were banded between 1983-2010
# GOAL: separate out birds banded as chicks with birds banded as adults, label them,
# then add them back together, so I know to split them up in analysis...

# NOTE: length of unique PengSeqs with only fledged = 30,259

#CRdat <- read.csv("PresenceAbsenceTOMwithageExport_9.20.21.csv") # with age data...ALL BANDED CHICKS, F and NF, nests and random
#CRdat <- read.csv("PresenceAbsenceTOMwithageExport_9.28.21.csv") # ONLY FLEDGED NEST CHICKS
#CRdat <- read.csv("PresenceAbsenceTOMwithageExport_nonfledged_9.28.21.csv") # ONLY UNFLEDGED NEST CHICKS
#CRdat <- read.csv("PresenceAbsenceTOMwithageExport_fledged_nestandrandom_9.28.21.csv") # FLEDGED NEST AND RANDOM CHICKS
CRdat <- read.csv("PresenceAbsenceTOMwithageExport_final_9.29.21.csv")

colnames(CRdat) <- c("sex", "bookyear", "pengid", "age")

# REMOVE DUPLICATES; ID = 58306...
# keep only distinct bookyears for each pengid; removes 94 rows..
CRdat <- CRdat %>% distinct(bookyear, pengid, .keep_all = T)

# extract chicks?
CRdat_chick <- CRdat %>% filter(age == "0") #48,538 chicks banded and followed...
# NOTE: only 9,724 here! (fledged nest chicks)
# ~1500 non-fledged nest chicks
# 45,359 fledged nest and random chicks

# and now merge later data...
CRdat_chick2 <- inner_join(CRdat_chick, CRdat, by = "pengid") %>%
  dplyr::select(-bookyear.x, -age.x, -sex.x, -age.y, -sex.y)

#NOTE: for now, remove duplicates that are goofy...
# CRdat_chick2 <- CRdat_chick2 %>%
#   #filter(pengid != 44487) %>% # FIX LATER FIX LATER FIX LATER!!!!!!
#   #filter(pengid != 47649) %>%
#   #filter(pengid != 57611) %>%
#   #filter(pengid != 57612) %>%
#   #filter(pengid != 65928) %>%
#   distinct(bookyear.y, pengid, .keep_all = T) # remove duplicates

CRdat2 <- CRdat_chick2 %>%
  mutate(blah = 1) %>% # adds observation for every year
  pivot_wider(names_from = bookyear.y, values_from = blah, names_sort = T) 

# change NA values to 0s
CRdat2[is.na(CRdat2)] <- 0

# now re-merge with sex data
CRdat3 <- inner_join(CRdat2, CRdat_chick, by = "pengid")

# remove duplicates
# identify duplicates
#View(CRdat3 %>% group_by(pengid) %>% summarize(n=n())) # see n = 2!
#blah <- CRdat3 %>% group_by(pengid) %>% filter(n()>1)

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

# check for duplicates
# NOTE: removed all juveniles and turned into adults...could change l8tr
#View(CRdat_adult2 %>% group_by(pengid) %>% summarize(n=n()))

# lastly, combine datasets (for now)
CRdat3_chick <- CRdat3 %>% dplyr::select(- bookyear) %>% mutate('1982' = 0)#, '2011' = 0)#; use for non-fledged
combo <- rbind(CRdat3_chick, CRdat_adult2) %>%
  dplyr::select(pengid, sex, age, '1982','1983','1984','1985','1986','1987','1988','1989','1990',
         '1991','1992','1993','1994','1995','1996','1997','1998','1999','2000','2001','2002',
         '2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014',
         '2015','2016','2017','2018','2019')

# NOTE: almost done, just need to fix the extra issues which will be done by Friday!

#### STATS
# banded as chicks: 695 Females, 1843 Males, 45,990 un-sexed (48528 total)
# banded as adults: 6297 Females, 7066 Males, 4548 un-sexed (17911 total)
# TOTAL COUNT: 66439 individual penguins

# WITH ONLY NEST CHECK FLEDGES
# banded as chicks -> fledges: 152 Females, 377 Males, 9195 un-sexed (9724 total)
# banded as adults: 6796 Females, 8517 Males, 5222 un-sexed (20535 total)

# WITH ONLY NEST CHECK AND RANDOM FLEDGES (FIXED)
# banded as chicks -> fledges: 677 Females, 1790 Males, 42892 un-sexed (45359 total)
# banded as adults: 6285 Females, 7119 Males, 4577 un-sexed (17981 total)
####

# LATER: split data into juvies, male adults, female adults!

##################################################################################################################

#### BUILD M-ARRAY DATA ####

# individuals marked as juveniles
CH.J <- combo %>% filter(age == 0) %>%
  dplyr::select(-pengid, -sex, -age) %>%
  as.matrix(.)

# individuals marked as adults
CH.A <- combo %>% filter(age == 2) %>%
  dplyr::select(-pengid, -sex, -age) %>%
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

#### TESTING SURVIVAL DATA ####
# NOTE: important to see differences b/t just banded vs. fledged chicks
# CURRENTLY: survival of just banded chicks

library(nimble)
library(MCMCvis)

peng.surv1 <- nimbleCode({
  
  # constraints
  for (t in 1:(nyears-1)){
    phi.juv[t] <-  mean.phijuv
    phi.ad[t] <-  mean.phiad
    p[t] <- mean.p
  }
  
  # detection constraints
  # remove 2011, cause there was no sampling? (verify later)
  # for (t in 1:29){p[t] <- mean.p}
  # p[30] <- 0.00001
  # for (t in 31:(nyears-1)){p[t] <- mean.p}
  
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
my.data <- list(marr.j = CH.J.marray,
                marr.a = CH.A.marray,
                rel.j = rowSums(CH.J.marray),
                rel.a = rowSums(CH.A.marray))

# my constants = 38 years of data!
my.constants <- list(nyears = dim(CH.J.marray)[2])

# inits
initial.values <- function(){list(
  mean.phijuv = runif(1,0,1),
  mean.phiad = runif(1,0,1),
  mean.p = runif(1,0,1))}

peng.surv1.out <- nimbleMCMC(code = peng.surv1,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values,
                             monitors = c("mean.phijuv", "mean.phiad", "mean.p"),
                             niter = 30000,
                             nburnin = 10000,
                             nchains = 3)
MCMCsummary(peng.surv1.out, round = 2) # AD = 84%; JUV = 8%
MCMCtrace(peng.surv1.out, pdf = F)
#MCMCplot(peng.surv1.out)

# try extracting data (from time-varying model)
phi.ad <- MCMCchains(object = peng.surv1.out,
                     params = "phi.ad",
                     mcmc.list = T)

# NOTE: currently years are labeled to represent non-breeding season
# so "1992" is the nonbreeding season of the "1991" cohort; and is really survival from 1991 -> 1992
MCMCplot(peng.surv1.out, params = "phi.ad",
         horiz = F,
         labels = c('1983','1984','1985','1986','1987','1988','1989','1990',
                  '1991','1992','1993','1994','1995','1996','1997','1998','1999','2000','2001','2002',
                  '2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014',
                  '2015','2016','2017','2018',"2019"),
         ylim = c(0,1),
         guide_lines=T)

phi.juv <- MCMCchains(object = peng.surv1.out,
                     params = "phi.juv",
                     mcmc.list = T)

MCMCplot(peng.surv1.out, params = "phi.juv",
         horiz = F,
         labels = c('1983','1984','1985','1986','1987','1988','1989','1990',
                    '1991','1992','1993','1994','1995','1996','1997','1998','1999','2000','2001','2002',
                    '2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014',
                    '2015','2016','2017','2018','2019'),
         ylim = c(0,0.6),
         guide_lines=T)  

# extract survival data b/c this is really annoying...
# adult survival!
phi.ad <- MCMCchains(peng.surv1.out, params = "phi.ad", mcmc.list = F)
phi.ad.med <- as.vector(apply(phi.ad,2,median))
phi.ad.CI <- apply(phi.ad,2,function (x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
phi.ad.data <- rbind(phi.ad.med,phi.ad.CI)
phi.ad.data <-tibble::rownames_to_column(data.frame(t(phi.ad.data)))
phi.ad.data$rowname <- c(seq(1983,2019,1))

# extract juvenile survival
phi.juv <- MCMCchains(peng.surv1.out, params = "phi.juv", mcmc.list = F)
phi.juv.med <- as.vector(apply(phi.juv,2,median))
phi.juv.CI <- apply(phi.juv,2,function (x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
phi.juv.data <- rbind(phi.juv.med,phi.juv.CI)
phi.juv.data <-tibble::rownames_to_column(data.frame(t(phi.juv.data)))
phi.juv.data$rowname <- c(seq(1983,2019,1))

# remove 2011 if you feel like it...just for visualization!!
phi.ad.data <- phi.ad.data[-30,]
phi.juv.data <- phi.juv.data[-30,]

# save data
survival.data <- data.frame(cbind(phi.ad.data, phi.juv.data))
survival.data$breedyear <- c(seq(1982,2010,1),seq(2012,2018,1)) # this is the breeding year, e.g., 2018 = 2018-2019
# the old year represents the year of the nonbreeding season, e.g., 2018 = 2017-2018
write.csv(survival.data,"survival_data_11.10.21.csv")

#plot
a <- ggplot(phi.ad.data, aes(x=rowname, y = phi.ad.med)) +
  geom_point(size = 3.5)+
  #geom_pointrange(aes(ymin = X25., ymax = X75.), size = 3.5, fatten = 0.5)+
  #geom_pointrange(aes(ymin=X2.5.,ymax = X97.5.), size = 1.15)+
  geom_line(color = "black", size = 1.05)+
  geom_ribbon(aes(ymin=X2.5., ymax = X97.5.), linetype = 2, alpha = 0.3, fill = "red")+
  ylim(0,1)+
  xlab("Year")+
  ylab("Adult Survival")+
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

#plot
b <- ggplot(phi.juv.data, aes(x=rowname, y = phi.juv.med)) +
  geom_point(size = 3.5)+
  #geom_pointrange(aes(ymin = X25., ymax = X75.), size = 3.5, fatten = 0.5)+
  #geom_pointrange(aes(ymin=X2.5.,ymax = X97.5.), size = 1.15)+
  geom_line(color = "black", size = 1.05)+
  geom_ribbon(aes(ymin=X2.5., ymax = X97.5.), linetype = 2, alpha = 0.3, fill = "red")+
  ylim(0,1)+
  xlab("Year")+
  ylab("Juvenile Survival")+
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

a+b

ggsave("overallsurvival_withonlyfledges_9.28.21.png", dpi = 600, width = 16)

# NOTES: does not seem to be a big difference between non-fledged and fledged survival? WHY?


#### SURVIVAL LABELING NOTES ####
# ABOVE TIME IS LIKELY WRONG. e.g, "2017" is the likelihood of surviving from 2017 -> 2018
# so there can't be any 2019, because we don't have 2020 data
# but there can be 1982, b/c it's the prob of surviving from 1982-1983. 
# inconsistency with Tasha's data is likely do to a mistake on her end?
# THIS IS TRUE = survival is always during nonbreeding season...
# so breeding year 1992's nonbreeding season is 1993, which is why estimates are weird!
################################

#### MISSING DATA NOTES ####
# NOTE: it might also be worthwhile to exclude 2011, because no data was collected then...or something else?
# one option - set p(2011) = 0 to indicate we would have never seen them!
################################

##################################################################################################################

#### FORM SEX-SPECIFIC SURVIVAL DATA (57979 inds)
# will have to censor the 4577 adults that were never sexed...
# and 784 recaptured banded juveniles that were never sexed... 
# will have to do a m-array for juveniles, adult females, adult males

#### BUILD M-ARRAY DATA ####
# NOTE: we are starting with only male and female adults, but maybe include M and F juveniles later?
# NOTE: randomly assign sex to the juvenile males and females...

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

# retain recaptured, sexed individuals (n = 2467)
combo4.5 <- combo3 %>%
  filter(age == 0, sex != "NULL") # this is recaptured, sexed individuals

# randomly assign sex to uncaptured, unsexed individuals (n = 42108)
combo5 <- combo3 %>%
  filter(age == 0, sum <= 1, sex =="NULL") %>% # uncaptured, unsexed individuals
  mutate(sex = sample(c("M","F"), n(), prob = c(0.5,0.5), replace = T)) # randomly assign sex

# retain uncaptured, sexed individuals (n = 30)
combo5.5 <- combo3 %>%
  filter(age == 0, sum <=1, sex != "NULL")

combo6 <- combo3 %>%
  filter(age == 2) %>%
  bind_rows(combo4) %>%
  bind_rows(combo4.5) %>%
  bind_rows(combo5) %>%
  bind_rows(combo5.5)

# individuals marked as juveniles, sexed as males
CH.JM <- combo6 %>% filter(age == 0) %>%
  filter(sex == "M") %>%
  select(-pengid, -sex, -age, -sum) %>%
  as.matrix(.)
# individuals marked as juveniles, sexed as females
CH.JF <- combo6 %>% filter(age == 0) %>%
  filter(sex == "F") %>%
  select(-pengid, -sex, -age, -sum) %>%
  as.matrix(.)
# individuals marked as juveniles, not sexed (NO RESIGHTS)
# CH.JN <- combo6 %>% filter(age == 0) %>%
#   filter(sex == "NULL") %>%
#   select(-pengid, -sex, -age, -sum) %>%
#   as.matrix(.)
# individuals marked as adults, sexed as males
CH.AM <- combo6 %>% filter(age == 2) %>%
  filter(sex == "M") %>%
  select(-pengid, -sex, -age, -sum) %>%
  as.matrix(.)
# individuals marked as adults, sexed as females
CH.AF <- combo6 %>% filter(age == 2) %>%
  filter(sex == "F") %>%
  select(-pengid, -sex, -age, -sum) %>%
  as.matrix(.)

# reformat into pure numerics
CH.JM2 <- matrix(as.numeric(CH.JM), nrow(CH.JM), ncol(CH.JM))
CH.JF2 <- matrix(as.numeric(CH.JF), nrow(CH.JF), ncol(CH.JF))
#CH.JN2 <- matrix(as.numeric(CH.JN), nrow(CH.JN), ncol(CH.JN))
CH.AM2 <- matrix(as.numeric(CH.AM), nrow(CH.AM), ncol(CH.AM))
CH.AF2 <- matrix(as.numeric(CH.AF), nrow(CH.AF), ncol(CH.AF))

# set-up
cap.JM <- apply(CH.JM2, 1, sum)
ind.JM <- which(cap.JM >= 2)
CH.JM.R <- CH.JM2[ind.JM,]
CH.JM.N <- CH.JM2[-ind.JM,]

cap.JF <- apply(CH.JF2, 1, sum)
ind.JF <- which(cap.JF >= 2)
CH.JF.R <- CH.JF2[ind.JF,]
CH.JF.N <- CH.JF2[-ind.JF,]

#CH.JN.N <- CH.JN2

# remove first capture of re-captured individuals (b/c they are adults with re-capture)
first <- numeric()
for (i in 1:dim(CH.JM.R)[1]){
  first[i] <- min(which(CH.JM.R[i,]==1))
}
CH.JM.R1 <- CH.JM.R
for (i in 1:dim(CH.JM.R)[1]){
  CH.JM.R1[i,first[i]] <- 0
}
first <- numeric()
for (i in 1:dim(CH.JF.R)[1]){
  first[i] <- min(which(CH.JF.R[i,]==1))
}
CH.JF.R1 <- CH.JF.R
for (i in 1:dim(CH.JF.R)[1]){
  CH.JF.R1[i,first[i]] <- 0
}

# add grown-up male juveniles to male adults and create m-array
CH.AM.m <- rbind(CH.AM2, CH.JM.R1)
CH.AM.marray <- marray(CH.AM.m)

# add grown-up female juveniles to female adults and create m-array
CH.AF.m <- rbind(CH.AF2, CH.JF.R1)
CH.AF.marray <- marray(CH.AF.m)

# create CH matrix for juveniles, ignoring subsequent recaptures (those who were recaptured)
second <- numeric()
for (i in 1:dim(CH.JM.R1)[1]){
  second[i] <- min(which(CH.JM.R1[i,]==1))
}
CH.JM.R2 <- matrix(0, nrow = dim(CH.JM.R)[1], ncol = dim(CH.JM.R)[2])
for (i in 1:dim(CH.JM.R)[1]){
  CH.JM.R2[i,first[i]] <- 1
  CH.JM.R2[i,second[i]] <- 1
}
second <- numeric()
for (i in 1:dim(CH.JF.R1)[1]){
  second[i] <- min(which(CH.JF.R1[i,]==1))
}
CH.JF.R2 <- matrix(0, nrow = dim(CH.JF.R)[1], ncol = dim(CH.JF.R)[2])
for (i in 1:dim(CH.JF.R)[1]){
  CH.JF.R2[i,first[i]] <- 1
  CH.JF.R2[i,second[i]] <- 1
}

# create m-array for these
CH.JM.R.marray <- marray(CH.JM.R2)
CH.JF.R.marray <- marray(CH.JF.R2)

# the last column ought to show # of juveniles not recaptured again,
# and should all be zeros, since all were released as adults
CH.JM.R.marray[,dim(CH.JM2)[2]] <- 0
CH.JF.R.marray[,dim(CH.JF2)[2]] <- 0

# create m-array for juvs never recaptured and add to previous array
CH.JM.N.marray <- marray(CH.JM.N)
CH.JM.marray <- CH.JM.R.marray + CH.JM.N.marray

CH.JF.N.marray <- marray(CH.JF.N)
CH.JF.marray <- CH.JF.R.marray + CH.JF.N.marray

##################################################################################################################

#### TESTING SURVIVAL DATA ####
# NOTE: important to see differences b/t just banded vs. fledged chicks
# CURRENTLY: survival of just banded chicks

library(nimble)
library(MCMCvis)

peng.surv1 <- nimbleCode({
  
  # constraints
  for (t in 1:(nyears-1)){
    phi.juvM[t] ~ dunif(0,1)#mean.phijuvM
    phi.juvF[t] ~ dunif(0,1)#mean.phijuvF
    phi.adM[t] ~ dunif(0,1)#mean.phiadM
    phi.adF[t] ~ dunif(0,1)#mean.phiadF
    
    phi.juv[t] <- (phi.juvM[t] + phi.juvF[t])/2 # try to calculate mean?
    #p[t] <- mean.p
  }
  
  # detection constraints - adults males
  for (t in 1:29){paM[t] <- mean.paM}
  paM[30] <- 0.00001 # OK - setting to 0 is bad, maybe this works?
  for (t in 31:(nyears-1)){paM[t] <- mean.paM}
  
  # detection constraints - adults females
  for (t in 1:29){paF[t] <- mean.paF}
  paF[30] <- 0.00001 # OK - setting to 0 is bad, maybe this works?
  for (t in 31:(nyears-1)){paF[t] <- mean.paF}
  
  # detection constraints - juveniles males
  # remove 2011, cause there was no sampling? (verify later)
  for (t in 1:29){pjM[t] <- mean.pjM}
  pjM[30] <- 0.00001 # OK - setting to 0 is bad, maybe this works?
  for (t in 31:(nyears-1)){pjM[t] <- mean.pjM}
  
  # detection constraints - juveniles females
  # remove 2011, cause there was no sampling? (verify later)
  for (t in 1:29){pjF[t] <- mean.pjF}
  pjF[30] <- 0.00001 # OK - setting to 0 is bad, maybe this works?
  for (t in 31:(nyears-1)){pjF[t] <- mean.pjF}
  
  # priors
  #mean.phijuvM ~ dunif(0,1) # prior for mean male juv. survival
  #mean.phijuvF ~ dunif(0,1) # prior for mean female juv. survival
  #mean.phiadM ~ dunif(0,1) # prior for mean ad. male survival
  #mean.phiadF ~ dunif(0,1) # prior for mean ad. female survival
  mean.paM ~ dunif(0,1) # prior for mean recapture, adults
  mean.paF ~ dunif(0,1)
  mean.pjM ~ dunif(0,1)
  mean.pjF ~ dunif(0,1)
  
  # define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr.jM[t,1:nyears] ~ dmulti(pr.jM[t,1:nyears], rel.jM[t])
    marr.jF[t,1:nyears] ~ dmulti(pr.jF[t,1:nyears], rel.jF[t])
    marr.aM[t,1:nyears] ~ dmulti(pr.aM[t,1:nyears], rel.aM[t])
    marr.aF[t,1:nyears] ~ dmulti(pr.aF[t,1:nyears], rel.aF[t])
  }
  
  # m-array cell probabilities for juveniles and adults
  for (t in 1:(nyears-1)){
    qaM[t] <- 1 - paM[t] # probability of recapture, adults
    qaF[t] <- 1 - paF[t]
    qjM[t] <- 1 - pjM[t] # probability of recapture, juveniles
    qjF[t] <- 1 - pjF[t]
    
    # main diagonal
    pr.jM[t,t] <- phi.juvM[t] * pjM[t]
    pr.jF[t,t] <- phi.juvF[t] * pjF[t]
    pr.aM[t,t] <- phi.adM[t] * paM[t]
    pr.aF[t,t] <- phi.adF[t] * paF[t]
    
    # above main diagonal
    for (j in (t+1):(nyears-1)){
      #pr.jM[t,j] <- phi.juvM[t] * prod(phi.adM[(t+1):j]) * prod(q[t:(j-1)]) * p[j] 
      #pr.jF[t,j] <- phi.juvF[t] * prod(phi.adF[(t+1):j]) * prod(q[t:(j-1)]) * p[j] 
      pr.jM[t,j] <- phi.juvM[t] * prod(phi.adM[(t+1):j]) * qjM[t] * prod(qaM[t:(j-1)]) * paM[j]/qaM[t]
      pr.jF[t,j] <- phi.juvF[t] * prod(phi.adF[(t+1):j]) * qjF[t] * prod(qaF[t:(j-1)]) * paF[j]/qaF[t]
      pr.aM[t,j] <- prod(phi.adM[t:j]) * prod(qaM[t:(j-1)]) * paM[j]
      pr.aF[t,j] <- prod(phi.adF[t:j]) * prod(qaF[t:(j-1)]) * paF[j]
    }
    
    # below main diagonal
    for (j in 1:(t-1)){
      pr.jM[t,j] <- 0
      pr.jF[t,j] <- 0
      pr.aM[t,j] <- 0
      pr.aF[t,j] <- 0
    } #j
  } #t
  
  for (t in 1:(nyears-1)){
    # last column: probability of non-recapture
    pr.jM[t,nyears] <- 1 - sum(pr.jM[t, 1:(nyears-1)])
    pr.jF[t,nyears] <- 1 - sum(pr.jF[t, 1:(nyears-1)])
    pr.aM[t,nyears] <- 1 - sum(pr.aM[t, 1:(nyears-1)])
    pr.aF[t,nyears] <- 1 - sum(pr.aF[t, 1:(nyears-1)])
  }
})

# my data
my.data <- list(marr.jM = CH.JM.marray,
                marr.jF = CH.JF.marray,
                marr.aM = CH.AM.marray,
                marr.aF = CH.AF.marray,
                rel.jM = rowSums(CH.JM.marray),
                rel.jF = rowSums(CH.JF.marray),
                rel.aM = rowSums(CH.AM.marray),
                rel.aF = rowSums(CH.AF.marray))

# my constants = 38 years of data!
my.constants <- list(nyears = dim(CH.J.marray)[2])

# inits
initial.values <- function(){list(
  #mean.phijuvM = runif(1,0,1),
  #mean.phijuvF = runif(1,0,1),
  #mean.phiadM = runif(1,0,1),
  #mean.phiadF = runif(1,0,1),
  phi.juvM = runif(my.constants$nyears-1,0,1),
  phi.juvF = runif(my.constants$nyears-1,0,1),
  phi.adM = runif(my.constants$nyears-1,0,1),
  phi.adF = runif(my.constants$nyears-1,0,1),
  mean.paM = runif(1,0,1),
  mean.paF = runif(1,0,1),
  mean.pjM = runif(1,0,1),
  mean.pjF = runif(1,0,1))}

peng.surv1.out <- nimbleMCMC(code = peng.surv1,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values,
                             #monitors = c("mean.phijuvM", "mean.phijuvF", "mean.phiadM", "mean.phiadF", "mean.p"),
                             monitors = c("phi.juv","phi.juvM", "phi.juvF", "phi.adM", "phi.adF", 
                                          "mean.paM", "mean.paF","mean.pjM", "mean.pjF"),
                             niter = 300000,
                             nburnin = 100000,
                             nchains = 3)
MCMCsummary(peng.surv1.out, round = 2) 

# extract survival data b/c this is really annoying...
# adult survival! = male
phi.adM <- MCMCchains(peng.surv1.out, params = "phi.adM", mcmc.list = F)
phi.adM.med <- as.vector(apply(phi.adM,2,median))
phi.adM.CI <- apply(phi.adM,2,function (x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
phi.adM.data <- rbind(phi.adM.med,phi.adM.CI)
phi.adM.data <-tibble::rownames_to_column(data.frame(t(phi.adM.data)))
phi.adM.data$rowname <- c(seq(1983,2019,1))

# adult survival! = female
phi.adF <- MCMCchains(peng.surv1.out, params = "phi.adF", mcmc.list = F)
phi.adF.med <- as.vector(apply(phi.adF,2,median))
phi.adF.CI <- apply(phi.adF,2,function (x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
phi.adF.data <- rbind(phi.adF.med,phi.adF.CI)
phi.adF.data <-tibble::rownames_to_column(data.frame(t(phi.adF.data)))
phi.adF.data$rowname <- c(seq(1983,2019,1))

# extract juvenile survival = male
phi.juvM <- MCMCchains(peng.surv1.out, params = "phi.juvM", mcmc.list = F)
phi.juvM.med <- as.vector(apply(phi.juvM,2,median))
phi.juvM.CI <- apply(phi.juvM,2,function (x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
phi.juvM.data <- rbind(phi.juvM.med,phi.juvM.CI)
phi.juvM.data <-tibble::rownames_to_column(data.frame(t(phi.juvM.data)))
phi.juvM.data$rowname <- c(seq(1983,2019,1))

# extract juvenile survival = female
phi.juvF <- MCMCchains(peng.surv1.out, params = "phi.juvF", mcmc.list = F)
phi.juvF.med <- as.vector(apply(phi.juvF,2,median))
phi.juvF.CI <- apply(phi.juvF,2,function (x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
phi.juvF.data <- rbind(phi.juvF.med,phi.juvF.CI)
phi.juvF.data <-tibble::rownames_to_column(data.frame(t(phi.juvF.data)))
phi.juvF.data$rowname <- c(seq(1983,2019,1))

# extract juvenile survival = averaged
phi.juv <- MCMCchains(peng.surv1.out, params = "phi.juv", mcmc.list = F)
phi.juv.med <- as.vector(apply(phi.juv,2,median))
phi.juv.CI <- apply(phi.juv,2,function (x) quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
phi.juv.data <- rbind(phi.juv.med,phi.juv.CI)
phi.juv.data <-tibble::rownames_to_column(data.frame(t(phi.juv.data)))
phi.juv.data$rowname <- c(seq(1983,2019,1))

# remove 2011 if you feel like it...just for visualization!!
phi.adM.data <- phi.adM.data[-30,]
phi.juvM.data <- phi.juvM.data[-30,]
phi.adF.data <- phi.adF.data[-30,]
phi.juvF.data <- phi.juvF.data[-30,]
phi.juv.data <- phi.juv.data[-30,]

#plot
a <- ggplot() +
  geom_line(aes(x=phi.adM.data$rowname, y = phi.adM.data$phi.adM.med), color = "red", size = 1.5)+
  geom_line(aes(x=phi.adF.data$rowname, y = phi.adF.data$phi.adF.med), color = "blue", size = 1.5)+
  geom_point(aes(x=phi.adM.data$rowname, y = phi.adM.data$phi.adM.med),size = 3.5, color = "red")+
  geom_point(aes(x=phi.adF.data$rowname, y = phi.adF.data$phi.adF.med),size = 3.5, color = "blue")+
  #geom_pointrange(aes(ymin = X25., ymax = X75.), size = 3.5, fatten = 0.5)+
  #geom_pointrange(aes(ymin=X2.5.,ymax = X97.5.), size = 1.15)+
  geom_ribbon(aes(x = phi.adM.data$rowname, ymin=phi.adM.data$X2.5., ymax = phi.adM.data$X97.5.), linetype = 2, alpha = 0.3, fill = "red")+
  geom_ribbon(aes(x = phi.adF.data$rowname, ymin=phi.adF.data$X2.5., ymax = phi.adF.data$X97.5.), linetype = 2, alpha = 0.3, fill = "blue")+
  ylim(0.5,1)+
  xlab("Year")+
  ylab("Adult Survival")+
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

# calculate and plot a sex-specific survival differential!!!!
diff <- data.frame(diff = phi.adM.data$phi.adM.med - phi.adF.data$phi.adF.med,
                   year = phi.adM.data$rowname)

c <- ggplot(diff, aes(x = year, y = diff)) + geom_point(size=3) + geom_line(size=1.25) + geom_smooth(method = "lm")+
  geom_hline(yintercept = 0, linetype = "dashed") +
    xlab("Year")+
    ylab("Sex-biased Adult Survival Differential")+
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

#plot
ggplot() +
  geom_line(aes(x=phi.juvM.data$rowname, y = phi.juvM.data$phi.juvM.med), color = "red", size = 1.5)+
  geom_line(aes(x=phi.juvF.data$rowname, y = phi.juvF.data$phi.juvF.med), color = "blue", size = 1.5)+
  geom_point(aes(x=phi.juvM.data$rowname, y = phi.juvM.data$phi.juvM.med),size = 3.5, color = "red")+
  geom_point(aes(x=phi.juvF.data$rowname, y = phi.juvF.data$phi.juvF.med),size = 3.5, color = "blue")+
  #geom_pointrange(aes(ymin = X25., ymax = X75.), size = 3.5, fatten = 0.5)+
  #geom_pointrange(aes(ymin=X2.5.,ymax = X97.5.), size = 1.15)+
  geom_ribbon(aes(x = phi.juvM.data$rowname, ymin=phi.juvM.data$X2.5., ymax = phi.juvM.data$X97.5.), linetype = 2, alpha = 0.3, fill = "red")+
  geom_ribbon(aes(x = phi.juvF.data$rowname, ymin=phi.juvF.data$X2.5., ymax = phi.juvF.data$X97.5.), linetype = 2, alpha = 0.3, fill = "blue")+
  ylim(0,1)+
  xlab("Year")+
  ylab("Juvenile Survival")+
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

#plot
ggplot(phi.juv.data, aes(x=rowname, y = phi.juv.med)) +
  geom_point(size = 3.5)+
  #geom_pointrange(aes(ymin = X25., ymax = X75.), size = 3.5, fatten = 0.5)+
  #geom_pointrange(aes(ymin=X2.5.,ymax = X97.5.), size = 1.15)+
  geom_line(color = "black", size = 1.05)+
  geom_ribbon(aes(ymin=X2.5., ymax = X97.5.), linetype = 2, alpha = 0.3, fill = "red")+
  ylim(0,1)+
  xlab("Year")+
  ylab("Juvenile Survival")+
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

(a | b)/c

ggsave("allsurvival_differential_10.7.21.png", dpi = 600, width = 16)


