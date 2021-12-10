# set directory
rm(list = ls())
setwd("~/Post-Doc/Data/Punta Tombo Data")

library(tidyverse)
library(lubridate)

# GOAL: manipulate other dem data for inclusion into model; also confirm potential bugs in the data

##################################################################################################################

#### Population Data ####

# NOTE (9/23/21): TWO ISSUES: A) Convert Population using 2012 Census; B) Create Conversion Factor from =>2016 dates

pop <- read.csv("StakeSurveyTOMACT22_9.21.21.csv", header = F)
pop[1,1] <- 1987 # fix weird stuff
colnames(pop) <- c("bookyear","00N01E", "00N02E", "00N03E", "00N04E", "00N05E", "00N06E", "00N07E", "00N08E",
                   "00N09E", "00N12E", "00N14E", "00N15E", "00N16E", "00N17E", "00N18E", "01N14E", "01S14E",
                   "02N14E", "02S14E", "03S14E", "04S14E", "05S14E", "total")
pop$density <- pop$total/22 # number of 100 meter survey areas = 22

ggplot(pop, aes(as.numeric(bookyear), density)) + geom_point(lwd=2) + geom_line(lwd=1.25) +
  geom_smooth(method = "lm")

# FIGURE OUT CONVERSION FACTOR

# load detailed survey data?
detail <- read.csv("StakeSurveyDetailed_9.23.21.csv")
colnames(detail) <- c("colony","bookyear","stakedate","staketime","stakenumber","substakenumber","bearing",
                      "stakeexception","barepct","grasspct","bushpct","bush1","bush2","bush3","livebush","deadbush",
                      "burrow","scrape","nestuse","nesttype","nestcover","nestcount","comment",
                      "stakegroup","stakeline","sortorder","entryby","entrydate","modby","moddate",
                      "stakenumseq","stakesurveyseq","stakenestseq","fieldentry","unknownbush")
detail$stakedate <- as.Date(detail$stakedate)

detail2 <- detail %>%
  mutate(year = year(stakedate)) %>%
  mutate(month = month(stakedate)) %>%
  mutate(day = day(stakedate)) %>%
  filter(month == "10") 

# sum over unique bookyear + stakesurveysequence
detail3 <- detail2 %>% 
  filter(bearing == "NULL") %>%
  filter(nestuse == "ACT" | nestuse == "ACT   ") %>%
  filter(stakeexception == "NULL" | stakeexception == "P" | stakeexception == "R") %>% # R is replicate!!!
  group_by(bookyear, stakenumber, stakedate) %>%
  summarize(total = sum(nestcount))

# OK - now let's keep just the early data
# extra surveys were done in 2017 and 2019. Started doing at the end of October in 2016
early <- detail3 %>%
  mutate(year = year(stakedate)) %>%
  mutate(month = month(stakedate)) %>%
  mutate(day = day(stakedate)) %>%
  filter(year == "1993"| year == "2017" | year == "2019") %>%
  filter(day < 15)

early_sum <- early %>%
  group_by(bookyear) %>%
  summarize(total = sum(total)/22)
# OK: so we had a sum of 5.77 (2017) and 5.36 (2019); 7.86 (1993)

# late data = this should match the density estimates from before
# extra surveys were done in 2017 and 2019. Started doing at the end of October in 2016
late <- detail3 %>%
  mutate(year = year(stakedate)) %>%
  mutate(month = month(stakedate)) %>%
  mutate(day = day(stakedate)) %>%
  filter(year == "1993" | year == "2017" | year == "2019") %>%
  filter(day > 15)

late_sum <- late %>%
  group_by(bookyear) %>%
  summarize(total = sum(total)/22)
# OK: so we had a sum of 4.14 (2017) and 4 (2019) = MATCHES; 7.64 (1993)

mean(c(early_sum$total/late_sum$total)) # conversion is about 1.36x higher with early counts

pop$blah <- c(pop$density[1:28]/1.36, pop$density[29:33])

ggplot(pop, aes()) +
  geom_point(aes(as.numeric(bookyear), density), lwd = 2) +
  geom_point(aes(as.numeric(bookyear), blah), lwd = 2, color = "red", alpha = 0.7)+
  geom_line(aes(as.numeric(bookyear), density),lwd=1.25)+
  geom_line(aes(as.numeric(bookyear), blah), lwd = 1.25, color = "red")+
  #geom_smooth(method="lm", aes(as.numeric(bookyear), blah))+
  ylim(c(2,10))+
  xlab("Year")+
  ylab("Mean Active Nests per 100 m2")

ggsave("nest_density_adjusted.png",dpi=600)

# convert to population estimate...
# estimates from Rebstock et al. (2016) Pop Ecol.
# 1987 = 314,000; 2012 (sep) = 209,300; 2014 = 201,000
# 2012: 209,300/35240 = 5.94 density; pretty close to fit line! use 35240!!!

pop$popest <- pop$density*35240 # 35000 is a guestimate based on the paper, clarify with Ginger!
pop$popest2 <- pop$blah*35240 # total habitat reported in the paper...(Table 1, 2012 survey)

ggplot(pop, aes(as.numeric(bookyear), popest)) + geom_point(lwd=2) + geom_line(lwd=1.25) +
  geom_smooth(method = "lm")
ggplot(pop, aes(as.numeric(bookyear), popest2)) + geom_point(lwd=2) + geom_line(lwd=1.25) +
  geom_smooth(method = "lm")+
  xlab("Year")+
  ylab("Total Active Nests")

ggsave("nest_total_adjusted.png", dpi = 600)

# ALTERNATE CONVERSION FACTOR?
blah <- early_sum
blah$late <- late_sum$total
colnames(blah) <- c("bookyear", "early", "late")
blah$convert <- blah$early/blah$late

ggplot(blah, aes(as.numeric(bookyear), convert)) + geom_point(lwd=2) + geom_smooth(method="lm")
summary(lm(convert ~ as.numeric(bookyear), data = blah))
# slope = 0.013344
boop <- (0.013344*seq(1,22,1))+1

pop$blah <- c(pop$density[1:6]/1.029, pop$density[7:28]/boop, pop$density[29:33])

ggplot(pop, aes()) +
  geom_point(aes(as.numeric(bookyear), density), lwd = 2) +
  geom_point(aes(as.numeric(bookyear), blah), lwd = 2, color = "red", alpha = 0.7)+
  geom_line(aes(as.numeric(bookyear), density),lwd=1.25)+
  geom_line(aes(as.numeric(bookyear), blah), lwd = 1.25, color = "red")+
  #geom_smooth(method="lm", aes(as.numeric(bookyear), blah))+
  ylim(c(2,10))+
  xlab("Year")+
  ylab("Mean Active Nests per 100 m2")

##################################################################################################################

#### REPRODUCTION DATA ####

# OK, here we have a sample they tracked every day
# which includes # of chicks hatched and fledged
# NOTE: might not match the # of chicks banded (which is when survival data starts)
# NOTE: so you could either change survival data (to when chicks hatch) OR change repro. data
# NOTE (9/23/21): Ginger suggested to stick with the RS data of fledged chicks/nest with eggs

repro <- read.csv("ReproSuccess_9.21.21.csv")
colnames(repro) <- c("bookyear", "numnests", "numeggs", "clutchsize", "numhatch", "meanhatch", "numfledged", "rs")

# number of chicks hatched per nest = 1.3 +/- 0.2 (sd)
ggplot(repro, aes(as.numeric(bookyear), meanhatch)) + geom_point(lwd=2) + geom_line(lwd=1.25) +
  geom_smooth(method = "lm")

# number of chicks fledged per nest = 0.52 +/- 0.24
ggplot(repro, aes(as.numeric(bookyear), rs)) + geom_point(lwd=2) + geom_line(lwd=1.25) +
  geom_smooth(method = "lm")





