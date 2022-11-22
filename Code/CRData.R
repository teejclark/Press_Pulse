# set directory
rm(list = ls())

library(tidyverse)
library(lubridate)
library(patchwork)

# GOAL: manipulate CR data for inclusion into model

##################################################################################################################

CRdat <- read.csv("PresenceAbsenceTOMwithageExport_final_9.29.21.csv")

colnames(CRdat) <- c("sex", "bookyear", "pengid", "age")

# REMOVE DUPLICATES; ID = 58306...
# keep only distinct bookyears for each pengid; removes 94 rows..
CRdat <- CRdat %>% distinct(bookyear, pengid, .keep_all = T)

# extract chicks?
CRdat_chick <- CRdat %>% filter(age == "0") 

# and now merge later data...
CRdat_chick2 <- inner_join(CRdat_chick, CRdat, by = "pengid") %>%
  dplyr::select(-bookyear.x, -age.x, -sex.x, -age.y, -sex.y)

CRdat2 <- CRdat_chick2 %>%
  mutate(blah = 1) %>% # adds observation for every year
  pivot_wider(names_from = bookyear.y, values_from = blah, names_sort = T) 

# change NA values to 0s
CRdat2[is.na(CRdat2)] <- 0

# now re-merge with sex data
CRdat3 <- inner_join(CRdat2, CRdat_chick, by = "pengid")

#### ADULT DATA
CRdat_adult <- CRdat %>% filter(age != "0")
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
  dplyr::select(pengid, sex, age, '1982','1983','1984','1985','1986','1987','1988','1989','1990',
         '1991','1992','1993','1994','1995','1996','1997','1998','1999','2000','2001','2002',
         '2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014',
         '2015','2016','2017','2018','2019')

