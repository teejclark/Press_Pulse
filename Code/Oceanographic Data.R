# set directory
rm(list = ls())
setwd("D:/Punta Tombo Oceanographic Data")

library(tidyverse)
library(lubridate)
library(raster)
library(stringr)
library(ggplot2)

# GOAL: manipulate data for inclusion into model; also confirm potential bugs in the data

##################################################################################################################

#### EXTRACT CORTAD FOR BREEDING ####
# Sea Surface Temps and Anomalies from 1982 to now
# could be used for both: a) breeding areas and b) migration areas
# NOTE: IMPORTANT: talk to Ginger before fully deciding on this vs. other ERDDAP databases
# stealing code from: https://github.com/brisneve/THE/blob/master/R/extract.timeseries.cordat.R

# CORTAD FOR BREEDING AREA
# CORE BREEDING AREA: 65° - 64° W, 43° - 44° S
# MAX BREEDING AREA: 65° - 63° W, 43° - 44° S (takes into account ~200 km, which is 95% of max foraging distances)
# NOTE: date is up in the air, but for now, let's go from hatching (~ Nov 20) to fledging (~ Feb 1)

# set coords
#coords <- extent(-65, -63, -44, -43) # xmin, xmax, ymin, ymax
# GINGER recommends going from -45 to -40 and -65 to -61
coords <- extent(-65.5, -61, -45, -40) # changed from -65 to -65.5

# ref datasets
#sst <- "cortadv6_WeeklySST_2.nc"
ssta <- "cortadv6_SSTA_2.nc"
#tsa <- "cortadv6_TSA.nc"

# first crop the data
#sst.brick <- brick(sst)
ssta.brick <- brick(ssta)
#tsa.brick <- brick(tsa)
#x <- crop(sst.brick, coords)
y <- crop(ssta.brick, coords)
#z <- crop(tsa.brick, coords)

# average over space and time-period of foraging
# NOTE: temp is in kelvin, if using raw SST, be sure to subract by 273.15
numyears <- seq(1982,2020,1)

# extract specific months - 11, 12, 01, 02 - THIS IS REALLY UGLY
all_yrs <- str_subset(names(y), c("\\.01\\.|\\.02\\.|\\.10\\.|\\.11\\.|\\.12\\."))

#x.mean <- rep(NA,length(numyears))
y.mean <- rep(NA,length(numyears))
#z.mean <- rep(NA,length(numyears))

blah <- data.frame(Oct = paste(numyears, c(".10"), sep = ""),
                   Nov = paste(numyears, c(".11"), sep = ""),
                   Dec = paste(numyears, c(".12"), sep = ""),
                   Jan = paste(numyears+1, c(".01"), sep = ""),
                   Feb = paste(numyears+1, c(".02"), sep = ""))

# run for loop
# for now, just calculating mean value across months, then across space
bloop.y <- list()

for (i in 1:length(numyears)){
  oct <- str_subset(all_yrs, as.character(blah$Oct[i]))
  nov <- str_subset(all_yrs, as.character(blah$Nov[i]))
  dec <- str_subset(all_yrs, as.character(blah$Dec[i]))
  jan <- str_subset(all_yrs, as.character(blah$Jan[i]))
  feb <- str_subset(all_yrs, as.character(blah$Feb[i]))
  month_stack <- c(oct, nov, dec, jan, feb)
  #bloop.x <- subset(x, month_stack)
  bloop.y[[i]] <- subset(y, month_stack)
  #bloop.z <- subset(z, month_stack)
  #x.mean[i] <- cellStats(mean(bloop.x), mean)-273.15 # SST in Kelvin
  #y.mean[i] <- cellStats(mean(bloop.y), mean)
  #z.mean[i] <- cellStats(mean(bloop.z), mean)
}

bloop.y <- stack(bloop.y)

# calculate median(?) of each raster
y.med <- rep(NA, length(names(bloop.y)))

for (i in 1:length(y.med)){
  y.med[i] <- cellStats(bloop.y[[i]], median)
}

sst.dat <- data.frame(date = as.Date(str_replace_all(names(bloop.y), "X", ""), "%Y.%m.%d"),
                      SSTA = y.med,
                      SSTA_s = scale(y.med))

ggplot(sst.dat, aes(date, SSTA_s)) + geom_point() + geom_smooth(method = "lm")

# average over data
sst.dat$month <- month(sst.dat$date)
sst.dat$year <- year(sst.dat$date)
sst.dat$bookyear <- ifelse(sst.dat$month == "1" | sst.dat$month == "2",
                           sst.dat$year-1, sst.dat$year)
sst.dat <- sst.dat[1:822,] # data in bookyear 2020 doesn't have full data
sst.dat.av <- sst.dat %>% group_by(bookyear) %>% summarize(av = mean(SSTA))
sst.dat.av$av_s <- scale(sst.dat.av$av)


ggplot(sst.dat.av, aes(bookyear, av_s)) + geom_point() + geom_smooth(method = "lm")

write.csv(sst.dat.av, "ssta_summarized.12.2.21.csv")

# density plot
sst.dat$split <- c(rep("1982-2000",402),rep("2001-2020",420))

cols <- c("#CA9502","#FDC835")
ggplot(sst.dat, aes(x=SSTA, color = split, linetype = split))+
  geom_density(alpha=0.4, size = 4, adjust = 2.5)+
  scale_color_manual(values= cols)+
  geom_vline(xintercept = 0.9485, linetype = "dashed", size = 1.5)+
  xlim(-3, 3)+
  xlab("SSTA - Breeding")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=32),
    axis.text = element_text(size=38, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank(),
    legend.position = "none")

ggsave("sstabreeding_densityplot_11.16.21.png", dpi = 600, width = 10, height = 7)

old <- sst.dat %>% filter(split == "1982-2000")
new <- sst.dat %>% filter(split == "2001-2020")

# fit norm distributions
fit_old <- fitdist(old$SSTA, "norm") # -0.1508, 0.801
fit_new <- fitdist(new$SSTA, "norm") # 0.1113, 0.665

# sst.dat <- data.frame(year = seq(1982,2020,1),
#                       SST = x.mean,
#                       SSTA = y.mean,
#                       TSA = z.mean)

# ggplot(sst.dat, aes(year, TSA)) + geom_point() + geom_smooth(method = "lm")
# # all vals show that ocean temperature is increasing
# 
# ggplot(sst.dat, aes(year, SST)) + geom_point() + geom_line() + ggtitle("SST - Breeding Region")
# ggplot(sst.dat, aes(year, SSTA)) + geom_point() + geom_line() + ggtitle("SSTA - Breeding Region")
# ggplot(sst.dat, aes(year, TSA)) + geom_point() + geom_line() + ggtitle("TSA - Breeding Region")

##################################################################################################################

#### PCA FOR BREEDING ####
# for now, work with sea surface temperature anomalies from Nov-Feb, can change later...
# NOTE: will verify this by re-doing PPI

# extract data for only certain months!
sst_stack <- list()
ssta_stack <- list()
#tsa_stack <- list()

for (i in 1:length(numyears)){
  
  nov <- str_subset(all_yrs, as.character(blah$Nov[i]))
  dec <- str_subset(all_yrs, as.character(blah$Dec[i]))
  jan <- str_subset(all_yrs, as.character(blah$Jan[i]))
  feb <- str_subset(all_yrs, as.character(blah$Feb[i]))
  month_stack <- c(nov, dec, jan, feb)
  #sst_stack[[i]] <- subset(x, month_stack)
  ssta_stack[[i]] <- subset(y, month_stack)
  #tsa_stack[[i]] <- subset(z, month_stack)
}

#sst_stack <- stack(sst_stack)
ssta_stack <- stack(ssta_stack)
#tsa_stack <- stack(tsa_stack)

# run pca - temporal patterns!
# NOTE: https://stackoverflow.com/questions/41022927/principal-component-analysis-pca-of-time-series-data-spatial-and-temporal-pat
# here, we are using the full dataset and all available weeks!
# "temporal patterns" = loadings of principal components (time-series)

# create data.frame, with rows as space, time as columns (use y)
weeks <- seq(1, nlayers(y), 1)
ssta_dat <- data.frame(matrix(NA, nrow = ncell(y), ncol = nlayers (y)+2))
ssta_dat[,1] <- as.data.frame(y[[1]], xy = T)$x
ssta_dat[,2] <- as.data.frame(y[[1]], xy = T)$y
colnames(ssta_dat) <- c("x","y",names(y))

for (i in 1:length(weeks)){
  ssta_dat[,i+2] <- as.vector(y[[i]]) # gives an error, but seems to work?
}

# extract certain months prior to PCA - Oct - Feb
blah <- as.data.frame(t(ssta_dat[3:ncol(ssta_dat)]))
blah$time <- as.Date(str_replace_all(names(ssta_dat[3:ncol(ssta_dat)]), "X", ""), "%Y.%m.%d")
blah$month <- month(blah$time)
blah2 <- blah %>%
  filter(month == "10" | month == "11" | month == "12" | month == "1" | month == "2")

ssta_datb <- as.data.frame(t(blah2[1:(ncol(blah2)-2)]))

# Take out chunks of NAs due to land???
#ssta_dat2 <- ssta_dat[3:ncol(ssta_dat)]
#ssta_dat2 <- ssta_dat2[rowSums(is.na(ssta_dat2)) != ncol(sst a_dat2),] # took out 50 NA rows
ssta_dat2 <- ssta_datb[rowSums(is.na(ssta_datb)) != ncol(ssta_datb),]

# run pca, exclude lat and lon
PCA_ssta <- prcomp(t(ssta_dat2), scale = F, center = F)

summaryPCA_ssta <- summary(PCA_ssta)
summaryPCA_ssta$importance[,c(1,2,3)] # 64.5,3.7 (Oct-Feb)

# extract loadings - should be temporal patterns
loading.ssta <- data.frame(time = as.Date(str_replace_all(names(ssta_dat2), "X", ""), "%Y.%m.%d"),
                          loadingPC1 = scale(as.vector(PCA_ssta$x[,'PC1'])),
                          loadingPC2 = scale(as.vector(PCA_ssta$x[,'PC2'])))

# plot temporal patterns over time
ggplot(loading.ssta, aes(x = time, y = loadingPC1)) + 
  #geom_point() + 
  geom_line() + 
  geom_smooth(method = "lm")+
  ggtitle("PCA of SSTA - Breeding")

#ggsave("ssta_breeding_temporal_overall_101821.png", dpi = 600)

# density plot
loading.ssta$split <- c(rep("1982-2000",991),rep("2001-2020",1043))

old <- loading.ssta %>% filter(split == "1982-2000")
new <- loading.ssta %>% filter(split == "2001-2020")

# fit gamma distributions
fit_old <- fitdist(old$loadingPC1, "norm") # -0.3098, 0.905
fit_new <- fitdist(new$loadingPC1, "norm") # 0.294, 0.996

cols <- c("#CA9502","#FDC835")
ggplot(loading.ssta, aes(x=loadingPC1, color = split, linetype = split))+
  geom_density(alpha=0.4, size = 4, adjust = 2.5)+
  scale_color_manual(values= cols)+
  geom_vline(xintercept = 1.346948, linetype = "dashed", size = 1.5)+
  ylim(0.0, 0.4)+
  xlab("SSTA - Breeding")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=32),
    axis.text = element_text(size=38, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank(),
    legend.position = "none")

ggsave("sstabreeding_densityplot_11.16.21.png", dpi = 600, width = 10, height = 7)

# match up averaged spatial patterns
# early = Nov-Dec (10-12)
# late = Jan-Feb (1-2)

loading.ssta$month <- month(loading.ssta$time)
loading.ssta$year <- year(loading.ssta$time)
loading.ssta$bookyear <- ifelse(loading.ssta$month == "1" | loading.ssta$month == "2",
                                loading.ssta$year-1, loading.ssta$year)
loading.ssta <- loading.ssta[9:830,]

loading.ssta.early <- loading.ssta %>%
  filter(month == "10" | month == "11" | month == "12") %>%
  group_by(bookyear) %>%
  summarize(mean = mean(loadingPC1),
            sd = sd(loadingPC1))

loading.ssta.late <- loading.ssta %>%
  filter(month == "1" | month == "2") %>%
  group_by(bookyear) %>%
  summarize(mean = mean(loadingPC1),
            sd = sd(loadingPC1))

loading.ssta.overall <- loading.ssta %>%
  filter(month == "10" | month == "11" | month == "12" | month == "1" | month == "2") %>%
  group_by(bookyear) %>%
  summarize(mean = mean(loadingPC1),
            sd = sd(loadingPC1))

loading.ssta.both <- data.frame(cbind(loading.ssta.early,loading.ssta.late, loading.ssta.overall))

#write.csv(loading.ssta.both, "ssta.PC2.breeding.11.15.21.csv")

# plot temporal patterns
ggplot(loading.ssta.both, aes(x = bookyear, y = mean)) +
  geom_line(lwd = 1.15)+
  geom_line(aes(x = bookyear, y = mean.1), lwd = 1.15, color = "red")+
  #geom_line(aes(x = year, y = mean.2), lwd = 1.15, color = "blue")+
  geom_smooth(method = "lm")+
  ggtitle("PCA of SSTA - Breeding; Black = Early, Red = Late")

#ggsave("ssta_breeding_temporal_average_101821.png", dpi = 600)

# dem vars
# repro (would affect the year of the bookyear, e.g., PCA of 1994 -> affect reproduction in the breeding season (1994))
ssta.repro <- loading.ssta.both[2:38,]
ssta.repro <- ssta.repro[-29,]
ssta.repro$repro <- repro$rs
ggplot(ssta.repro, aes(mean, repro)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.repro, aes(sd.2, repro)) + geom_point() + geom_smooth(method = "lm")

# survival (don't need to adjust year, b/c the data is set to nonbreeding year)
survival.data <- read.csv("survival_data_11.10.21.csv")
ssta.surv <- loading.ssta.both[2:37,]
ssta.surv <- ssta.surv[-29,]
ssta.surv$asurv <- survival.data$phi.ad.med[2:36]
ssta.surv$jsurv <- survival.data$phi.juv.med[2:36]
ggplot(ssta.surv, aes(mean.2, asurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(mean.2, jsurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(sd.2, asurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(sd.2, jsurv)) + geom_point() + geom_smooth(method = "lm")


# run pca - spatial patterns!
blah <- na.omit(ssta_dat)
PC1.ssta <- data.frame(x = blah$x, y = blah$y, val = scale(PCA_ssta$rotation[,'PC1']))

PC1.ssta.raster <- rasterFromXYZ(PC1.ssta)

# plot
palette <- colorRampPalette(colors = c("blue","turquoise","green","yellow","orange","red"))
plot(PC1.ssta.raster, col = palette(100), zlim = c(-2,2))


##################################################################################################################

#### EXTRACT CORTAD FOR MIGRATION ####
# https://cran.r-project.org/web/packages/marmap/vignettes/marmap-DataAnalysis.pdf

library(marmap)

# extract bathymetric dataset!!!
bat <- getNOAA.bathy(lon1 = -65, lon2 = -40, lat1 = -40, lat2 = -23, res = 4)

# plot
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))
# Plot
plot(bat, image = TRUE, land = TRUE, lwd = 0.1, bpal = list(c(0, max(bat), greys), c(min(bat), 0, blues)))
plot(bat, lwd = 0.8, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline

# 3d plot?
lattice::wireframe(unclass(bat), shade = T, aspect = c(1/2,0.1))

# extract area
area <- get.area(bat, level.inf = -1000)

plot(bat, lwd = 0.2)
col.area <- rgb(0.7,0,0,0.3)
plotArea(area, col = col.area)

# convert to raster
bat.raster <- marmap::as.raster(bat)

# get just -1000 to 0 elevation
bat.raster[bat.raster > 0,] <- NA
bat.raster[bat.raster < -1000,] <- NA
plot(bat.raster)

writeRaster(bat.raster, "bathraster.tif")

# extract CORTAD data
# set coords
coords <- extent(-65, -40, -40, -23) # xmin, xmax, ymin, ymax

# ref datasets
#sst <- "cortadv6_WeeklySST.nc"
ssta <- "cortadv6_SSTA.nc"
#tsa <- "cortadv6_TSA.nc"

# first crop the data, then crop by bathymetry
#sst.brick <- brick(sst)
ssta.brick <- brick(ssta)
#tsa.brick <- brick(tsa)
#x <- crop(sst.brick, coords)
y <- crop(ssta.brick, coords)
#z <- crop(tsa.brick, coords)

bat.raster <- raster("bathraster.tif")
bat.raster <- setExtent(bat.raster, y)
bat.raster <- projectRaster(bat.raster, y)
yy <- mask(y, bat.raster)

# plume weak
plot(yy$X2012.05.29, zlim = c(-2.3,4))
# plume strong
plot(yy$X2007.08.07, zlim = c(-6.4,1.4))

# run pca - temporal patterns!
# NOTE: https://stackoverflow.com/questions/41022927/principal-component-analysis-pca-of-time-series-data-spatial-and-temporal-pat
# here, we are using the full dataset and all available weeks!
# "temporal patterns" = loadings of principal components (time-series)

# create data.frame, with rows as space, time as columns (use y)
weeks <- seq(1, nlayers(yy), 1)
ssta_dat <- data.frame(matrix(NA, nrow = ncell(yy), ncol = nlayers (yy)+2))
ssta_dat[,1] <- as.data.frame(yy[[1]], xy = T)$x
ssta_dat[,2] <- as.data.frame(yy[[1]], xy = T)$y
colnames(ssta_dat) <- c("x","y",names(yy))

for (i in 1:length(weeks)){
  ssta_dat[,i+2] <- as.vector(yy[[i]]) # gives an error, but seems to work?
}

# remove NA rows
ssta_dat2 <- ssta_dat[3:ncol(ssta_dat)]
ssta_dat2 <- ssta_dat2[rowSums(is.na(ssta_dat2)) != ncol(ssta_dat2),] # took out like 200,000 NA rows, makes sense!

# run pca, exclude lat and lon - takes like 5 mins to run
PCA_ssta <- prcomp(t(ssta_dat2), center = F, scale. = F)

summaryPCA_ssta <- summary(PCA_ssta)
summaryPCA_ssta$importance[,c(1,2,3)] #30.4, 14.0, 6.3 of variation!

# extract loadings - should be temporal patterns
# loading.ssta <- data.frame(time = as.Date(str_replace_all(names(ssta_dat2), "X", ""), "%Y.%m.%d"),
#                            loadingPC1 = scale(as.vector(PCA_ssta$rotation[,'PC1'])),
#                            loadingPC2 = scale(as.vector(PCA_ssta$rotation[,'PC2'])))

loading.ssta <- data.frame(time = as.Date(str_replace_all(names(ssta_dat2), "X", ""), "%Y.%m.%d"),
                           loadingPC1 = scale(as.vector(PCA_ssta$x[,'PC1'])),
                           loadingPC2 = scale(as.vector(PCA_ssta$x[,'PC2'])))

# plot temporal patterns over time
ggplot(loading.ssta, aes(x = time, y = loadingPC1)) + 
  #geom_point() + 
  geom_line() + 
  geom_smooth(method = "lm")+
  ggtitle("PCA of SSTA - Migration")

#ggsave("ssta_migration_temporal_overall_101821.png", dpi = 600)

# density plot
loading.ssta$split <- c(rep("1982-2000",991),rep("2001-2020",1043))

cols <- c("#2F8936","#67CB6F")
ggplot(loading.ssta, aes(x=loadingPC1, color = split, linetype = split))+
  geom_density(alpha=0.4, size = 4, adjust = 2.5)+
  scale_color_manual(values= cols)+
  geom_vline(xintercept = 1.30719, linetype = "dashed", size = 1.5)+
  ylim(0.0, 0.4)+
  xlab("Plume Index")+
  ylab("Density")+
  theme_bw()+
  theme(
    plot.title=element_text(size=20,hjust=0.5),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=32),
    axis.text = element_text(size=38, color="black"),
    legend.text = element_text(size=20,color="black"),
    legend.title = element_blank(),
    legend.position = "none")

ggsave("sstamigration_densityplot_11.16.21.png", dpi = 600, width = 10, height = 7)
  
# match up averaged spatial patterns
# early = May-June (5-6)
# late = July-Aug (7-8)

loading.ssta$month <- month(loading.ssta$time)
loading.ssta$year <- year(loading.ssta$time)

loading.ssta.early <- loading.ssta %>%
  filter(month == "5" | month == "6") %>%
  group_by(year) %>%
  summarize(mean = mean(loadingPC1),
            sd = sd(loadingPC1))

loading.ssta.late <- loading.ssta %>%
  filter(month == "7" | month == "8") %>%
  group_by(year) %>%
  summarize(mean = mean(loadingPC1),
            sd = sd(loadingPC1))

loading.ssta.overall <- loading.ssta %>%
  filter(month == "5" | month == "6" | month == "7" | month =="8") %>%
  group_by(year) %>%
  summarize(mean = mean(loadingPC1),
            sd = sd(loadingPC1))

loading.ssta.both <- data.frame(cbind(loading.ssta.early,loading.ssta.late, loading.ssta.overall))

# plot temporal patterns in means
ggplot(loading.ssta.both, aes(x = year, y = mean)) +
  geom_line(lwd = 1.15)+
  geom_line(aes(x = year, y = mean.1), lwd = 1.15, color = "red")+
  #geom_line(aes(x = year, y = mean.2), lwd = 1.15, color = "blue")+
  geom_smooth(method = "lm")+
  ggtitle("PCA of SSTA - Migration; Black = Early, Red = Late")

# plot temporal patterns in variances
ggplot(loading.ssta.both, aes(x = year, y = sd)) +
  geom_line(lwd = 1.15)+
  geom_line(aes(x = year, y = sd.1), lwd = 1.15, color = "red")+
  #geom_line(aes(x = year, y = sd.2), lwd = 1.15, color = "blue")+
  geom_smooth(method = "lm")+
  ggtitle("PCA of SSTA VAR - Migration; Black = Early, Red = Late")

# extra plot
ggplot(loading.ssta.both, aes(x = year, y = mean.2)) + geom_line()+
geom_text(aes(label = year), hjust = 0, vjust = 0, size = 5)+
  xlab("Year")+
  ylab("PCA of SSTA")+
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

#ggsave("ssta_migration_temporal_average_101821.png", dpi = 600)

# dem vars
# repro (would affect the year of the bookyear, e.g., PCA of 1994 -> affect reproduction in the breeding season (1994))
ssta.repro <- loading.ssta.both[2:38,]
ssta.repro <- ssta.repro[-30,]
ssta.repro$repro <- repro$rs
ggplot(ssta.repro, aes(mean.2, repro)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.repro, aes(sd.1, repro)) + geom_point() + geom_smooth(method = "lm")

# survival (don't need to adjust year, b/c the data is set to nonbreeding year)
survival.data <- read.csv("survival_data_11.10.21.csv")
ssta.surv <- loading.ssta.both[2:38,]
ssta.surv <- ssta.surv[-30,]
ssta.surv$asurv <- survival.data$phi.ad.med
ssta.surv$jsurv <- survival.data$phi.juv.med
ggplot(ssta.surv, aes(mean.2, asurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(mean.2, jsurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(sd, asurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(sd, jsurv)) + geom_point() + geom_smooth(method = "lm")

ggplot(ssta.surv, aes(mean.2, asurv)) + geom_point()+ geom_smooth(method = "lm")+
geom_text(aes(label = year), hjust = 0, vjust = 0, size = 5)+
  xlab("PCA of SSTA")+
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

# spatial patterns
blah <- na.omit(ssta_dat)
#PC1.ssta <- data.frame(x = blah$x, y = blah$y, val = PCA_ssta$x[,'PC1'])
PC1.ssta <- data.frame(x = blah$x, y = blah$y, val = scale(PCA_ssta$rotation[,'PC1']))

PC1.ssta.raster <- rasterFromXYZ(PC1.ssta)

# plot
palette <- colorRampPalette(colors = c("blue","turquoise","green","yellow","orange","red"))
plot(PC1.ssta.raster, col = palette(100), zlim = c(-2,2)) # pretty close!


##################################################################################################################

#### EXTRACT ESA CCI OCEAN COLOR FOR BREEDING ####
# following vignette in: https://rmendels.github.io/UsingrerddapXtracto.html
# NOTE: can't extract version 5.0 for some reason...

library(rerddapXtracto)
library(gganimate)
library(ggplot2)
library(plotdap)
library(raster)

xpos <- c(-65.5, -61) # edited from -65, -61
ypos <- c(-45, -40)
tpos <- c("1997-09-06", "last") # for some reason, can't start with 1997-09-04???
info <- rerddap::info("pmlEsaCCI42OceanColor8Day")

chlor <- rxtracto_3D(info, parameter= "chlor_a", xcoord = xpos, ycoord = ypos, tcoord = tpos)


# helper function to reshape data for ggplot - editing it to extract time as well
mapframe <- function(longitude,latitude,value){
  dims <- dim(value)
  value <- array(value, dims[1]*dims[2])
  longitude <- longitude
  valueframe <- expand.grid(x=longitude, y = latitude)
  valueframe$value <- value
  return(valueframe)
}

# extract vals
weeks <- seq(1+2, dim(chlor$chlor_a)[3]+2, 1)
valueframe <- array(NA, dim(chlor$chlor_a)[1]*dim(chlor$chlor_a)[2])
chlorframe <- expand.grid(x=chlor$longitude, y =chlor$latitude)
chlorframe[,weeks] <- valueframe
colnames(chlorframe) <- c("x","y",as.character(chlor$time))
chlorstack <- list()

for (i in 1:length(weeks)){ # this takes awhile...
  chlorframe[,i+2] <- as.vector(chlor$chlor_a[,,i])
  e <- extent(chlorframe[,1:2]) # could be like min 288
  r <- raster(e, 121, 97)
  chlorstack[[i]] <- rasterize(chlorframe[,1:2], r, chlorframe[,i+2], fun = mean)
  #chlorstack[[i]] <- rasterize(trial5[,1:2], r, trial5[,i+2], fun = mean)
}
chlorraster <- stack(chlorstack) 
crs(chlorraster) <- "+proj=longlat +datum=WGS84 +no_defs"
plot(log(chlorraster[[333]]))

# interpolate data
chlorraster3 <- approxNA(chlorraster, rule = 2)
names(chlorraster3) <- as.character(chlor$time)

# run pca - temporal patterns!
# create data.frame, with rows as space, time as columns (use y)
weeks <- seq(1, nlayers(chlorraster3), 1)
ssta_dat <- data.frame(matrix(NA, nrow = ncell(chlorraster3), ncol = nlayers (chlorraster3)+2))
ssta_dat[,1] <- as.data.frame(chlorraster3[[1]], xy = T)$x
ssta_dat[,2] <- as.data.frame(chlorraster3[[1]], xy = T)$y
colnames(ssta_dat) <- c("x","y",names(chlorraster3))

for (i in 1:length(weeks)){
  ssta_dat[,i+2] <- as.vector(chlorraster3[[i]])
}

# remove NA rows
ssta_dat2 <- ssta_dat[3:ncol(ssta_dat)]
ssta_dat2 <- ssta_dat2[rowSums(is.na(ssta_dat2)) != ncol(ssta_dat2),]
ssta_dat2 <- log(ssta_dat2)
#ssta_dat3 <- ssta_dat2[,colSums(is.na(ssta_dat2)) != nrow(ssta_dat2)]

# take seasonality out of the values (next code needs to be run first...)
# years as rows, weeks as columns
# NOTE: really does not work! find chlorophyll anomaly data instead...
blah <- data.frame(t(ssta_dat2))
blah$date <- as.Date(chlor$time)
blah$year <- year(blah$date)
blah$month <- month(blah$date)
blah$day <- day(blah$date)

#try to create some labeling system for unique days?s
a <- table(blah$month, blah$year)
a[a == 0] <- 1 # 0s mess up things, remove 1s later
blargh <- list()
blahblah <- as.vector(a)
for (i in 1:length(blahblah)){
  blargh[[i]] <- seq(1, blahblah[i], 1)
}
b <- unlist(blargh)
b2 <- b[-c(1:8)] # remove first ones
b2 <- b2[-c(1040:1048)] # remove last ones
blah$order <- b2 # YAY
blah$yikes <- (blah$month*4)+blah$order

# try something out?
trial1 <- blah %>%
  group_by(month, order) %>%
  summarize_at(vars(-date, -year, -day, -yikes),
               funs(mean(., na.rm = T))) %>%
  inner_join(blah, by = c("month", "order"))

# cut out months (Oct-Feb)
trial1 <- trial1 %>%
  filter(month == "10" | month == "11" | month == "12" | month == "1" | month == "2")

# split datasets and subtract
trial2a <- trial1[,3:9124]
trial2b <- trial1[,9125:18246]
trial3 <- trial2b - trial2a
trial3$date <- trial1$date
trial3 <- arrange(trial3, date)
trial4 <- data.frame(t(trial3[,1:9122]))
colnames(trial4) <- trial3$date

# run pca, exclude lat and lon - takes like 5 mins to run
# NOTE: remember that I logged the values before PCA-ing!
#PCA_ssta <- prcomp(t(log(ssta_dat2)), center = F, scale. = F)
PCA_ssta <- prcomp(t(trial4), center = F, scale. = F)

summaryPCA_ssta <- summary(PCA_ssta)
summaryPCA_ssta$importance[,c(1,2)] #22.5, 7.2 of variation!

# extract loadings - should be temporal patterns
# loading.ssta <- data.frame(time = as.Date(str_replace_all(names(ssta_dat2), "X", ""), "%Y.%m.%d"),
#                            loadingPC1 = scale(as.vector(PCA_ssta$rotation[,'PC1'])),
#                            loadingPC2 = scale(as.vector(PCA_ssta$rotation[,'PC2'])))

# NOTE 11.11.21: had to edit this, not sure why...
#loading.ssta <- data.frame(time = as.Date(str_replace_all(names(trial4), "X", ""), "%Y.%m.%d"),
#                           loadingPC1 = scale(as.vector(PCA_ssta$x[,'PC1'])),
#                           loadingPC2 = scale(as.vector(PCA_ssta$x[,'PC2'])))

loading.ssta <- data.frame(time = as.Date(names(trial4)),
                           loadingPC1 = scale(as.vector(PCA_ssta$x[,'PC1'])),
                           loadingPC2 = scale(as.vector(PCA_ssta$x[,'PC2'])))

# plot temporal patterns over time
ggplot(loading.ssta, aes(x = time, y = loadingPC1)) + 
  #geom_point() + 
  geom_line() + 
  geom_smooth(method = "lm")+
  ggtitle("PCA of SSTA - Migration")

ggsave("ssta_migration_temporal_overall_101821.png", dpi = 600)

# match up averaged spatial patterns

loading.ssta$month <- month(loading.ssta$time)
loading.ssta$year <- year(loading.ssta$time)
loading.ssta$bookyear <- ifelse(loading.ssta$month == "1" | loading.ssta$month == "2",
                                loading.ssta$year-1, loading.ssta$year)

loading.ssta.early <- loading.ssta %>%
  filter(month == "10" | month == "11" | month == "12") %>%
  group_by(bookyear) %>%
  summarize(mean = mean(loadingPC2),
            sd = sd(loadingPC2))

loading.ssta.late <- loading.ssta %>%
  filter(month == "1" | month == "2") %>%
  group_by(bookyear) %>%
  summarize(mean = mean(loadingPC2),
            sd = sd(loadingPC2))

loading.ssta.overall <- loading.ssta %>%
  filter(month == "10" | month == "11" | month == "12" | month =="1" | month == "2") %>%
  group_by(bookyear) %>%
  summarize(mean = mean(loadingPC2),
            sd = sd(loadingPC2))

loading.ssta.both <- data.frame(cbind(loading.ssta.early,loading.ssta.late, loading.ssta.overall))

write.csv(loading.ssta.both, "chlor.PC2.breeding.11.15.21.csv")

# plot temporal patterns in means
ggplot(loading.ssta.both, aes(x = bookyear, y = mean)) +
  geom_line(lwd = 1.15)+
  geom_line(aes(x = bookyear, y = mean.1), lwd = 1.15, color = "red")+
  #geom_line(aes(x = bookyear, y = mean.2), lwd = 1.15, color = "blue")+
  geom_smooth(method = "lm")+
  ggtitle("PCA of Chlor-a Breeding; Black = Early, Red = Late")

# plot temporal patterns in variances
ggplot(loading.ssta.both, aes(x = bookyear, y = sd)) +
  geom_line(lwd = 1.15)+
  geom_line(aes(x = bookyear, y = sd.1), lwd = 1.15, color = "red")+
  #geom_line(aes(x = bookyear, y = sd.2), lwd = 1.15, color = "blue")+
  geom_smooth(method = "lm")+
  ggtitle("PCA of Chlor-a VAR - Breeding; Black = Early, Red = Late")

# extra plot
ggplot(loading.ssta.both, aes(x = bookyear, y = mean.2)) + geom_line()+
  geom_text(aes(label = bookyear), hjust = 0, vjust = 0, size = 5)+
  xlab("Year")+
  ylab("PCA of Chlor-a")+
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

#ggsave("ssta_migration_temporal_average_101821.png", dpi = 600)

# dem vars
# repro (would affect the year of the bookyear, e.g., PCA of 1994 -> affect reproduction in the breeding season (1994))
ssta.repro <- loading.ssta.both
ssta.repro <- ssta.repro[-15,]
ssta.repro$repro <- repro$rs[15:36]
ggplot(ssta.repro, aes(mean, repro)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.repro, aes(sd, repro)) + geom_point() + geom_smooth(method = "lm")

# survival (don't need to adjust year, b/c the data is set to nonbreeding year)
survival.data <- read.csv("survival_data_11.10.21.csv")
ssta.surv <- loading.ssta.both[1:22,]
ssta.surv <- ssta.surv[-15,]
ssta.surv$asurv <- survival.data$phi.ad.med[16:36]
ssta.surv$jsurv <- survival.data$phi.juv.med[16:36]
ggplot(ssta.surv, aes(mean.2, asurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(mean.2, jsurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(sd.2, asurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(sd.2, jsurv)) + geom_point() + geom_smooth(method = "lm")

ggplot(ssta.surv, aes(mean.2, asurv)) + geom_point()+ geom_smooth(method = "lm")+
  geom_text(aes(label = bookyear), hjust = 0, vjust = 0, size = 5)+
  xlab("PCA of SSTA")+
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


# spatial patterns
blah <- na.omit(ssta_dat)
PC1.ssta <- data.frame(x = blah$x, y = blah$y, val = scale(PCA_ssta$rotation[,'PC1']))

PC1.ssta.raster <- rasterFromXYZ(PC1.ssta)

# plot
palette <- colorRampPalette(colors = c("blue","turquoise","green","yellow","orange","red"))
plot(PC1.ssta.raster, col = palette(100), zlim = c(-3,3)) 


##################################################################################################################

#### EXTRACT ESA CCI OCEAN COLOR FOR MIGRATION ####

library(rerddapXtracto)
library(gganimate)
library(ggplot2)
library(plotdap)

xpos <- c(-65, -40) 
ypos <- c(-40, -23)
tpos <- c("1997-09-06", "last") # for some reason, can't start with 1997-09-04???
info <- rerddap::info("pmlEsaCCI42OceanColor8Day")

chlor <- rxtracto_3D(info, parameter= "chlor_a", xcoord = xpos, ycoord = ypos, tcoord = tpos)

# helper function to extract data
mapframe <- function(longitude,latitude,value){
  dims <- dim(value)
  value <- array(value, dims[1]*dims[2])
  longitude <- longitude
  valueframe <- expand.grid(x=longitude, y = latitude)
  valueframe$value <- value
  return(valueframe)
}

# extract vals
weeks <- seq(1+2, dim(chlor$chlor_a)[3]+2, 1)
valueframe <- array(NA, dim(chlor$chlor_a)[1]*dim(chlor$chlor_a)[2])
chlorframe <- expand.grid(x=chlor$longitude, y =chlor$latitude)
chlorframe[,weeks] <- valueframe
colnames(chlorframe) <- c("x","y",as.character(chlor$time))

chlorstack <- list()

for (i in 1:length(weeks)){ # this takes awhile...
  chlorframe[,i+2] <- as.vector(chlor$chlor_a[,,i])
  e <- extent(chlorframe[,1:2]) # could be like min 288
  r <- raster(e, 400, 409) # 601, 409
  chlorstack[[i]] <- rasterize(chlorframe[,1:2], r, chlorframe[,i+2], fun = mean)
}
chlorraster <- stack(chlorstack) 
crs(chlorraster) <- "+proj=longlat +datum=WGS84 +no_defs"
#writeRaster(chlorraster, "chlorraster.tif", format = "GTiff") # DOESNT WORK!!!
plot(log(chlorraster[[333]]))

#chlorraster <- raster("chlorraster.tif")

# clip by bathymetry
bat.raster <- raster("bathraster.tif")
bat.raster <- setExtent(bat.raster, chlorraster)
bat.raster <- projectRaster(bat.raster, chlorraster)
chlorraster2 <- mask(chlorraster, bat.raster)
names(chlorraster2) <- as.character(chlor$time)
plot(log(chlorraster2[[333]]))

# plume weak? = 2012.05.29
plot(log(chlorraster2$X2012.05.24),zlim=c(-2,4))
# plume strong? = 2007.08.07
plot(log(chlorraster2$X2007.08.05),zlim=c(-2,4))

# looking at strong vs. weak?
blah <- rep(NA, nlayers(chlorraster2))
for (i in 1:nlayers(chlorraster2)){
  blah[i] <- cellStats(log(chlorraster2[[i]]), "sum")
}

# min vals
plot(log(chlorraster2[[114]]))
# max vals???? keep messing with this
plot(log(chlorraster2[[418]]))

# think I need to interpolate the missing data...
# https://gis.stackexchange.com/questions/278476/filling-na-gaps-in-raster-with-r
chlorraster3 <- approxNA(chlorraster2, rule = 2)

# run pca - temporal patterns!
# need to do some interpolation/gap filling on some of the data!!!!
# NOTE: https://stackoverflow.com/questions/41022927/principal-component-analysis-pca-of-time-series-data-spatial-and-temporal-pat
# here, we are using the full dataset and all available weeks!
# "temporal patterns" = loadings of principal components (time-series)

# create data.frame, with rows as space, time as columns (use y)
weeks <- seq(1, nlayers(chlorraster3), 1)
ssta_dat <- data.frame(matrix(NA, nrow = ncell(chlorraster3), ncol = nlayers (chlorraster2)+2))
ssta_dat[,1] <- as.data.frame(chlorraster3[[1]], xy = T)$x
ssta_dat[,2] <- as.data.frame(chlorraster3[[1]], xy = T)$y
colnames(ssta_dat) <- c("x","y",names(chlorraster2))

for (i in 1:length(weeks)){
  ssta_dat[,i+2] <- as.vector(chlorraster3[[i]])
}

# remove NA rows
ssta_dat2 <- ssta_dat[3:ncol(ssta_dat)]
ssta_dat2 <- ssta_dat2[rowSums(is.na(ssta_dat2)) != ncol(ssta_dat2),] 
ssta_dat2 <- log(ssta_dat2)
#ssta_dat3 <- ssta_dat2[,colSums(is.na(ssta_dat2)) != nrow(ssta_dat2)]

# take seasonality out of the values (next code needs to be run first...)
# years as rows, weeks as columns
# NOTE: really does not work! find chlorophyll anomaly data instead...
blah <- data.frame(t(ssta_dat2))
blah$date <- as.Date(chlor$time)
blah$year <- year(blah$date)
blah$month <- month(blah$date)
blah$day <- day(blah$date)

#try to create some labeling system for unique days?s
a <- table(blah$month, blah$year)
a[a == 0] <- 1 # 0s mess up things, remove 1s later
blargh <- list()
blahblah <- as.vector(a)
for (i in 1:length(blahblah)){
  blargh[[i]] <- seq(1, blahblah[i], 1)
}
b <- unlist(blargh)
b2 <- b[-c(1:8)] # remove first ones
b2 <- b2[-c(1040:1048)] # remove last ones
blah$order <- b2 # YAY
blah$yikes <- (blah$month*4)+blah$order

# try something out?
trial1 <- blah %>%
  group_by(month, order) %>%
  summarize_at(vars(-date, -year, -day, -yikes),
               funs(mean(., na.rm = T))) %>%
  inner_join(blah, by = c("month", "order"))

# split datasets and subtract
trial2a <- trial1[,3:21730]
trial2b <- trial1[,21731:43458]
trial3 <- trial2b - trial2a
trial3$date <- trial1$date
trial3 <- arrange(trial3, date)
trial4 <- data.frame(t(trial3[,1:21728]))
colnames(trial4) <- trial3$date

# run pca, exclude lat and lon - takes like 5 mins to run
# NOTE: remember that I logged the values before PCA-ing
#PCA_ssta <- prcomp(t(log(ssta_dat2)), center = F, scale. = F)
PCA_ssta <- prcomp(t(trial4), center = F, scale. = F) # logged earlier

summaryPCA_ssta <- summary(PCA_ssta)
summaryPCA_ssta$importance[,c(1,2)] #14.2, 8.9

# extract loadings - should be temporal patterns
# loading.ssta <- data.frame(time = as.Date(str_replace_all(names(ssta_dat2), "X", ""), "%Y.%m.%d"),
#                            loadingPC1 = scale(as.vector(PCA_ssta$rotation[,'PC1'])),
#                            loadingPC2 = scale(as.vector(PCA_ssta$rotation[,'PC2'])))

loading.ssta <- data.frame(time = as.Date(str_replace_all(names(ssta_dat2), "X", ""), "%Y.%m.%d"),
                           loadingPC1 = scale(as.vector(PCA_ssta$x[,'PC1'])),
                           loadingPC2 = scale(as.vector(PCA_ssta$x[,'PC2'])))

# plot temporal patterns over time
ggplot(loading.ssta, aes(x = time, y = loadingPC1)) + 
  #geom_point() + 
  geom_line() + 
  geom_smooth(method = "lm")+
  ggtitle("PCA of Chlorophyll-a - Migration")

ggsave("ssta_migration_temporal_overall_101821.png", dpi = 600)

# match up averaged spatial patterns
# early = May-June (5-6)
# late = July-Aug (7-8)

loading.ssta$month <- month(loading.ssta$time)
loading.ssta$year <- year(loading.ssta$time)

loading.ssta.early <- loading.ssta %>%
  filter(month == "5" | month == "6") %>%
  group_by(year) %>%
  summarize(mean = mean(loadingPC1),
            sd = sd(loadingPC1))

loading.ssta.late <- loading.ssta %>%
  filter(month == "7" | month == "8") %>%
  group_by(year) %>%
  summarize(mean = mean(loadingPC1),
            sd = sd(loadingPC1))

loading.ssta.overall <- loading.ssta %>%
  filter(month == "5" | month == "6" | month == "7" | month =="8") %>%
  group_by(year) %>%
  summarize(mean = mean(loadingPC1),
            sd = sd(loadingPC1))

loading.ssta.both <- data.frame(cbind(loading.ssta.early,loading.ssta.late,loading.ssta.overall))

#write.csv(loading.ssta.both, "chlor.PC2.migration.11.16.21.csv")

# plot temporal patterns
ggplot(loading.ssta.both, aes(x = year, y = mean)) +
  geom_line(lwd = 1.15)+
  geom_line(aes(x = year, y = mean.1), lwd = 1.15, color = "red")+
  #geom_line(aes(x = year, y = mean.2), lwd = 1.15, color = "blue")+
  #geom_smooth(method = "lm")+
  ggtitle("PCA of Chlorophyll-a - Migration; Black = Early, Red = Late")

#ggsave("ssta_migration_temporal_average_101821.png", dpi = 600)

# dem vars
# repro (would affect the year of the bookyear, e.g., PCA of 1994 -> affect reproduction in the breeding season (1994))
ssta.repro <- loading.ssta.both
ssta.repro <- ssta.repro[-14,]
ssta.repro$repro <- repro$rs[16:36]
ggplot(ssta.repro, aes(mean.1, repro)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.repro, aes(sd.2, repro)) + geom_point() + geom_smooth(method = "lm")

# survival (don't need to adjust year, b/c the data is set to nonbreeding year)
survival.data <- read.csv("survival_data_11.10.21.csv")
ssta.surv <- loading.ssta.both
ssta.surv <- ssta.surv[-15,]
ssta.surv$asurv <- survival.data$phi.ad.med[16:36]
ssta.surv$jsurv <- survival.data$phi.juv.med[16:36]
ggplot(ssta.surv, aes(mean.2, asurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(mean.2, jsurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(sd, asurv)) + geom_point() + geom_smooth(method = "lm")
ggplot(ssta.surv, aes(sd.2, jsurv)) + geom_point() + geom_smooth(method = "lm")

# spatial patterns
blah <- na.omit(ssta_dat)
PC1.ssta <- data.frame(x = blah$x, y = blah$y, val = scale(PCA_ssta$rotation[,'PC1']))

PC1.ssta.raster <- rasterFromXYZ(PC1.ssta)

# plot
palette <- colorRampPalette(colors = c("blue","turquoise","green","yellow","orange","red"))
plot(PC1.ssta.raster, col = palette(100), zlim = c(-3,2)) 





