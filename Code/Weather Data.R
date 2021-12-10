# set directory
rm(list = ls())
setwd("~/Post-Doc/Data/Punta Tombo Data")

library(tidyverse)
library(lubridate)

# GOAL: manipulate data for inclusion into model; also confirm potential bugs in the data

##################################################################################################################

#### 1) Weather Data = Precipitation/Temperature/Wind
# NOTE: check for other data? ocean surface temps?

weather <- read.csv("Site Weather.csv")

# replace NAs in precipitation
weather1 <- weather %>% 
  mutate(Precipitation = replace(Precipitation, which(is.na(Precipitation)), 0.0))

# convert character dates to actual dates
weather1$WeathDate <- as.Date(weather1$WeathDate)

#### TEMP

# temperature a la Katie's measurements
# 1) % of days with max temp >25 per breeding season
num_days_25 <- weather1 %>% group_by(BookYear) %>%
  filter(!is.na(MaxTemp)) %>%
  summarize(n_days = n(),
            n_gt25 = sum(MaxTemp > 25), # 25
            p_gt25 = n_gt25/n_days)
ggplot(num_days_25, aes(BookYear, p_gt25)) + geom_point() + geom_line() + #geom_smooth(method = "lm")+
  geom_text(aes(label = BookYear), hjust = 0, vjust = 0, size = 5)+
  xlab("Year")+
  ylab("% of Days with Max Temp > 25C")+
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
# OK - this makes it seem that temps have been decreasing quite a lot!

# plot against dem vars
num_days_25$repro <- repro$rs
survivalb <- survival.data[-1,]
heatb <- num_days_25[-36,]
heatb$asurv <- survivalb$phi.ad.med
heatb$jsurv <- survivalb$phi.juv.med

ggplot(num_days_25, aes(p_gt25, repro)) + geom_point() + geom_smooth(method = "lm") #nothing
ggplot(heatb, aes(p_gt25, asurv)) + geom_point() + geom_smooth(method = "lm") + 
  geom_text(aes(label = BookYear), hjust = 0, vjust = 0, size = 5)+
  xlab("% of Days with Max Temp > 25C")+
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
ggplot(heatb, aes(p_gt25, jsurv)) + geom_point() + geom_smooth(method = "lm")

# 2) Max temp per breeding season
max_temp_yr <- weather1 %>% group_by(BookYear) %>%
  filter(!is.na(MaxTemp)) %>%
  summarize(max = max(MaxTemp))
ggplot(max_temp_yr, aes(BookYear, max)) + geom_point() + geom_smooth(method = "lm")
# MORE EQUIVOCAL....

# 3) Average max temp per breeding season per year
av_temp_yr <- weather1 %>% group_by(BookYear) %>%
  filter(!is.na(MaxTemp)) %>%
  summarize(av = mean(MaxTemp))
ggplot(av_temp_yr, aes(BookYear, av)) + geom_point() + geom_smooth(method = "lm")
# this is also decreasing! Weird

# 4) Average temp with diff between max and min?
weather1$AvTemp <- rowMeans(weather1[,c("MaxTemp", "MinTemp")], na.rm = T)
av_temp <- weather1 %>% group_by(BookYear) %>%
  filter(!is.na(MaxTemp)) %>%
  filter(!is.na(MinTemp)) %>%
  summarize(diff = mean(MaxTemp) - mean(MinTemp),
            avmax = mean(MaxTemp),
            avmin = mean(MinTemp),
            avtot = mean(avmax, avmin),
            n_days = n(),
            blah = sum(AvTemp > 20),
            blah_av = blah/n_days)
ggplot(av_temp, aes(BookYear, blah_av)) + geom_point() + geom_smooth(method = "lm") # AVERAGE TEMP IS INCREASING

#### RAIN
rain <- weather1 %>% group_by(BookYear) %>%
  filter(!is.na(Precipitation)) %>%
  filter(!is.na(MinTemp)) %>%
  summarize(n_days = n(),
            min_temp = min(MinTemp),
            min_temp_av = mean(MinTemp),
            days_rain = sum(Precipitation > 0),
            days_rain_60 = sum(Precipitation > 30)/n_days,
            days_rain_yr = days_rain/n_days, # how many days did it rain/yr
            rain_yr = sum(Precipitation),
            mean_yr = mean(Precipitation),
            total_rain_yr = rain_yr/n_days) # total rain/
ggplot(rain, aes(BookYear, days_rain_yr)) + geom_point() + geom_smooth(method = "lm") # BIG INCREASE IN RAIN
ggplot(rain, aes(BookYear, rain_yr)) + geom_point() + geom_smooth(method = "lm") # NOT MORE RAIN

# calculate years when at least 60mm of rain fell between 15 October and 15 December
rain_60 <-  weather1 %>%
  filter(WeathDate >= as.Date(paste(year(WeathDate), 10, 15, sep = "-")),
         WeathDate <= as.Date(paste(year(WeathDate), 12, 15, sep = "-"))) %>%
  group_by(BookYear) %>%
  filter(!is.na(Precipitation)) %>%
  summarize(amt_rain = sum(Precipitation))
ggplot(rain_60, aes(BookYear, amt_rain)) + geom_point() + geom_line()+
  geom_text(aes(label = BookYear), hjust = 0, vjust = 0, size = 5)+
  xlab("Year")+
  ylab("Total Rain (15Oct and 15Dec)")+
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

# compare to reproductive success
rain_60$repro <- repro$rs
ggplot(rain_60, aes(amt_rain, repro)) + geom_point() + geom_smooth(method = "lm") + #+ geom_vline(xintercept = 60)
  geom_text(aes(label = BookYear), hjust = 0, vjust = 0, size = 5)+
  xlab("Total Rain (15Oct and 15Dec)")+
  ylab("Chicks Fledged Per Nest")+
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

# compare to survival
survivalb <- survival.data[-1,]
rain_60b <- rain_60[-36,]
rain_60b$asurv <- survivalb$phi.ad.med
rain_60b$jsurv <- survivalb$phi.juv.med
ggplot(rain_60b, aes(amt_rain, asurv)) + geom_point() # not a big correlation
ggplot(rain_60b, aes(amt_rain, jsurv)) + geom_point() # slight correlation

#### WIND
wind <- weather1 %>% group_by(BookYear) %>%
  filter(!is.na(MaxWind)) %>%
  summarize(n_days = n(),
            max_wind = mean(MaxWind),
            max_wind_day = max_wind/n_days,
            sum_wind = sum(MaxWind > 10),
            sum_wind_day = sum_wind/n_days)
ggplot(wind, aes(BookYear, sum_wind_day)) + geom_point() + geom_smooth(method = "lm") # NO REAL PATTERN
ggplot(wind, aes(BookYear, max_wind)) + geom_point() + geom_smooth(method = "lm")

##################################################################################################################

#### Messing with textfile from Trelew for ClimPact

weather <- read.table("Trelew Weather Station.txt", fill = T)
weather <- weather[-2,]
colnames(weather) <- weather[1,]
weather <- weather[-1,]
weather[weather == -9999] <- -99.9
colnames(weather) <- c("a","b","c","d","e","f","g","h","i","j")

weather2 <- data.frame(cbind(weather$e, weather$f, weather$i, weather$j))
colnames(weather2) <- c("Date", "PR", "TX", "TN")

# edit dates...must be three numerical columns: year, month, day
weather2$Date <- as.Date(weather2$Date, format = "%Y%m%d")
weather2$Year <- year(weather2$Date)
weather2$Month <- month(weather2$Date)
weather2$Day <- day(weather2$Date)

# moving around things
weather2 <- weather2[,-1]
weather3 <- data.frame(cbind(as.numeric(weather2$Year),as.numeric(weather2$Month),as.numeric(weather2$Day),
                             as.numeric(weather2$PR), as.numeric(weather2$TX), as.numeric(weather2$TN)))
#colnames(weather3) <- c("Year","Month","Day","PR","TX","TN")
write.table(weather3, "Trelew Edited Weather.txt", sep = ",", col.names=F, row.names = F)

#### UH?????
weather <- read.table("Trelew Weather Station.txt", fill = T)
weather <- weather[-2,]
colnames(weather) <- weather[1,]
weather <- weather[-1,]
weather[weather == -9999] <- NA
colnames(weather) <- c("a","b","c","d","e","f","g","h","i","j")

weather2 <- data.frame(cbind(weather$e, weather$f, weather$i, weather$j))
colnames(weather2) <- c("Date", "PR", "TX", "TN")
# edit dates...must be three numerical columns: year, month, day
weather2$Date <- as.Date(weather2$Date, format = "%Y%m%d")
weather2$PR <- as.numeric(weather2$PR)
weather2$TX <- as.numeric(weather2$TX)
weather2$TN <- as.numeric(weather2$TN)

# convert to celcius and mm
weather2$PR_mm <- weather2$PR*25.4 
weather2$TX_C <- (weather2$TX - 32) * (5/9)
weather2$TN_C <- (weather2$TN - 32) * (5/9)

ggplot(weather2, aes(Date, TX_C)) + geom_point() + geom_smooth(method = "lm")
ggplot(weather2, aes(Date, PR_mm)) + geom_point() + geom_smooth(method = "lm")

# predict weather in the future
a <- lm(TX_C ~ Date, data = weather2)
time <- seq.Date(from = as.Date("2021-09-13"), to = as.Date("2050-12-31"), by = "day")
new <- data.frame(Date = time)
cool <- data.frame(predict.lm(a, new, interval = "confidence"))
cool$date <- time

blah <- data.frame(weather2$Date, weather2$TX_C)
colnames(blah) <- c("Date","TX_C")
blahblah <- data.frame(cool$date, cool$fit)
colnames(blahblah) <- c("Date","TX_C")
yeah <- rbind(blah,blahblah)

ggplot(yeah, aes(Date, TX_C)) + geom_point()

##################################################################################################################

#### FUCKING around with Trelew Weather Station Data = DOESNT WORK?
library(climate)

# find weather station
nearest_stations_ogimet(country = "Argentina",
                        date = Sys.Date(),
                        add_map = T,
                        point = c(-1,53),
                        no_of_stations = 100) #ID: 87828

# download daily data from OGIMET repository
weather_dat <- meteo_ogimet(interval = "daily",
                            coords = F,
                            station = 87828,
                            date = c(Sys.Date()-20000,Sys.Date()))
View(weather_dat)
ggplot(weather_dat, aes(Date, TemperatureCMax)) + geom_point() + geom_smooth(method = "lm")

#### OTHER TRY

weather_dat <- meteo_noaa_hourly(station = "878280-99999",
                                 year = 2019)

weather_dat <- meteo_imgw_daily(year = 2019,
                                station = "878280-99999")

##################################################################################################################

# NEXT THING IS TO COLLECT DATA FROM PUNTA TOMBO VIA TERRACLIMATE
# THIS SHOULD BE THE OFFICIAL DATA FOR THE REGION
# VERIFY THAT IT MATCHES WITH THE AIRPORT DATA

library(AOI)
library(climateR)
library(sf)

# POTENTIAL VARS: "prcp" (precip); "palmer" (drought index); tmax; tmin; wind (wind speed)
# SITES HAVE A ~ 4km resolution. Better than weather station, which is ~100 km away!

#### PRECIPITATION
blah <- data.frame(Longitude = -65.2238, Latitude = -44.0458)
site <- st_as_sf(x = blah, coords = c("Longitude", "Latitude"),
                 crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

sites_stack <- getTerraClim(AOI = site, param = "prcp",
                            startDate = "1958-01-01",
                            endDate = "2019-12-30")

# manually add the first of each month for it to be a date, but this is fake.
dates <- as.Date(paste(sites_stack$date, "-01", sep =""))
weather <- data.frame(date = dates,
                      year = year(dates),
                      month = month(dates),
                      ppt = sites_stack$prcp)
ggplot(weather, aes(date, ppt)) + geom_point() + geom_smooth(method = "lm")

# split graph up
weather$split <- c(rep("1958-1988",372),rep("1989-2019",372))

old <- weather %>% filter(split == "1958-1988")
new <- weather %>% filter(split == "1989-2019")

# fit gamma distributions
fit_old <- fitdist(old$ppt, "gamma", method = "mme") # 1.395, 0.101
fit_new <- fitdist(new$ppt, "gamma", method = "mme") # 1.183, 0.0764

cols <- c("#26567E","#3F88C5")
ggplot(weather, aes(x=ppt, color = split, linetype = split))+
  geom_density(alpha=0.4, size = 4, adjust = 2.5)+
  scale_color_manual(values = cols)+
  geom_vline(xintercept = 32.92, linetype = "dashed", size = 1.5)+
  xlab("Average Monthly Precipitation")+
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

ggsave("precipitation_densityplot_11.16.21.png", dpi = 600, width = 10, height = 7)

# calculate yearly mean and variance 

# filter out only Oct - Dec
storm <- weather %>%
  filter(month == "10" | month == "11" | month == "12") %>%
  group_by(year) %>%
  summarize(total = sum(ppt))
ggplot(storm, aes(year, total)) + geom_point() + geom_abline(intercept = 60, slope = 0)

# split graph up
storm$split <- c(rep("old",31),rep("new",31))

ggplot(storm, aes(x=total, color = split))+
  geom_density(alpha=0.4, size = 2, adjust = 2.5)+
  scale_fill_manual(values= cols)+
  geom_vline(xintercept = 32.92, linetype = "dashed")+
  xlab("Precipitation")+
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

# dem relationships
stormb <- storm[26:62,]
stormb <- stormb[-29,]
stormb$repro <- repro$rs

ggplot(stormb, aes(total, repro)) + geom_point() + geom_smooth(method = "lm")+
  geom_text(aes(label = year))

ggplot(rain_60, aes(amt_rain, repro)) + geom_point() + geom_smooth(method = "lm")
  geom_text(aes(label = BookYear), hjust = 0, vjust = 0, size = 5)+
  xlab("Total Rain (15Oct and 15Dec)")+
  ylab("Chicks Fledged Per Nest")+
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


#### TEMPERATURE MAXIMUM - not as strong as weather data, but still an increase
# 0.6 degrees in 50 years vs. 1.3 degrees in 50 years (weather stn data)
# point is, the trend is still consistent!!!
blah <- data.frame(Longitude = -65.2238, Latitude = -44.0458)
site <- st_as_sf(x = blah, coords = c("Longitude", "Latitude"),
                 crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

sites_stack <- getTerraClim(AOI = site, param = "tmax",
                            startDate = "1958-01-01",
                            endDate = "2019-12-30")

# manually add the first of each month for it to be a date, but this is fake.
dates <- as.Date(paste(sites_stack$date, "-01", sep =""))
weather <- data.frame(date = dates,
                      year = year(dates),
                      month = month(dates),
                      tmax = sites_stack$tmax)
ggplot(weather, aes(date, tmax)) + geom_point() + geom_smooth(method = "lm")

weather$split <- c(rep("1958-1988",372),rep("1989-2019",372))

old <- weather %>% filter(split == "1958-1988")
new <- weather %>% filter(split == "1989-2019")

# fit gamma distributions
fit_old <- fitdist(old$tmax, "norm") # 18.38, 5.089
fit_new <- fitdist(new$tmax, "norm") # 18.83, 5.14

cols <- c("#B80000","#FF4747")
ggplot(weather, aes(x=tmax, color = split, linetype = split))+
  geom_density(alpha=0.4, size = 4, adjust = 1.5)+
  scale_color_manual(values= cols)+
  geom_vline(xintercept = 25.207, linetype = "dashed", size = 1.5)+
  xlab("Average Monthly Max Temp")+
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

ggsave("tempmax_densityplot_11.16.21.png", dpi = 600, width = 10, height = 7)

# calculate average max temps from Oct-Feb
heat <- weather %>%
  filter(month == "10" | month == "11" | month == "12" | month == "1"| month == "2") %>%
  group_by(year) %>%
  summarize(mean = mean(tmax))
ggplot(heat, aes(year, mean)) + geom_point() + geom_smooth(method = "lm")

heatb <- heat[26:62,]
heatb <- heatb[-29,]
heatb$repro <- repro$rs
ggplot(heatb, aes(mean, repro)) + geom_point() + geom_smooth(method = "lm")

survivalb$heat <- heatb$mean[-36]
ggplot(survivalb, aes(heat, phi.ad.med)) + geom_point()
ggplot(survivalb, aes(heat, phi.juv.med)) + geom_point()


#### TEMPERATURE MINIMUM
blah <- data.frame(Longitude = -65.2238, Latitude = -44.0458)
site <- st_as_sf(x = blah, coords = c("Longitude", "Latitude"),
                 crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

sites_stack <- getTerraClim(AOI = site, param = "tmin",
                            startDate = "1958-01-01",
                            endDate = "2019-12-30")

# manually add the first of each month for it to be a date, but this is fake.
dates <- as.Date(paste(sites_stack$date, "-01", sep =""))
weather <- data.frame(date = dates,
                      tmin = sites_stack$tmin)
ggplot(weather, aes(date, tmin)) + geom_point() + geom_smooth(method = "lm")

#### WIND
blah <- data.frame(Longitude = -65.2238, Latitude = -44.0458)
site <- st_as_sf(x = blah, coords = c("Longitude", "Latitude"),
                 crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

sites_stack <- getTerraClim(AOI = site, param = "wind",
                            startDate = "1958-01-01",
                            endDate = "2019-12-30")

# manually add the first of each month for it to be a date, but this is fake.
dates <- as.Date(paste(sites_stack$date, "-01", sep =""))
weather <- data.frame(date = dates,
                      wind = sites_stack$wind)
ggplot(weather, aes(date, wind)) + geom_point() + geom_smooth(method = "lm")









