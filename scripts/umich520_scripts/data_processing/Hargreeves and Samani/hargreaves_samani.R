###############################################################################
###############################################################################
###
###     HARGREAVES AND SAMANI METHOD FOR ESTIMATION OF EVAPORATION AT MONO LAKE 
###
###############################################################################
###############################################################################


## Load relevant packages

library(lubridate)
library(dplyr)
library(tidyr)


## Define Constants
lambda = 2.45               # latent heat of vaporization (MJ/kg)
roww = 1000                  # density of water (kg/m^3)
L = 0.663225                    # latitude (radians)
Gsc = 118.11                # Solar constant (MJ/m^2/day)



## Check the working directory (set the directory)

getwd()
setwd("G:/.shortcut-targets-by-id/1z6VpmA5-jrDJz7t8GhrIMasEXUCDJwxo/520 Mono Lake/data/air temps/")
getwd()


## First read in the temperature data

temps = read.csv(file = "temps.csv", header = T)
head(temps)

## Isolate the average temperatures

# temps = temps[,c('date','max_temp','min_temp','ave_temp')]
# head(temps)


## Rename the column headers

colnames(temps) = c("date","maxtemp","mintemp","avetemp")
head(temps)


## Convert temps to Celcius

temps$mintemp = (temps$mintemp-32) * 5/9
temps$avetemp = (temps$avetemp-32) * 5/9
temps$maxtemp = (temps$maxtemp-32) * 5/9
head(temps)


## Define dates as dates

class (temps$date)
temps$date = lubridate::ymd(temps$date)
class (temps$date)


# ## Create columns for month and year
# 
# temps$month = months(temps$date)
# temps$year = format(temps$date,format = "%Y")
# head(temps)

## Create an evap data frame

evap = temps
head(evap)

# ## Extra-terrestrial solar radiation for 38 degrees north obtained by linear
# ## interpolation of data provided by Allen and Pruitt (1986) - hydrolearn.org
# ##
# ## Calculate evaporation for each row based on which month the day falls into
# 
# jan = evap[evap$month == 'January',]
# jan$evap = 0.0023*(jan$avetemp+17.8)*sqrt(jan$maxtemp-jan$mintemp)*6.88
# 
# 
# feb = evap[evap$month == 'February',]
# feb$evap = 0.0023*(feb$avetemp+17.8)*sqrt(feb$maxtemp-feb$mintemp)*8.94
# 
# mar = evap[evap$month == 'March',]
# mar$evap = 0.0023*(mar$avetemp+17.8)*sqrt(mar$maxtemp-mar$mintemp)*11.66
# 
# apr = evap[evap$month == 'April',]
# apr$evap = 0.0023*(apr$avetemp+17.8)*sqrt(apr$maxtemp-apr$mintemp)*14.4
# 
# may = evap[evap$month == 'May',]
# may$evap = 0.0023*(may$avetemp+17.8)*sqrt(may$maxtemp-may$mintemp)*16.34
# 
# jun = evap[evap$month == 'June',]
# jun$evap = 0.0023*(jun$avetemp+17.8)*sqrt(jun$maxtemp-jun$mintemp)*17.26
# 
# jul = evap[evap$month == 'July',]
# jul$evap = 0.0023*(jul$avetemp+17.8)*sqrt(jul$maxtemp-jul$mintemp)*16.74
# 
# aug = evap[evap$month == 'August',]
# aug$evap = 0.0023*(aug$avetemp+17.8)*sqrt(aug$maxtemp-aug$mintemp)*15.26
# 
# sep = evap[evap$month == 'September',]
# sep$evap = 0.0023*(sep$avetemp+17.8)*sqrt(sep$maxtemp-sep$mintemp)*12.78
# 
# oct = evap[evap$month == 'October',]
# oct$evap = 0.0023*(oct$avetemp+17.8)*sqrt(oct$maxtemp-oct$mintemp)*10.04
# 
# nov = evap[evap$month == 'November',]
# nov$evap = 0.0023*(nov$avetemp+17.8)*sqrt(nov$maxtemp-nov$mintemp)*7.48
# 
# dec = evap[evap$month == 'December',]
# dec$evap = 0.0023*(dec$avetemp+17.8)*sqrt(dec$maxtemp-dec$mintemp)*6.22
# 


## Calculate daily evap
evap$j = yday(evap$date)
evap$dr = 1+0.033*cos(2*pi*evap$j/365)
evap$delta = 0.409*sin(2*pi*evap$j/365-1.39)
evap$ws = acos(-1*tan(L)*tan(evap$delta))
evap$ra = 1/pi*Gsc*evap$dr*(evap$ws*sin(L)*sin(evap$delta)+cos(L)*cos(evap$delta)*sin(evap$ws))
evap$ra.term = 1000*evap$ra/lambda/roww
evap$PET = 0.0023*(evap$avetemp+17.8)*sqrt(evap$maxtemp-evap$mintemp)*1000*evap$ra/lambda/roww


## Store values for evaporation estimation
# 
# evap = rbind(jan,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec)
# head(evap)


## Aggregate (sum) evaporation for each month in the time series

# monevap = aggregate (evap ~ month + year, 
#                        data = evap, sum)


## Order the time series in chronological order

# monevap = monevap %>%
#   mutate(
#     month = factor(month, levels = month.name)
#   ) %>%
#   arrange(year,month)
# 
# head(monevap)


## Make a date column

# monevap$m = substring(monevap$month,1,3)
# head(monevap$m)
# 
# monevap$d = rep(1,nrow(monevap))
# 
# monevap$ye = monevap$year
# head(monevap$ye)
# head(monevap)
# 
# monevap = unite(monevap,date,c(m,d,ye))
# class(monevap$date)
# head(monevap)
# 
# monevap$date = lubridate::mdy(monevap$date)
# class(monevap$date)
# head(monevap)
# 
# monevap = monevap[,c('date','evap')]
# head(monevap)
# 

## Convert evaporation from mm to inches

evap$PET = evap$PET*0.0393701
head(evap)


## Write output to .csv

setwd("G:/.shortcut-targets-by-id/1z6VpmA5-jrDJz7t8GhrIMasEXUCDJwxo/520 Mono Lake/data/t-based_evap")
getwd()
data =  subset(evap, select = -c(mintemp,avetemp,maxtemp,j,dr,delta,ws,ra,ra.term))
write.csv(data, file = "hargreaves_samani_pet_estimate.csv",row.names = F)

