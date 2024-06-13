###############################################################################
###############################################################################
###
###     HAMON METHOD FOR ESTIMATION OF EVAPORATION AT MONO LAKE 
###
###############################################################################
###############################################################################


## Load relevant packages

library(lubridate)
library(dplyr)
library(tidyr)

## Check the working directory (set the directory)

getwd()
setwd("G:/.shortcut-targets-by-id/1z6VpmA5-jrDJz7t8GhrIMasEXUCDJwxo/520 Mono Lake/data/air temps/")
getwd()


## First read in the temperature data

temps = read.csv(file = "temps.csv", header = T)
head(temps)


## Isolate the average temperatures
temps = temps[,c('date','avg_temp')]
head(temps)


## Rename the column headers

colnames(temps) = c("date","avetemp")
head(temps)


## Convert temps to Celcius

temps$avetemp = (temps$avetemp-32)*5/9

## Define dates as dates

class (temps$date)
temps$date = lubridate::ymd(temps$date)
class (temps$date)


## Create columns for month and year
# temps$month = months(temps$date)
# temps$year = format(temps$date,format = "%Y")
# head(temps)


## Aggregate (average) temperatures for each month in the time series

#montemps = aggregate (avetemp ~ month + year, 
#                      data = temps, mean)
#head(montemps)


## Order the time series in chronological order

# montemps = montemps %>%
#   mutate(
#     month = factor(month, levels = month.name)
#   ) %>%
#   arrange(year,month)
# 
# head(montemps)


## Store values for evaporation estimation

evap = temps
head(evap)
k = 0.0065
evap$eo = 610.8*exp((17.27*evap$avetemp)/(evap$avetemp+237.3))
evap$hsat = 18.016*evap$eo/8.314/(273.16+evap$avetemp)

# ## Manually assign the median day of each month in the time series
# 
# jstart = c(136,167,197,228,259,289,320,350)
# j=c(15,45,75,106,136,167,197,228,259,289,320,350)
# j91 = c(15,45,75,106,167,197,228,259,289,320,350)
# j93=c(15,45,75,106,136,167,197,228,289,320,350)
# j97=c(15,45,75,106,136,167,197,228,259,289,350)
# j98=c(15,45,75,106,136,197,228,259,289,320,350)
# j99=c(15,45,106,136,167,197,228,259,289,320,350)
# jend = c(15,45,75,106,136,167,197,228,259)
# evap$j = c(jstart,j,j,j91,j,j93,j,j,j,j97,j98,j99,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,jend)

# ## Manually assign the number of days in each month
# 
# dstart = c(31,30,31,31,30,31,30,31)
# d=c(31,28.25,31,30,31,30,31,31,30,31,30,31)
# d91 = c(31,28.25,31,30,30,31,31,30,31,30,31)
# d93=c(31,28.25,31,30,31,30,31,31,31,30,31)
# d97=c(31,28.25,31,30,31,30,31,31,30,31,31)
# d98=c(31,28.25,31,30,31,31,31,30,31,30,31)
# d99=c(31,28.25,30,31,30,31,31,30,31,30,31)
# dend = c(31,28.25,31,30,31,30,31,31,30)
# evap$mond = c(dstart,d,d,d91,d,d93,d,d,d,d97,d98,d99,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,dend)


evap$theta = 0.2163108+2*atan(0.9671396*tan(0.0086*(yday(evap$date)-186)))
evap$phi = asin(0.39795*cos(evap$theta))
lat = 38
p=0.8333
ld=1/12*(24-24/pi*acos((sin(p*pi/180)+sin(lat*pi/180)*sin(evap$phi))/cos(lat*pi/180)/cos(evap$phi)))
evap$PET = k*ld*evap$hsat    # units of (in/day)

# evapmon = aggregate(PET ~ month + year, data = evap, sum)


# evapmon = evapmon %>%
#   mutate(
#     month = factor(month, levels = month.name)
#   ) %>%
#   arrange(year,month)

## Look at the data

# head(evapmon)

## Make a date column
# 
# evapmon$m = substring(evapmon$month,1,3)
# head(evapmon$m)
# 
# evap$d = rep(1,nrow(evap))
# 
# evap$ye = evap$year
# head(evap$ye)
# head(evap)
# 
# evap = unite(evap,date,c(m,d,ye))
# class(evap$date)
# head(evap)
# 
# evap$date = lubridate::mdy(evap$date)
# class(evap$date)


## Write output to .csv

setwd("G:/.shortcut-targets-by-id/1z6VpmA5-jrDJz7t8GhrIMasEXUCDJwxo/520 Mono Lake/data/t-based_evap/")
getwd()
data =  subset(evap, select = -c(avetemp,eo,hsat,theta,phi))
write.csv(data, file = "hamon_pet_estimate.csv",row.names = F)

