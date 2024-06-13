## the purpose of this script is to take hourly temperature data from CQ251
## and produce a daily time series of max, min, and ave temps

# Load required libraries
library(dplyr)
library(lubridate)

getwd()
setwd("G:/.shortcut-targets-by-id/1z6VpmA5-jrDJz7t8GhrIMasEXUCDJwxo/520 Mono Lake/data")
getwd()

temps = read.csv(file = "CQ251_Sep_Nov_23_Meso.csv", header = T,skip=7)

df = temps[,2:3]
colnames(df) = c("Time","Temp")

# Convert UTC time to PST
df$Time <- as.POSIXct(df$Time, format = "%m/%d/%Y %H:%M", tz = "UTC")
df$Time <- with_tz(df$Time, tzone = "America/Los_Angeles")

# Aggregate data into daily max, min, and average temperatures
df <- df %>%
  mutate(date = as.Date(Time)) %>%
  group_by(date) %>%
  summarise(
    max_temp = max(Temp, na.rm = TRUE),
    min_temp = min(Temp, na.rm = TRUE),
    avg_temp = mean(Temp, na.rm = TRUE)
  ) %>%
  ungroup()


## Write output to .csv

setwd("G:/.shortcut-targets-by-id/1z6VpmA5-jrDJz7t8GhrIMasEXUCDJwxo/520 Mono Lake/data/")
getwd()
data =  df
write.csv(data, file = "temps.csv",row.names = F)
