# Load necessary libraries
library(readr)        # For reading CSV files
library(readxl)       # For reading Excel files
library(lubridate)    # For handling date-time data
library(dplyr)        # For data manipulation
library(ggplot2)      # For plotting
library(RColorBrewer) # For plotting

# Set working directory
setwd("C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/eva pan")

# Open pdf for plot
pdf(file = "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/eva pan/plot/topo_factor_barplot.pdf", height = 6, width = 9, paper = "special")

# Clear workspace
rm(list=ls())
# options(error=NULL)

# Set color palette
coul <- brewer.pal(5, "Set2") 

# Function to import water level and corresponding date and time from CSV
water_level_data_import  <- function(filename) {
  
  # Read the whole water level file, skip the first row in excel
  A <- read.csv(file = filename, header = FALSE, skip = 2, sep = ',')
  
  h <- A[[3]] # Assuming the water level data is in the 3rd column
  date_time <- as.POSIXct(A[[2]], format="%m/%d/%Y %H:%M") # Assuming the time data is in the 2nd column
  return(data.frame(date_time=date_time, h=h))
}

# Import data from two CSV files
data1 <- water_level_data_import('ML_EP_20230907.csv')
data2 <- water_level_data_import('2023-11-2MLevap (2).csv')

# Merge the two water level data and time data
h <- c(data1$h, data2$h)
date_time <- c(data1$date_time, data2$date_time)

# Function to calculate evaporation rate
eva_rate_calc <- function(h, date_time) {
  eva_rate <- c()
  eva_date_time <- c()
  
  for (i in 1:(length(h) - 2)) {
    window <- h[i:(i + 2)]
    duration <- as.numeric(difftime(date_time[i + 2], date_time[i], units="hours"))
    if (all(diff(window) <= 0)) { # Check if the water level is in descending order
      rate <- (window[1] - window[3]) / duration
      eva_rate <- c(eva_rate, rate)
      eva_date_time <- c(eva_date_time, as.POSIXct(mean(c(date_time[i], date_time[i + 2]))))
    }
  }
  return(data.frame(eva_date_time=eva_date_time, eva_rate=eva_rate))
}

# Calculate evaporation rate and average evaporation rate
eva_data <- eva_rate_calc(h, date_time)
eva_rate <- eva_data$eva_rate
eva_date_time <- eva_data$eva_date_time
eva_rate_avg <- mean(eva_rate[eva_rate < 0.75]) # mm/hr
eva_rate_avg

# Function to import weather data and plot it
weather_factor <- function(filename) {
  sheets <- excel_sheets(filename)
  num_plot <- length(sheets)
  
  for (sheet in sheets) {
    col_names <- array(read_excel(filename, sheet = sheet, n_max = 1, col_names = FALSE))
    data <- data.frame(read_excel(filename, sheet = sheet, skip = 2, col_names = FALSE))
    colnames(data) <- col_names
    t <- data[[1]] # Assuming the time data is in the 1nd column
    values <- data[, -1]
    
    if (sheet != "Precipitation,in") {
      matplot(t, values, type="l", lty=1, lwd=1, col=1:ncol(values),ylab = sheet)
    } else {
      values[is.na(values)] <- 0
      barplot(-1 * as.matrix(t(values)), beside=TRUE, names.arg=t,ylab = sheet, border = coul, col = coul)
    }
    
    legend("topright", legend=col_names, col=1:ncol(values), lty=1, cex=0.8)
  }
}

# Plot weather factors
weather_factor('Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx')

dev.off()

# Import precipitation data
P_d <- read_excel('Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx', sheet='Precipitation,in', skip = 2)[[2]]
t <- ymd(read_excel('Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx', sheet='Precipitation,in', skip = 2)[[1]])
P_d[is.na(P_d)] <- 0

# Create data frames for plotting
df <- data.frame(eva_date_time=as.POSIXct(eva_date_time, format = "%Y-%m-%d", tz = "UTC"), eva_rate=eva_rate)
df_precip <- data.frame(date_time=as.POSIXct(t, format = "%Y-%m-%d", tz = "UTC"), P_d=-1 * P_d * 25.4 / 24)

# Set plot directory
plot_dir <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/eva pan/plot"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Plot evaporation rate and precipitation
evarate_prep_barplot <- ggplot() +
  geom_bar(data=df, aes(x=eva_date_time, y=eva_rate), stat="identity", fill="blue", color="blue") +
  geom_bar(data=df_precip, aes(x=date_time, y=P_d), stat="identity", fill="red", color="red") +
  geom_hline(yintercept=eva_rate_avg, linetype="dashed", color="blue") +
  labs(x="Date Time", y="mm/hr", title="Evaporation Rate and Precipitation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%b %Y") # Adjust date_breaks as needed for your data

# Save ggplot
ggsave(filename = "evarate_prep_barplot.pdf", plot = evarate_prep_barplot, path = plot_dir, device = "pdf", width = 9, height = 6)
