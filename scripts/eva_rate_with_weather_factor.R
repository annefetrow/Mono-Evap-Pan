# Load necessary libraries
library(readr)        # For reading CSV files
library(readxl)       # For reading Excel files
library(lubridate)    # For handling date-time data
library(dplyr)        # For data manipulation
library(ggplot2)      # For plotting
library(RColorBrewer) # For plotting

# Define paths for input data and output plots
data_directory <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/eva pan"
plot_directory <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/plots"

# Ensure the plot directory exists
if (!dir.exists(plot_directory)) {
  dir.create(plot_directory)
}

# Source function files
source("C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/scripts/eva_climateFactor_functions.R")
source("C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/scripts/eva_calc_functions.R")

# Set color palette
coul <- brewer.pal(5, "Set2")

# Import and preprocess water level data from multiple sources
data1 <- water_level_data_import(file.path(data_directory, 'ML_EP_20230907.csv'))
data2 <- water_level_data_import(file.path(data_directory, '2023-11-2MLevap (2).csv'))
h <- c(data1$h, data2$h)
date_time <- as.POSIXct(c(data1$date_time, data2$date_time))

# Calculate evaporation rate
eva_data <- eva_rate_calc(h, date_time)
eva_rate <- eva_data$eva_rate
eva_date_time <- as.POSIXct(eva_data$eva_date_time)
eva_data$date_only <- as.Date(eva_data$eva_date_time) # Convert the eva_date_time column to Date objects to remove the time component


# Compute average evaporation rate
eva_rate_avg <- mean(eva_rate[eva_rate < 0.75], na.rm = TRUE) # Ignore evaporation that is larger than 0.75 mm/hr
gg_daily <- daily_eva_rate_calc(eva_data)
gg_monthly <- monthly_eva_rate_clac(eva_data)

# Create data frames for plotting evaporation and precipitation
df_eva <- data.frame(eva_date_time, eva_rate)
df_precip <- precip(data_directory, 'Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx')

# Generate and save plot for evaporation rate and precipitation
eva_prep_plot <- ggplot() +
  geom_bar(data = df_eva, aes(x = eva_date_time, y = eva_rate), stat = "identity", fill = "blue", color = "blue") +
  geom_bar(data = df_precip, aes(x = date_time, y = P_d), stat = "identity", fill = "red", color = "red") +
  geom_hline(yintercept = eva_rate_avg, linetype = "dashed", color = "blue") +
  labs(x = "Date Time", y = "mm/hr", title = "Evaporation Rate and Precipitation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%b %Y")  # Adjust date_breaks as needed

# Save plot to file
ggsave(filename = "evarate_prep_barplot.pdf", plot = eva_prep_plot, path = plot_directory, device = "pdf", width = 9, height = 6)
ggsave(filename = "daily_evap_rate.pdf", plot = gg_daily, path = plot_directory, device = "pdf", width = 9, height = 6)
ggsave(filename = "monthly_evap_rate.pdf", plot = gg_monthly, path = plot_directory, device = "pdf", width = 9, height = 6)

# Plot weather factors and save to PDF
weather_factors_filename <- 'Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx'
plot_weather_factors(file.path(data_directory, weather_factors_filename), plot_directory)

