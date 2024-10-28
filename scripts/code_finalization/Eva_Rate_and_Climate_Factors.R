library(ggplot2)
library(dplyr)
library(readr)
library(lubridate)

# Function to read and combine CSV files from a folder
read_and_combine_csv <- function(folder_path, date_col, value_col) {
  files <- list.files(folder_path, full.names = TRUE, pattern = "\\.csv$")
  data_list <- lapply(files, read_csv)
  combined_data <- bind_rows(data_list)
  # Assuming the first column is date and the second is the value
  combined_data <- combined_data %>% select(date_col, value_col)
  combined_data$date_col <- ymd_hms(combined_data[[date_col]])  # Convert to datetime
  combined_data[[value_col]] <- as.numeric(combined_data[[value_col]])  # Ensure it's numeric
  return(list(combined_data[[value_col]], combined_data$date_col))
}

# Function to calculate evaporation rate (stub, replace with actual implementation)
eva_rate_cal <- function(water_level, date_time) {
  # Placeholder: replace this with your actual evaporation rate calculation logic
  return(water_level * 0.1)  # Example calculation
}

# Function to calculate daily average and export to CSV
calculate_and_export_daily_avg <- function(data, date_time, output_folder, filename, variable_name) {
  daily_avg <- data.frame(Date = as.Date(date_time), Value = data) %>%
    group_by(Date) %>%
    summarize(Daily_Average = mean(Value, na.rm = TRUE))
  
  write_csv(daily_avg, file.path(output_folder, filename))
  return(daily_avg)
}

# Initialization
# Water Level
folder_path <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/Pan Water Level"
all_h, all_date_time <- read_and_combine_csv(folder_path, 5, 7)

# Plot Water Level
ggplot(data = data.frame(DateTime = all_date_time, WaterLevel = all_h), aes(x = DateTime, y = WaterLevel)) +
  geom_line(color = "blue") +
  labs(x = "Date and Time", y = "Water Level (mm)", title = "Water Level over Time") +
  theme_minimal()

# Evaporation Rate
eva_rate <- eva_rate_cal(all_h, all_date_time)
eva_rate_avg <- mean(eva_rate[eva_rate < 3], na.rm = TRUE)  # Exclude outliers

# Bar plot for Evaporation Rate
ggplot(data = data.frame(DateTime = all_date_time, EvaporationRate = eva_rate), aes(x = DateTime, y = EvaporationRate)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_hline(yintercept = eva_rate_avg, linetype = "dashed", color = "blue") +
  annotate("text", x = tail(all_date_time, 1) - days(1), y = max(eva_rate, na.rm = TRUE) + 0.2,
           label = paste("Average:", round(eva_rate_avg, 2), "mm/hr"),
           hjust = 1, vjust = 1, fontface = "bold", color = "blue") +
  labs(x = "Date", y = "Evaporation Rate (mm/hr)", title = "Evaporation Rate Over Time") +
  theme_minimal()

# Calculate and export daily evaporation rates
calculate_and_export_daily_avg(eva_rate, eva_date_time, "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/2024", "2024_daily_average_evaporation_rates.csv", "Evaporation Rate mm/hr")

# Water Temperature
folder_path <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/Water Temp"
T_water_all, T_water_all_date_time <- read_and_combine_csv(folder_path, 1, 1)

# Plot Water Temperature
ggplot(data = data.frame(DateTime = T_water_all_date_time, Temperature = T_water_all), aes(x = DateTime, y = Temperature)) +
  geom_line(color = "red") +
  labs(x = "Date", y = "Temperature (C)", title = "Water Temperature over Time") +
  theme_minimal()

# Calculate and export daily water temperature
calculate_and_export_daily_avg(T_water_all, T_water_all_date_time, "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/2024", "2024_daily_average_water_temperature.csv", "Water Temp C")

# Air Temperature
folder_path <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/Air Temp"
T_air_all, T_air_all_date_time <- read_and_combine_csv(folder_path, 1, 1)

# Plot Air Temperature
ggplot(data = data.frame(DateTime = T_air_all_date_time, Temperature = T_air_all), aes(x = DateTime, y = Temperature)) +
  geom_line(color = "green") +
  labs(x = "Date", y = "Temp (C)", title = "Air Temperature over Time") +
  theme_minimal()

# Calculate and export daily air temperature
calculate_and_export_daily_avg(T_air_all, T_air_all_date_time, "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/2024", "2024_daily_average_air_temperature.csv", "Air Temp C")

# Relative Humidity
folder_path <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/RH"
RH_all, RH_all_date_time <- read_and_combine_csv(folder_path, 1, 1)

# Plot Relative Humidity
ggplot(data = data.frame(DateTime = RH_all_date_time, RH = RH_all), aes(x = DateTime, y = RH)) +
  geom_line(color = "black") +
  labs(x = "Date", y = "RH (%)", title = "Relative Humidity over Time") +
  theme_minimal()

# Calculate and export daily RH
calculate_and_export_daily_avg(RH_all, RH_all_date_time, "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/2024", "2024_daily_average_RH.csv", "RH %")
