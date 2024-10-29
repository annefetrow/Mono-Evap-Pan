# Load necessary libraries
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(vroom)

data_import <- function(file_name) {
  # data_import reads the data and timestamps from a CSV file.
  #
  # Args:
  #   file_name: Path to the CSV file
  #
  # Returns:
  #   A list containing:
  #     - data: Water level or other data (third column)
  #     - timestamp: Timestamps corresponding to the data (second column)
  
  # Read the CSV file
  data_table <- read.csv(file = file_name, skip = 1, head = TRUE, sep=",")
  
  # Keep only the 2nd and 3rd columns
  data_table <- data_table[, c(2, 3)]
  
  # Optionally rename columns for easier access (if desired)
  colnames(data_table) <- c("Timestamp", "Value")  # Renaming columns (optional)
  data_table$Timestamp <- mdy_hms(data_table$Timestamp) # Convert data type to datetime
  
  # Display the resulting data
  print(data_table)
  
  # Return the data and timestamps as a list
  return(list(value = data_table$Value, timestamp = data_table$Timestamp))
}

read_and_combine_csv <- function(folder_path, start_file, end_file) {
  # Get a list of all CSV files in the folder
  csv_files <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)
  
  # Check if the specified range is valid
  num_files <- length(csv_files)
  if (start_file < 1 || end_file > num_files || start_file > end_file) {
    stop("Invalid file range specified. Please check start_file and end_file indices.")
  }
  
  # Initialize empty data frames to hold combined data
  all_data <- data.frame()
  
  # Loop through each file and read the data
  for (k in start_file:end_file) {
    # Read the CSV file
    data <- data_import(csv_files[k])
    
    # Combine the data
    all_data <- bind_rows(all_data, data)
  }
  
  # Remove duplicates and NaN values
  unique_data <- all_data %>%
    distinct() %>%
    filter(!is.na(value) & !is.na(timestamp))
  
  return(unique_data)
}

eva_rate_cal <- function(h, date_time) {
  eva_rate <- c()
  eva_date_time <- as.POSIXct(character())
  
  for (i in 1:(length(h) - 2)) {
    window <- h[i:(i + 2)]
    duration <- as.numeric(difftime(date_time[i + 2], date_time[i], units = "hours"))
    
    if (is.unsorted(-window) == FALSE) {  # Check if window is sorted in descending order
      rate <- (window[1] - window[3]) / duration
      eva_rate <- c(eva_rate, rate)
      mean_timestamp <- mean(date_time[i:(i + 2)])
      eva_date_time <- c(eva_date_time, mean_timestamp)  # Mean timestamp for the window
    }
  }
  
  # Return eva_rate and eva_date_time as a list
  return(list(eva_rate = eva_rate, eva_date_time = eva_date_time))
}

# Set default line width globally for plots
theme_set(theme_minimal(base_size = 14))

# Water Level
# Path to the folder containing the CSV files for water level
folderPath <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/Pan Water Level'
all_data <- read_and_combine_csv(folderPath, 5, 7)

# Plot water level over time
ggplot(all_data, aes(x = all_data$timestamp, y = all_data$value)) +
  geom_line(color = 'blue') +
  labs(x = 'Date and Time', y = 'Water Level (mm)', title = 'Water Level over Time') +
  theme_minimal()

# Evaporation Rate Calculation Function
eva_rate <- eva_rate_cal(all_data$value, all_data$timestamp)
eva_rate <- list(
  eva_rate = eva_rate$eva_rate[eva_rate$eva_rate < 3],
  eva_date_time = eva_rate$eva_date_time[eva_rate$eva_rate < 3]
)
# eva_rate_avg <- mean(eva_rate<3, na.rm = TRUE) # exclude outliers

# # Plot evaporation rate over time
# ggplot() +
#   geom_bar(aes(x = eva_date_time, y = eva_rate), stat = "identity", fill = "blue") +
#   geom_hline(yintercept = eva_rate_avg, linetype = "dashed", color = "blue") +
#   annotate("text", x = max(eva_date_time) - days(1), y = max(eva_rate) + 0.2, 
#            label = paste("Average:", round(eva_rate_avg, 2), "mm/hr"), hjust = 1) +
#   labs(x = 'Date', y = 'Evaporation Rate (mm/hr)', title = 'Evaporation Rate Over Time') +
#   theme_minimal()
# 
# # Function to calculate daily averages and export to CSV
# calculate_and_export_daily_avg <- function(values, date_times, output_path, file_name, col_name) {
#   daily_avg <- tibble(Date = as.Date(date_times), Value = values) %>%
#     group_by(Date) %>%
#     summarize(Daily_Avg = mean(Value, na.rm = TRUE)) %>%
#     rename(!!col_name := Daily_Avg)
#   write_csv(daily_avg, file.path(output_path, file_name))
#   daily_avg
# }
# 
# # Calculate daily evaporation rate and export
# output_folder <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/2024'
# calculate_and_export_daily_avg(eva_rate, eva_date_time, output_folder, '2024_daily_average_evaporation_rates.csv', 'Evaporation Rate (mm/hr)')
# 
# # Water Temperature
# folderPath <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/Water Temp'
# water_temp_data <- read_and_combine_csv(folderPath, skip = 1, col_to_use = c(1, 2))
# 
# ggplot(water_temp_data, aes(x = Date_Time, y = Value)) +
#   geom_line(color = 'red') +
#   labs(x = 'Date', y = 'Temperature (C)', title = 'Water Temperature over Time') +
#   theme_minimal()
# 
# calculate_and_export_daily_avg(water_temp_data$Value, water_temp_data$Date_Time, output_folder, '2024_daily_average_water_temperature.csv', 'Water Temp (C)')
# 
# # Air Temperature
# folderPath <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/Air Temp'
# air_temp_data <- read_and_combine_csv(folderPath, skip = 1, col_to_use = c(1, 2))
# 
# ggplot(air_temp_data, aes(x = Date_Time, y = Value)) +
#   geom_line(color = 'green') +
#   labs(x = 'Date', y = 'Temperature (C)', title = 'Air Temperature over Time') +
#   theme_minimal()
# 
# calculate_and_export_daily_avg(air_temp_data$Value, air_temp_data$Date_Time, output_folder, '2024_daily_average_air_temperature.csv', 'Air Temp (C)')
# 
# # Relative Humidity
# folderPath <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/RH'
# rh_data <- read_and_combine_csv(folderPath, skip = 1, col_to_use = c(1, 2))
# 
# ggplot(rh_data, aes(x = Date_Time, y = Value)) +
#   geom_line(color = 'black') +
#   labs(x = 'Date', y = 'RH (%)', title = 'Relative Humidity over Time') +
#   theme_minimal()
# 
# calculate_and_export_daily_avg(rh_data$Value, rh_data$Date_Time, output_folder, '2024_daily_average_RH.csv', 'RH (%)')
