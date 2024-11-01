# Initialization
cat("\014")
rm(list = ls())
graphics.off()

# Load necessary libraries
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(vroom)
library(gridExtra)

# Function to parse datetime flexibly for mdy and ymd formats
parse_datetime_flexible <- function(ori_datetime_str) {
  # Attempt initial parsing with mdy_hms format
  parsed_datetime <- mdy_hms(ori_datetime_str, quiet = TRUE)
  
  # If parsing with mdy_hms fails, try ymd_hms
  if (any(is.na(parsed_datetime))) {
    parsed_datetime <- ymd_hms(ori_datetime_str, quiet = TRUE)
  }
  
  # If both mdy_hms and ymd_hms failed, append ":00" to entries that resulted in NA
  if (any(is.na(parsed_datetime))) {
    datetime_str <- ori_datetime_str
    datetime_str[is.na(parsed_datetime)] <- paste0(ori_datetime_str[is.na(parsed_datetime)], ":00")
    
    # Try parsing again with mdy_hms and then ymd_hms on modified datetime strings
    parsed_datetime <- mdy_hms(datetime_str, quiet = TRUE)
    if (any(is.na(parsed_datetime))) {
      parsed_datetime <- ymd_hms(datetime_str, quiet = TRUE)
    }
  }
  
  # Return only YY:MM:DD HH:MM format
  return(format(parsed_datetime, "%Y-%m-%d %H:%M"))
}


data_import <- function(file_name, data_column) {
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
  data_table <- data_table[, c(2, data_column)]
  
  # Optionally rename columns for easier access (if desired)
  colnames(data_table) <- c("Timestamp", "Value")  # Renaming columns (optional)
  data_table$Timestamp <- parse_datetime_flexible(data_table$Timestamp) # Convert data type to datetime
  
  # Display the resulting data
  print(data_table)
  
  # Return the data and timestamps as a list
  return(list(value = data_table$Value, timestamp = data_table$Timestamp))
}

read_and_combine_csv <- function(folder_path, start_file, end_file, data_column) {
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
    data <- data_import(csv_files[k], data_column)
    data$timestamp <- parse_datetime_flexible(data$timestamp)
    
    # Combine the data
    all_data <- bind_rows(all_data, data)
  }
  
  # Remove duplicates and NaN values
  unique_data <- all_data %>%
    distinct() %>%
    filter(!is.na(value) & !is.na(timestamp))
  
  unique_data$timestamp <- as.POSIXct(unique_data$timestamp, format = "%Y-%m-%d %H:%M")
  
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
all_data <- read_and_combine_csv(folderPath, 5, 7, 3)

# Plot water level over time
ggplot(all_data, aes(x = timestamp, y = value)) +
  geom_line(color = 'blue') +
  labs(x = 'Date and Time', y = 'Water Level (mm)', title = 'Water Level over Time') +
  theme_minimal()

# Evaporation Rate Calculation Function
eva_rate <- eva_rate_cal(all_data$value, all_data$timestamp)
eva_rate <- list(
  eva_rate = eva_rate$eva_rate[eva_rate$eva_rate < 3],
  eva_date_time = eva_rate$eva_date_time[eva_rate$eva_rate < 3]
)
eva_rate_avg <- mean(eva_rate$eva_rate, na.rm = TRUE) # exclude outliers

# Plot evaporation rate over time
ggplot() +
  geom_bar(aes(x = eva_rate$eva_date_time, y = eva_rate$eva_rate), stat = "identity", fill = "blue") +
  geom_hline(yintercept = eva_rate_avg, linetype = "dashed", color = "blue") +
  annotate("text", x = max(eva_rate$eva_date_time) - days(1), y = max(eva_rate$eva_rate) + 0.2,
           label = paste("Average:", round(eva_rate_avg, 2), "mm/hr"), hjust = 1) +
  labs(x = 'Date', y = 'Eva Rate (mm/hr)', title = 'Evaporation Rate Over Time') +
  theme_minimal()

# Function to calculate daily averages and export to CSV
calculate_and_export_daily_avg <- function(values, date_times, output_path, file_name, col_name) {
  daily_avg <- tibble(Date = as.Date(date_times), Value = values) %>%
    group_by(Date) %>%
    summarize(Daily_Avg = mean(Value, na.rm = TRUE)) %>%
    rename(!!col_name := Daily_Avg)
  write_csv(daily_avg, file.path(output_path, file_name))
  daily_avg
}

# Calculate daily evaporation rate and export
output_folder <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/2024'
calculate_and_export_daily_avg(eva_rate$eva_rate, eva_rate$eva_date_time, output_folder, '2024_daily_average_evaporation_rates.csv', 'Evaporation Rate (mm/hr)')

# Water Temperature
folderPath <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/Water Temp'
water_temp_data <- read_and_combine_csv(folderPath, 1, 2, 3)

ggplot(water_temp_data, aes(x = timestamp, y = value)) +
  geom_line(color = 'red') +
  labs(x = 'Date', y = 'Temperature (C)', title = 'Water Temperature over Time') +
  theme_minimal()

calculate_and_export_daily_avg(water_temp_data$value, water_temp_data$timestamp, output_folder, '2024_daily_average_water_temperature.csv', 'Water Temp (C)')

# Air Temperature
folderPath <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/Air Temp'
air_temp_data <- read_and_combine_csv(folderPath, 1, 2, 3)

ggplot(air_temp_data, aes(x = timestamp, y = value)) +
  geom_line(color = 'green') +
  labs(x = 'Date', y = 'Temperature (C)', title = 'Air Temperature over Time') +
  theme_minimal()

calculate_and_export_daily_avg(air_temp_data$value, air_temp_data$timestamp, output_folder, '2024_daily_average_air_temperature.csv', 'Air Temp (C)')

# Relative Humidity
folderPath <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/RH'
rh_data <- read_and_combine_csv(folderPath, 1, 2, 4)

ggplot(rh_data, aes(x = timestamp, y = value)) +
  geom_line(color = 'black') +
  labs(x = 'Date', y = 'RH (%)', title = 'Relative Humidity over Time') +
  theme_minimal()

calculate_and_export_daily_avg(rh_data$value, rh_data$timestamp, output_folder, '2024_daily_average_RH.csv', 'RH (%)')

plotByMonth <- function(T_water_all_date_time, T_water_all, T_air_all_date_time, T_air_all, RH_all_date_time, RH_all) {

  # Convert date-time columns to Date format if they are not already
  T_water_all_date <- as.Date(T_water_all_date_time)
  T_air_all_date <- as.Date(T_air_all_date_time)
  RH_all_date <- as.Date(RH_all_date_time)

  # Get unique months and years from the datetime arrays
  months <- format(T_water_all_date, "%m")
  years <- format(T_water_all_date, "%Y")
  uniqueMonths <- unique(data.frame(year = years, month = months))

  # List to hold plots
  plots <- list()

  # Loop through each unique month and year
  for (i in 1:nrow(uniqueMonths)) {
    currentMonth <- uniqueMonths$month[i]
    currentYear <- uniqueMonths$year[i]

    # Filter data for the current month and year
    idxWater <- which(format(T_water_all_date, "%Y-%m") == paste(currentYear, currentMonth, sep = "-"))
    idxAir <- which(format(T_air_all_date, "%Y-%m") == paste(currentYear, currentMonth, sep = "-"))
    idxRH <- which(format(RH_all_date, "%Y-%m") == paste(currentYear, currentMonth, sep = "-"))

    # Create a data frame for plotting
    plot_data <- data.frame(
      Date_time = c(T_water_all_date_time[idxWater], T_air_all_date_time[idxAir], RH_all_date_time[idxRH]),
      Value = c(T_water_all[idxWater], T_air_all[idxAir], RH_all[idxRH]),
      Variable = factor(
        rep(c("Water Temp", "Air Temp", "RH"),
            times = c(length(idxWater), length(idxAir), length(idxRH)))
      )
    )

    # Make sure the format is correct for Date_time
    plot_data$Date_time <- as.POSIXct(plot_data$Date_time, format = "%Y-%m-%d %H:%M:%S")

    # Define your color mapping
    color_mapping <- c("Water Temp" = "red", "Air Temp" = "green", "RH" = "black")

    p <- ggplot(plot_data, aes(x = Date_time, y = Value, color = Variable)) +
      geom_line() +
      scale_y_continuous(
        name = "Temp (°C)",
        sec.axis = sec_axis(~ ., name = "RH (%)")
      ) +
      labs(title = paste("Month:", currentYear, "-", currentMonth)) +
      scale_color_manual(values = color_mapping) +  # Set manual colors
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = if (i == 1) "top" else "none"
      )

    # Add plot to the list
    plots[[i]] <- p
  }

  # Arrange all plots in a single figure with multiple rows
  grid.arrange(grobs = plots, ncol = 1, top = "Water Temperature, Air Temperature, and Relative Humidity by Month")
}


# Plot Water Temperature, Air Temperature, and Relative Humidity
plotByMonth(water_temp_data$timestamp, water_temp_data$value, air_temp_data$timestamp, air_temp_data$value, rh_data$timestamp, rh_data$value)
