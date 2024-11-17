# Things to input:
# For each 

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

# Function to trim data to a common date range
trim_to_common_date_boundaries <- function(data1, time1, data2, time2) {
  # Find the maximum of the start dates and the minimum of the end dates
  common_start <- max(min(time1), min(time2))
  common_end <- min(max(time1), max(time2))
  
  # Create the common date range
  common_date_range <- c(common_start, common_end)
  
  # Trim data1 and data2 to this common date range
  time1_trimmed <- time1[time1 >= common_start & time1 <= common_end]
  data1_trimmed <- data1[time1 >= common_start & time1 <= common_end]
  
  time2_trimmed <- time2[time2 >= common_start & time2 <= common_end]
  data2_trimmed <- data2[time2 >= common_start & time2 <= common_end]
  
  # Return the common date range and trimmed data
  list(common_date_range = common_date_range,
       data1_trimmed = data1_trimmed,
       data2_trimmed = data2_trimmed,
       time1_trimmed = time1_trimmed,
       time2_trimmed = time2_trimmed)
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
all_data <- read_and_combine_csv(folderPath, 5, 8, 3)

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
output_folder <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output'
calculate_and_export_daily_avg(eva_rate$eva_rate, eva_rate$eva_date_time, output_folder, 'daily_average_evaporation_rates.csv', 'Evaporation_Rate_mm_hr')

# Water Temperature
folderPath <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/Water Temp'
water_temp_data <- read_and_combine_csv(folderPath, 1, 2, 3)

ggplot(water_temp_data, aes(x = timestamp, y = value)) +
  geom_line(color = 'red') +
  labs(x = 'Date', y = 'Temperature (C)', title = 'Water Temperature over Time') +
  theme_minimal()

calculate_and_export_daily_avg(water_temp_data$value, water_temp_data$timestamp, output_folder, 'daily_average_water_temperature.csv', 'WaterTemp_C')

# Air Temperature
folderPath <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/Air Temp and RH'
air_temp_data <- read_and_combine_csv(folderPath, 1, 2, 3)

ggplot(air_temp_data, aes(x = timestamp, y = value)) +
  geom_line(color = 'green') +
  labs(x = 'Date', y = 'Temperature (C)', title = 'Air Temperature over Time') +
  theme_minimal()

calculate_and_export_daily_avg(air_temp_data$value, air_temp_data$timestamp, output_folder, 'daily_average_air_temperature.csv', 'AirTemp_C')

# Relative Humidity
folderPath <- 'C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/Raw data – Evaporation Pan/Air Temp and RH'
rh_data <- read_and_combine_csv(folderPath, 1, 2, 4)

ggplot(rh_data, aes(x = timestamp, y = value)) +
  geom_line(color = 'black') +
  labs(x = 'Date', y = 'RH (%)', title = 'Relative Humidity over Time') +
  theme_minimal()

calculate_and_export_daily_avg(rh_data$value, rh_data$timestamp, output_folder, 'daily_average_RH.csv', 'RH_Percent')

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

# Plot function for monthly evaporation rate with air and water temperatures
plotMonthlyEvapWaterAir <- function(eva_rate, water_temp_data, air_temp_data) {
  
  # Trim the data to the same date range using the helper function
  trimmed_data <- trim_to_common_date_boundaries(eva_rate$eva_rate, eva_rate$eva_date_time, water_temp_data$value, water_temp_data$timestamp)
  common_date_range <- trimmed_data$common_date_range
  trimmed_eva_time <- trimmed_data$time1_trimmed
  trimmed_eva_rate <- trimmed_data$data1_trimmed
  trimmed_temp_time <- trimmed_data$time2_trimmed
  trimmed_temp <- trimmed_data$data2_trimmed
  
  # Filter air temperature data to the common date range
  filter_idx_air <- air_temp_data$timestamp >= common_date_range[1] & air_temp_data$timestamp <= common_date_range[2]
  trimmed_T_air_time <- air_temp_data$timestamp[filter_idx_air]
  trimmed_T_air <- air_temp_data$value[filter_idx_air]
  
  # Extract unique months in the common date range
  unique_months <- unique(month(trimmed_eva_time))
  num_months <- length(unique_months)
  
  # List to store ggplot objects for each subplot
  plots <- vector("list", num_months)
  
  # Loop over each month to create subplots
  for (i in 1:num_months) { 
    # Filter data for the current month
    month_idx <- month(trimmed_eva_time) == unique_months[i]
    month_data <- data.frame(
      eva_time = trimmed_eva_time[month_idx],
      eva_rate = trimmed_eva_rate[month_idx]
    )
    
    temp_month_idx <- month(trimmed_temp_time) == unique_months[i]
    temp_data <- data.frame(
      temp_time = trimmed_temp_time[temp_month_idx],
      temp_value = trimmed_temp[temp_month_idx]
    )
    
    air_month_idx <- month(trimmed_T_air_time) == unique_months[i]
    air_data <- data.frame(
      air_time = trimmed_T_air_time[air_month_idx],
      air_value = trimmed_T_air[air_month_idx]
    )
    
    # Convert Date_time to POSIXct format for month_data, temp_data, and air_data
    month_data$eva_time <- as.POSIXct(month_data$eva_time, format = "%Y-%m-%d %H:%M:%S")
    temp_data$temp_time <- as.POSIXct(temp_data$temp_time, format = "%Y-%m-%d %H:%M:%S")
    air_data$air_time <- as.POSIXct(air_data$air_time, format = "%Y-%m-%d %H:%M:%S")
    
    # Define the color mapping for the legend
    color_mapping <- c("Water Temp" = "red", "Air Temp" = "green")
    
    # Plot for the current month
    p <- ggplot() +
      # Bar plot for Evaporation Rate
      geom_bar(data = month_data, aes(x = eva_time, y = eva_rate, fill = "Evaporation Rate"), stat = "identity", alpha = 0.6) +
      # Line plot for Water Temperature (dashed)
      geom_line(data = temp_data, aes(x = temp_time, y = temp_value, color = "Water Temp"), linetype = "dashed") +
      # Line plot for Air Temperature (dashed)
      geom_line(data = air_data, aes(x = air_time, y = air_value, color = "Air Temp"), linetype = "dashed") +
      # Customize the y-axis labels and secondary axis for Water Temp and Air Temp
      scale_y_continuous(
        name = "Evaporation Rate (mm/hr)",
        sec.axis = sec_axis(~ ., name = "Temperature (°C)")
      ) +
      # Set title and x-axis label
      labs(title = paste("Evaporation Rate, Water Temp, and Air Temp - Month", unique_months[i]),
           x = "Date") +
      # Apply color mappings to the legend
      scale_color_manual(values = color_mapping) +
      scale_fill_manual(values = c("Evaporation Rate" = "blue")) +  # Fill color for Evaporation Rate bar
      # Apply minimal theme with customized axis labels and title
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = if (i == 1) "top" else "none",  # Display the legend at the top
        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
      )
    
    # Store plot in list
    plots[[i]] <- p
  }
  
  # Arrange all plots in a single figure
  grid.arrange(grobs = plots, ncol = 1)
}

plotMonthlyEvapWaterAir(eva_rate, water_temp_data, air_temp_data)

