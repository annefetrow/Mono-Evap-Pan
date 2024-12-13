# Initialization
cat("\014")
rm(list = ls())
graphics.off()

# Assume the climate station is measured at 4 m height
h <- 4

# Load necessary libraries
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(vroom)
library(gridExtra)
library(readxl)
library(tidyr)

# Define a global base directory, please change it when download the code to your personal desktop
global_save_dir <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/plots"

# Define folder path and date range
folder_path <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output"
start_date <- as.Date("2024-07-05")
end_date <- as.Date("2024-08-13")

read_and_crop_data <- function(folder_path, start_date, end_date) {
  # This function reads all the .csv files in the specified folder,
  # crops them to the specified date range, and combines them into a single data frame.
  # Inputs:
  # - folder_path: the folder where the data files are located
  # - start_date: the earliest date to include in the output (Date class)
  # - end_date: the latest date to include in the output (Date class)
  # Output:
  # - combined_data: a data frame with the combined data, where the first column is 'Date'
  
  # List all .csv files in the folder
  file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Initialize an empty list to store each file's data
  data_list <- list()
  
  for (file_path in file_list) {
    # Read each file
    data <- read.csv(file_path)
    
    # Ensure the Date column is in Date format
    data$Date <- as.POSIXct(data$Date, format = "%Y-%m-%d")
    
    # Extract the variable name (second column name)
    variable_name <- colnames(data)[2]
    
    # Crop the data to the specified date range
    cropped_data <- data[data$Date >= start_date & data$Date <= end_date, ]
    
    # Rename the second column with the variable name
    colnames(cropped_data)[2] <- variable_name
    
    # Add the cropped data to the list
    data_list[[variable_name]] <- cropped_data
  }
  
  # Combine all data frames in the list by 'Date' column using a full join
  combined_data <- Reduce(function(x, y) merge(x, y, by = "Date", all = TRUE), data_list)
  
  # Ensure that the combined data is sorted by date
  combined_data <- combined_data[order(combined_data$Date), ]
  
  return(combined_data)
}

# Call the function to read and crop data
combined_data <- read_and_crop_data(folder_path, start_date, end_date)

# Assuming combined_data is a data frame with a Date column
# Initialize the Salinity_g_per_kg column with zeros
combined_data$Salinity_g_per_kg <- 0  # Set initial values to 0

# Define the date boundaries for setting different salinity values
start_date <- as.Date("2024-07-20")
mid_date <- as.Date("2024-08-15")

# Set salinity values based on date conditions
combined_data$Salinity_g_per_kg <- ifelse(
  combined_data$Date < start_date, 75,  # Salinity for earlier dates
  ifelse(
    combined_data$Date < mid_date, 81,  # Salinity for middle dates
    85  # Salinity for later dates
  )
)


match_table_sizes <- function(table1, table2) {
  # Ensure that the Date columns are in Date format
  table1$Date <- as.Date(table1$Date)
  table2$Date <- as.Date(table2$Date)
  
  # Find the common date range (intersection of dates)
  common_start_date <- max(min(table1$Date), min(table2$Date))
  common_end_date <- min(max(table1$Date), max(table2$Date))
  
  # Crop table1 to the common date range
  cropped_table1 <- table1 %>%
    filter(Date >= common_start_date & Date <= common_end_date)
  
  # Crop table2 to the common date range
  cropped_table2 <- table2 %>%
    filter(Date >= common_start_date & Date <= common_end_date)
  
  # Display messages if either table was cropped
  if (nrow(table1) != nrow(cropped_table1)) {
    message("Table 1 has been cropped to match Table 2.")
  }
  if (nrow(table2) != nrow(cropped_table2)) {
    message("Table 2 has been cropped to match Table 1.")
  }
  
  # Return the cropped tables as a list
  return(list(cropped_table1 = cropped_table1, cropped_table2 = cropped_table2))
}


read_and_average_air_temperature <- function(file_path) {
  # Read the data from the specified Excel file
  data <- read_excel(file_path)
  
  # Extract the date and temperature columns
  # Assuming the first column contains date and time, and the last column contains temperature
  date_time <- data[[1]]  # First column
  lake_air_temperature <- data[[ncol(data)]]  # Last column
  
  # Convert date_time to POSIXct format
  date_time <- ymd_h(date_time)  # Format as "yyyy-MM-dd HH" using lubridate
  
  # Filter out NA values from lake air temperature
  valid_data <- !is.na(lake_air_temperature)
  date_time <- date_time[valid_data]
  lake_air_temperature <- lake_air_temperature[valid_data]
  
  # Group by date (day only) and calculate daily average temperature
  daily_air_temp_data <- data.frame(Date = as.Date(date_time), 
                                    Lake_Air_Temperature_C = lake_air_temperature) %>%
    group_by(Date) %>%
    summarize(Lake_Air_Temperature_C = mean(Lake_Air_Temperature_C, na.rm = TRUE)) %>%
    ungroup()
  
  return(daily_air_temp_data)
}

export_plot_to_png <- function(plot, file_name, save_dir = global_save_dir, width = 1200, base_height = 200, res = 150, num_rows = 1, combined = FALSE) {
  
  # Adjust height only if it's a combined plot and there are multiple rows
  if (combined && num_rows > 1) {
    # Add extra space for the x-axis labels on the first plot (top plot)
    extra_space_for_labels <- 300  # Adjust this value as needed
  } else {
    extra_space_for_labels <- 0  # No extra space needed for non-combined or single-row plots
    base_height = 400
  }
  
  # Calculate the total height
  height <- base_height * num_rows + extra_space_for_labels
  
  # Create the directory if it doesn't exist
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Construct the full file path
  file_path <- file.path(save_dir, file_name)
  
  # Open a PNG device with calculated height
  png(filename = file_path, width = width, height = height, res = res)
  
  # Check if the input is a single plot or a grid object
  if (inherits(plot, "ggplot")) {
    # For a single ggplot, use print
    print(plot)
  } else {
    # For grid-arranged objects, use grid.draw
    grid.draw(plot)
  }
  
  # Close the device
  dev.off()
  
  cat("Plot exported to:", file_path, "\n")
}

# Load air temperature data
file_path <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/2024_station_data/Air Temperature C.xlsx"
air_temp_data <- read_and_average_air_temperature(file_path)

# Ensure combined_data and air_temp_data have Date columns in datetime format
combined_data$Date <- as_date(combined_data$Date)
air_temp_data$Date <- as_date(air_temp_data$Date)

# Align the date ranges between combined_data and air_temp_data
common_start_date <- max(min(combined_data$Date), min(air_temp_data$Date))
common_end_date <- min(max(combined_data$Date), max(air_temp_data$Date))

combined_data <- combined_data %>% filter(Date >= common_start_date & Date <= common_end_date)
air_temp_data <- air_temp_data %>% filter(Date >= common_start_date & Date <= common_end_date)

# Add the Lake Air Temperature from air_temp_data to combined_data
combined_data <- combined_data %>%
  left_join(air_temp_data %>% select(Date, Lake_Air_Temperature_C), by = "Date")

# Define function to read, process, and average solar radiation data
read_and_average_solar_radiation <- function(file_path) {
  # Read the data from the specified Excel file
  data <- read_excel(file_path)
  
  # Extract the date and solar radiation columns
  # Assuming the first column contains date and time, and the last column contains solar radiation
  date_time <- data[[1]]  # First column
  solar_radiation <- data[[ncol(data)]]  # Last column
  
  # Convert date_time to POSIXct format
  date_time <- ymd_h(date_time)  # Format as "yyyy-MM-dd HH" using lubridate
  
  # Filter out NA values from solar radiation
  valid_data <- !is.na(solar_radiation)
  date_time <- date_time[valid_data]
  solar_radiation <- solar_radiation[valid_data]
  
  # Group by date (day only) and calculate daily average solar radiation
  daily_solar_radiation_data <- data.frame(Date = as.Date(date_time), 
                                           Solar_Radiation_W_m2 = solar_radiation) %>%
    group_by(Date) %>%
    summarize(Solar_Radiation_W_m2 = mean(Solar_Radiation_W_m2, na.rm = TRUE)) %>%
    ungroup()
  
  return(daily_solar_radiation_data)
}

# Example usage with solar radiation data
file_path <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/2024_station_data/Radiation.xlsx"

# Assume combined_data is pre-loaded or defined
solar_radiation_data <- read_and_average_solar_radiation(file_path)

# Ensure combined_data and solar_radiation_data have Date columns in datetime format
combined_data$Date <- as_date(combined_data$Date)
solar_radiation_data$Date <- as_date(solar_radiation_data$Date)

# Align the date ranges between combined_data and solar_radiation_data
common_start_date <- max(min(combined_data$Date), min(solar_radiation_data$Date))
common_end_date <- min(max(combined_data$Date), max(solar_radiation_data$Date))

combined_data <- combined_data %>% filter(Date >= common_start_date & Date <= common_end_date)
solar_radiation_data <- solar_radiation_data %>% filter(Date >= common_start_date & Date <= common_end_date)

# Add the Solar Radiation data to combined_data
combined_data <- combined_data %>%
  left_join(solar_radiation_data %>% select(Date, Solar_Radiation_W_m2), by = "Date")

# Define function to read, process, and average wind speed data
read_and_average_wind_speed <- function(file_path) {
  # Read the data from the specified Excel file
  data <- read_excel(file_path)
  
  # Extract the date and wind speed columns
  # Assuming the first column contains date and time, and the last column contains wind speed
  date_time <- data[[1]]  # First column
  wind_speed <- data[[ncol(data)]]  # Last column
  
  # Convert date_time to POSIXct format
  date_time <- ymd_h(date_time)  # Format as "yyyy-MM-dd HH" using lubridate
  
  # Filter out NA values from wind speed
  valid_data <- !is.na(wind_speed)
  date_time <- date_time[valid_data]
  wind_speed <- wind_speed[valid_data]
  
  # Group by date (day only) and calculate daily average wind speed
  daily_wind_speed_data <- data.frame(Date = as.Date(date_time), 
                                      Wind_Speed_m_s = wind_speed) %>%
    group_by(Date) %>%
    summarize(Wind_Speed_m_s = mean(Wind_Speed_m_s, na.rm = TRUE)) %>%
    ungroup()
  
  return(daily_wind_speed_data)
}

export_plot_to_png <- function(plot, file_name, save_dir = global_save_dir, width = 1200, base_height = 200, res = 150, num_rows = 1, combined = FALSE) {
  
  # Adjust height only if it's a combined plot and there are multiple rows
  if (combined && num_rows > 1) {
    # Add extra space for the x-axis labels on the first plot (top plot)
    extra_space_for_labels <- 300  # Adjust this value as needed
  } else {
    extra_space_for_labels <- 0  # No extra space needed for non-combined or single-row plots
    base_height = 400
  }
  
  # Calculate the total height
  height <- base_height * num_rows + extra_space_for_labels
  
  # Create the directory if it doesn't exist
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Construct the full file path
  file_path <- file.path(save_dir, file_name)
  
  # Open a PNG device with calculated height
  png(filename = file_path, width = width, height = height, res = res)
  
  # Check if the input is a single plot or a grid object
  if (inherits(plot, "ggplot")) {
    # For a single ggplot, use print
    print(plot)
  } else {
    # For grid-arranged objects, use grid.draw
    grid.draw(plot)
  }
  
  # Close the device
  dev.off()
  
  cat("Plot exported to:", file_path, "\n")
}

# Example usage with wind speed data
file_path <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/2024_station_data/Wind Speed m_s.xlsx"

# Assume combined_data is pre-loaded or defined
wind_speed_data <- read_and_average_wind_speed(file_path)

# Ensure combined_data and wind_speed_data have Date columns in datetime format
combined_data$Date <- as_date(combined_data$Date)
wind_speed_data$Date <- as_date(wind_speed_data$Date)

# Align the date ranges between combined_data and wind_speed_data
common_start_date <- max(min(combined_data$Date), min(wind_speed_data$Date))
common_end_date <- min(max(combined_data$Date), max(wind_speed_data$Date))

combined_data <- combined_data %>% filter(Date >= common_start_date & Date <= common_end_date)
wind_speed_data <- wind_speed_data %>% filter(Date >= common_start_date & Date <= common_end_date)

# Add the Wind Speed data to combined_data
combined_data <- combined_data %>%
  left_join(wind_speed_data %>% select(Date, Wind_Speed_m_s), by = "Date")

# Monthly Lake Temperature Data
# Create a data frame with month names and average temperatures from Anne's field measurements and online database
MonthNames <- c("January", "February", "March", "April", "May", "June", 
                "July", "August", "September", "October", "November", "December")

# Average temperature should be updated when new data comes in
AverageTemperatures_C <- (c(30, 32, 35, 45, 55, 60, 69.6, 71.8, 62.275, 52.75, 49.8, 35) - 32) * 5 / 9

# Create a data frame for the monthly average temperatures
monthly_avg_temp_df <- data.frame(MonthNames, AverageTemperatures_C)

# Assuming combined_data is already a data frame with a 'Date' column in Date format
# Preallocate the lake temperature column in combined_data
combined_data$Lake_Temperature_C <- rep(NA, nrow(combined_data))

# Loop through combined_data and assign lake temperature based on the month
combined_data$Lake_Temperature_C <- sapply(combined_data$Date, function(date) {
  # Get the month index for the current date in combined_data
  current_month <- as.integer(format(date, "%m"))
  
  # Assign the corresponding average lake temperature based on the month
  monthly_avg_temp_df$AverageTemperatures_C[current_month]
})

# Display the first few rows of the updated combined_data
head(combined_data[, c("Date", "Lake_Temperature_C")])

# Function for saturated water vapor pressure
# p_sat in kPa, T in °C
saturated_water_vapor_pressure <- function(T) {
  p_sat <- 0.611 * exp(17.27 * T / (T + 237.3))
  
  # Antoine Equation (optional)
  # A <- 8.07131
  # B <- 1730.63
  # C <- 233.426
  # p_sat <- 10^((A - B / (T + C))) * 0.133322  # kPa
  
  return(p_sat)
}

# Function for water surface vapor pressure
# p_sur in kPa, T in °C, RH as a fraction (e.g., 0.6 for 60% RH)
water_surface_vapor_pressure <- function(T, RH) {
  p_sur <- 0.611 * exp(17.27 * T / (T + 237.3)) * RH
  
  # Antoine Equation (optional)
  # A <- 8.07131
  # B <- 1730.63
  # C <- 233.426
  # p_sur <- 10^((A - B / (T + C))) * 0.133322 * RH  # kPa
  
  return(p_sur)
}

# Function for vapor pressure at a given height
# p0 is initial pressure in kPa, T in °C, height in m
vapor_pressure_height <- function(p0, T, height) {
  M <- 18.015e-3  # kg/mol
  gravity <- 9.81  # m/s^2
  R <- 8.314  # J/(mol*K)
  
  # Convert T from °C to K for calculation
  T_K <- T + 273.15
  
  ph <- p0 * exp(-1 * (M * gravity * height) / (R * T_K))
  return(ph)
}

# Function for vapor pressure with salinity effect
# T in °C, S in g/kg (salinity)
salinity_vapor_pressure <- function(T, S) {
  pw <- saturated_water_vapor_pressure(T)
  
  a1 <- -2.1609e-4
  a2 <- -3.5015e-7
  
  p <- pw * 10^(a1 * S + a2 * S^2)
  return(p)
}

# Calculate Saturation Vapor Pressure
e_p_kPa <- saturated_water_vapor_pressure(combined_data$WaterTemp_C)
e_p_0_kPa <- water_surface_vapor_pressure(combined_data$AirTemp_C, combined_data$RH_Percent / 100)
e_p_4_kPa <- vapor_pressure_height(e_p_0_kPa, combined_data$AirTemp_C + 273.15, h)

e_l_kPa <- salinity_vapor_pressure(combined_data$Lake_Temperature_C, combined_data$Salinity_g_per_kg)
e_l_0_kPa <- water_surface_vapor_pressure(combined_data$Lake_Air_Temperature_C, combined_data$RH_Percent / 100)
e_l_4_kPa <- vapor_pressure_height(e_l_0_kPa, combined_data$Lake_Air_Temperature_C + 273.15, h)

# Calculate Lake Evaporation
E_lake <- ((e_l_kPa - e_l_4_kPa) / (e_p_kPa - e_p_4_kPa)) * combined_data$Evaporation_Rate_mm_hr
coefficient <- ((e_l_kPa - e_l_4_kPa) / (e_p_kPa - e_p_4_kPa))

# Create a data frame for ggplot
coefficient_data <- data.frame(
  Date = combined_data$Date,
  Coefficient = coefficient
)

# Plot using ggplot
p <- ggplot(coefficient_data, aes(x = Date, y = Coefficient)) +
  geom_line(color = "blue") +
  labs(title = "Conversion Coefficient Over Time",
       x = "Date",
       y = "Conversion Coefficient") +
  theme_minimal()

export_plot_to_png(p, file_name = "Webb_Coefficient.png", num_rows = 1, combined = FALSE)

# Create a data frame for plotting
vapor_pressure_data <- data.frame(
  Date = combined_data$Date,  # Assuming 'Date' is in combined_data
  e_p_4_kPa = e_p_4_kPa,
  e_p_kPa = e_p_kPa,
  e_l_kPa = e_l_kPa,
  e_l_4_kPa = e_l_4_kPa
)

# Pivot to long format for ggplot
vapor_pressure_long <- vapor_pressure_data %>%
  pivot_longer(cols = c(e_p_4_kPa, e_p_kPa, e_l_kPa, e_l_4_kPa), 
               names_to = "VaporPressure", 
               values_to = "Value")

# Ensure VaporPressure column is a factor with consistent levels
vapor_pressure_long$VaporPressure <- factor(
  vapor_pressure_long$VaporPressure, 
  levels = c("e_p_4_kPa", "e_p_kPa", "e_l_kPa", "e_l_4_kPa") # Matches scale_color_manual
)

# Create the plot
p <- ggplot(vapor_pressure_long, aes(x = Date, y = Value, color = VaporPressure)) +
  geom_line() +
  labs(title = "Vapor Pressures Over Time",
       x = "Date",
       y = "Vapor Pressure (kPa)") +
  scale_color_manual(values = c("e_p_4_kPa" = "blue", 
                                "e_p_kPa" = "red", 
                                "e_l_kPa" = "green", 
                                "e_l_4_kPa" = "purple"),
                     labels = c("Pan Air, 4m Above", "Pan Water", "Lake Water", "Lake Air, 4m Above")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "top")

export_plot_to_png(p, file_name = "Webb_Vapor_Pressure.png", num_rows = 1, combined = FALSE)

# Pivot the temperature data to long format
temperature_data <- combined_data %>%
  select(Date, Lake_Temperature_C, WaterTemp_C, AirTemp_C, Lake_Air_Temperature_C) %>%
  pivot_longer(cols = c(Lake_Temperature_C, WaterTemp_C, AirTemp_C, Lake_Air_Temperature_C), 
               names_to = "TemperatureType", 
               values_to = "Temperature")
# Pivot the temperature data to long format
temperature_data <- combined_data %>%
  select(Date, Lake_Temperature_C, WaterTemp_C, AirTemp_C, Lake_Air_Temperature_C) %>%
  pivot_longer(cols = c(Lake_Temperature_C, WaterTemp_C, AirTemp_C, Lake_Air_Temperature_C), 
               names_to = "TemperatureType", 
               values_to = "Temperature")

# Set the order for the legend using factor levels
temperature_data$TemperatureType <- factor(
  temperature_data$TemperatureType, 
  levels = c("AirTemp_C", "WaterTemp_C", "Lake_Temperature_C", "Lake_Air_Temperature_C") # Desired order
)

# Create the ggplot
p <- ggplot(temperature_data, aes(x = Date, y = Temperature, color = TemperatureType)) +
  geom_line() +
  labs(title = "Temperature Over Time",
       x = "Date",
       y = "Temperature (°C)") +
  scale_color_manual(
    values = c("AirTemp_C" = "blue", 
               "WaterTemp_C" = "red", 
               "Lake_Temperature_C" = "green", 
               "Lake_Air_Temperature_C" = "purple"),
    labels = c("Pan Air, 4m Above", "Pan Water", "Lake Water", "Lake Air, 4m Above")
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

export_plot_to_png(p, file_name = "Webb_Temperature.png", num_rows = 1, combined = FALSE)

# Create a data frame with calculated coefficient
output_coeff <- data.frame(
  Date = coefficient_data$Date,
  Coeff = coefficient_data$Coefficient
)

# Write the data frame to a CSV file
write.csv(output_coeff, "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/eva_coeff_Webb_RStudio.csv", row.names = FALSE)