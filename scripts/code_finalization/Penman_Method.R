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
library(readxl)

# Define global constants
RHO_W <- 1000        # Density of water in kg/m^3
RHO_L <- 1060        # Density of liquid in kg/m^3
SIGMA <- 5.67e-8     # Stefan-Boltzmann constant in W/(m^2 * K^4)
P <- 101325          # Pressure in Pa
Cp_w <- 3.7794e-3    # Specific heat of water in MJ/(kg * K)
Cp_a <- 1.005e-3     # Specific heat of air in MJ/(kg * K)
LAMBDA <- 2.45       # Latent heat in MJ/kg



# Define folder path and date range
folder_path <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/2024"
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

