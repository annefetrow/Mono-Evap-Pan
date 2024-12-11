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

# Define global constants in an environment
globalVars <- new.env()
globalVars$RHO_W <- 1000       # Density of water in kg/m^3
globalVars$RHO_L <- 1060       # Density of liquid in kg/m^3
globalVars$SIGMA <- 5.67e-8    # Stefan-Boltzmann constant in W/(m^2 * K^4)
globalVars$P <- 101325         # Pressure in Pa
globalVars$Cp_w <- 3.7794e-3   # Specific heat of water in MJ/(kg * K)
globalVars$Cp_a <- 1.005e-3    # Specific heat of air in MJ/(kg * K)
globalVars$LAMBDA <- 2.45      # Latent heat in MJ/kg

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

# Function to calculate E
penman_calc <- function(delta, gamma, Rn, Ea) {
  G <- 0 # MJ/m2d
  E <- 1 / globalVars$LAMBDA * (delta * (Rn - G) + gamma * Ea) / ((delta + gamma) * globalVars$RHO_L)
  return(E)
}

# Function to calculate Rn
Rn_calc <- function(Ta, RH, R1) {
  r <- 0.05
  a <- 0.4
  b <- 0.274
  n_N <- 0.8
  ed_0_kPa <- water_surface_vapor_pressure(Ta, RH / 100)
  ed_4_kPa <- vapor_pressure_height(ed_0_kPa, Ta + 273.15, 4)
  
  # Calculate RB
  RB <- globalVars$SIGMA * (Ta + 273.15)^4 * (0.56 - 0.09 * sqrt(ed_4_kPa)) * (0.1 + 0.9 * n_N)
  
  # Calculate Rn
  Rn <- R1 * (1 - r) - RB
  return(Rn)
}

# Function for saturated water vapor pressure
saturated_water_vapor_pressure <- function(T) {
  p_sat <- 0.611 * exp(17.27 * T / (T + 237.3))
  return(p_sat)
}

# Function for water surface vapor pressure
water_surface_vapor_pressure <- function(T, RH) {
  p_sur <- 0.611 * exp(17.27 * T / (T + 237.3)) * RH
  return(p_sur)
}

# Function for vapor pressure at a certain height
vapor_pressure_height <- function(p0, T, height) {
  M <- 18.015 * 10^(-3)  # kg/mol
  gravity <- 9.81
  R <- 8.314
  ph <- p0 * exp(-1 * (M * gravity * height) / (R * T))
  return(ph)
}

# Function for salinity-adjusted vapor pressure
salinity_vapor_pressure <- function(T, S) {
  pw <- saturated_water_vapor_pressure(T)
  a1 <- -2.1609 * 10^(-4)
  a2 <- -3.5015 * 10^(-7)
  p <- pw * 10^(a1 * S + a2 * S^2)
  return(p)
}

# Function to calculate delta
delta_calc <- function(S, T) {
  dT <- 0.1
  delta <- (salinity_vapor_pressure(T + dT, S) - salinity_vapor_pressure(T - dT, S)) / (2 * dT)
  return(delta)
}

# Function for gamma calculation
salinity_sigma_calc <- function() {
  mu <- 0.622
  gamma <- (globalVars$Cp_a * globalVars$P / 1000) / (mu * globalVars$LAMBDA)
  return(gamma)
}

# Calculate delta
delta <- delta_calc(combined_data$Salinity_g_per_kg, combined_data$Lake_Air_Temperature_C)  # kPa/K

# Specific Heat
gamma <- salinity_sigma_calc()  # kPa/K

# Radiation calculation
Rn <- Rn_calc(combined_data$Lake_Air_Temperature_C, combined_data$RH_Percent, combined_data$Solar_Radiation_W_m2)  # W/m^2

# Check Ea with theoretical Ea
Ea_theo <- 6.43 * (0.18 + 0.55 * wind_speed_data$Wind_Speed_m_s) * 
  (saturated_water_vapor_pressure(combined_data$WaterTemp_C) - 
     water_surface_vapor_pressure(combined_data$AirTemp_C, combined_data$RH_Percent / 100))

# Calculate En
Ea_pan <- combined_data$Evaporation_Rate_mm_hr * globalVars$RHO_W * globalVars$LAMBDA * 24 / 1000  # MJ/m^2/day
E_pan <- penman_calc(delta, gamma, Rn * 0.0864, Ea_pan) * 1000  # mm/day
E_theo <- penman_calc(delta, gamma, Rn * 0.0864, Ea_theo) * 1000  # mm/day

ggplot(data = combined_data, aes(x = Date)) +
  geom_line(aes(y = E_pan, color = "Lake, from Pan Eva")) +
  geom_line(aes(y = E_theo, color = "Lake, from Theoretical Eva")) +
  geom_line(aes(y = Evaporation_Rate_mm_hr * 24, color = "Pan Eva")) +
  labs(y = "Eva Rate, mm/d", x = "Date") +
  scale_color_manual(values = c("Lake, from Pan Eva" = "blue",
                                "Lake, from Theoretical Eva" = "red",
                                "Pan Eva" = "green")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "top",       # Move legend to the top
        legend.direction = "horizontal")  # Arrange legend items horizontally


# Calculate Conversion Coefficients
C_pan <- E_pan / (combined_data$Evaporation_Rate_mm_hr * 24)
C_theo <- E_theo / (combined_data$Evaporation_Rate_mm_hr * 24)

# Add conversion coefficients to combined_data for plotting
combined_data$C_pan <- C_pan
combined_data$C_theo <- C_theo

ggplot(data = combined_data, aes(x = Date)) +
  geom_line(aes(y = C_pan, color = "From Pan Eva")) +
  geom_line(aes(y = C_theo, color = "From Theoretical Eva")) +
  labs(y = "Conversion Coefficient", x = "Date") +
  scale_color_manual(values = c("From Pan Eva" = "blue", "From Theoretical Eva" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "top",       # Move legend to the top
        legend.direction = "horizontal") +  # Arrange legend items horizontally
  labs(color = "Evaporation Source")  # Adds legend label


# Create a data frame with the calculated data
output_table <- data.frame(
  Date = combined_data$Date,
  Lake_From_Pan_Eva_mm_d = E_pan,
  Lake_From_Theoretical_Eva_mm_d = E_theo,
  Pan_Eva_mm_d = combined_data$Evaporation_Rate_mm_hr * 24
)

# Write the data frame to a CSV file
write.csv(output_table, "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/eva_estimate_Penman_RStudio.csv", row.names = FALSE)
