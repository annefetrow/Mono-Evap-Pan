# Function to import water level and corresponding date and time from CSV
water_level_data_import <- function(filename) {
  water_level_file <- read.csv(file = filename, header = FALSE, skip = 2, sep = ',')
  h <- water_level_file[[3]]          # Water level data (column 3)
  date_time <- as.POSIXct(water_level_file[[2]], format = "%m/%d/%Y %H:%M")  # Date and time data (column 2)
  data.frame(date_time = date_time, h = h)
}

# Function to calculate total evaporation rate
eva_rate_calc <- function(h, date_time) {
  eva_rate <- c()
  eva_date_time <- c()
  
  for (i in 1:(length(h) - 2)) {
    window <- h[i:(i + 2)]
    duration <- as.numeric(difftime(date_time[i + 2], date_time[i], units = "hours"))
    if (all(diff(window) <= 0)) { # Check for descending water level
      rate <- (window[1] - window[3]) / duration
      eva_rate <- c(eva_rate, rate)
      eva_date_time <- c(eva_date_time, as.POSIXct(mean(c(date_time[i], date_time[i + 2]))))
    }
  }
  data.frame(eva_date_time = as.POSIXct(eva_date_time), eva_rate = eva_rate)
}

# Function to calculate daily evaporation rate
daily_eva_rate_calc <- function(eva_data) {
  
  # Get the unique dates to loop through
  unique_dates <- unique(eva_data$date_only)
  
  # Prepare a data frame to store the results
  results <- data.frame(Date = as.Date(character()), AverageEvaporationRateMmperHr = numeric())
  
  for (specific_date in unique_dates) {
    # Identify the rows with the matching date
    daily_data <- eva_data %>% filter(date_only == specific_date)
    
    # Extract the evaporation rates for the current date while excluding zeros
    daily_eva_rate <- daily_data %>% 
      filter(eva_rate != 0) %>%
      pull(eva_rate)
    
    # Calculate the average evaporation rate for the current date
    average_value <- mean(daily_eva_rate, na.rm = TRUE)  # na.rm = TRUE to remove NA values
    
    # Save the results in the data frame
    results <- rbind(results, data.frame(Date = as.Date(specific_date), AverageEvaporationRateMmperHr = average_value))
  }
  
  # Write the results to a CSV file
  write.csv(results, "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/daily_average_evaporation_rates.csv", row.names = FALSE)
  
  # Plot the results
  gg_daily <- ggplot(results, aes(x = Date, y = AverageEvaporationRateMmperHr)) +
    geom_col() +
    labs(title = "Daily Evaporation Rate", x = "Date", y = "Daily Evaporation Rate (mm/hr)") +
    theme_minimal()
  
  return(gg_daily)
  
}

# Function to calculate monthly evaporation rate
monthly_eva_rate_clac <- function(eva_data) {
  
  # Now prepare the data frame for grouping by year and month
  eva_data <- eva_data %>%
    mutate(
      year = year(date_only), 
      month = month(date_only),
      year_month = floor_date(date_only, "month")
    )
  
  # Calculate monthly average evaporation rate
  monthly_avg_evap <- eva_data %>%
    group_by(year_month) %>%
    summarise(
      monthly_avg_evap = mean(eva_rate, na.rm = TRUE)  # Use na.rm to ignore NA values
    ) %>%
    ungroup()
  
  # Save the monthly averages to a CSV file
  write.csv(monthly_avg_evap, "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/monthly_evaporation_averages.csv", row.names = FALSE)
  
  # Plot the results
  gg_monthly <- ggplot(monthly_avg_evap, aes(x = year_month, y = monthly_avg_evap)) +
    geom_col() +
    labs(title = "Monthly Evaporation Rate", x = "Date", y = "Monthly Evaporation Rate (mm/hr)") +
    theme_minimal()
  
  return(gg_monthly)
  
}