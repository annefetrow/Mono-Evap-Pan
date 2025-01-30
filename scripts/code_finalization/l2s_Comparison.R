# Initialize
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
library(grid)
library(cowplot)

# Define a global base directory, please change it when download the code to your personal desktop
global_save_dir <- "C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/plots"

export_plot_to_png <- function(plot, file_name, save_dir = global_save_dir, width = 1200, base_height = 200, res = 150, num_rows = 1, combined = FALSE) {
  
  # Adjust height only if it's a combined plot and there are multiple rows
  if (combined && num_rows > 1) {
    # Add extra space for the x-axis labels on the first plot (top plot)
    extra_space_for_labels <- 300  # Adjust this value as needed
  } else {
    extra_space_for_labels <- 0  # No extra space needed for non-combined or single-row plots
    base_height = 500
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

# Load the CSV file
data <- read.csv('C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/data/l2s/Joker_L2SWBM_jSum_Formatted.csv')

# Convert the 'Date' column to Date type
data$Date <- as.Date(data$Date, format = "%Y-%m-%d")

# Extract the relevant columns (evaporation rate, upper boundary, lower boundary)
l2s_evaporation_rate <- 304.8 / 30.44 * data[, 17]
l2s_upper_boundary <- 304.8 / 30.44 * data[, 19]
l2s_lower_boundary <- 304.8 / 30.44 * data[, 23]

# Create a new data frame with the necessary columns
tbl <- data.frame(Date = data$Date, l2s_evaporation_rate, l2s_upper_boundary, l2s_lower_boundary)

# Extract the month from the Date column
tbl$Month <- month(tbl$Date)

# Calculate the monthly averages across all years using dplyr
l2s_monthly_avg <- tbl %>%
  group_by(Month) %>%
  summarise(
    monthly_avg_evaporation_rate = mean(l2s_evaporation_rate, na.rm = TRUE),
    monthly_avg_upper_boundary = mean(l2s_upper_boundary, na.rm = TRUE),
    monthly_avg_lower_boundary = mean(l2s_lower_boundary, na.rm = TRUE)
  )

# Create Date column for the monthly averages (using an arbitrary year, e.g., 2000)
l2s_monthly_avg$Date <- as.Date(paste(2000, l2s_monthly_avg$Month, 1, sep = "-"))

# View the resulting data frame
print(l2s_monthly_avg)


# Load the CSV file with stringsAsFactors = FALSE to avoid factor conversion
data <- read.csv('C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/eva_estimate_Penman_RStudio.csv', stringsAsFactors = FALSE)

# Check the first few rows of the Date column to ensure it's in the correct format
head(data$Date)

# Convert the 'Date' column to Date type (YYYY-MM-DD format)
data$Date <- as.Date(data$Date, format = "%Y-%m-%d")  # YYYY-MM-DD format

# Extract the relevant columns
dates <- data$Date
lake_from_pan_eva <- data$Lake_From_Pan_Eva_mm_d
lake_from_theoretical_eva <- data$Lake_From_Theoretical_Eva_mm_d
pan_eva <- data$Pan_Eva_mm_d

# Create a data frame with the necessary columns
tbl <- data.frame(dates, lake_from_pan_eva, lake_from_theoretical_eva, pan_eva)

# Extract the year and month from the dates
tbl$Year <- year(tbl$dates)
tbl$Month <- month(tbl$dates)

# Calculate the monthly averages across all years using dplyr
monthly_avg <- tbl %>%
  group_by(Month) %>%
  summarise(
    monthly_avg_lake_from_pan_eva = mean(lake_from_pan_eva, na.rm = TRUE),
    monthly_avg_lake_from_theoretical_eva = mean(lake_from_theoretical_eva, na.rm = TRUE),
    monthly_avg_pan_eva = mean(pan_eva, na.rm = TRUE)
  )

# Create Date column for the monthly averages (using an arbitrary year, e.g., 2000)
monthly_avg$YearMonth <- as.Date(paste(2000, monthly_avg$Month, 1, sep = "-"))

# View the resulting data frame
print(monthly_avg)


# Load the new CSV file
data_webb <- read.csv('C:/Users/24468/Desktop/Research/SEAS-HYDRO/Mono Lake/Mono-Evap-Pan/output/eva_coeff_Webb_RStudio.csv', stringsAsFactors = FALSE)

# Convert the 'Date' column to Date type
data_webb$Date <- ymd(data_webb$Date)

# Extract Year-Month as a string (e.g., "2024-July")
data_webb$YearMonth <- paste(year(data_webb$Date), month(data_webb$Date, label = TRUE, abbr = FALSE), sep = "-")

# Replace the year part with "2000" while keeping the month intact
data_webb$YearMonth <- gsub("^\\d{4}", "2000", data_webb$YearMonth)

# Convert to Date (assuming 'YearMonth' is now in "2000-July" format)
data_webb$YearMonth <- as.Date(paste0(data_webb$YearMonth, "-01"), format = "%Y-%B-%d")

# Sort the table by YearMonth
data_webb <- data_webb %>% arrange(YearMonth)

# Calculate the monthly averages for the 'Coeff' column
monthly_avg_webb <- data_webb %>%
  group_by(YearMonth) %>%
  summarise(mean_Coeff = mean(Coeff, na.rm = TRUE))

# Combine the two data frames to get the result
monthly_avg <- left_join(monthly_avg, monthly_avg_webb, by = "YearMonth")

# Add the new result as a column (monthly evaporation based on coefficients)
monthly_avg$monthly_avg_lake_from_webb <- monthly_avg$monthly_avg_pan_eva * monthly_avg_webb$mean_Coeff

# View the resulting data frame
print(monthly_avg)

# Start ggplot without specifying data initially
p <- ggplot()

# Plot the l2s data as lines
p <- p +
  # Evaporation Rate (Solid Black Line)
  geom_line(data = l2s_monthly_avg, 
            aes(x = Date, y = monthly_avg_evaporation_rate, color = "Evaporation Rate"), 
            size = 1) +
  
  # Upper Boundary (Dashed Black Line)
  geom_line(data = l2s_monthly_avg, 
            aes(x = Date, y = monthly_avg_upper_boundary, color = "Upper Boundary"), 
            linetype = "dashed", size = 1) +
  
  # Lower Boundary (Dashed Black Line)
  geom_line(data = l2s_monthly_avg, 
            aes(x = Date, y = monthly_avg_lower_boundary, color = "Lower Boundary"), 
            linetype = "dashed", size = 1) +
  
  # Plot the Penman Method Scatter Points
  # Lake Theoretical Eva, Penman (Red Filled Circles)
  geom_point(data = monthly_avg, 
             aes(x = YearMonth, y = monthly_avg_lake_from_theoretical_eva, color = "Penman, Theoretical"), 
             size = 3, shape = 21, fill = "red") +
  
  # Actual Pan Eva (Green Filled Circles)
  geom_point(data = monthly_avg, 
             aes(x = YearMonth, y = monthly_avg_pan_eva, color = "Actual Pan Eva"), 
             size = 3, shape = 21, fill = "green") +
  
  # Plot the Webb Method Scatter Points
  # Lake Eva, Webb (Blue Filled Circles)
  geom_point(data = monthly_avg, 
             aes(x = YearMonth, y = monthly_avg_lake_from_webb, color = "Webb"), 
             size = 3, shape = 21, fill = "blue") +
  
  # Lake Eva Incorporating Pan Eva, Penman (Yellow Filled Circles)
  geom_point(data = monthly_avg, 
             aes(x = YearMonth, y = monthly_avg_lake_from_pan_eva, color = "Penman, Pan Eva"), 
             size = 3, shape = 21, fill = "yellow") +
  
  # Labels and Title
  labs(
    y = 'Evaporation Rate (mm/d)',
    x = 'Month',
    title = 'Monthly Average Lake Evaporation Rate from Multiple Methods with l2s Model',
    color = 'Legend'
  ) +
  
  # Customize X-axis to show each month with labels
  scale_x_date(
    date_breaks = "1 month", 
    date_labels = "%B",
    limits = c(as.Date("2000-01-01"), as.Date("2000-12-01"))
  ) +
  
  # Define manual colors for the legend
  scale_color_manual(
    values = c(
      "Evaporation Rate" = "black",
      "Upper Boundary" = "black",
      "Lower Boundary" = "black",
      "Penman, Theoretical" = "red",
      "Actual Pan Eva" = "green",
      "Webb" = "blue",
      "Penman, Pan Eva" = "yellow"
    )
  ) +
  
  # Apply a minimal theme and rotate x-axis labels for better readability
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"  # Place the legend at the top
  ) +
  
  # Add grid lines
  theme(panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90"))

# Print the plot
print(p)

export_plot_to_png(p, file_name = "l2s_comparison.png", num_rows = 1, combined = FALSE)