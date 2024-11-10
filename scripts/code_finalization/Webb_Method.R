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
    data$Date <- as.Date(data$Date, format = "%Y-%m-%d")
    
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

# Usage:
# combined_data <- read_and_crop_data("path/to/your/folder", as.Date("2023-01-01"), as.Date("2023-12-31"))
