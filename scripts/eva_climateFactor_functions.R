# Function to create plots for topographic factors
plot_weather_factors <- function(filename, plot_dir) {
  sheets_list <- excel_sheets(filename)
  
  pdf(file = file.path(plot_dir, "weather_factors_plots.pdf"), height = 6, width = 9, paper = "special")
  
  for (sheet in sheets_list) {
    weather_data <- read_excel(filename, sheet = sheet, skip = 2)
    col_names <- names(weather_data)
    weather_data <- na.omit(weather_data)  # Remove rows with NA
    
    if (sheet != "Precipitation,in") {
      matplot(as.numeric(weather_data[[1]]), weather_data[, -1], type = "l", lty = 1, lwd = 1, 
              col = 1:ncol(weather_data[, -1]), ylab = sheet, main = paste("Weather Factor:", sheet))
    } else {
      barplot(-1 * as.matrix(t(weather_data[[2]])), beside = TRUE, names.arg = weather_data[[1]], ylab = sheet, 
              border = coul, col = coul, main = paste("Weather Factor:", sheet))
    }
    
    legend("topright", legend = col_names[-1], col = 1:length(col_names[-1]), lty = 1, cex = 0.8)
  }
  
  dev.off()
}

# Function to obtain precipitation data
precip <- function(data_directory, filename) {
  
  # Import precipitation data
  precip_data <- read_excel(file.path(data_directory, filename), 
                            sheet = 'Precipitation,in', skip = 2)
  
  # Preprocess precipitation data
  P_d <- precip_data[[2]]
  t <- ymd(precip_data[[1]])
  P_d[is.na(P_d)] <- 0
  
  data.frame(date_time = as.POSIXct(t), P_d = -1 * P_d * 25.4 / 24)  # Conversion from inches to millimeters
}