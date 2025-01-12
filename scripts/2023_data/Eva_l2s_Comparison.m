%% Initialization
clearvars; close all; clc;

%% Import l2s Data
% Load the CSV file
data = readtable('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\l2s\Joker_L2SWBM_jSum_Formatted.csv');
% Convert the 'Date' column to datetime
data.Date = datetime(data.Date, 'InputFormat', 'yyyy-MM-dd');

% Extract the relevant columns (date, evaporation rate, upper boundary, lower boundary)
l2s_dates = data.Date;
l2s_evaporation_rate = 304.8./30.44.*data{:, 17};
l2s_upper_boundary = 304.8./30.44.*data{:, 19};
l2s_lower_boundary = 304.8./30.44.*data{:, 23};

% Create a table with the necessary columns
tbl = table(l2s_dates, l2s_evaporation_rate, l2s_upper_boundary, l2s_lower_boundary);

% Extract the month from the dates
tbl.Month = month(tbl.l2s_dates);

% Calculate the monthly averages across all years
l2s_monthly_avg = varfun(@mean, tbl, 'InputVariables', {'l2s_evaporation_rate', 'l2s_upper_boundary', 'l2s_lower_boundary'}, ...
                     'GroupingVariables', 'Month');

% Create datetime values for the monthly averages (using an arbitrary year, e.g., 2000)
l2s_monthly_avg.Date = datetime(2000, l2s_monthly_avg.Month, 1);

%% Import Penman Method Results
% Load the new CSV file
data = readtable('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\eva_estimate_Penman_RStudio.csv');

% Convert the 'Date' column to datetime
data.Date = datetime(data.Date, 'InputFormat', 'dd-MMM-yyyy');

% Extract the relevant columns
dates = data.Date;
lake_from_pan_eva = data.Lake_From_Pan_Eva_mm_d;
lake_from_theoretical_eva = data.Lake_From_Theoretical_Eva_mm_d;
pan_eva = data.Pan_Eva_mm_d;

% Create a table with the necessary columns
tbl = table(dates, lake_from_pan_eva, lake_from_theoretical_eva, pan_eva);

% Extract the year and month from the dates
tbl.Year = year(tbl.dates);
tbl.Month = month(tbl.dates);

% Calculate the monthly averages across all years
monthly_avg = varfun(@mean, tbl, 'InputVariables', {'lake_from_pan_eva', 'lake_from_theoretical_eva', 'pan_eva'}, ...
                     'GroupingVariables', 'Month');

% Create datetime values for the monthly averages (using an arbitrary year, e.g., 2000)
monthly_avg.Date = datetime(2000, monthly_avg.Month, 1);

%% Import Webb Methods Results

% Load the new CSV file
data_webb = readtable('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\eva_coeff_Webb_RStudio.csv');

% Convert the date column to datetime format (assuming the first column is the date)
data_webb.Date = datetime(data_webb.Date, 'InputFormat', 'yyyy-MM-dd');  % Adjust format if needed

% Extract year and month from the Date column
data_webb.YearMonth = strcat(string(year(data_webb.Date)), '-', string(month(data_webb.Date, 'long')));

% Convert YearMonth to datetime to ensure correct sorting
data_webb.YearMonth = datetime(data_webb.YearMonth, 'InputFormat', 'yyyy-MMMM');  % Adjust format as necessary

% Sort the table by YearMonth
data_webb = sortrows(data_webb, 'YearMonth');

% Group the data by Year-Month and calculate the average coefficient for each month
monthly_avg_webb = varfun(@mean, data_webb, 'InputVariables', 'Coeff', ...
    'GroupingVariables', 'YearMonth');

% Add the new result as a column in the monthly_avg table
monthly_avg.mean_lake_from_webb = monthly_avg.mean_pan_eva .* monthly_avg_webb.mean_Coeff;

%% Plot all data
% Plot the data from l2s
figure;
plot(l2s_monthly_avg.Date, l2s_monthly_avg.mean_l2s_evaporation_rate, 'k', 'DisplayName', 'Evaporation Rate');
hold on;
plot(l2s_monthly_avg.Date, l2s_monthly_avg.mean_l2s_upper_boundary, 'k--', 'DisplayName', 'Upper Boundary');
plot(l2s_monthly_avg.Date, l2s_monthly_avg.mean_l2s_lower_boundary, 'k--', 'DisplayName', 'Lower Boundary');

% Plot the data as scatter plots from Penman
scatter(monthly_avg.Date, monthly_avg.mean_lake_from_theoretical_eva, 50, 'r', 'o', 'filled', 'DisplayName', 'Lake Theoretical Eva, Penman');
scatter(monthly_avg.Date, monthly_avg.mean_pan_eva, 50, 'g', 'o', 'filled', 'DisplayName', 'Actual Pan Eva');

% Plot the data as scatter plots from Webb
scatter(monthly_avg.Date, monthly_avg.mean_lake_from_webb, 50, 'b', 'o', 'filled', 'DisplayName', 'Lake Eva, Webb');

% Plot pan evaporation
scatter(monthly_avg.Date, monthly_avg.mean_lake_from_pan_eva, 50, 'y', 'o', 'filled', 'DisplayName', 'Lake Eva Incorporating Pan Eva, Penman');

hold off;

% Plot Settings
ylabel('Evaporation Rate (mm/d)');
xlabel('Month');
xticks(datetime(2000,1:12,1)); % Set x-axis to display each month
xticklabels(month(datetime(2000,1:12,1), 'name')); % Display month names
legend('show');
title('Monthly Average Lake Evaporation Rate from Multiple Methods, Plotted along with l2s Model');
grid on; % Add grid lines for better visualization