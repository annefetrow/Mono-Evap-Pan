%% Initialization
clearvars; close all; clc;
set(groot, 'defaultLineLineWidth', 1.5);  % Sets the default line width to 1.5

%% Constants
h = 4; % Height of the climate station in meters
salinities = [0, 50, 100]; % Define the salinity levels to investigate

%% Data Import
folder_path = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\2024';
start_date = datetime('2024-07-05');
end_date = datetime('2024-08-22');

% Call the function to read and crop the data
combined_data = read_and_crop_data(folder_path, start_date, end_date);

% Air temperature data
file_path = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\2024_station_data\Air Temperature C.xlsx';
air_temp_data = read_and_average_air_temperature(file_path);

% Align the data tables
[combined_data, air_temp_data] = match_table_sizes(combined_data, air_temp_data);

% Add the Lake Air Temperature from air_temp_data to combined_data
combined_data.Lake_Air_Temperature_C = air_temp_data.Lake_Air_Temperature_C;

%% Monthly Lake Temperature Data
MonthNames = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'}';
AverageTemperatures_C = (([30, 32, 35, 45, 55, 60, 69.6, 71.8, 62.275, 52.75, 49.8, 35] - 32).*5./9)';

monthly_avg_temp_table = table(MonthNames, AverageTemperatures_C);

% Assign average lake temperature based on month
combined_data.Lake_Temperature_C = zeros(height(combined_data), 1);
for i = 1:height(combined_data)
    current_month = month(combined_data.Date(i));
    combined_data.Lake_Temperature_C(i) = monthly_avg_temp_table.AverageTemperatures_C(current_month);
end


% Define line styles for different salinities
line_styles = {'-', '--', '-.'}; % Solid for 0, dashed for 50, dot-dash for 100
colors = lines(length(salinities)); % Get distinct colors for the plots

% Create separate figures for each parameter
figure(1); hold on;
figure(2); hold on;
figure(3); hold on;

for j = 1:length(salinities)
    % Set current salinity
    current_salinity = salinities(j);
    combined_data.Salinity_g_per_kg = current_salinity * ones(height(combined_data), 1);

    %% Calculate Saturation Vapor Pressure
    e_p_kPa = saturated_water_vapor_pressure(combined_data.Daily_Avg_WaterTempC);
    e_p_0_kPa = water_surface_vapor_pressure(combined_data.Daily_Avg_AirTempC, combined_data.Daily_Avg_RH_./100);
    e_p_4_kPa = vapor_pressure_height(e_p_0_kPa, combined_data.Daily_Avg_AirTempC + 273.15, h);
    
    e_l_kPa = salinity_vapor_pressure(saturated_water_vapor_pressure(combined_data.Lake_Temperature_C), combined_data.Salinity_g_per_kg, combined_data.Lake_Temperature_C);
    e_l_0_kPa = water_surface_vapor_pressure(combined_data.Lake_Air_Temperature_C, combined_data.Daily_Avg_RH_./100);
    e_l_4_kPa = vapor_pressure_height(e_l_0_kPa, combined_data.Lake_Air_Temperature_C + 273.15, h);

    %% Calculate Lake Evaporation and Coefficient
    E_lake = ((e_l_kPa - e_l_4_kPa) ./ (e_p_kPa - e_p_4_kPa)) .* combined_data.Daily_Avg_EvaporationRateMm_hr;
    coefficient = ((e_l_kPa - e_l_4_kPa) ./ (e_p_kPa - e_p_4_kPa));

    %% Plot Conversion Coefficient (Figure 1)
    figure(1);
    plot(combined_data.Date, coefficient, 'LineStyle', line_styles{j}, 'Color', colors(j, :), 'DisplayName', sprintf('Salinity = %d g/kg', current_salinity));
    ylabel('Conversion Coefficient');
    xlabel('Date');
    % Only show the legend once
    if j == length(salinities)
        legend show;
    end

    %% Plot Vapor Pressures (Figure 2)
    figure(2);
    plot(combined_data.Date, e_p_4_kPa, 'LineStyle', line_styles{j}, 'Color', colors(j, :), 'DisplayName', sprintf('Salinity = %d g/kg', current_salinity));
    ylabel('Vapor Pressure (kPa)');
    xlabel('Date');
    % Only show the legend once
    if j == length(salinities)
        legend('e_{p,4}', 'e_p', 'e_l', 'e_{l,4}');
    end

    %% Plot Temperatures (Figure 3)
    figure(3);
    plot(combined_data.Date, combined_data.Lake_Temperature_C, 'LineStyle', line_styles{j}, 'Color', colors(j, :), 'DisplayName', sprintf('Salinity = %d g/kg', current_salinity));
    ylabel('Temperature (C)');
    xlabel('Date');
    % Only show the legend once
    if j == length(salinities)
        legend('Lake Water', 'Pan Water', 'Pan Air', 'Lake Air');
    end
end

hold off;


%% Function: Read data and put all into a table
function combined_data = read_and_crop_data(folder_path, start_date, end_date)
    % This function reads all the files in the specified folder,
    % crops them to the specified date range, and combines them into a single table.
    % Inputs:
    % - folder_path: the folder where the data files are located
    % - start_date: the earliest date to include in the output (datetime)
    % - end_date: the latest date to include in the output (datetime)
    % Output:
    % - combined_data: a table with the combined data, where the first column is 'Date'

    % List all .csv files in the folder
    file_list = dir(fullfile(folder_path, '*.csv'));
    
    % Initialize an empty table for the combined data
    combined_data = table();
    
    for i = 1:length(file_list)
        % Read each file
        file_path = fullfile(folder_path, file_list(i).name);
        data = readtable(file_path);

        % Extract the variable name from the file (second column)
        variable_name = data.Properties.VariableNames{2};

        % Convert the date column to datetime
        data.Date = datetime(data.Date, 'Format', 'yyyy-MM-dd');
        
        % Crop the data to the specified date range
        cropped_data = data(data.Date >= start_date & data.Date <= end_date, :);

        % Merge the cropped data into the combined table
        if isempty(combined_data)
            % If it's the first file, initialize the combined table with the date and the variable data
            combined_data = cropped_data;
        else
            % Otherwise, add the variable data to the existing combined table
            combined_data = outerjoin(combined_data, cropped_data, 'Keys', 'Date', ...
                'MergeKeys', true, 'Type', 'left', 'LeftVariables', combined_data.Properties.VariableNames, ...
                'RightVariables', {variable_name});
        end
    end

    % Ensure that the combined data is sorted by date
    combined_data = sortrows(combined_data, 'Date');
end

%% Function: Lake Air Temperature
function daily_air_temp_data = read_and_average_air_temperature(file_path)
    % This function reads the lake air temperature data from the specified Excel file.
    % It calculates the daily average of the lake air temperature, ignoring NaN values.
    % Inputs:
    % - file_path: full path to the "Air Temperature C.xlsx" file
    % Output:
    % - daily_air_temp_data: a table with 'Date' and 'Average_Air_Temperature_C' columns

    % Read the entire file
    data = readtable(file_path);
    
    % Extract the first column as DateTime (yyyy-MM-dd HH format)
    date_time = data{:, 1};  % Assuming the first column contains date and time
    
    % Extract the last column as lake air temperature
    lake_air_temperature = data{:, end};  % Assuming the last column is the lake air temperature in °C

    % Convert the date-time column to datetime format if necessary
    date_time = datetime(date_time, 'InputFormat', 'yyyy-MM-dd HH');
    
    % Convert to just date for grouping
    dates = dateshift(date_time, 'start', 'day');
    
    % Remove NaN values from lake air temperature data
    valid_data_idx = ~isnan(lake_air_temperature);
    dates = dates(valid_data_idx);
    lake_air_temperature = lake_air_temperature(valid_data_idx);
    
    % Calculate daily averages
    [unique_dates, ~, idx] = unique(dates);
    daily_avg_air_temp = accumarray(idx, lake_air_temperature, [], @mean);

    % Create a table with the daily averages
    daily_air_temp_data = table(unique_dates, daily_avg_air_temp, ...
        'VariableNames', {'Date', 'Lake_Air_Temperature_C'});
end


%% Function: crop matriex to be the same size
function [cropped_combined_data, cropped_air_temp_data] = match_table_sizes(combined_data, air_temp_data)
    % This function crops one of the two tables (combined_data or air_temp_data)
    % so that their date ranges match.
    % Inputs:
    % - combined_data: table with the first column as Date
    % - air_temp_data: table with the first column as Date
    % Outputs:
    % - cropped_combined_data: cropped combined_data table
    % - cropped_air_temp_data: cropped air_temp_data table

    % Extract the date columns from both tables
    combined_dates = combined_data.Date;
    air_temp_dates = air_temp_data.Date;

    % Find the common date range (intersection of dates)
    common_start_date = max(min(combined_dates), min(air_temp_dates));
    common_end_date = min(max(combined_dates), max(air_temp_dates));

    % Crop the combined_data to the common date range
    cropped_combined_data = combined_data(combined_dates >= common_start_date & combined_dates <= common_end_date, :);

    % Crop the air_temp_data to the common date range
    cropped_air_temp_data = air_temp_data(air_temp_dates >= common_start_date & air_temp_dates <= common_end_date, :);

    % Display a message if the tables were cropped
    if height(combined_data) ~= height(cropped_combined_data)
        disp('Combined data table has been cropped to match air temperature data.');
    end
    if height(air_temp_data) ~= height(cropped_air_temp_data)
        disp('Air temperature data has been cropped to match combined data.');
    end
end


%% Function
function p_sat = saturated_water_vapor_pressure(T) % p_s in kPa, T in C
    p_sat = 0.611 .* exp(17.27.*T ./ (T + 237.3));

    % Antoine Equation

    % A = 8.07131;
    % B = 1730.63;
    % C = 233.426;
    % 
    % p_s = 10.^((A - B./(T + C))) * 0.133322; % kPa

end

function p_sur = water_surface_vapor_pressure(T,RH) % pw in kPa, T in C

    p_sur = 0.611 .* exp(17.27.*T ./ (T + 237.3)) .* RH;

    % Antoine Equation

    % A = 8.07131;
    % B = 1730.63;
    % C = 233.426;
    % 
    % p_sur = 10.^((A - B./(T + C))) * 0.133322.*RH; % kPa

end

function ph = vapor_pressure_height(p0, T, height)
    M = 18.015*10^(-3);% kg/mol
    gravity = 9.81;
    R = 8.314;
    ph = p0 .* exp(-1 .* (M.*gravity.*height) ./ (R.*T));
end


function p = salinity_vapor_pressure(pw, S, T) %pw in kPa, T in K
    a1 = -2.1609*10^(-4);
    a2 = -3.5015*10^(-7);

    p = pw .* 10.^(a1.*S + a2.*S.^2);
end