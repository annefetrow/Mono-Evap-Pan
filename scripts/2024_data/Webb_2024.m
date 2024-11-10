%% Initialization
clearvars; close all; clc;
set(groot, 'defaultLineLineWidth', 1.5);  % Sets the default line width to 1.5

%% Constant
h = 4; % assume the climate station is measured at 4 m height

%% Data Import
folder_path = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\2024';
start_date = datetime('2024-07-05');
end_date = datetime('2024-08-22');

% Call the function to read and crop the data
combined_data = read_and_crop_data(folder_path, start_date, end_date);

% data for air temperature
file_path = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\\data\2024_station_data\Air Temperature C.xlsx';
air_temp_data = read_and_average_air_temperature(file_path);

% Assuming combined_data and air_temp_data have been loaded
[combined_data, air_temp_data] = match_table_sizes(combined_data, air_temp_data);

% Add the Lake Air Temperature from air_temp_data to combined_data
combined_data.Lake_Air_Temperature_C = air_temp_data.Lake_Air_Temperature_C;

%% Dynamic Salinity Input
% Create a column for salinity in your table based on your data or rules.
combined_data.Salinity_g_per_kg = zeros(height(combined_data), 1); % Assuming salinity is initially zero

% Define a rule to set different salinity values based on date (e.g., arbitrary example)
for i = 1:height(combined_data)
    current_date = combined_data.Date(i);
    if current_date < datetime('2024-07-20')
        combined_data.Salinity_g_per_kg(i) = 75; % Salinity for earlier dates
    elseif current_date < datetime('2024-08-15')
        combined_data.Salinity_g_per_kg(i) = 81; % Salinity for middle dates
    else
        combined_data.Salinity_g_per_kg(i) = 85; % Salinity for later dates
    end
end

%% Monthly Lake Temperature Data
% Create a table with month names and average temperatures from Anne's field measurement and online database
MonthNames = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'}';
AverageTemperatures_C = (([30, 32, 35, 45, 55, 60, 69.6, 71.8, 62.275, 52.75, 49.8, 35] - 32).*5./9)';

% Create a table:
monthly_avg_temp_table = table(MonthNames, AverageTemperatures_C);

% Preallocate the lake temperature column in combined_data
combined_data.Lake_Temperature_C = zeros(height(combined_data), 1);

% Loop through combined_data and assign lake temperature based on the month
for i = 1:height(combined_data)
    % Get the month index for the current date in combined_data
    current_month = month(combined_data.Date(i));
    
    % Assign the corresponding average lake temperature based on the month
    combined_data.Lake_Temperature_C(i) = monthly_avg_temp_table.AverageTemperatures_C(current_month);
end

% Display the first few rows of the updated combined_data
disp(combined_data(1:10, {'Date', 'Lake_Temperature_C'}));

%% Calculate Saturation Vapor Pressure
e_p_kPa = saturated_water_vapor_pressure(combined_data.WaterTemp_C);
e_p_0_kPa = water_surface_vapor_pressure(combined_data.AirTemp_C,combined_data.RH_Percent./100);
e_p_4_kPa = vapor_pressure_height(e_p_0_kPa, combined_data.AirTemp_C+273.15, h);

e_l_kPa = salinity_vapor_pressure(combined_data.Lake_Temperature_C, combined_data.Salinity_g_per_kg);
e_l_0_kPa = water_surface_vapor_pressure(combined_data.Lake_Air_Temperature_C,combined_data.RH_Percent./100);
e_l_4_kPa = vapor_pressure_height(e_l_0_kPa, combined_data.Lake_Air_Temperature_C+273.15, h);


%% Calculate Lake Evaporation
E_lake = ((e_l_kPa - e_l_4_kPa) ./ (e_p_kPa - e_p_4_kPa)) .* combined_data.Evaporation_Rate_mm_hr;
coefficient = ((e_l_kPa - e_l_4_kPa) ./ (e_p_kPa - e_p_4_kPa));

%% Plot
plot(combined_data.Date,coefficient);
ylabel('Conversion Coefficient');
xlabel('Date');

figure;
plot(combined_data.Date,e_p_4_kPa);
hold on;
plot(combined_data.Date,e_p_kPa);
plot(combined_data.Date,e_l_kPa);
plot(combined_data.Date,e_l_4_kPa);
hold off;
legend('e_{p,4}','e_p','e_l','e_{l,4}');
ylabel('Vapor Pressure (kPa)');
xlabel('Date');

figure;
plot(combined_data.Date,combined_data.Lake_Temperature_C);
hold on;
plot(combined_data.Date,combined_data.WaterTemp_C);
plot(combined_data.Date,combined_data.AirTemp_C);
plot(combined_data.Date,combined_data.Lake_Air_Temperature_C);
hold off;
xlabel('Date');
ylabel('Temperature (C)');
legend('Lake Water','Pan Water','Pan Air','Lake Air');

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


%% Function: Crop tables to match date ranges
function [cropped_table1, cropped_table2] = match_table_sizes(table1, table2)
    % This function crops two tables so that their date ranges match.
    % Inputs:
    % - table1: first table with the first column as Date
    % - table2: second table with the first column as Date
    % Outputs:
    % - cropped_table1: cropped table1 to match the date range of table2
    % - cropped_table2: cropped table2 to match the date range of table1

    % Extract the date columns from both tables
    dates_table1 = table1.Date;
    dates_table2 = table2.Date;

    % Find the common date range (intersection of dates)
    common_start_date = max(min(dates_table1), min(dates_table2));
    common_end_date = min(max(dates_table1), max(dates_table2));

    % Crop table1 to the common date range
    cropped_table1 = table1(dates_table1 >= common_start_date & dates_table1 <= common_end_date, :);

    % Crop table2 to the common date range
    cropped_table2 = table2(dates_table2 >= common_start_date & dates_table2 <= common_end_date, :);

    % Display a message if either table was cropped
    if height(table1) ~= height(cropped_table1)
        disp('Table 1 has been cropped to match Table 2.');
    end
    if height(table2) ~= height(cropped_table2)
        disp('Table 2 has been cropped to match Table 1.');
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

function ph = vapor_pressure_height(p0, T, height) % T in C, height in m
    M = 18.015*10^(-3);% kg/mol
    gravity = 9.81;
    R = 8.314;
    ph = p0 .* exp(-1 .* (M.*gravity.*height) ./ (R.*T));
end


function p = salinity_vapor_pressure(T, S) % p in kPa, T in C
    pw = saturated_water_vapor_pressure(T);

    a1 = -2.1609*10^(-4);
    a2 = -3.5015*10^(-7);

    p = pw .* 10.^(a1.*S + a2.*S.^2);
end