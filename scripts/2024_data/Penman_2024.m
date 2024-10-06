%% Initialization
clearvars; close all; clc;
global RHO_W RHO_L SIGMA P Cp_w Cp_a LAMBDA
set(groot, 'defaultLineLineWidth', 1.5);  % Sets the default line width to 1.5

RHO_W = 1000;
RHO_L = 1060; % kg/m3
SIGMA = 5.67 * 10^(-8); % 1/(s m2 K4)
P = 101325; % Pa
Cp_w = 3.7794*10^(-3); % MJ/kg K
Cp_a = 1.005*10^(-3); % MJ/kg K
LAMBDA = 2.45; % MJ/kg

%% Data Import
folder_path = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\2024';
start_date = datetime('2024-07-05');
end_date = datetime('2024-08-13');

% Call the function to read and crop the data
combined_data = read_and_crop_data(folder_path, start_date, end_date);

% data for air temperature
file_path = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\\data\2024_station_data\Air Temperature C.xlsx';
air_temp_data = read_and_average_air_temperature(file_path);

% Assuming combined_data and air_temp_data have been loaded
[combined_data, air_temp_data] = match_table_sizes(combined_data, air_temp_data);

% Add the Lake Air Temperature from air_temp_data to combined_data
combined_data.Lake_Air_Temperature_C = air_temp_data.Lake_Air_Temperature_C;

% data for radiation
solar_radiation_data = read_and_average_solar_radiation('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\2024_station_data\Radiation.xlsx');

% Assuming combined_data and solar radiation have been loaded
[combined_data, solar_radiation_data] = match_table_sizes(combined_data, solar_radiation_data);

% Add the Lake Air Temperature from air_temp_data to combined_data
combined_data.Solar_Radiation_W_m2 = solar_radiation_data.Average_Solar_Radiation_W_m2;

% data for wind speed
wind_speed_data = read_and_average_wind_speed('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\2024_station_data\Wind Speed m_s.xlsx');

% Assuming combined_data and solar radiation have been loaded
[combined_data, wind_speed_data] = match_table_sizes(combined_data, wind_speed_data);

% Add the Lake Air Temperature from air_temp_data to combined_data
combined_data.Wind_Speed_m_s = wind_speed_data.Average_Wind_Speed_m_s;

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

%% Calculation
% Calculate delta
delta = delta_calc(combined_data.Salinity_g_per_kg, combined_data.Lake_Air_Temperature_C); % kPa/K

%% Specific Heat
gamma = salinity_sigma_calc; % kPa/K

%% Radation
Rn = Rn_calc(combined_data.Lake_Air_Temperature_C, combined_data.Daily_Avg_RH_, combined_data.Solar_Radiation_W_m2); % W/m2

%% Check Ea with theoretical Ea
Ea_theo = 6.43 .* (0.18 + 0.55 .* wind_speed_data.Average_Wind_Speed_m_s) .* (saturated_water_vapor_pressure(combined_data.Daily_Avg_WaterTempC) - water_surface_vapor_pressure(combined_data.Daily_Avg_AirTempC,combined_data.Daily_Avg_RH_./100));

%% Calculate En
Ea_pan = combined_data.Daily_Avg_EvaporationRateMm_hr.*RHO_W.*2.45.*24./1000; % MJ/m2d
E_pan = penman_calc(delta, gamma, Rn.*0.0864, Ea_pan) .* 1000; % mm/d
E_theo = penman_calc(delta, gamma, Rn.*0.0864, Ea_theo) .* 1000; % mm/d

plot(combined_data.Date,E_pan);
ylabel('mm/d');
xlabel('Date');
hold on;
plot(combined_data.Date,E_theo);
plot(combined_data.Date,combined_data.Daily_Avg_EvaporationRateMm_hr.*24);

legend('Lake, from Pan Eva','Lake, from Theoretical Eva','Pan Eva');

%% Calculate Conversion Coefficient
C_pan = E_pan ./ (combined_data.Daily_Avg_EvaporationRateMm_hr.*24);
C_theo = E_theo ./ (combined_data.Daily_Avg_EvaporationRateMm_hr.*24);

figure;
plot(combined_data.Date,C_pan);
ylabel('Conversion Coefficient');
xlabel('Date');
hold on;
plot(combined_data.Date,C_theo);

legend('From Pan Eva','From Theoretical Eva');

%% Save the Results
% Create a table with the calculated data
output_table = table(combined_data.Date, E_pan, E_theo, combined_data.Daily_Avg_EvaporationRateMm_hr.*24, ...
    'VariableNames', {'Date', 'Lake_From_Pan_Eva_mm_d', 'Lake_From_Theoretical_Eva_mm_d', 'Pan_Eva_mm_d'});

% Write the table to a CSV file
writetable(output_table, 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\eva_estimate_Penman.csv');

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

%% Function: Read solar radiation data
function daily_solar_radiation_data = read_and_average_solar_radiation(file_path)
    % This function reads the solar radiation data from the specified Excel file.
    % It calculates the daily average of the solar radiation.
    % Inputs:
    % - file_path: full path to the "Solar Radiation.xlsx" file
    % Output:
    % - daily_solar_radiation_data: a table with 'Date' and 'Average_Solar_Radiation_W_m2' columns

    % Read the entire file
    data = readtable(file_path);
    
    % Extract the first column as DateTime (yyyy-MM-dd HH format)
    date_time = data{:, 1};  % Assuming the first column contains date and time
    
    % Extract the last column as solar radiation
    solar_radiation = data{:, end};  % Assuming the last column is the solar radiation in W/m²

    % Convert the date-time column to datetime format if necessary
    date_time = datetime(date_time, 'InputFormat', 'yyyy-MM-dd HH');
    
    % Convert to just date for grouping
    dates = dateshift(date_time, 'start', 'day');
    
    % Calculate daily averages
    [unique_dates, ~, idx] = unique(dates);
    daily_avg_solar_radiation = accumarray(idx, solar_radiation, [], @mean);

    % Create a table with the daily averages
    daily_solar_radiation_data = table(unique_dates, daily_avg_solar_radiation, ...
        'VariableNames', {'Date', 'Average_Solar_Radiation_W_m2'});
end

%% Function: Read Windspeed Data
function daily_wind_speed_data = read_and_average_wind_speed(file_path)
    % This function reads the wind speed data from the specified Excel file.
    % It calculates the daily average of the wind speed, ignoring NaN values.
    % Inputs:
    % - file_path: full path to the "Wind Speed.xlsx" file
    % Output:
    % - daily_wind_speed_data: a table with 'Date' and 'Average_Wind_Speed_m_s' columns

    % Read the entire file
    data = readtable(file_path);
    
    % Extract the first column as DateTime (yyyy-MM-dd HH format)
    date_time = data{:, 1};  % Assuming the first column contains date and time
    
    % Extract the last column as wind speed
    wind_speed = data{:, end};  % Assuming the last column is the wind speed in m/s

    % Convert the date-time column to datetime format if necessary
    date_time = datetime(date_time, 'InputFormat', 'yyyy-MM-dd HH');
    
    % Convert to just date for grouping
    dates = dateshift(date_time, 'start', 'day');
    
    % Remove NaN values from wind speed data
    valid_data_idx = ~isnan(wind_speed);
    dates = dates(valid_data_idx);
    wind_speed = wind_speed(valid_data_idx);
    
    % Calculate daily averages
    [unique_dates, ~, idx] = unique(dates);
    daily_avg_wind_speed = accumarray(idx, wind_speed, [], @mean);

    % Create a table with the daily averages
    daily_wind_speed_data = table(unique_dates, daily_avg_wind_speed, ...
        'VariableNames', {'Date', 'Average_Wind_Speed_m_s'});
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
function E = penman_calc(delta, gamma, Rn, Ea) % delta in kPa/C, gamma in kPa/C, Ea in MJ/m^2, Rn in MJ/m2d
    
    % Global Variables
    global RHO_L LAMBDA

    % Define Constants
    G = 0; % MJ/m2d

    % Calculate E
    E = 1 ./ LAMBDA .* (delta .* (Rn - G) + gamma .* Ea) ./ ((delta + gamma) .* RHO_L);

end

function Rn = Rn_calc(Ta, RH, R1)
    
    % Global Variables
    global SIGMA

    % Define Constants
    r = 0.05;
    a = 0.4; % assumption, vary with latitude
    b = 0.274; % assumption, vary with latitude
    n_N = 0.8; % assumption
    ed_0_kPa = water_surface_vapor_pressure(Ta,RH./100);
    ed_4_kPa = vapor_pressure_height(ed_0_kPa, Ta+273.15, 4); % kPa, assume measured at 4m height

    % % Calculate R1
    % R1 = RA .* (a + b .* n_N);

    % Calculate R2
    RB = SIGMA .* (Ta+273.15).^4 .*(0.56 - 0.09 .* sqrt(ed_4_kPa)) .* (0.1 + 0.9 .* n_N);

    % Calculate Rn
    Rn = R1 .* (1-r) - RB;

end

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

function delta = delta_calc(S, T) % T in C
    
    dT = 0.1;
    delta = (salinity_vapor_pressure(T+dT, S) - salinity_vapor_pressure(T - dT, S)) ./ (2 .* dT);

end

function gamma = salinity_sigma_calc()
    
    % Global Variables
    global P Cp_a LAMBDA

    % Define Constants
    mu = 0.622;

    % Gamma Calculation
    gamma = (Cp_a .* P./1000) ./ (mu .* LAMBDA);

end

function Cp = specific_heat()

    global S
        
    a1 =  5.328;
    a2 =  -9.76 * 10^(-2);
    a3 =  4.04 * 10^(-4);
    
    b1 = -6.913 * 10^(-3);
    b2 = 7.351 * 10^(-4);
    b3 = -3.15 * 10^(-6);
    
    c1 = 9.6 * 10^(-6);
    c2 = -1.927 * 10^(-6);
    c3 = 8.23 * 10^(-9);
    
    d1 = 2.5 * 10^(-9); 
    d2 = 1.666 * 10^(-9);
    d3 = -7.125 * 10^(-12);
    
    T = 20 + 273; %K
    Cp = (a1 + a2*S + a3*S^2) + (b1 + b2*S + b3*S^2)*T + (c1 + c2*S + c3*S^2)*T^2 +(d1 + d2*S + d3*S^2)*T^3;

end

function resultsTable = read_weather_factor(filePath)
    
    % Define the date range
    startDate = datetime('16-Aug-2023', 'InputFormat', 'dd-MMM-yyyy');
    endDate = datetime('26-Oct-2023', 'InputFormat', 'dd-MMM-yyyy');
    
    % Get information about the Excel file
    [~, sheetNames] = xlsfinfo(filePath);
    
    % Initialize a table to store the results
    resultsTable = table;
    
    % Loop through each sheet
    for i = 1:length(sheetNames)
        % Read the data from the current sheet
        sheetName = sheetNames{i};
        dataTable = readtable(filePath, 'Sheet', sheetName);
        
        % Assuming the first column is the date and the rest are data columns
        dates = datetime(dataTable{:,1}, 'InputFormat', 'dd-MMM-yyyy');
        
        % Filter the data based on the date range
        dateFilter = dates >= startDate & dates <= endDate;
        filteredData = dataTable(dateFilter, :);
        
        % Calculate the average of the data columns (excluding NaN values)
        avgValues = mean(filteredData{:, 2:end}, 2, 'omitnan');
        
        % Store the dates and average values in the results table
        if isempty(resultsTable)
            % Add dates as the first column in the table
            resultsTable.Date = dates(dateFilter);
        end
        resultsTable.(sheetName) = avgValues(:);

        % Display the results table
        disp(resultsTable);
    end
end

function resultsTable = read_eva_rate(filePath)

    % Define the date range
    startDate = datetime('16-Aug-2023', 'InputFormat', 'dd-MMM-yyyy');
    endDate = datetime('26-Oct-2023', 'InputFormat', 'dd-MMM-yyyy');
    
    % Get information about the Excel file
    [~, sheetNames] = xlsfinfo('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\daily_average_evaporation_rates.csv');
    
    
    % Initialize a table to store the results
    resultsTable = table;
    
        dataTable = readtable('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\daily_average_evaporation_rates.csv');
        
        % Assuming the first column is the date and the rest are data columns
        dates = datetime(dataTable{:,1}, 'InputFormat', 'dd-MMM-yyyy');
        
        % Filter the data based on the date range
        dateFilter = dates >= startDate & dates <= endDate;
        filteredData = dataTable(dateFilter, :);
        
        % Store the dates and average values in the results table
        resultsTable.Date = dates(dateFilter);
        resultsTable.(sheetNames{1}) = mean(filteredData{:, 2:end}, 2, 'omitnan');

        % Display the results table
        disp(resultsTable);
end

function T = Ferenheit_to_Kelvin(T)
    T = (T - 32) .* 5./9 + 273.15;
end
