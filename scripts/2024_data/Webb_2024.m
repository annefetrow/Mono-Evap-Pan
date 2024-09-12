%% Initialization
clearvars; close all; clc;
set(groot, 'defaultLineLineWidth', 1.5);  % Sets the default line width to 1.5

%% Constant
S = 81;
height = 4;

%% Data Import
folder_path = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\2024';
start_date = datetime('2024-07-05');
end_date = datetime('2024-08-22');

% Call the function to read and crop the data
combined_data = read_and_crop_data(folder_path, start_date, end_date);

% data for air temperature
file_path = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\\data\2024_station_data\Air Temperature C.xlsx';
air_temp_data = read_air_temperature(file_path);

%% Align Matrix
% Assuming combined_data and air_temp_data have been loaded
[combined_data, air_temp_data] = match_table_sizes(combined_data, air_temp_data);

%% Monthly Lake Temperature Data
% Create a table with month names and average temperatures, data from
% Anne's field measurement and online database
MonthNames = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'}';
AverageTemperatures_C = (([30, 32, 35, 45, 55, 60, 69.6, 71.8, 62.275, 52.75, 49.8, 35] - 32).*5./9)';

% Create a table:
monthly_avg_temp_table = table(MonthNames, AverageTemperatures_C);
monthly_avg_temp_table.SalinityPressure_kPa = salinity_vapor_pressure(saturated_pressure(AverageTemperatures_C), S, (monthly_avg_temp_table.AverageTemperatures_C + 273.15));
monthly_avg_temp_table.SurfacePressure_kPa = vapor_pressure(AverageTemperatures_C+273.15);

%% Calculate Saturation Vapor Pressure
e_p_kPa = saturated_pressure(combined_data.Daily_Avg_WaterTempC);
e_p_0_kPa = vapor_pressure(combined_data.Daily_Avg_WaterTempC+273.15);
e_p_4_kPa = vapor_pressure_height(e_p_0_kPa, combined_data.Daily_Avg_AirTempC+273.15, height);

e_l_0_kPa = vapor_pressure(air_temp_data.Lake_Air_Temperature_C+273.15);
e_l_kPa = zeros(length(air_temp_data.Lake_Air_Temperature_C),1);
T_lake = zeros(length(air_temp_data.Lake_Air_Temperature_C),1);

for i = 1:length(e_l_kPa)
    monthIndex = month(combined_data.Date(i));
    e_l_kPa(i) = monthly_avg_temp_table.SalinityPressure_kPa(monthIndex);
    T_lake(i) = AverageTemperatures_C(monthIndex);
end

e_l_4_kPa = vapor_pressure_height(e_l_0_kPa, air_temp_data.Lake_Air_Temperature_C+273.15, height);

%% Calculate Lake Evaporation
E_lake = 1.5 .* ((e_l_kPa - e_l_4_kPa) ./ (e_p_kPa - e_p_4_kPa)) .* combined_data.Daily_Avg_EvaporationRateMm_hr;
coefficient = 1.5 .* ((e_l_kPa - e_l_4_kPa) ./ (e_p_kPa - e_p_4_kPa));

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
plot(combined_data.Date,T_lake);
hold on;
plot(combined_data.Date,combined_data.Daily_Avg_WaterTempC);
plot(combined_data.Date,combined_data.Daily_Avg_AirTempC);
plot(combined_data.Date,air_temp_data.Lake_Air_Temperature_C);
hold off;
xlabel('Date');
ylabel('Temperature');
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
function air_temp_data = read_air_temperature(file_path)
    % This function reads the air temperature data from the specified Excel file.
    % It extracts the first column as the date and the last column as lake air temperature.
    % Inputs:
    % - file_path: full path to the "Air Temperature C.xlsx" file
    % Output:
    % - air_temp_data: a table with 'Date' and 'Lake_Air_Temperature_C' columns

    % Read the entire file
    data = readtable(file_path);
    
    % Extract the first column as Date
    dates = data{:, 1};  % Assuming the first column is the date
    
    % Extract the last column as lake air temperature
    lake_air_temp = (data{:, end} - 32) * 5/9;  % Assuming the last column is the lake air temperature

    % Convert the date column to datetime format if necessary
    dates = datetime(dates, 'Format', 'MM-dd-yyyy');
    
    % Create a table with the extracted data
    air_temp_data = table(dates, lake_air_temp, ...
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
function p_s = saturated_pressure(T) % p_s in kPa, T in C
    p_s = 0.611 .* exp(17.27.*T ./ (T + 237.3));
end

function pw = vapor_pressure(T) % pw in kPa, T in K

    A = 5.04221;
    B = 1838.675;
    C = -31.737;
    
    pw = 10.^((A - B./(C+T))); % bar

    % Convert unit from bar to kPa
    pw = 10^2 .* pw;
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