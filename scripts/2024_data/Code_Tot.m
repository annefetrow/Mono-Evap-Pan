%% Initialization
clearvars; close all; clc;
syms x

set(gca, 'LineWidth', 1.5); % Set line width for axes
set(0, 'DefaultLineLineWidth', 1.5); % Set default line width for all lines

%% Water Level
% Path to the folder containing the CSV files
folderPath = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\Raw data – Evaporation Pan\Pan Water Level';

% Read and combine data from all CSV files
[all_h, all_date_time] = read_and_combine_csv(folderPath,5,7); % most recent 3 files

plot(all_date_time, all_h, 'b-');
xlabel('Date and Time');
ylabel('Water Level (mm)');
title('Water Level over Time');
grid on;

%% Evaporation Rate
[eva_rate,eva_date_time] = eva_rate_cal(all_h,all_date_time);

eva_rate_avg = mean(eva_rate(eva_rate<3)); %mm/hr

figure;
bar(eva_date_time,eva_rate,'b',EdgeColor='b');
hold on;
plot(eva_date_time,ones(1,length(eva_date_time))*eva_rate_avg,'--b');

% Add the average evaporation rate text in the top right corner
x_position = eva_date_time(end) - 1; % Right-most x-axis position
y_position = max(eva_rate) + 0.2; % Slightly above the highest evaporation rate
text(x_position, y_position, ['Average: ', num2str(eva_rate_avg, '%.2f'), ' mm/hr'], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', 'Color', 'b', 'FontSize', 10);

xlabel('Date Time');
ylabel('mm/hr');
title('Eavporation Rate Over Time');
hold off;

% Calculate daily evaporation rate and export to a file
daily_eva_rate_avg = calculate_and_export_daily_avg(eva_rate, eva_date_time, 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\2024', '2024_daily_average_evaporation_rates.csv','Evaporation Rate mm/hr');

%% Water Temperature

% Path to the folder containing the CSV files
folderPath = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\Raw data – Evaporation Pan\Water Temp';

% Read and combine data from all CSV files
[T_all, T_all_date_time] = read_and_combine_csv(folderPath,1,1); % most recent 1 file

% Plot temperature with time
figure;
plot(T_all_date_time,T_all, 'r-');
xlabel('Date and Time');
ylabel('Temperature (C)');
title('Water Temperature over Time');
grid on;
axis tight;

% Calculate daily evaporation rate and export to a file
daily_water_temp_avg = calculate_and_export_daily_avg(T_all, T_all_date_time, 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\2024', '2024_daily_average_water_temperature.csv','Water Temp C');

%% Air Temperature

% Path to the folder containing the CSV files
folderPath = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\Raw data – Evaporation Pan\Air Temp';

% Read and combine data from all CSV files
[T_air_all, T_air_all_date_time] = read_and_combine_csv(folderPath,1,1); % most recent 1 file

% Plot temperature with time
figure;
plot(T_air_all_date_time,T_air_all, 'g-');
xlabel('Date and Time');
ylabel('Temperature (C)');
title('Air Temperature over Time');
grid on;
axis tight;

% Calculate daily evaporation rate and export to a file
daily_air_temp_avg = calculate_and_export_daily_avg(T_air_all, T_air_all_date_time, 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\2024', '2024_daily_average_air_temperature.csv','Air Temp C');


%% Relative Humidity

% Path to the folder containing the CSV files
folderPath = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\Raw data – Evaporation Pan\RH';

% Read and combine data from all CSV files
[RH_all, RH_all_date_time] = read_and_combine_csv(folderPath,1,1); % most recent 1 file

% Plot temperature with time
figure;
plot(RH_all_date_time,RH_all, 'k-');
xlabel('Date and Time');
ylabel('Relative Humidity (%)');
title('Relative Humidity over Time');
grid on;
axis tight;

% Calculate daily evaporation rate and export to a file
daily_RH_avg = calculate_and_export_daily_avg(RH_all, RH_all_date_time, 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\2024', '2024_daily_average_RH.csv','RH %');


%% Plot Water Temperature, Air Temperature, and Relative Humidity
% Create figure
figure;

% Plot water temperature and air temperature on the left y-axis
yyaxis left;
plot(T_all_date_time, T_all, '-r', 'DisplayName', 'Water Temperature');
hold on;
plot(T_air_all_date_time, T_air_all, '-g', 'DisplayName', 'Air Temperature');
ylabel('Temperature (°C)');
ylim([min([T_all; T_air_all]) - 1, max([T_all; T_air_all]) + 1]);
legend('show');

% Plot relative humidity on the right y-axis
yyaxis right;
plot(RH_all_date_time, RH_all, '--k', 'DisplayName', 'Relative Humidity');
ylabel('Relative Humidity (%)');
ylim([min(RH_all) - 10, max(RH_all) + 10]);

% Add labels and title
xlabel('Date and Time');
title('Water Temperature, Air Temperature, and Relative Humidity Over Time');
grid on;

% Ensure the legend is visible
legend('show');
axis tight;

%% Combine Water Temperature, Evaporation Rate, and Air Temperature
% Trim the data to have the same date range
[common_date_range, trimmed_eva_rate, trimmed_temp, trimmed_eva_time, trimmed_temp_time] = trim_to_common_date_boundaries(eva_rate, eva_date_time, T_all, T_all_date_time);

% Filter air temperature data to the common date range
filter_idx_air = T_air_all_date_time >= common_date_range(1) & T_air_all_date_time <= common_date_range(2);
trimmed_T_air_time = T_air_all_date_time(filter_idx_air);
trimmed_T_air = T_air_all(filter_idx_air);

% Create figure
figure;

% Plot the evaporation rate
yyaxis left;
bar(trimmed_eva_time, trimmed_eva_rate, 'b', 'EdgeColor', 'b');
ylabel('Evaporation Rate (mm/hr)');
ylim([0, max(trimmed_eva_rate) + 1]); % Adjust limits for better visualization

% Plot the water temperature
yyaxis right;
plot(trimmed_temp_time, trimmed_temp, '--r');
ylabel('Water Temperature (°C)');
ylim([min(trimmed_temp) - 1, max(trimmed_temp) + 1]); % Adjust limits for better visualization

% Plot the air temperature
hold on;
yyaxis right; % Set to the left y-axis
plot(trimmed_T_air_time, trimmed_T_air, '--g'); % Use a different line style for clarity
ylabel('Water Temperature / Air Temperature (C)');
legend('Evaporation Rate', 'Water Temperature', 'Air Temperature');

% Add labels and title
xlabel('Date');
title('Evaporation Rate, Water Temperature, and Air Temperature Over Time');
grid on;

%% Water Temperature and Evaporation Rate Correspondance
plot_evaporation_data(trimmed_eva_time, trimmed_eva_rate, T_air_all_date_time, T_air_all, T_all_date_time, T_all, RH_all_date_time, RH_all, datetime('2024-08-15'));

%% Webb's Method

%% Function: Read and Combine Data from CSV Files
function [unique_data, unique_timestamp] = read_and_combine_csv(folderPath, startfile, endfile)
    % read_and_combine_csv reads CSV files from a specified folder
    % and combines the data and timestamps while removing duplicates and NaNs.
    %
    % Inputs:
    %   folderPath - Path to the folder containing the CSV files
    %   startfile  - Index of the first file to read (1-based index)
    %   endfile    - Index of the last file to read
    %
    % Outputs:
    %   unique_data      - Combined data with duplicates and NaNs removed
    %   unique_timestamp - Combined timestamps with duplicates and NaNs removed

    % Get a list of all CSV files in the folder
    filePattern = fullfile(folderPath, '*.csv');
    csvFiles = dir(filePattern);

    % Initialize arrays to hold the combined data
    all_data = [];
    all_timestamp = [];

    % Loop through each file and read the data
    for k = startfile:endfile
        baseFileName = csvFiles(k).name;
        fullFileName = fullfile(folderPath, baseFileName);

        % Import the data from the CSV file
        [data, timestamp] = data_import(fullFileName);

        % Append the data to the combined arrays
        all_data = [all_data; data];
        all_timestamp = [all_timestamp; timestamp];
    end

    % Remove duplicates and NaN values from the combined data and timestamps
    [unique_data, unique_timestamp] = remove_duplicates_and_nans(all_data, all_timestamp);

    % Optionally, sort the data by timestamp if needed
    [unique_timestamp, sortIdx] = sort(unique_timestamp);
    unique_data = unique_data(sortIdx);
end

%% Function: Data Import with Corrected Year
function [data, timestamp] = data_import(fileName)
    % data_import reads the data and timestamps from a CSV file.
    % It ensures the correct interpretation of timestamps, especially the year.
    %
    % Inputs:
    %   fileName - Path to the CSV file
    %
    % Outputs:
    %   data - Water level or other data (third column)
    %   timestamp - Timestamps corresponding to the data (second column)

    % Import the table from the CSV file
    opts = detectImportOptions(fileName);
    
    % Read the table from the CSV file
    dataTable = readtable(fileName, opts);

    % Extract the second column for the date (timestamp) and the third column for the data
    dateStrings = dataTable{:, 2};  % Assuming the date is in the second column
    data = dataTable{:, 3};         % Assuming the water level data is in the third column

    % Convert the date strings to datetime
    timestamp = datetime(dateStrings, 'InputFormat', 'MM/dd/yyyy HH:mm');  % Adjust 'MM/dd/yyyy HH:mm' as needed
    
    % Fix the year if it's being read as '0024' instead of '2024'
    % This example assumes all dates should be in the 2000s
    wrongYearIdx = year(timestamp) < 100;  % Find rows where the year is less than 100 (e.g., 0024)
    % Create a new datetime array with the corrected year, preserving the hour and minute
    corrected_years = year(timestamp(wrongYearIdx)) + 2000;
    timestamp(wrongYearIdx) = datetime(corrected_years, month(timestamp(wrongYearIdx)), day(timestamp(wrongYearIdx)), ...
                                        hour(timestamp(wrongYearIdx)), minute(timestamp(wrongYearIdx)), second(timestamp(wrongYearIdx)));
end

%% Function: Remove Duplicates and NaNs
function [unique_data, unique_timestamp] = remove_duplicates_and_nans(all_data, all_timestamp)
    % Removes duplicate timestamp entries and NaN values from the data.
    %
    % Inputs:
    %   all_data - Combined data array
    %   all_timestamp - Combined timestamp array
    %
    % Outputs:
    %   unique_data - Data with duplicates and NaNs removed
    %   unique_timestamp - Timestamps with duplicates and NaNs removed

    % Remove NaN values
    validIdx = ~isnan(all_data);
    all_data = all_data(validIdx);
    all_timestamp = all_timestamp(validIdx);

    % Remove duplicates based on timestamps, keeping the first occurrence
    [unique_timestamp, uniqueIdx] = unique(all_timestamp, 'stable');
    unique_data = all_data(uniqueIdx);
end


%% Function: Trim Data to Common Date Boundaries
function [common_date_range, data1_trimmed, data2_trimmed, time1_trimmed, time2_trimmed] = trim_to_common_date_boundaries(data1, time1, data2, time2)
    % trim_to_common_date_boundaries trims two data sets to have the same start
    % and end date, without aligning the data points in between. It ensures both 
    % datasets share the same date range.
    %
    % Inputs:
    %   data1 - First dataset (e.g., evaporation rate)
    %   time1 - Corresponding datetime array for data1
    %   data2 - Second dataset (e.g., water temperature)
    %   time2 - Corresponding datetime array for data2
    %
    % Outputs:
    %   common_date_range - The common start and end dates
    %   data1_trimmed - Trimmed data1 within the common date range
    %   data2_trimmed - Trimmed data2 within the common date range
    %   time1_trimmed - Trimmed time1 within the common date range
    %   time2_trimmed - Trimmed time2 within the common date range

    % Find the maximum of the start dates and the minimum of the end dates
    common_start = max([min(time1), min(time2)]);
    common_end = min([max(time1), max(time2)]);

    % Create the common date range
    common_date_range = [common_start, common_end];

    % Trim data1 and data2 to this common date range
    time1_trimmed = time1(time1 >= common_start & time1 <= common_end);
    data1_trimmed = data1(time1 >= common_start & time1 <= common_end);

    time2_trimmed = time2(time2 >= common_start & time2 <= common_end);
    data2_trimmed = data2(time2 >= common_start & time2 <= common_end);
end

%% Function: Nearest Neighbor
function [nearest_values, nearest_indices] = nearest_neighbor(target_times, reference_times, reference_values)
    % nearest_neighbor finds the nearest reference values for each target time
    % based on the nearest neighbor principle.
    %
    % Inputs:
    %   target_times - Array of target datetime values
    %   reference_times - Array of reference datetime values
    %   reference_values - Array of values corresponding to reference_times
    %
    % Outputs:
    %   nearest_values - Nearest reference values for each target time
    %   nearest_indices - Indices of the nearest reference times

    % Preallocate arrays for results
    nearest_values = NaN(size(target_times));
    nearest_indices = NaN(size(target_times));
    
    % Find the nearest reference time for each target time
    for i = 1:length(target_times)
        % Calculate the time differences
        time_diffs = abs(reference_times - target_times(i));
        
        % Find the index of the minimum time difference
        [~, idx] = min(time_diffs);
        
        % Store the nearest value and index
        nearest_values(i) = reference_values(idx);
        nearest_indices(i) = idx;
    end
end

%% Function: Correspondance Between Evaporation Rate and Other Climate Factors
function plot_evaporation_data(eva_time, eva_rate, T_air_time, T_air, T_water_time, T_water, RH_time, RH, cutoff_date)
    % plot_evaporation_data creates a figure with three subplots:
    % 1. Evaporation rate vs. air temperature
    % 2. Evaporation rate vs. water temperature
    % 3. Evaporation rate vs. relative humidity

    % Filter data based on cutoff_date
    filter_idx_eva = eva_time <= cutoff_date;
    filtered_eva_time = eva_time(filter_idx_eva);
    filtered_eva_rate = eva_rate(filter_idx_eva);

    % Air Temperature
    filter_idx_air = T_air_time <= cutoff_date;
    filtered_T_air_time = T_air_time(filter_idx_air);
    filtered_T_air = T_air(filter_idx_air);

    % Water Temperature
    filter_idx_water = T_water_time <= cutoff_date;
    filtered_T_water_time = T_water_time(filter_idx_water);
    filtered_T_water = T_water(filter_idx_water);

    % Relative Humidity
    filter_idx_RH = RH_time <= cutoff_date;
    filtered_RH_time = RH_time(filter_idx_RH);
    filtered_RH = RH(filter_idx_RH);

    % Find nearest evaporation rate for each temperature measurement
    nearest_eva_rate_air = nearest_neighbor(filtered_T_air_time, filtered_eva_time, filtered_eva_rate);
    nearest_eva_rate_water = nearest_neighbor(filtered_T_water_time, filtered_eva_time, filtered_eva_rate);
    nearest_eva_rate_RH = nearest_neighbor(filtered_RH_time, filtered_eva_time, filtered_eva_rate);

    % Determine common x-axis limits
    x_limits_air = [min(filtered_T_air), max(filtered_T_air)];
    x_limits_water = [min(filtered_T_water), max(filtered_T_water)];
    x_limits_RH = [min(filtered_RH), max(filtered_RH)];
    
    % Find the overall x-axis limits
    common_x_limits = [min([x_limits_air(1), x_limits_water(1), x_limits_RH(1)]), ...
                       max([x_limits_air(2), x_limits_water(2), x_limits_RH(2)])];

    % Create figure and subplots
    figure;

    % Subplot 1: Evaporation Rate vs Air Temperature
    subplot(3, 1, 1);
    plot(filtered_T_air, nearest_eva_rate_air, '.');
    xlabel('Air Temperature (°C)');
    ylabel('Evaporation Rate (mm/hr)');
    title('Evaporation Rate vs. Air Temperature');
    xlim(common_x_limits); % Align x-axis
    grid on;

    % Subplot 2: Evaporation Rate vs Water Temperature
    subplot(3, 1, 2);
    plot(filtered_T_water, nearest_eva_rate_water, '.');
    xlabel('Water Temperature (°C)');
    ylabel('Evaporation Rate (mm/hr)');
    title('Evaporation Rate vs. Water Temperature');
    xlim(common_x_limits); % Align x-axis
    grid on;

    % Subplot 3: Evaporation Rate vs Relative Humidity
    subplot(3, 1, 3);
    plot(filtered_RH, nearest_eva_rate_RH, '.');
    xlabel('Relative Humidity (%)');
    ylabel('Evaporation Rate (mm/hr)');
    title('Evaporation Rate vs. Relative Humidity');
    xlim(common_x_limits); % Align x-axis
    grid on;
end

%% Function: Calculate Daily Average and Export
function daily_avg = calculate_and_export_daily_avg(data, date_time, folder_path, file_name, variable_name)
    % This function calculates the daily average of any input data and exports it to a file.
    % It returns the daily average values.
    % Inputs:
    % - data: the data array (e.g., evaporation rate, temperature, etc.)
    % - date_time: corresponding datetime array
    % - folder_path: folder path to save the output file
    % - file_name: name of the output file (e.g., 'Daily_Averages.csv')
    % - variable_name: name of the variable being averaged (e.g., 'Evaporation Rate')

    % Extract only the date part
    date_only = dateshift(date_time, 'start', 'day');
    
    % Calculate daily average
    [unique_dates, ~, idx] = unique(date_only);
    daily_avg = accumarray(idx, data, [], @mean); % Average for each day

    % Ensure unique_dates and daily_avg are column vectors
    if isrow(unique_dates)
        unique_dates = unique_dates'; % Convert to column if needed
    end
    if isrow(daily_avg)
        daily_avg = daily_avg'; % Convert to column if needed
    end

    % Create a table with the results
    daily_avg_table = table(unique_dates, daily_avg, ...
        'VariableNames', {'Date', ['Daily_Avg_', variable_name]});

    % Define full file path
    full_file_path = fullfile(folder_path, file_name);

    % Write the table to the file
    writetable(daily_avg_table, full_file_path);
end



%% Function: Get Date and Time from Excel
function t = get_time(filename,sheetname)
    opts = detectImportOptions(filename);
    A = readtable(filename,opts,'Sheet',sheetname);
    t = datetime(A{:,1});
end

%% Function: Get Value
function value = get_value(filename,sheetname)
    opts = detectImportOptions(filename);
    A = readtable(filename,opts,'Sheet',sheetname);
    value = table2array(A(:,2:end));
end

% %% Function: Import Water Level and Corresponding Date and Time
% function [h, t] = data_import(filename)
%     opts = detectImportOptions(filename);
%     A = readtable(filename,opts); % import whole table from excel
% 
%     % extract water level and time to array
%     h = table2array(A(:,3)); % water level date at row 3
%     t = datetime(A{:,2}); % time data at row 2
% end

%% Function: Evaporation Rate Calculation
function [eva_rate,eva_date_time] = eva_rate_cal(h,date_time)
    eva_rate = [];
    eva_date_time = [];
    
    for i = 1:(length(h)-2)
        window = [h(i:i+2)];
        duration = hours(time(between(date_time(i),date_time(i+2),'time')));
        if issorted(window,'descend')
            rate = (window(1) - window(3))/duration;
            eva_rate = [eva_rate rate];
            eva_date_time = [eva_date_time  mean(date_time(i:i+2))]; % time corresponding to the rate was 
        end
    end
end

%% Function: Weather Factor Plot
function weather_factor(filename,eva_date_time,eva_rate)
    sheets = sheetnames(filename);
    num_plot = length(sheets);
    figure;
    
    for i = 1:num_plot
        A = readtable(filename,'Sheet',sheets(i),'PreserveVariableNames',true);
        header = A.Properties.VariableNames;
        station = header(2:end);
        
        value = table2array(A(:,2:end));
        t = datetime(A{:,1});
        subplot(num_plot,1,i);
        if sheets(i) ~= "Precipitation,in"
            plot(t,value,'LineWidth',1);
            legend(station);
            % legend('USC00044881','KCALEEVI12','BTN (CDEC)','USS0019L13S','USR0000CBEN','USR0000CBR4','USR0000CCRE','NOHRSC');
            ylabel(sheets(i));
        else
            bar(t,-1*value);
            hold on;
            bar(eva_date_time,eva_rate,'b',EdgeColor='b');
            hold off;
            legend([station,{'Evaporation'}]);
            % legend('USC00044881','KCALEEVI12','BTN (CDEC)','USS0019L13S','USR0000CBEN','USR0000CBR4','USR0000CCRE','NOHRSC');
            ylabel(sheets(i));
        end
    end
end

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