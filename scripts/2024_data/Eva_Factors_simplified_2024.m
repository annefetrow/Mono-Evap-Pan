%% Initialization
clearvars; close all; clc;
syms x

set(gca, 'LineWidth', 1.5); % Set line width for axes
set(0, 'DefaultLineLineWidth', 1.5); % Set default line width for all lines

%% Water Level
% Path to the folder containing the CSV files for water level
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

eva_rate_avg = mean(eva_rate(eva_rate<3)); %mm/hr, exclude out liers that have extremely high eva rate

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

xlabel('Date');
ylabel('Eva Rate (mm/hr)');
title('Eavporation Rate Over Time');
hold off;

% Calculate daily evaporation rate and export to a file
daily_eva_rate_avg = calculate_and_export_daily_avg(eva_rate, eva_date_time, 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\2024', '2024_daily_average_evaporation_rates.csv','Evaporation Rate mm/hr');

%% Water Temperature

% Path to the folder containing the CSV files
folderPath = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\Raw data – Evaporation Pan\Water Temp';

% Read and combine data from all CSV files
[T_water_all, T_water_all_date_time] = read_and_combine_csv(folderPath,1,1); % most recent 1 file

% Plot temperature with time
figure;
plot(T_water_all_date_time,T_water_all, 'r-');
xlabel('Date');
ylabel('Temperature (C)');
title('Water Temperature over Time');
grid on;
axis tight;

% Calculate daily water temperature and export to a file
daily_water_temp_avg = calculate_and_export_daily_avg(T_water_all, T_water_all_date_time, 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\2024', '2024_daily_average_water_temperature.csv','Water Temp C');

%% Air Temperature

% Path to the folder containing the CSV files
folderPath = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\Raw data – Evaporation Pan\Air Temp';

% Read and combine data from all CSV files
[T_air_all, T_air_all_date_time] = read_and_combine_csv(folderPath,1,1); % most recent 1 file

% Plot temperature with time
figure;
plot(T_air_all_date_time,T_air_all, 'g-');
xlabel('Date');
ylabel('Temp (C)');
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

% Plot RH with time
figure;
plot(RH_all_date_time,RH_all, 'k-');
xlabel('Date');
ylabel('RH (%)');
title('Relative Humidity over Time');
grid on;
axis tight;

% Calculate daily RH and export to a file
daily_RH_avg = calculate_and_export_daily_avg(RH_all, RH_all_date_time, 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\2024', '2024_daily_average_RH.csv','RH %');

%% Plot Water Temperature, Air Temperature, and Relative Humidity
plotByMonth(T_water_all_date_time, T_water_all, T_air_all_date_time, T_air_all, RH_all_date_time, RH_all);

%% Combine Water Temperature, Evaporation Rate, and Air Temperature
[trimmed_eva_time, trimmed_eva_rate] = plotMonthlyEvapWaterAir(eva_rate, eva_date_time, T_water_all, T_water_all_date_time, T_air_all, T_air_all_date_time);

%% Water Temperature and Evaporation Rate Correspondance
plot_evaporation_data(trimmed_eva_time, trimmed_eva_rate, T_air_all_date_time, T_air_all, T_water_all_date_time, T_water_all, RH_all_date_time, RH_all, datetime('2024-08-15'));

%% Save the plot to folder plots
% Specify the full path where you want to save the PDF file
fileLocation = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\plots\eva_rate_and_climate_factors.pdf';

% Get the list of all open figures
figures = findall(groot, 'Type', 'figure');

% Loop through each figure and append to the same PDF
for i = 1:length(figures)
    % Specify to append to the PDF after the first page
    if i == 1
        exportgraphics(figures(i), fileLocation, 'ContentType', 'vector'); % First page
    else
        exportgraphics(figures(i), fileLocation, 'ContentType', 'vector', 'Append', true); % Append subsequent pages
    end
end


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

%% Function: Plot all climate factors together by different months
function plotByMonth(T_water_all_date_time, T_water_all, T_air_all_date_time, T_air_all, RH_all_date_time, RH_all)
    % Get unique months and years from the datetime arrays
    months = month(T_water_all_date_time);
    years = year(T_water_all_date_time);
    uniqueMonths = unique([years, months], 'rows');

    % Create figure
    figure;
    nMonths = size(uniqueMonths, 1);

    for i = 1:nMonths
        % Extract the data for the current month
        currentMonth = uniqueMonths(i, 2);
        currentYear = uniqueMonths(i, 1);

        % Select data corresponding to the current month and year
        idxWater = (month(T_water_all_date_time) == currentMonth) & (year(T_water_all_date_time) == currentYear);
        idxAir = (month(T_air_all_date_time) == currentMonth) & (year(T_air_all_date_time) == currentYear);
        idxRH = (month(RH_all_date_time) == currentMonth) & (year(RH_all_date_time) == currentYear);

        % Create subplot for the current month (single column)
        subplot(nMonths, 1, i); % One column, multiple rows

        % Plot water and air temperature on the left y-axis
        yyaxis left;
        plot(T_water_all_date_time(idxWater), T_water_all(idxWater), '-r', 'DisplayName', 'Water Temp');
        hold on;
        plot(T_air_all_date_time(idxAir), T_air_all(idxAir), '-g', 'DisplayName', 'Air Temp');
        ylabel('Temp (°C)');
        ylim([min([T_water_all(idxWater); T_air_all(idxAir)]) - 1, max([T_water_all(idxWater); T_air_all(idxAir)]) + 1]);
        
        % Plot relative humidity on the right y-axis
        yyaxis right;
        plot(RH_all_date_time(idxRH), RH_all(idxRH), '--k', 'DisplayName', 'RH');
        ylabel('RH (%)');
        ylim([min(RH_all(idxRH)) - 10, max(RH_all(idxRH)) + 10]);

        % Add labels, title, and grid
        xlabel('Date');
        title(sprintf('Month: %d-%d', currentMonth, currentYear));
        grid on;

        % Show the legend for the first subplot
        if i == 1
            legend('show');
        end
    end

    % Adjust the layout of the figure
    sgtitle('Water Temperature, Air Temperature, and Relative Humidity by Month');
end

%% Function: Plot all evaporation rate with air and water temperature by different months
function [trimmed_eva_time, trimmed_eva_rate] = plotMonthlyEvapWaterAir(eva_rate, eva_date_time, T_water_all, T_water_all_date_time, T_air_all, T_air_all_date_time)
    % Trim the data to the same date range
    [common_date_range, trimmed_eva_rate, trimmed_temp, trimmed_eva_time, trimmed_temp_time] = ...
        trim_to_common_date_boundaries(eva_rate, eva_date_time, T_water_all, T_water_all_date_time);
    
    % Filter air temperature data to the common date range
    filter_idx_air = T_air_all_date_time >= common_date_range(1) & T_air_all_date_time <= common_date_range(2);
    trimmed_T_air_time = T_air_all_date_time(filter_idx_air);
    trimmed_T_air = T_air_all(filter_idx_air);

    % Extract unique months in the common date range
    unique_months = unique(month(trimmed_eva_time));
    num_months = numel(unique_months);

    % Create figure with subplots
    figure;
    for i = 1:num_months
        % Extract data for the current month
        month_idx = month(trimmed_eva_time) == unique_months(i);
        eva_time_month = trimmed_eva_time(month_idx);
        eva_rate_month = trimmed_eva_rate(month_idx);

        temp_time_month = trimmed_temp_time(month(trimmed_temp_time) == unique_months(i));
        temp_month = trimmed_temp(month(trimmed_temp_time) == unique_months(i));

        air_time_month = trimmed_T_air_time(month(trimmed_T_air_time) == unique_months(i));
        air_month = trimmed_T_air(month(trimmed_T_air_time) == unique_months(i));

        % Create a subplot for each month
        subplot(num_months, 1, i);
        
        % Plot evaporation rate (left y-axis)
        yyaxis left;
        bar(eva_time_month, eva_rate_month, 'b', 'EdgeColor', 'b');
        ylabel('Eva Rate (mm/hr)');
        ylim([0, max(eva_rate_month) + 1]);
        
        % Plot water temperature (right y-axis)
        yyaxis right;
        plot(temp_time_month, temp_month, '--r');
        hold on;
        plot(air_time_month, air_month, '--g'); % Air temperature
        ylabel('Water Temp / Air Temp (°C)');
        ylim([min([temp_month; air_month]) - 1, max([temp_month; air_month]) + 1]);
        
        % Add labels and title
        title(sprintf('Evaporation Rate, Water Temp, and Air Temp - Month %d', unique_months(i)));
        if i == num_months
            xlabel('Date');
        end
        grid on;
    end

    % Add a legend to the last subplot
    legend('Evaporation Rate', 'Water Temperature', 'Air Temperature');
end
