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

%% Combine Water Temperature and Evaporation Rate
% Combine Water Temperature and Evaporation Rate
% Trim the evaporation rate and water temperature data to have the same date range
[common_date_range, trimmed_eva_rate, trimmed_temp, trimmed_eva_time, trimmed_temp_time] = trim_to_common_date_boundaries(eva_rate, eva_date_time, T_all, T_all_date_time);

% Plot both evaporation rate and water temperature on the same graph
figure;

% Plot the evaporation rate
yyaxis left;
bar(trimmed_eva_time, trimmed_eva_rate, 'b',EdgeColor='b');
ylabel('Evaporation Rate (mm/hr)');
ylim([0, max(trimmed_eva_rate) + 1]); % Adjust limits for better visualization

% Plot the water temperature
yyaxis right;
plot(trimmed_temp_time, trimmed_temp, '-r');
ylabel('Water Temperature (°C)');
ylim([min(trimmed_temp) - 1, max(trimmed_temp) + 1]); % Adjust limits for better visualization

% Add labels and title
xlabel('Date');
title('Evaporation Rate and Water Temperature Over Time (by Date)');
grid on;

%% Water Temperature and Evaporation Rate Correspondance

% % Cut after Aug 15
% % Define the cutoff date
% cutoff_date = datetime('2024-08-15'); % Adjust year if needed
% 
% % Filter out data after the cutoff date
% filter_idx_temp = trimmed_temp_time <= cutoff_date;
% filter_idx_eva = trimmed_eva_time <= cutoff_date;
% 
% % Apply the filter
% trimmed_temp_time = trimmed_temp_time(filter_idx_temp);
% trimmed_temp = trimmed_temp(filter_idx_temp);
% 
% trimmed_eva_time = trimmed_eva_time(filter_idx_eva);
% trimmed_eva_rate = trimmed_eva_rate(filter_idx_eva);

% Find nearest evaporation rate for each temperature measurement
[nearest_eva_rate, ~] = nearest_neighbor(trimmed_temp_time, trimmed_eva_time, trimmed_eva_rate);

% Plot temperature vs evaporation rate
figure;

% Plot temperature against evaporation rate
plot(trimmed_temp, nearest_eva_rate, '.');
xlabel('Water Temperature (°C)');
ylabel('Evaporation Rate (mm/hr)');
title('Evaporation Rate vs. Water Temperature');
grid on;


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