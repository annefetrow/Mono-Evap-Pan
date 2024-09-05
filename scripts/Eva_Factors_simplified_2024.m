%% Initialization
clearvars; close all; clc;
syms x

set(gca, 'LineWidth', 2); % Set line width for axes
set(0, 'DefaultLineLineWidth', 2); % Set default line width for all lines

%% Data Import
[h_1,date_time_1] = data_import('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\ML_EP_20230907.csv');
[h_2,date_time_2] = data_import('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\2023-11-2MLevap (2).csv');

% merge the two water level data and time data
h = [h_1;h_2];
date_time = [date_time_1;date_time_2];

% Path to the folder containing the CSV files
folderPath = 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\Raw data – Evaporation Pan\Pan Water Level';

% Read and combine data from all CSV files
[all_h, all_date_time] = read_and_combine_csv(folderPath,5,7); % most recent 3 files

%% Plot Water Level
figure;
plot(all_date_time, all_h, '-');
xlabel('Date and Time');
ylabel('Water Level (mm)');
title('Water Level over Time');
grid on;

%% Evaporation Rate
[eva_rate,eva_date_time] = eva_rate_cal(all_h,all_date_time);

eva_rate_avg = mean(eva_rate(eva_rate<0.75))%mm/hr

%% Plot Data (Evaporation)
figure;
bar(eva_date_time,eva_rate,'b',EdgeColor='b');
%bar(date_time,h,'b',EdgeColor='b');
hold on;
plot(eva_date_time,ones(1,length(eva_date_time))*eva_rate_avg,'--b');
xlabel('Date Time');
ylabel('mm/hr');
title('Eavporation Rate Over Time');
hold off;

% %% Weather Factor
% weather_factor('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx',eva_date_time,eva_rate);
% 
% %% Plot Data (Evaporation+Precipitation)
% figure;
% P_d = get_value('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx','Precipitation,in');
% t = get_time('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx','Precipitation,in');
% bar(eva_date_time,eva_rate,'b',EdgeColor='b');
% %bar(date_time,h,'b',EdgeColor='b');
% hold on;
% bar(t,-1*P_d*25.4/24,'r',EdgeColor='none');
% plot(eva_date_time,ones(1,length(eva_date_time))*eva_rate_avg,'--b');
% xlabel('Date Time');
% ylabel('mm/hr');

%% Function: Read Water Level Data
function [unique_h, unique_date_time] = read_and_combine_csv(folderPath,startfile,endfile)
    % read_and_combine_csv reads CSV files from a specified folder
    % and combines their water level and date_time data.
    %
    % Inputs:
    %   folderPath - Path to the folder containing the CSV files
    %
    % Outputs:
    %   all_h - Combined water level data
    %   all_date_time - Combined date_time data
    
    % Get a list of all CSV files in the folder
    filePattern = fullfile(folderPath, '*.csv');
    csvFiles = dir(filePattern);
    
    % Initialize arrays to hold the combined data
    all_h = [];
    all_date_time = [];
    
    % Loop through each file and read the data
    for k = startfile:endfile
        baseFileName = csvFiles(k).name;
        fullFileName = fullfile(folderPath, baseFileName);
        
        % Import the data from the CSV file
        [h, date_time] = data_import(fullFileName);
        
        % Append the data to the combined arrays
        all_h = [all_h; h];
        all_date_time = [all_date_time; date_time];

        % Remove repetitive entry
        [unique_h, unique_date_time] = remove_duplicates_and_nans(all_h, all_date_time);
    end

    % Optionally, sort the data by date_time if needed
    [unique_date_time, sortIdx] = sort(unique_date_time);
    unique_h = unique_h(sortIdx);
end

%% Function: Remove Repetitive Data Entry
function [unique_h, unique_date_time] = remove_duplicates_and_nans(all_h, all_date_time)
    % remove_duplicates_and_nans removes duplicate date_time entries
    % and NaN values from the data while keeping the corresponding water level data.
    %
    % Inputs:
    %   all_h - Combined water level data
    %   all_date_time - Combined date_time data
    %
    % Outputs:
    %   unique_h - Water level data with duplicates and NaNs removed
    %   unique_date_time - date_time data with duplicates and NaNs removed

     % Remove rows where water level data (all_h) is NaN
    validIdx = ~isnan(all_h);
    all_h = all_h(validIdx);
    all_date_time = all_date_time(validIdx);
    
    % Convert date_time to a numeric format for easier processing if necessary
    % all_date_time = datenum(all_date_time); % Uncomment if date_time is a datetime object
    
    % Remove duplicate date_time entries, keeping the first occurrence
    [unique_date_time, unique_idx] = unique(all_date_time, 'stable');
    unique_h = all_h(unique_idx);
    
    % Optionally, convert date_time back to datetime if you converted it earlier
    % unique_date_time = datetime(unique_date_time, 'ConvertFrom', 'datenum'); % Uncomment if needed
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

%% Function: Import Water Level and Corresponding Date and Time
function [h, t] = data_import(filename)
    opts = detectImportOptions(filename);
    A = readtable(filename,opts); % import whole table from excel
    
    % extract water level and time to array
    h = table2array(A(:,3)); % water level date at row 3
    t = datetime(A{:,2}); % time data at row 2
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