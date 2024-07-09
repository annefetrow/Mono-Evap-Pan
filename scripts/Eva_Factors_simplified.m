%% Initialization
clearvars; close all; clc;
syms x

%% Data Import
[h_1,date_time_1] = data_import('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\ML_EP_20230907.csv');
[h_2,date_time_2] = data_import('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\2023-11-2MLevap (2).csv');

% merge the two water level data and time data
h = [h_1;h_2];
date_time = [date_time_1;date_time_2];

%% Evaporation Rate
[eva_rate,eva_date_time] = eva_rate_cal(h,date_time);

eva_rate_avg = mean(eva_rate(eva_rate<0.75))%mm/hr

%% Weather Factor
weather_factor('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx',eva_date_time,eva_rate);

%% Plot Data (Evaporation+Precipitation)
figure;
P_d = get_value('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx','Precipitation,in');
t = get_time('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx','Precipitation,in');
bar(eva_date_time,eva_rate,'b',EdgeColor='b');
%bar(date_time,h,'b',EdgeColor='b');
hold on;
bar(t,-1*P_d*25.4/24,'r',EdgeColor='none');
plot(eva_date_time,ones(1,length(eva_date_time))*eva_rate_avg,'--b');
xlabel('Date Time');
ylabel('mm/hr');

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