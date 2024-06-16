%% Initialization
clearvars; close all; clc;

%% Data Import
A = readtable('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\daily_average_evaporation_rates.csv');
B = readtable('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx','Sheet','Temperature,F','PreserveVariableNames',true);
header = B.Properties.VariableNames;
station = header(2:end); % need 1,6,7 (7 for pan

T_tot = table2array(B(:,2:end));
t_cli = datetime(B{:,1});
       
t_eva = datetime(A{:,1});
evaRate_mm_hr = table2array(A(:,2));

%% Calculate Temperature for Lake and Pan
T_lake = (T_tot(:,1) + T_tot(:,6) + T_tot(:,7))/3;
T_pan = T_tot(:,7);