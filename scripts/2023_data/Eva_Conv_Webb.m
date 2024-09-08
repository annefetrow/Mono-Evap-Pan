%% Initialization
clearvars; close all; clc;
set(groot, 'defaultLineLineWidth', 1);  % Sets the default line width to 2

%% Constant
S = 81;
height = 4;

%% Data Import
A = readtable('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\daily_average_evaporation_rates.csv');
B = readtable('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx','Sheet','Temperature,F','PreserveVariableNames',true);
header = B.Properties.VariableNames;
station = header(2:end); % need 1,7 (7 for pan

T_tot = table2array(B(:,2:end));
t_cli = datetime(B{:,1});
       
t_eva = datetime(A{:,1});
evaRate_mm_hr = table2array(A(:,2));

% Find min and max datetime values in the smaller array
minDatetime = min(t_eva);
maxDatetime = max(t_eva);

% Use logical indexing to crop the larger array 
% to the date & time range of the smaller array
mask = (t_cli >= minDatetime & t_cli <= maxDatetime);

t_cli = t_cli(mask);

%% Monthly Lake Temperature Data
% Create a table with month names and average temperatures:
MonthNames = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'}';
AverageTemperatures_F = [30, 32, 35, 45, 55, 60, 69.6, 71.8, 62.275, 52.75, 49.8, 35]';

% Create a table:
monthly_avg_temp_table = table(MonthNames, AverageTemperatures_F);
monthly_avg_temp_table.SalinityPressure_kPa = salinity_vapor_pressure(vapor_pressure(AverageTemperatures_F)./7.501, S, (monthly_avg_temp_table.AverageTemperatures_F - 32) .* 5./9 + 273.15);
%monthly_avg_temp_table.SalinityPressure_kPa = vapor_pressure(AverageTemperatures_F)./7.501;

%% Calculate Temperature for 4m above Ground and Pan
T_4 = (T_tot(:,1) + T_tot(:,7))/3; % assume lake and pan have the same temperature above ground and it's the average of 3 stations
T_pan = T_tot(:,7); % assume pan's temperature is the same as air temperature
T_4 = T_4(mask);
T_4 = T_4(T_4 ~= 0) - 2.5*(6770 - 6383.8)/1000;
T_pan = T_pan(mask);
T_pan = 0.8*T_pan(T_pan ~= 0);

%% Calculate Saturation Vapor Pressure
e_p_kPa = vapor_pressure(T_pan);
e_p_4_kPa = vapor_pressure_height(e_p_kPa, T_pan, height);

e_l_kPa = zeros(length(T_4),1);
T_lake = zeros(length(T_4),1);

for i = 1:length(e_l_kPa)
    monthIndex = month(t_cli(i));
    e_l_kPa(i) = monthly_avg_temp_table.SalinityPressure_kPa(monthIndex);
    T_lake(i) = AverageTemperatures_F(monthIndex);
end

e_l_4_kPa = vapor_pressure_height(e_l_kPa, T_4, height);

%% Calculate Lake Evaporation
% E_lake = 1.5 .* ((e_l_kPa - e_l_4_kPa) ./ (e_p_kPa - e_p_4_kPa)) .* evaRate_mm_hr;
coefficient = 1.5 .* ((e_l_kPa - e_l_4_kPa) ./ (e_p_kPa - e_p_4_kPa));

%% Plot
plot(t_cli,coefficient);
ylabel('Conversion Coefficient');
xlabel('Date');

figure;
plot(t_cli,e_p_4_kPa);
hold on;
plot(t_cli,e_p_kPa);
plot(t_cli,e_l_kPa);
plot(t_cli,e_l_4_kPa);
hold off;
legend('e_{p,4}','e_p','e_l','e_{l,4}');
ylabel('Vapor Pressure (kPa)');
xlabel('Date');

figure;
plot(t_cli,T_lake);
hold on;
plot(t_cli,T_pan);
plot(t_cli,T_4);
hold off;
xlabel('Date');
ylabel('Temperature');
legend('Lake','Pan','4 m Above Ground');

%% Function
function p_s = saturated_pressure(T) %kPa
    T_F = 5./9.*(T-32);
    p_s = 0.611 .* exp(17.27.*T_F ./ (T_F + 237.3));
end

function pw = vapor_pressure(T) % mmHg
    A = 8.07131;
    B = 1730.63;
    C = 233.426;
    
    pw = 10.^((A - B./(C+T)));
end

function ph = vapor_pressure_height(p0, T, height)
    M = 18;% g/mol
    gravity = 9.8;
    R = 8.314;
    ph = p0 .* exp(-1 .* (M.*gravity.*height) ./ (R.*T));
end


function p = salinity_vapor_pressure(pw, S, T) %T in K
    a1 = -2.1609*10^(-4);
    a2 = -3.5015*10^(-7);

    p = pw .* 10.^(a1.*S + a2.*S.^2);
end