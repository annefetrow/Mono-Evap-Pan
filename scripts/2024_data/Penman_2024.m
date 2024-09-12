%% Initialization
clearvars; close all; clc;
global RHO_W RHO_L SIGMA S P Cp_w Cp_a LAMBDA

RHO_W = 1000;
RHO_L = 1060; % kg/m3
SIGMA = 5.67 * 10^(-8); % 1/(s m2 K4)
S = 81; % g/kg
P = 0.1 .* 10^6; % Pa
Cp_w = 3.7794*10^(-3); % MJ/kg K
Cp_a = 1.005*10^(-3); % MJ/kg K
LAMBDA = 2.25; % MJ/kg

%% Import Data
% Import weather factor data
table_weather_factor = read_weather_factor('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\data\eva pan\Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx');

% Import evaporation data
table_eva_rate = read_eva_rate('C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\daily_average_evaporation_rates.csv'); % mm/hr

% Ensure that the datetime columns are properly formatted
table_weather_factor.Date = datetime(table_weather_factor.Date);
table_eva_rate.Date = datetime(table_eva_rate.Date);

% Perform the inner join based on the datetime column
combined_table = innerjoin(table_weather_factor, table_eva_rate, 'Keys', 'Date');

%% Calculation
% Calculate delta
delta = delta_calc(Ferenheit_to_Kelvin(table2array(combined_table(:,2)))); % Pa per K

%% Specific Heat
gamma = salinity_sigma_calc; % kPa/K

%% Radiation
Rn = Rn_calc(Ferenheit_to_Kelvin(table2array(combined_table(:,2))), 1000); % W/m2

%% Check Ea with theoretical Ea
Ea_theo = 6.43 .* (0.18 + 0.55 .* 0.44704 .* combined_table.("Windspeed,mph")) .* (saturated_pressure(combined_table.("Temperature,F")) - 10^(-3) .* vapor_pressure(Ferenheit_to_Kelvin(combined_table.("Temperature,F"))));
% figure;
% plot(combined_table.Date, Ea_theo);
% hold on;
% plot(combined_table.Date,Ea);
% legend('Theoretical','Experimental');

%% Calculate En
Ea = table2array(combined_table(:,7)).*RHO_W.*2.45.*24./1000; % MJ/m2d
E_pan = penman_calc(delta./1000, gamma, Rn.*0.0864, Ea) .* 1000; % mm/d
E_theo = penman_calc(delta./1000, gamma, Rn.*0.0864, Ea_theo) .* 1000; % mm/d

plot(combined_table.Date,E_pan);
ylabel('mm/d');
xlabel('Date');
hold on;
plot(combined_table.Date,E_theo);
plot(combined_table.Date,combined_table.daily_average_evaporation_rates.*24);

legend('Lake, from Pan Eva','Lake, from Theoretical Eva','Pan Eva');

%% Save the Results
% Create a table with the calculated data
output_table = table(combined_table.Date, E_pan, E_theo, combined_table.daily_average_evaporation_rates.*24, ...
    'VariableNames', {'Date', 'Lake_From_Pan_Eva_mm_d', 'Lake_From_Theoretical_Eva_mm_d', 'Pan_Eva_mm_d'});

% Write the table to a CSV file
writetable(output_table, 'C:\Users\24468\Desktop\Research\SEAS-HYDRO\Mono Lake\Mono-Evap-Pan\output\eva_estimate_Penman.csv');

%% Function
function E = penman_calc(delta, gamma, Rn, Ea) % delta in kPa/C, gamma in kPa/C, Ea in MJ/m^2, Rn in MJ/m2d
    
    % Global Variables
    global RHO_L LAMBDA

    % Define Constants
    G = 0; % MJ/m2d

    % Calculate E
    E = 1 ./ LAMBDA .* (delta .* (Rn - G) + gamma .* Ea) ./ ((delta + gamma) .* RHO_L);

end

function Rn = Rn_calc(Ta, RA)
    
    % Global Variables
    global SIGMA

    % Define Constants
    r = 0.05;
    a = 0.4; % assumption, vary with latitude
    b = 0.274; % assumption, vary with latitude
    n_N = 0.8; % assumption
    %ed = vapor_pressure_height(vapor_pressure(Ta), Ta, 4).*10^(-3); % Pa, assume measured at 4m height
    ed = vapor_pressure_height(vapor_pressure(Ta), Ta, 4).*10^(-3); % Pa, assume measured at 4m height


    % Calculate R1
    R1 = RA .* (a + b .* n_N);

    % Calculate R2
    RB = SIGMA .* Ta.^4 .*(0.56 - 0.09 .* sqrt(ed)) .* (0.1 + 0.9 .* n_N);

    % Calculate Rn
    Rn = R1 .* (1-r) - RB;

end

function p_s = saturated_pressure(T) %kPa
    T_C = 5./9.*(T-32);
    p_s = 0.611 .* exp(17.27.*T_C ./ (T_C + 237.3));
end

function pw = vapor_pressure(T) % pw in Pa, T in K
    A = 5.04221;
    B = 1838.675;
    C = -31.737;
    
    pw = 10.^((A - B./(C+T))); % bar

    % Convert unit from bar to Pa
    pw = 10^5 .* pw;
end

function ph = vapor_pressure_height(p0, T, height)
    M = 18;% g/mol
    gravity = 9.8;
    R = 8.314;
    ph = p0 .* exp(-1 .* (M.*gravity.*height) ./ (R.*T));
end

function p = salinity_vapor_pressure(T) % p in Pa, T in K

    % Global Variables
    global S

    % Calculate saturated vapor pressure for pure water
    pw = vapor_pressure(T); % pw in Pa

    % Convert to salty water
    a1 = -2.1609*10^(-4);
    a2 = -3.5015*10^(-7);

    p = pw .* 10.^(a1.*S + a2.*S.^2); % Pa
end

function delta = delta_calc(T) % T in K
    
    dT = 0.1;
    delta = (salinity_vapor_pressure(T + dT) - salinity_vapor_pressure(T - dT)) ./ (2 .* dT);

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
