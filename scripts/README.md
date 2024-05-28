# Scripts

---

This directory should hold all script files that are sourced/loaded from your RMarkdown or Jupyter notebooks.

## `R_from_MATLAB.R`

Translated from `Eva_Factors_simplified.m` file that Yolanda processes all data. 

1. Create bar plot for evaporation rate and precipitation (precipitation from climate stations)
2. Calculate average evaporation rate
3. Create line plot for different topographical factors, data stored in  `Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx` in the `data` folder

## Evaporation Calculation Method

During Aug, 2023 to Nov, 2023, eva pan data are logged every 2 hours. Refill was automatically done for every 2 hours to reach a certain water level and no precipitation data is available.

Therefore, Yolanda used a *mid-point approximation method* to evaluate evaporation rate.

### Assumptions

1. If the water level has an increase within 6 hour window, there's precipitation or refill happening
2. Evaporation rate during the 6 hour are averaged to get a mean value instead of instantaneous results

### Method

1. Loop through the whole water level data set by evaluating the water level change for every 6 hours
2. If the water level increases, `eva rate = 0`
3. If the water level decreases, `eva_rate = 6 hour water level change / 6 hour`
4. Bar plot evaporation rate along with precipitation results obtained from climate stations