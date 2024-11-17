# Date Operation

---

This directory should holds all of our raw data files, no matter the format (e.g. Excel, CSV, image data, etc.). All data files should be treated as **read only**.

If you want to prevent any particular data file or data subfolder from being automatically added to GitHub, simply list their names in the `.gitignore` file that is in this directory.

**Note: Current Code are designed uder a specific data format. So if your code is not running please check the data format first.**

## Raw data- Evaporation Pan

The folder stores all the raw data from the evaporation pan and on-site relevant gauges, downloaded monthly by Maureen. Note that the common data renewal place is [Google Drive](https://drive.google.com/drive/folders/1xlh3-rdhs1PhuJu4rxLKHyPeRMTOjEsy?usp=drive_link) and it's further downloaded to GitHub for data processing purposes. Therefore, please make sure you download the latest Google folder for complete data set.

### Air Temp

Stores the temperature of air above the evaporation pan. The file data type should be `.csv` and it can be converted manually from `.xlsx`.

| Data Structure |                      |                          |
|---------------------------|----------------------|--------------------------|
| #                         | Date Time, GMT-07:00 | EV pan Level, Millimeter |
| 1                         | MM/DD/YY HH:MM       | num                      |
| 2                         | MM/DD/YY HH:MM       | num                      |

## Evaporation Pan Water Level

*File location: MONO-EVAP-PAN -> data -> eva pan*

The evaporation pan is refilled at 2:00 am everyday and the data is logged every 10 minutes.

Evapration pan data are stored in `.csv` file, with the following format:

| Plot Title: Eva Pan Level |                      |                          |
|---------------------------|----------------------|--------------------------|
| #                         | Date Time, GMT-07:00 | EV pan Level, Millimeter |
| 1                         | MM/DD/YY HH:MM       | num                      |
| 2                         | MM/DD/YY HH:MM       | num                      |

## Air Temperature

TBD

## Relative Humidity (w/ Radiation Shield)

TBD

## Topographic Factors from Climate Stations

Yolanda researched on nearby climate station on possible topographic factors that might influence evaporation, including temperature, precipitation, windspeed, relative humidity and solar radiation. 

All factors are put in one `Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx` file in multiple sheets. Different sheets are named as such: `Temperature,F` for the y-axis label of the plot.

`Mono Lake_Evaporation_Factors_Yolanda_Rawdata.xlsx` contains all hourly and daily data. However, note that the final plot exclude all hourly data and use  `Mono Lake_Evaporation_Factors_Daily_Yolanda.xlsx` instead.

Within each sheet, the data formatting are as follows

|  |  Station 1  |  Station 2  |
|---------------------------|----------------------|--------------------------|
| Date | Unit | Unit |
| MM/DD/YY | num | num |
| MM/DD/YY | num | num |

## Climate Station Location

The location of the station list is *MONO-EVAP-PAN -> data -> arcgis_climate_station_loc -> Eva_Station_Satellite_In_Use.xlsx*. It graphically shows where all stations are located.

Information includes 

* Agency name
* Station name and url address
* Station ID
* Elevation (ft)
* Latitude (N)
* Longitude (W)

