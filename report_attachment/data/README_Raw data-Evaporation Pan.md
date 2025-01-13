# Folder Information

---

The directory stores all the raw data from the evaporation pan and on-site relevant gauges, downloaded monthly by Maureen. Note that the common data renewal place is [Google Drive](https://drive.google.com/drive/folders/1xlh3-rdhs1PhuJu4rxLKHyPeRMTOjEsy?usp=drive_link) and it's further downloaded to GitHub for data processing purposes. Therefore, please make sure you download the latest Google folder for complete data set.

## Water Temp

Stores the temperature of pan water data.

| Data Sturcture |                      |                          |
|---------------------------|----------------------|--------------------------|
| #                         | Date-Time (PDT) | Temperature (C) |
| 1                         | MM/DD/YY HH:MM:SS       | `double`                       |
| 2                         | MM/DD/YY HH:MM:SS       | `double`                       |

## Mono Lake Surface Temp

Stores the surface temperature of Mono Lake, measured by Anne during her field trip. However, the data only extends from 6/1/2023 to 5/22/2024. 

| Data Sturcture |                      |                          |
|---------------------------|----------------------|--------------------------|
| #                         | datetime | Tsurf |
| 1                         | MM/DD/YY HH:MM       | `double`                       |
| 2                         | MM/DD/YY HH:MM       | `double`                       |

## Pan Water Level

Stores the water level change with time. The data is logged every 10 minutes and pan refill happens when (TBD by Anne).

| Plot Title: Eva Pan Level |                      |                          |
|---------------------------|----------------------|--------------------------|
| #                         | Date Time, GMT-07:00 | EV pan Level, Millimeter |
| 1                         | MM/DD/YY HH:MM       | `double`                      |
| 2                         | MM/DD/YY HH:MM       | `double`|

## Precipitation Gauge

Stores the precipitation happens on the pan. Note that the precipitation records the accumulative water level, but not the precipitation rate with time and therefore it's necessary to check climate station data for the time span of rain fall and intensity.


| Data Structure |                      |                          | | | | |
|---------------------------|----------------------|--------------------------|----------------------|----------------------|----------------------|----------------------|
| #                         |  Date  | Time  | Time since last measurement | Precip (cm) | Recorder | Notes
| 1                         | MM/DD/YY AM/PM     | HH:MM:SS    | |  `double`       | | |
| 2                         | MM/DD/YY   |     HH:MM:SS AM/PM | |`double`| | |


## Air Temp and RH

Stores the temperature of air above the evaporation pan. The file data type should be `.csv` and it can be converted manually from `.xlsx`.

| Data Structure |                      |                          | | |
|---------------------------|----------------------|--------------------------|----------------------|----------------------|
| #                         | Date-Time (PDT)  | Temperature (C)  | RH (%) | Dew Point (C) |
| 1                         | MM/DD/YY HH:MM:SS       | `double`    |`double` |  `double`       |
| 2                         | MM/DD/YY HH:MM:SS       |`double` |`double` |`double`|