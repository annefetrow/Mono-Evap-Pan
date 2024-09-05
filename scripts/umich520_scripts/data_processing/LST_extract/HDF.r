## The purpose of this script is to unpack HDF files from remote sensing
## products for land surface temperature (LST) and other parameters related
## to the Mono Lake evaporation project for CEE 520 (April 2024)

install.packages("terra")
library(terra)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5")


file <- H5Fopen("C:/cygwin64/home/Caleb/MOD11A1.A2023250.h08v05.061.2023252003330.hdf", "H5F_ACC_RDONLY")

h5ls(file)

data <- h5read(file, "datasetName")

H5Fclose(file)
