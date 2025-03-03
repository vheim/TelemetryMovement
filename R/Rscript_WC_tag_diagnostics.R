# TOP OF CODE ----
### ====================================================================================================
### Project:    NA
### Analysis:   Run a series of diagnostics and visualise tag performance for WC satellite tags
### Script:     Rscript_WC_tag_diagnostics.R
### Author:     Vital Heim
### Version:    1.0
### ====================================================================================================

### ....................................................................................................
### Content: this script contains the code to run and visualise a series of diagnostics (incl. plots)
###          for satellite tags from Wildlife Computers.
###          At the moment the written code focuses on data from Smart Position and Temperature tags but
###          but the code will be expanded on once additional tag types are deployed/got their diagnostics
###          checked.
###
###          For SPOT tags we check the following:
###          Battery voltage over time (-Status.csv, column; BattVoltage)
###          Transmission attempts over time (-Status.csv, column; Transmits)
###          Temperature at time of status message (-Status.csv, column; Temperature)
###          WetDry information: current, dailymin and dailymax (-Status.csv, columns; WetDry, MinWetDry, MaxWetDry)
###          Best Rx level (-Argos.csv, column; Power)
###          Corrupted vs uncorrupted receptions (-Argos.csv, column; Corrupt)
###          Queued data summaries (-Status.csv, column; Xmit Queue)
###          Fastloc Success vs. Failure (not yet implemented)
###
### Please see TODO list at end of script for open issues.
### ....................................................................................................

### ....................................................................................................
### [A] Ready environment, load packages ----
### ....................................................................................................

# A1: clear memory ----

rm(list = ls())

# A2: load necessary packages ----

## if first time
# install.packages("tidyverse")
# install.packages("magrittr")
# install.packages("patchwork")
# install.packages("ggplot2")

## load packages
library(tidyverse)
library(magrittr)
library(patchwork)
library(ggplot2)

# A3: Specify dataloc and saveloc ----

## you need to manually copy-paste the newest datafolders into the input folder
spotloc <- file.path("//Sharktank","Science","Projects_current", "2025_WC_tag_diagnostics","Data_input","SPOT_input") # Adjust this
saveloc <- file.path("//Sharktank","Science","Projects_current", "2025_WC_tag_diagnostics","Data_output") # Adjust this

# A4: define needed functions, universal variables etc.----

## avoid error if timeout reached
options(timeout = 3000) # manually increase time out threshold (needed when downloading basemap)

### ....................................................................................................
### [B] Data import ----
### ....................................................................................................

# B1: Import data ----

## Status data
spotdata <- list.files(path = spotloc,
                       pattern = paste0("\\-Status.csv"),
                       recursive = TRUE,
                       full.names = TRUE ) %>%
  purrr::map_dfr(~read_csv(.x) %>%
                   mutate(DeployID = as.character(DeployID),
                          LocationQuality = as.character(LocationQuality)))
### TODO have to adjust all columns so they are the same across
### DEAL WITH warnings()


### ....................................................................................................
### [C] Data housekeeping and preparation ----
### ....................................................................................................


### ....................................................................................................
### [D] Plot animated data with fixed plot window for individual animals ----
### ....................................................................................................


### ....................................................................................................
### [E] Plot animated data with dynamic axes for multiple individuals ----
### ....................................................................................................


# END OF SCRIPT ----
#### TODO LIST ####
# TODO1: ----
# add location of data per csv file in descripiton box
