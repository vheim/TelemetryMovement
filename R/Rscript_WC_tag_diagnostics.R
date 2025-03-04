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
###          Battery voltage over time (-Status.csv, column; BattVoltage) - done
###          Transmission attempts over time (-Status.csv, column; Transmits) - done, double check with WC why outliers in 264015 email
###          Temperature at time of status message (-Status.csv, column; Temperature) - done
###          WetDry information: current, dailymin and dailymax (-Status.csv, columns; WetDry, MinWetDry, MaxWetDry) - done
###          Best Rx level (-All.csv, column; Best level) - done
###          Corrupted vs uncorrupted receptions (-Argos.csv, column; Corrupt) - TBD
###          Queued data summaries (-Status.csv, column; Xmit Queue) - done
###          Fastloc Success vs. Failure (not yet implemented) - TBD
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
# install.packages("lubridate")
# install.packages("plyr")
# install.packages("ggplot2")

## load packages
library(tidyverse)
library(magrittr)
library(lubridate)
library(plyr)
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
all_status = dir(spotloc, recursive=T, full.names=T, pattern="\\-Status.csv$") # import files in folders in path directory all at once
mystatus = lapply(all_status, read.csv,sep=",",dec=".",stringsAsFactor=F,header=T) # import all .csv files containing TAT-Hiso data, but skip header lines
mystatus <- do.call("rbind.fill",mystatus) #rbind.fill() is a dplyr function that drops columnnames as we have dfs with different colnames

## For the fun of it, here is a way to import all dataframes as one using the purrr package
# spotdata <- list.files(path = spotloc,
#                        pattern = paste0("\\-Status.csv"),
#                        recursive = TRUE,
#                        full.names = TRUE ) %>%
#   purrr::map_dfr(~read_csv(.x) %>%
#                    mutate(DeployID = as.character(DeployID),
#                           LocationQuality = as.character(LocationQuality))) ### DEAL WITH warnings()


## Argos data
# all_argos = dir(spotloc, recursive=T, full.names=T, pattern="\\-Argos.csv$") # import files in folders in path directory all at once
# myargos = lapply(all_argos, read.csv,sep=",",dec=".",stringsAsFactor=F,header=T) # import all .csv files containing TAT-Hiso data, but skip header lines
# myargos <- do.call("rbind.fill",myargos) #rbind.fill() is a dplyr function that drops columnnames as we have dfs with different colnames

## All data
all_all = dir(spotloc, recursive=T, full.names=T, pattern="\\-All.csv$") # import files in folders in path directory all at once
myall = lapply(all_all, read.csv,sep=",",dec=".",stringsAsFactor=F,header=T) # import all .csv files containing TAT-Hiso data, but skip header lines
myall <- do.call("rbind.fill",myall) #rbind.fill() is a dplyr function that drops columnnames as we have dfs with different colnames

### ....................................................................................................
### [C] Data housekeeping and preparation ----
### ....................................................................................................

# C1: Status info - keep necessary columns ----

## We need
## Ptt info
## Datetime info
## Battery voltage over time - column; BattVoltage
## Transmission attempts over time - column; Transmits
## Temperature at time of status message - column; Temperature
## WetDry information: current, dailymin and dailymax - columns; WetDry, MinWetDry, MaxWetDry

status <- mystatus %>%
  dplyr::mutate(
    Ptt = as.character(Ptt),
    DatetimeUTC = as.POSIXct(Received,format="%H:%M:%S %d-%b-%Y", tz="UTC", usetz = T),
    DatetimeEST = lubridate::with_tz(DatetimeUTC, tzone = "US/Eastern"),
    Date_local = format(DatetimeEST, "%d-%b-%y"),
    Transmits = as.numeric(Transmits),
    WetDry = as.numeric(WetDry),
    MinWetDry = as.numeric(MinWetDry),
    MaxWetDry = as.numeric(MaxWetDry)
  ) %>%
  dplyr::select(
    Ptt,
    DatetimeUTC,
    DatetimeEST,
    Date_local,
    BattVoltage,
    Transmits,
    Temperature,
    WetDry,
    MinWetDry,
    MaxWetDry
  )


# C2: Argos info - keep necessary columns ----

## We need
## Ptt info
## Datetime info
## Best Rx level - column; Power
## Corrupted vs uncorrupted receptions - column; Corrupt

# argos <- myargos %>%
#   dplyr::mutate(
#     Ptt = as.character(Ptt),
#     DatetimeUTC = as.POSIXct(Date,format="%H:%M:%S %d-%b-%Y", tz="UTC", usetz = T),
#     DatetimeEST = lubridate::with_tz(DatetimeUTC, tzone = "US/Eastern"),
#     Date_local = format(DatetimeEST, "%d-%b-%y"),
#     Power = as.numeric(Power),
#   ) %>%
#   dplyr::select(
#     # select needed columns only
#   )

# C3: All info - keep necessary columns ----

## We need
## Ptt info
## Datetime info
## Best Rx level - column; Best level


all <- myall %>%
  dplyr::mutate(
    Ptt = as.character(Platform.ID.No.),
    DatetimeUTC = as.POSIXct(Loc..date,format="%m/%d/%Y %H:%M:%S", tz="UTC", usetz = T),
    DatetimeEST = lubridate::with_tz(DatetimeUTC, tzone = "US/Eastern"),
    Date_local = format(DatetimeEST, "%d-%b-%y"),
    BestRx = as.numeric(Best.level),
  ) %>%
  dplyr::select(
    Ptt,
    DatetimeUTC,
    DatetimeEST,
    Date_local,
    BestRx
  )

### ....................................................................................................
### [D] Plot diagnostic plots ----
### ....................................................................................................

# D0: prepare loop (not yet ready)

status_i <- status %>%
  dplyr::filter(
    Ptt == "264015"
  ) %>%
  dplyr::mutate(
    Index = row_number()
  )

par(mar = c(2.85,3.1,2,0.85), mgp = c(1.95, 0.6,0))

# D1: Battery voltage plot ----

## Plot
plot(status_i$DatetimeEST, status_i$BattVoltage,
     pch = 21,
     bg = alpha(ifelse(status_i$BattVoltage >= 3.2, "blue", "red"),0.75),
     col = "black",
     xaxt = "n", # remove x axis labels
     cex.axis = .7,
     las = 1,
     xlab = expression(bold("Date")),
     ylab = expression(bold("Voltage [V]")))

points(centroids2, col = "black", pch = 10, cex = 2, lwd = 4)

## label the axes
### X-axis
# mtext(expression(bold(paste(delta^{15}, "N (\u2030)"))), side = 2, line = 1.5, font = 2, las = 0, cex = 0.7)
axis(1, at = c(min(status_i$DatetimeEST), max(status_i$DatetimeEST)),
     labels = c(status_i$Date_local[which.min(status_i$DatetimeEST)], status_i$Date_local[which.max(status_i$DatetimeEST)]),
     cex.axis = 0.7,
     tck = -0.035,
     las = 1)
##TODO add a lable where votlage drops below 3.2V

## add a legend
legend("bottomleft",
       legend = c(expression(paste("\u2265", " 3.2 V")),  # ≥ symbol
                  # expression(paste("\u2264", " 3.2 V"))),  # ≤ symbol
                  expression(paste("<", " 3.2 V"))),
       pch = 21,
       pt.bg = c("blue","red"),
       horiz = T,
       cex = 1, bty = "n")

## label the plot
title(main = expression(bold("Battery Voltage")))


### ....................................................................................................
### [E] Plot animated data with dynamic axes for multiple individuals ----
### ....................................................................................................


# END OF SCRIPT ----
#### TODO LIST ####
# TODO1: ----
# check with WC where which diagnostic variable is taken from
# implement diagnostics for Corrupted and uncorrupted plot, queued data summaries
# check with Devon WC why transmit attempts has outliers in 264015
