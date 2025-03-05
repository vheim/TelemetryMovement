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
mystatus = lapply(all_status, read.csv,sep=",",dec=".",stringsAsFactor=F,header=T)
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
# myargos = lapply(all_argos, read.csv,sep=",",dec=".",stringsAsFactor=F,header=T)
# myargos <- do.call("rbind.fill",myargos) #rbind.fill() is a dplyr function that drops columnnames as we have dfs with different colnames

## All data
all_all = dir(spotloc, recursive=T, full.names=T, pattern="\\-All.csv$") # import files in folders in path directory all at once
myall = lapply(all_all, read.csv,sep=",",dec=".",stringsAsFactor=F,header=T)
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
    DateEST = as.Date(DatetimeEST),
    Date_local = format(DatetimeEST, "%d-%b-%y"),
    Transmits = as.numeric(Transmits),
    Temperature = if_else(Temperature < -5 | Temperature > 50, NA, Temperature), # Temperature readings below 0 or above 50 dont make no sense
    WetDry = as.numeric(WetDry),
    MinWetDry = as.numeric(MinWetDry),
    MaxWetDry = as.numeric(MaxWetDry)
  ) %>%
  dplyr::select(
    Ptt,
    DatetimeUTC,
    DatetimeEST,
    DateEST,
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
    DateEST = as.Date(DatetimeEST),
    Date_local = format(DatetimeEST, "%d-%b-%y"),
    BestRx = as.numeric(Best.level),
  ) %>%
  dplyr::select(
    Ptt,
    DatetimeUTC,
    DatetimeEST,
    DateEST,
    Date_local,
    BestRx
  )

### ....................................................................................................
### [D] Plot diagnostic plots ----
### ....................................................................................................

# D0: prepare loop (not yet ready)

## create a list with all IDs
fishlist <- unique(status$Ptt)

for (i in fishlist){
  ## subset the dfs
  # i <- 264015
  status_i <- status[which(status$Ptt == i),]
  all_i <- all[which(all$Ptt == i),]

  ## open plotting window
  png(filename = paste0(saveloc, "/SPOT_output/Deployment_diagnostics_Ptt_",i,"_", today(),".png"), width = 297, height = 210, units = "mm", bg = "white", res = 150, family = "", type = "cairo-png")

  ## define plot window parameters
  par(mar = c(2.85,3.1,2,0.85), mgp = c(1.95, 0.6,0),
      oma = c(0, 0, 3, 0))  # Outer margins for the overall title

  ## setup multipanel layout - this is for 5 plots and might need to be adjusted
  layout(matrix(c(1,1, # first plot (covers top spaces 1-2)
                  2,2, # second plot (covers top spaces 3-4)
                  3,3, # third plot (covers top spaces 5-6)
                  0, # blank (covers bottom space 1 (0 = blank plot)
                  4,4,# fourth plot (covers bottom spaces 2-3)
                  5,5, # 5th plot (covers bottom spaces 4-5)
                  0), # blank (covers bottom space 6)
                nrow = 2, byrow = TRUE))

# D1: Battery voltage plot ----

## Plot
plot(status_i$DatetimeEST, status_i$BattVoltage,
     pch = 21,
     bg = alpha(ifelse(status_i$BattVoltage >= 3.2, "blue", "red"),0.75),
     col = "black",
     xaxt = "n", # remove x axis labels
     cex.axis = .9,
     las = 1,
     xlab = expression(bold("Date")),
     ylab = expression(bold("Voltage [V]")))

## label the axes
### X-axis
# mtext(expression(bold(paste(delta^{15}, "N (\u2030)"))), side = 2, line = 1.5, font = 2, las = 0, cex = 0.7)
axis(1, at = c(min(status_i$DatetimeEST), status_i$DatetimeEST[which(status_i$BattVoltage < 3.2)[1]], max(status_i$DatetimeEST)),
     labels = c(status_i$Date_local[which.min(status_i$DatetimeEST)], status_i$Date_local[which(status_i$BattVoltage < 3.2)[1]], status_i$Date_local[which.max(status_i$DatetimeEST)]),
     cex.axis = 0.9,
     tck = -0.025,
     las = 1)

## add a legend
legend("bottomleft",
       legend = c(expression(paste("\u2265", " 3.2 V")),  # ≥ symbol
                  # expression(paste("\u2264", " 3.2 V"))),  # ≤ symbol
                  expression(paste("<", " 3.2 V"))),
       pch = 21,
       pt.bg = c("blue","red"),
       horiz = T,
       cex = 1.1, bty = "n")

## label the plot
title(main = expression(bold("Battery Voltage ")))

# D2: Transmission attempts plot ----

## Plot
plot(status_i$DatetimeEST, status_i$Transmits,
     pch = 21,
     bg = alpha("blue",0.75),
     col = "black",
     xaxt = "n", # remove x axis labels
     cex.axis = .9,
     las = 1,
     xlab = expression(bold("Date")),
     ylab = expression(bold("# Transmits")))

## label the axes
### X-axis
# mtext(expression(bold(paste(delta^{15}, "N (\u2030)"))), side = 2, line = 1.5, font = 2, las = 0, cex = 0.7)
axis(1, at = c(min(status_i$DatetimeEST), max(status_i$DatetimeEST)),
     labels = c(status_i$Date_local[which.min(status_i$DatetimeEST)], status_i$Date_local[which.max(status_i$DatetimeEST)]),
     cex.axis = 0.9,
     tck = -0.025,
     las = 1)

## add a legend
###NA

## as this is the center plot in the upper row we add main title for multipanel figure here
## Last step: add overall title to multipanel figure
mtext(bquote(bold(paste("Deployment diagnostics for Ptt ID ", .(i)))), line = 2.5)

## label the plot
title(main = expression(bold("Transmit Attempts ")))

# D3: Temp at time of status message plot ----

## Plot
plot(status_i$DatetimeEST, status_i$Temperature,
     pch = 21,
     bg = alpha("blue",0.75),
     col = "black",
     xaxt = "n", # remove x axis labels
     cex.axis = .9,
     las = 1,
     xlab = expression(bold("Date")),
     ylab = expression(bold("Temperature [°C]")))

## label the axes
### X-axis
# mtext(expression(bold(paste(delta^{15}, "N (\u2030)"))), side = 2, line = 1.5, font = 2, las = 0, cex = 0.7)
axis(1, at = c(min(status_i$DatetimeEST), max(status_i$DatetimeEST)),
     labels = c(status_i$Date_local[which.min(status_i$DatetimeEST)],  status_i$Date_local[which.max(status_i$DatetimeEST)]),
     cex.axis = 0.9,
     tck = -0.025,
     las = 1)

## add a legend
### NA

## label the plot
title(main = expression(bold("Temperature at time of status message")))

# D4: WetDry patterns plot----

## Plot WetDry values
plot(status_i$DatetimeEST, status_i$WetDry,
     # ylim = c(0,max(status_i$WetDry+10, na.rm = T)),
     ylim = c(0,255),
     pch = 3,
     lwd = 2,
     # bg = alpha(ifelse(status_i$BattVoltage >= 3.2, "blue", "red"),0.75),
     col = "forestgreen",
     xaxt = "n", # remove x axis labels
     cex.axis = .9,
     las = 1,
     xlab = expression(bold("Date")),
     ylab = expression(bold("WetDry")))

## add minimum values
points(status_i$DatetimeEST, status_i$MinWetDry,
       pch = 25, bg = "blue", col = "black")

## add maximum values
points(status_i$DatetimeEST, status_i$MaxWetDry,
       pch = 24, bg = "goldenrod", col = "black")

## label the axes
### X-axis
# mtext(expression(bold(paste(delta^{15}, "N (\u2030)"))), side = 2, line = 1.5, font = 2, las = 0, cex = 0.7)
axis(1, at = c(min(status_i$DatetimeEST), max(status_i$DatetimeEST)),
     labels = c(status_i$Date_local[which.min(status_i$DatetimeEST)], status_i$Date_local[which.max(status_i$DatetimeEST)]),
     cex.axis = 0.9,
     tck = -0.025,
     las = 1)

## add a legend
legend("bottomright",
       legend = c(expression(paste("Current")),
                  expression(paste("DailyMin")),
                  expression(paste("DailyMax"))),
       pch = c(3,25,24),
       pt.bg = c("forestgreen", "blue","goldenrod"),
       col = c("forestgreen", "black", "black"),
       pt.lwd = c(2,1,1),
       horiz = T,
       cex = 1.1, bty = "n")

## label the plot
title(main = expression(bold("WetDry Patterns")))

# D5: Best Rx Level plot ----

## Plot db of Signal values
plot(all_i$DateEST, all_i$BestRx,
     # ylim = c(0,max(status_i$WetDry+10, na.rm = T)),
     pch = 21,
     bg = alpha("blue", 0.75),
     col = "black",
     xaxt = "n", # remove x axis labels
     cex.axis = .9,
     las = 1,
     xlab = expression(bold("Date")),
     ylab = expression(bold("dB of signal")))

## add average
daily_avg <- all_i %>% # calculate mean first
  dplyr::group_by(DateEST) %>%
  dplyr::summarise(dailymean = mean(BestRx))

lines(daily_avg$DateEST, daily_avg$dailymean, # add mean as a line
       type = "l", col = "red", lwd = 2)

points(min(daily_avg$DateEST, na.rm = T), daily_avg$dailymean[which.min(daily_avg$DateEST)],
       pch = 21, bg = "red", col = "black", cex = 2) # add first BestRx level

points(max(daily_avg$DateEST, na.rm = T), daily_avg$dailymean[which.max(daily_avg$DateEST)],
       pch = 21, bg = "red", col = "black", cex = 2) # add first BestRx level


## label the axes
### X-axis
# mtext(expression(bold(paste(delta^{15}, "N (\u2030)"))), side = 2, line = 1.5, font = 2, las = 0, cex = 0.7)
axis(1, at = c(min(all_i$DateEST, na.rm = T), max(all_i$DateEST, na.rm = T)),
     labels = c(all_i$Date_local[which.min(all_i$DateEST)], all_i$Date_local[which.max(all_i$DateEST)]),
     cex.axis = 0.9,
     tck = -0.025,
     las = 1)

## add a legend
legend("topright",
       legend = c(expression(paste("BestRx level")),
                  expression(paste("Daily mean level"))),
       pch = c(21, 21),
       pt.bg = c("blue", "red"),
       col = c("black", "black"),
       pt.lwd = c(1,1),
       horiz = T,
       cex = 1.1, bty = "n")

## close and save plot window
dev.off()
}


# END OF SCRIPT ----
#### TODO LIST ####
# TODO1: ----
# check with WC where which diagnostic variable is taken from
# implement diagnostics for Corrupted and uncorrupted plot, queued data summaries
# check with Devon WC why transmit attempts has outliers in 264015
# check with Devon regarding Temperature column in status file with readings up to 3000°C??
