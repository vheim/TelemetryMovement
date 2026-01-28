# TOP OF CODE ----
### ====================================================================================================
### Project:    Large bodied hammerhead complex in the western North Atlantic
### Analysis:   Processing and cleaning satellite telemetry data of fin-mounted SPOT tags for further steps
### Script:     ~2025_NWA_hammerhead_complex/R/Rscript1_Filter_SPOT_data_argosfilter.R
### Author:     Vital Heim
### Version:    1.0
### ====================================================================================================

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Content: this R script contains the code to filter, clean and process raw data collected by fin-
###          mounted Smart Position and Temperature (SPOT) transmitters. We use the argosfilter package
###          by Freitas et al. 2008.
###          Here processed data can then be used in further scripts and steps to model and analyze
###          SPOT tag data.
###          The filtering step is parallelized allowing to filter multiple IDs at once.
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Note: Parts of this R script were strongly influenced by a tutorial on the R package crawl by
###       London J.M. and Johnson D. S. (see here: https://jmlondon.github.io/crawl-workshop/index.html),
###       which is a very useful resource and I can highly recommend for people working with movement data
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### ....................................................................................................
### [A] Setwd, paths and parameters ----
### ....................................................................................................

# A1: clear memory ----

rm(list = ls())

# A2: load necessary packages ----

## install
#install.packages("devtools")
require("devtools")

# devtools::install_version("argosfilter", version = "0.70")
# install.packages("trip")
# install.packages("tidyverse")
# install.packages("magrittr")
# install.packages("patchwork")
# install.packages("sf")
# install.packages("sp")
# install.packages("purrr")
# install.packages("furrr")
# install.packages("pander")
# install.packages("ggplot2")
# install.packages("ggspatial")
# install.packages("RColorBrewer")
# install.packages("xts")
# install.packages("prettymapr")

## load
library(tidyverse)
library(magrittr)
library(patchwork)
library(sf)
library(sp)
library(purrr)
library(furrr)
library(pander)
library(ggplot2)
library(ggspatial)
library(RColorBrewer)
library(prettymapr)
library(xts)
library(argosfilter) # to filter raw Argos data, i.e. SPOT tag data
#library(trip)

# A3: Specify needed functions

## Function to deal with near duplicate  Argos observations
make_unique <- function(x) {
  xts::make.time.unique(x$date,eps = 10) # eps = number of seconds to make unique
}

# A4: Specify data and saveloc ----

YOUR_IP <- "NA" # add your IP address, server name or similar if you connect via a shared drive

## Input data
dataloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2025_NWA_hammerhead_complex","Data_input", "Telemetry")

## Output data
saveloc <- file.path("/", YOUR_IP,"Science","Projects_current", "2025_NWA_hammerhead_complex","Data_output", "Telemetry") # Adjust this

# A5: Define universal variables (e.g. for plotting) ----

# NA

# A6: Define universal options, variables, etc. (e.g. for plotting) ----

options(scipen=999) # so that R doesn't act up for pit numbers  

### ...............................................................................................
### [B] Data import and initial check and filtering of missing information ----
### ................................................................................................

# B1: Import data ----

## Movement data
all_csv = dir(dataloc, recursive=T, full.names=T, pattern="\\-Locations.csv$") # import files in folders in path directory all at once
mydata = lapply(all_csv, read.csv,sep=",",dec=".",stringsAsFactor=F,header=T) # import all .csv files containing TAT-Hiso data, but skip header lines
mydets <- do.call("rbind",mydata)
# sort(unique(mydets$Ptt)) # 96 IDs
# [1]  23596 174515 174516 175427 177940 177941 177942 179468 179471 179472 180910 180911 180912
# [14] 180913 180914 180915 182746 183619 183620 183621 183622 183623 183624 198201 198202 198203
# [27] 198204 198205 198206 200367 200368 200369 201510 209020 210602 210603 210604 210605 210606
# [40] 222133 222134 222135 222136 222137 222138 233536 233537 233538 235278 235281 235283 244597
# [53] 244598 244599 244600 244601 244602 244603 244604 244605 244606 244607 244608 244609 244610
# [66] 244611 261302 261303 261304 261305 261306 261307 261308 261309 261742 261743 264015 264019
# [79] 264020 264021 264022 264023 264024 264025 264026 264029 284465 286887 286888 286889 286890
# [92] 286895 286896 286897 286898 286899

## Shark metadata
tags_all <- read.csv(file = file.path(dataloc, "Data_AllSpecies_SPOT_tags_metadata.csv"),sep=",",dec=".",header=T,na.strings=c(""," ",NA))
tags_all$datetime_deployment_local <- as.POSIXct(tags_all$datetime_deployment_local,format="%Y/%m/%d %H:%M:%S",tz="US/Eastern")
attr(tags_all$datetime_deployment, "tzone") <- "UTC"
tags <- dplyr::select(tags_all, ptt_id, group, species, sex, datetime_deployment, deployment_lat, deployment_lon, pcl, fl, stl)
colnames(tags) <- c("id", "group", "species", "sex", "date", "lat", "lon", "pcl", "fl", "stl")
tags$lc <- "G" # we add a location class criteria for later joining with the movement data. we define the LC as "G" for gps, so that we do not need to worry about smaj,smin,eor
tags$id <- as.character(tags$id)

### if needed, filter for target species and
# sort(unique(tags$group)); sort(unique(tags$species))
# [1] "Andros"                "Bimini"                "Florida Keys"         
# [4] "Jupiter"               "Jupiter, FL"           "Marquesas"            
# [7] "San Jose del Cabo, MX" "South Carolina"        "Tampa"                
# [1] "C. carcharias" "C.plumbeus"    "S.lewini"      "S.mokarran"    "S.zygaena"

tags <- tags %>%
  dplyr::filter(
    species %in% c("S.mokarran", "S.lewini","S.zygaena"),
    group %in% c("Florida Keys", "Jupiter", "Jupiter, FL", "Marquesas", "South Carolina", "Tampa")
  )

### make a df with just the target ptt-ids
target_ids <- tags[,1]

### make a df with just the tagging date
tagging_date <- tags[,c(1,5)] # df for filtering of osbervations pre-deployment

### make a df with just the tagging location information
tagging_location <- tags[,c(1,5,11,7,6)] # tagging location as observation for later ctcrw fitting
tagging_location$smaj <- 50
tagging_location$smin <- 50
tagging_location$eor <- 0
colnames(tagging_location) <- c("id", "date", "lc", "lon", "lat", "smaj", "smin", "eor")
# head(tagging_location)
# id                date lc       lon      lat smaj smin eor
# 1 179468 2024-02-17 18:27:00  G -79.99082 26.98042   50   50   0
# 2 179471 2021-04-21 17:18:00  G -80.72854 24.75576   50   50   0
# 3 179472 2019-08-22 12:30:00  G -80.09000 32.35550   50   50   0
# 4 180911 2019-05-01 04:00:00  G -80.38520 32.44750   50   50   0
# 5 180912 2019-05-07 12:30:00  G -80.18100 32.36770   50   50   0
# 6 180914 2019-05-09 12:30:00  G -79.53040 32.85320   50   50   0

### make a df with all the potential post-release mortalities
pprm <- tags_all[which(tags_all$status == "pprm"), 1]
# > pprm
# [1] 180911 180913 180915 201510 222135 200367 244599 174516 264020 ## this is all shark species all tagging locs!

# B2: Basic housekeeping ----

ddet <-mydets %>%
  dplyr::select( # select relevant columns, here: id, date, location class (lc), lon, lat,
    Ptt,
    Date,
    Type,
    Quality,
    Longitude,
    Latitude,
    Error.Semi.major.axis,
    Error.Semi.minor.axis,
    Error.Ellipse.orientation
  ) %>%
  filter( # locations that were user specified within the Wildlife Data Portal
    !Type %in% c("User"),
    !Quality %in% c("Z"),
    !Ptt %in% pprm, # remove post release mortality sharks or tags without data
    Ptt %in% target_ids # filter for target ids only
  ) %>%
  mutate( # define Date format
    Date = as.POSIXct(Date,format="%H:%M:%S %d-%b-%Y", tz="UTC", usetz = T),
    Ptt = as.character(Ptt) # define tag id as character class
  ) %>%
  dplyr::rename( # rename the columns so they fit the requirements for the fit functions
    id = Ptt,
    date = Date,
    lc = Quality,
    lon = Longitude,
    lat = Latitude,
    smaj = Error.Semi.major.axis,
    smin = Error.Semi.minor.axis,
    eor = Error.Ellipse.orientation
  ) %>%
  dplyr::select( # remove Type column as it is not needed later
    -Type
  )

## add tagging date for subsequent filtering
ddet$tagging.date <- NA

for (i in 1:nrow(ddet)){
  ddet[i, 9] <- as.character(tagging_date[which(tagging_date$id == ddet[i, 1]), 2])
}

## remove detections pre-tag deployment and and tagging location as data point
ddet %<>%
  filter( # remove occurences that happened before the release time
    date >= tagging.date
  ) %>%
  dplyr::select( #get rid of unneeded columns
    -tagging.date
  ) %>%
  bind_rows( # add tagging location as observation of class "gps"
    semi_join(tagging_location, ddet, by = "id")
  ) %>%
  arrange( # arrange by timestamp by individual so df can be used for fit_() functions
    id,
    date
  )

## if you need/want, remove detections after a certain date
## this is useful if you have active tags, but need to write a report or plan on submitting a manuscript soon

filter_needed <- "no" # change this
max_date <- as.POSIXct("2025-07-16 00:00:00",format="%Y-%m-%d %H:%M:%S", tz="UTC", usetz = T) # change this

if (filter_needed == "yes"){ # do NOT change this
  ddet <- ddet %>%
    dplyr::filter(
      !date >= as.POSIXct(max_date)
    )
} else (
  ddet <- ddet
)

# B3: Deal with near duplicate observations ----

# Argos data often contain near duplicate records. These are identified by location estimates
# with the same date-time but differing coordinate or error values. In theory, SSM/CTCRW models should be able
# to deal with these situations, but according to https://jmlondon.github.io/crawl-workshop/crawl-practical.html#duplicate-times
# it is more reliable to fix these occurences beforehand.
# The first option for fixing the records would be to eliminate one of the duplicate records.
# However, reliably identifying which record is more suited to be kept or to be discarded can be difficult and near impossible
# Alas, it might be better to adjust the timestamp for one of the observations
# and increasing the value by 10 second. To do this, we will rely on the
# xts::make.time.unique() function.

ddet_tc <- ddet %>% # create time corrected df
  dplyr::arrange(id,date) %>%
  dplyr::group_by(id) %>% tidyr::nest() %>%
  dplyr::mutate(unique_time = purrr::map(data, make_unique)) %>%
  tidyr::unnest(cols = c(data, unique_time)) %>%
  dplyr::select(-date) %>% rename(date = unique_time)

## check for dublicated time stamps (should be 0 now...)
dup_times <- ddet_tc %>% group_by(id) %>%
  filter(duplicated(date)); nrow(dup_times) # should be 0 after correcting for near duplicates

## data summary
ddet_tc %>% dplyr::group_by(id) %>%
  dplyr::summarise(num_locs = n(),
                   start_date = min(date),
                   end_date = max(date)) ## We see that there are also datapoints in there from the tag initiation, we will deal with this later

### write data summary to .csv
write.csv(ddet_tc %>% dplyr::group_by(id) %>%
            dplyr::summarise(num_locs = n(),
                             start_date = min(date),
                             end_date = max(date)),
          file.path(saveloc, "Data_NWA_hammerheads_SPOT_location_summary_first_to_last.csv"))

# B4: visualise raw data ----

sf_ddet <- sf::st_as_sf(ddet_tc, coords = c("lon","lat")) %>%
  sf::st_set_crs(4326)

sf_lines <- sf_ddet %>%
  dplyr::arrange(id, date) %>%
  sf::st_geometry() %>%
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_ddet$id))) %>%
  sf::st_cast("MULTILINESTRING") %>%
  sf::st_sf(deployid = as.factor(unique(sf_ddet$id)))

esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
                     'Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')

## Define the number of colors you want
nb.cols <- length(unique(ddet$id))
mycolors <- colorRampPalette(brewer.pal(nb.cols, "YlOrRd"))(nb.cols)

ggplot() +
  annotation_map_tile(type = esri_ocean,zoomin = 1,progress = "none") +
  layer_spatial(sf_ddet, size = 0.5) +
  layer_spatial(sf_lines, size = 0.75,aes(color = deployid)) +
  scale_x_continuous(expand = expansion(mult = c(.6, .6))) +
  scale_fill_manual(values = mycolors) +
  theme() +
  ggtitle("Observed Argos Location Paths - Raw data",
          subtitle = paste0("SPOT tagged Sphyrna spp. from the U.S. Atlantic (n = ", length(unique(ddet_tc$id)), ")"))

## save if needed
ggsave(file.path(saveloc,"Raw_argos_detections_pre_sda_filter.tiff"),
       width = 21, height = 15, units = "cm", device ="tiff", dpi=300)

### ....................................................................................................
### [C] Filter data based on speed, distance and turning angles using argosfilter::sdafilter() ----
### ....................................................................................................

## To filter the track based on speed and turning angles, we use the package
## "argosfilter" from Freitas et al. 2008

## citations("argosfilter")
# Freitas C (2022). _argosfilter: Argos Locations
# Filter_. R package version 0.70,
# <https://CRAN.R-project.org/package=argosfilter>.

## Information below is from the package vignette and correspondin paper:
# Freitas C., Lydersen C., Fedak M. A., Kovacs K. M. (2008). A simple new algorithm to filter marine
# mammal Argos locations. Marine Mammal Science 24(2): 315-325

## To filter data based on turning angles, speed, distance etc. one can use the sdafilter()
## function. The locations are filtered using the Freitas et al. 2008 algorithm.

## Locations are filtered using the algorithm described in Freitas et al. (2008). The algorithm first
## removes all locations with location class Z (-9), which are the points for which the location process
## failed. Then all locations requiring unrealistic swimming speeds are removed, using the MacConnell et al. (1992) algorithm, unless the point is located at less than 5 km from the previous
## location. This procedure enables retaining good quality locations for which high swimming speeds
## result from location being taken very close to each other in time. The default maximum speed
## threshold is 2 m/s. The last step is optional, and enables to remove unlikely spikes from the animal’s path. The angles of the spikes should be specified in ang, and their respective length in
## distlim. The default is c(15,25) for ang and c(2500,5000) for distlim, meaning that all spikes
## with angles smaller than 15 and 25 degrees will be removed if their extension is higher than 2500
## m and 5000 m respectively. No spikes are removed if ang=-1. ang and distlim vectors must have
## the same length

## The output will be a vector with the following elements:
## "removed" = location removed by filter
## "not" = location not removed
## "end_location" = location a tht end of the track where the algorithm could not be applied

## The function needs the following arguments
## lat = numeric vector of latitudes in decimal degrees
## lon = numeric vector of longitudes in decimal degrees
## dtime = a vector of class POSIXct with date and time for each location
## lc = a numeric or character vector of Argos location classes. Z classes can be entered as "Z", "z", or -9
## vmax = speed threshold in ms-1. Default is 2ms-1
## ang = anles of the spikes to be removed. Default is c(15,25), No spikes are removed if ang = -1.
## distlim = lengths of the above spikes, in meters. Default is c(2500,5000).

# This analysis can be run in series, however, we can also parallelise it so
# that multiple ids can be filtered at once. The parallel package can be used for taking advantage of multiple processors.
# However, we can also keep working with the purrr and nested column tibble data structure. The multidplyr package is
# therefore a viable option and can be installed via the devtools package and source code hosted on GitHub.

# C1: apply filter ----

future::plan(multisession)

prefilter_obs <- ddet_tc

ddet_af <- ddet_tc %>%
  dplyr::arrange(id, date) %>%
  dplyr::group_by(id) %>%
  tidyr::nest()

#tbl_locs %>% dplyr::summarise(n = n())

ddet_af <- ddet_af %>%
  dplyr::mutate(filtered = furrr::future_map(data, ~ argosfilter::sdafilter(
    lat = .x$lat,
    lon = .x$lon,
    dtime = .x$date,
    lc = .x$lc,
    vmax = 2.1, # based on Payne et al. 2017
    ang = c(15,25),
    distlim = c(5000,8000) # based on Vaudo et al. 2017
  ))) %>%
  tidyr::unnest(cols = c(data, filtered)) %>%
  dplyr::filter(filtered %in% c("not", "end_location")) %>%
  dplyr::select(-filtered) %>%
  dplyr::arrange(id,date) %>%
  dplyr::select( # not to select but to order the columns so they are in order for further use
    id,
    date,
    lc,
    lon,
    lat,
    smaj,
    smin,
    eor
  )

cat("You removed ", nrow(prefilter_obs)-nrow(ddet_af)," locations.")
# You removed  4956  locations.
prefilter_obs %>% group_by(id) %>% dplyr::summarise(n = n())
ddet_af %>% dplyr::summarise(n = n())

## write csv to show new data structure
write.csv(ddet_af %>% dplyr::summarise(n = n()),
          file.path(saveloc, "Data_NWA_hammerheads_SPOT_nr_locations_post_sda_filter.csv"))

## Visualise the filtered tracks

# esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
#                      'Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')

af_sf <- sf::st_as_sf(ddet_af, coords = c("lon","lat")) %>%
  sf::st_set_crs(4326)

af_lines <- af_sf %>%
  dplyr::arrange(id, date) %>%
  sf::st_geometry() %>%
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(af_sf$id))) %>%
  sf::st_cast("MULTILINESTRING") %>%
  sf::st_sf(id = as.factor(unique(af_sf$id)))

ggplot() +
  annotation_map_tile(type = esri_ocean,zoomin = 1,progress = "none") +
  layer_spatial(sf_ddet, size = 0.5) +
  layer_spatial(af_lines, size = 0.75,aes(color = id)) +
  scale_x_continuous(expand = expansion(mult = c(.6, .6))) +
  #scale_fill_manual() +
  theme() +
  ggtitle("Argos detections with argosfilter::sdafilter()",
          subtitle = paste0("SPOT tagged Sphyrna spp. from the U.S. Atlantic (n = ", length(unique(ddet$id)), ");
argsofilter::sdafilter() removes ", nrow(prefilter_obs)-nrow(ddet_af)," locations."))

## save if needed
ggsave(file.path(saveloc,"Argos_detections_with_argosfilter_sdafilter.tiff"),
       width = 21, height = 15, units = "cm", device ="tiff", dpi=300)

# C2: Save the data ----

start <- as.Date(min(ddet_af$date), format = "%Y-%m")
end <- as.Date(max(ddet_af$date), format = "%Y-%m")

## txt
#write.table(ddet_af, paste0(saveloc, "Argosfilter_filtered_Sphyrna_SPOT_tracks_multiID_", start, "_", end,".txt"), row.names=F,sep=",",dec=".")

## csv
write.table(ddet_af, file.path(saveloc, paste0("Argosfilter_filtered_Sphyrna_SPOT_tracks_multiID_", start, "_", end,".csv")), row.names=F,sep=",",dec=".")
# write.table(ddet_af, paste0("C:/Users/Vital Heim/switchdrive/Science/Projects_and_Manuscripts/Andros_Hammerheads/InputData/CTCRW/", "Argosfilter_filtered_Sphyrna_SPOT_tracks_multiID_", start, "_", end,".csv"), row.names=F,sep=",",dec=".")

## RDS
saveRDS(ddet_af, file.path(saveloc, paste0("Argosfilter_filtered_Sphyrna_SPOT_tracks_multiID_", start, "_", end,".R")))
saveRDS(ddet_af, file.path(dataloc, paste0("Argosfilter_filtered_Sphyrna_SPOT_tracks_multiID_", start, "_", end,".R")))
# saveRDS(ddet_af, paste0("C:/Users/Vital Heim/switchdrive/Science/Projects_and_Manuscripts/Andros_Hammerheads/InputData/CTCRW/", "Argosfilter_filtered_Sphyrna_SPOT_tracks_multiID_CTCRW_input.R"))

# END OF CODE ----