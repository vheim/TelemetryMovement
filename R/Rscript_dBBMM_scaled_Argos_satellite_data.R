# TOP OF CODE ----
### ====================================================================================================
### Project:    Large bodied hammerhead complex in the western North Atlantic
### Analysis:   Utilisation distribution calculations using dBBMMs
### Script:     ~2025_NWA_hammerhead_complex/R/Rscript5_dBBMM_scaled_Argos_satellite_data.R
### Author:     Vital Heim
### Version:    1.0
### ====================================================================================================

### ....................................................................................................
### Content: this code allows to fit dynamic Brownian Bridge Movement Models (Kranstauber et al. 2012)
### to data collected by Argos-satellite-linked geolocators, specifically SPOT tags. The data used here
### was previously filtered and cleaned using various packages. The output will be scaled individual-
### and group-level Utilisation Distributions (UD) raster files.
###
### This script applies the functions provided in the movegroup package by
### Dedman & van Zinnicq Bergmann
###
### ....................................................................................................

### ....................................................................................................
### [A] Ready environment, load packages ----
### ....................................................................................................

# A1: clear memory ----

rm(list = ls())

# A2: install and load necessary packages ----

## if first time or needed updates
# install.packages("tidyverse")
# install.packages("magrittr")
# install.packages("lubridate")
# install.packages("sp")
# install.packages("sf")
# install.packages("dplyr")
# install.packages("tidylog")
# install.packages("beepr")
# install.packages("raster")
# install.packages("terra")
# install.packages("reshape2")
# install.packages("maptools")
# install.packages("data.table")
# install.packages("purrr")
# install.packages("vctrs")
# install.packages("ggplot2")
library(remotes)
# remotes::install_github("SimonDedman/movegroup", force = T)

## load packages and source needed functions
#library(dBBMMhomeRange)
library(movegroup)
#library(vctrs)
library(tidyverse)
library(magrittr)
library(lubridate)
library(sp)
library(sf)
library(dplyr)
# library(tidylog)
library(beepr)
library(raster)
library(terra)
library(reshape2)
library(data.table)
library(purrr)
library(ggplot2)

# A3: Specify needed functions ----

# *3.1: reprojecting coordinates
reproject <- function(x, coorda, coordb, latloncrs, projectedcrs) {
  x <- sf::st_as_sf(x, coords = c(coorda, coordb)) |>
    sf::st_set_crs(latloncrs) |> # latlon degrees sf object
    st_transform(projectedcrs) |> # eastings northings units metres
    dplyr::select(-everything()) # remove all columns. Geometry is protected and retained
  return(x)
}

# *A3.2.: function to calculate location error ----
# moveLocErrorCalc <- function(x,
#                              loncol = "lon",
#                              latcol = "lat",
#                              latloncrs = 4326,
#                              projectedcrs = 32617,
#                              lon025 = "lon025",
#                              lon975 = "lon975",
#                              lat025 = "lat025",
#                              lat975 = "lat975"
# ) { # open moveLocErrorCalc function
# 
#   # build reproject function for later use
#   reproject <- function(x,
#                         loncol = loncol,
#                         latcol = latcol,
#                         latloncrs = latloncrs,
#                         projectedcrs = projectedcrs) {
#     x <- sf::st_as_sf(x, coords = c(loncol, latcol)) |>
#       sf::st_set_crs(latloncrs) |> # latlon degrees sf object
#       sf::st_transform(projectedcrs) |> # eastings northings units metres
#       dplyr::select(-tidyselect::everything()) # remove all columns. Geometry is protected and retained
#     return(x)
#   }
# 
#   tracksfmean <- reproject(x = x,
#                            loncol = loncol,
#                            latcol = latcol,
#                            latloncrs = latloncrs,
#                            projectedcrs = projectedcrs)
# 
#   meanMoveLocDist <- list(
#     data.frame(loncol = x[,loncol], latcol = x[,lat975]), # U # were originally c(loncol, lat975), check format is right,
#     data.frame(loncol = x[,lon975], latcol = x[,latcol]), # R # this block should create a list of 4 dfs with 2 columns, instead has created 4 vectors of length 2n
#     data.frame(loncol = x[,loncol], latcol = x[,lat025]), # D # now created static names since this is only internal use
#     data.frame(loncol = x[,lon025], latcol = x[,latcol]) # L # allows static name call in reproject function below
#   ) |>
#     rlang::set_names(c("U", "R", "D", "L")) |> # set names of list elements
#     lapply(function(x) reproject(x = x,
#                                  loncol = "loncol", # x[1]
#                                  latcol = "latcol", # x[2]
#                                  latloncrs = latloncrs,
#                                  projectedcrs = projectedcrs
#     )) |>
#     lapply(
#       function(vertextrack) { # distance from vertices to centre
#         sf::st_distance(
#           x = tracksfmean,
#           y = vertextrack,
#           by_element = TRUE
#         )
#       }
#     ) |>
#     purrr::map_df(~.x) |> # collapse list to df of 4 columns
#     rowMeans() # make row means
# 
#   rm(tracksfmean)
#   # return(tracksfmean)
#   return(meanMoveLocDist)
# } # close moveLocErrorCalc function

# A4: Specify data and saveloc ----

YOUR_IP <- "NA" # add your IP address, server name or similar if you connect via a shared drive

## Input data
dataloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2025_NWA_hammerhead_complex","Data_input", "Telemetry")

## Output data
saveloc <- file.path("/", YOUR_IP,"Science","Projects_current", "2025_NWA_hammerhead_complex","Data_output", "Telemetry", "dBBMMs") # Adjust this

### ....................................................................................................
### [B] Data import ----
### ....................................................................................................

# B1: Data import ----

## movement data - filtered and standardised
hammers <- readRDS(file = file.path(dataloc, paste0("dBBMMs/Data_aniMotum_CRW_output_segmented_rerouted_proj_WGS84_converted_with_coord_CIs_with_Argosfilter_data.rds"))) |>
  mutate(shark = as.numeric(str_sub(id, start = 1, end = str_locate(id, "\\_")[,1] - 1)))

## metadata for species
meta <- read_csv(file.path(dataloc, "Data_AllSpecies_SPOT_tags_metadata.csv"))|>
  dplyr::filter( # keep species of movement data
    species %in% c("S.mokarran", "S.lewini","S.zygaena"),
    group %in% c("Florida Keys", "Jupiter", "Jupiter, FL", "Marquesas", "South Carolina", "Tampa")
  ) |>
  dplyr::select( # keep needed columns only
    ptt_id,
    datetime_deployment_local,
    stl,
    sex,
    species
  ) |>
  dplyr::rename(
    shark = ptt_id,
    FishLengthCm = stl)

## combine movement data with metadata
hammers %<>% left_join(meta)  # , by = join_by(shark == id) # doesn't work naming columns, has gotten worse.

## VERY IMPORTANT: The movement data needs to have time stamps in ascending order,
## double check that this is true by ordering the dataframe
hammers %<>%
  dplyr::arrange(
    shark,
    date
  ) %>%
  dplyr::select( # and since we are at it, the old id column can be removed
    -id
  )

### ....................................................................................................
### [C] Data housekeeping and preparation, location error, and timeDiffLong estimation ----
### ....................................................................................................

# C1: housekeeping movement dataframe ----

## define your grouping variables
### seasons by month numbers (based on definitions in Kroetz et al. 2021)
swinter <- c(12,1,2) # december, january, february
sspring <- c(3,4,5) # march, april, may
ssummer <- c(6,7,8) # june, july, august
sautumn <- c(9,10,11) # september, october, november

## Now combine detections file with grouping variables
hammers <- hammers %>%
  mutate(
    shark = make.names(shark), #Prefixes "X" to numerical-named sharks to avoid issues later
    month = as.numeric(format(date, format = "%m", tz = "UTC")) # get rid of year for later season variable definitions
    ) %>%
  # dplyr::rename(
  #   Datetime = date
  #   ) %>%
  # dplyr::select(
  #   !id
  # ) %>%
  mutate(# add a grouping variable for the season
    season = with(.,case_when(
      (month %in% swinter) ~ "winter",
      (month %in% sspring)~ "spring",
      (month %in% ssummer) ~ "summer",
      (month %in% sautumn) ~ "autumn",
      TRUE ~ "seasonnonexistent"
    )
    )
  )

# C3: calculate location error estimates

## set variables for reproject function
# loncol = "lon"
# latcol = "lat"
# lon025 = "lon025"
# lon975 = "lon975"
# lat025 = "lat025"
# lat975 = "lat975"
# latloncrs = 4326
# projectedcrs = 3395

# tracksfmean <- reproject(x = hammers,
#                          coorda = loncol,
#                          coordb = latcol,
#                          latloncrs = latloncrs,
#                          projectedcrs = projectedcrs)
# 
# meanMoveLocDist <- list(
#   c(loncol, lat975), # U
#   c(lon975, latcol), # R
#   c(loncol, lat025), # D
#   c(lon025, latcol) # L
# ) |>
#   lapply(function(x) reproject(x = hammers,
#                                coorda = x[1],
#                                coordb = x[2],
#                                latloncrs = latloncrs,
#                                projectedcrs = projectedcrs
#   )) |>
#   set_names(c("U", "R", "D", "L")) |> # set names of list elements
#   lapply(
#     function(vertextrack) { # distance from vertices to centre
#       st_distance(
#         x = tracksfmean,
#         y = vertextrack,
#         by_element = TRUE
#       )
#     }
#   ) |>
#   purrr::map_df(~.x) |> # collapse list to df of 4 columns
#   rowMeans()# make row means
# # (make overall mean of them?)
# hammers$meanMoveLocDist <- meanMoveLocDist

## alternative method with function directly from movegroup
hammers$meanMoveLocDist <- moveLocErrorCalc(hammers,
                                             loncol = "lon",
                                             latcol = "lat",
                                             latloncrs = 4326,
                                             projectedcrs = 3395,
                                             lon025 = "lon025",
                                             lon975 = "lon975",
                                             lat025 = "lat025",
                                             lat975 = "lat975"
                                             )

# C3: calculatetimeDiffLong ----

hammers$diffmins <- c(as.numeric(NA), as.numeric(hammers$date[2:length(hammers$date)] - hammers$date[1:length(hammers$date) - 1]))
# gives error but still works
# get index of first row per id
firstrows <- hammers |>
  mutate(Index = 1:nrow(hammers)) |>
  group_by(shark) |>
  summarise(firstrowid = first(Index))
# use this to blank out the diffmins since it's the time difference between 1 shark ending & another staring which is meaningless
hammers[firstrows$firstrowid, "diffmins"] <- NA

## works, doesn't need tee pipe. Uses braces to evaluate that bit immediately
# hammers %>%
#   group_by(shark) %>%
#   group_walk(~ {
#     p <- ggplot(.x) + geom_histogram(aes(x = diffmins))
#     # ggsave(plot = p, filename = paste0(saveloc, "movegroup dBBMMs/timeDiffLong histograms/", lubridate::today(), "_diffMinsHistGG_", .y$id, ".png")) #SD
#     ggsave(plot = p, filename = file.path(saveloc, paste0("/timeDiffLong_histograms/", lubridate::today(), "_diffMinsHistGG_", .y$shark, ".png"))
# 
#     .x
#   }) %>%
#   summarise(meandiffmins = mean(diffmins, na.rm = TRUE),
#             sd3 = sd(diffmins, na.rm = TRUE) * 3)

### ..........................................................................................
### [D] Prepare raster CRS, define dBBMM parameters ----
### ....................................................................................................

# D1: define rasterResolution ----
# With default = 6: Error: cannot allocate vector of size 2286.7 Gb
2 * mean(hammers$meanMoveLocDist) # 111255.6
hist(hammers$meanMoveLocDist) # very left skewed
summary(hammers$meanMoveLocDist) # median 16808
# choose 10000

# D3: find min and max values for rasterCRS ----
# mean(hammers$lon, na.rm = F) # -80.45663
# mean(hammers$lat, na.rm = F) # 28.49443
# min(hammers$lon, na.rm = F) # -97.79646
# max(hammers$lon, na.rm = F) # -67.83846
# min(hammers$lat, na.rm = F) # 22.35267
# max(hammers$lat, na.rm = F) # 41.31454

# https://tmackinnon.com/2005/images/utmworld.gif
# Mean point is zone 17
# rasterCRS = CRS("+proj=utm +zone=17 +datum=WGS84")

### ....................................................................................................
### [E] Prefilter data for species, region and season ----
### ....................................................................................................

# We want to calculate dBBMMs for all species (n = 3) and all seasons (n = 4). We can create subset
# combinations and loop over them

# E1: prepare subsets and folders ----

## First we need to create a new column that contains the species without the "." as we will create directories
## and a "." might create issues in filenames
hammers$spp <- gsub("\\.", "", hammers$species) # create a character string to name files when saving, "." might cause issues downstream

# ## create directories for your species - not needed as dir.create is part of movegroup::movegroup
# spps <- sort(unique(hammers$spp))
# 
# for (i in spps){
#   dir.create(file.path(saveloc, i))
# }

# E2: define your subsets ----

# CHANGE THIS ONLY ...........................................................................
s = 2 # select for species: 1 = Slewini, 2 = Smokarran, 3 = Szygaena
t = 1 # select for season: 1 = autumn, 2 = spring, 3 = summer, 4 = winter
# ............................................................................................

# E2.sp: species selection - DONT change this ----

sp <- sort(unique(hammers$spp))
# > sp
# [1] "Slewini"   "Smokarran" "Szygaena" 

# Select species
sp.f <- sp[s]

# Filter for species
hammersubset <- filter(hammers, spp == sp.f)

# E2.tp: Season selection - DONT change this ----

tp <- sort(unique(hammers$season))
# > tp
# [1] "autumn" "spring" "summer" "winter"

# Select season
tp.f <- tp[t]

# Filter for studyperiod
hammersubset <- filter(hammersubset, season == tp.f)

## check how many relocations
table(hammersubset$shark)
check <- hammersubset[which(hammersubset$lat975 > 50),]
if(nrow(check) == 0){
  hammersubset <- hammersubset
} else {
  hammersubset <- hammersubset[which(hammersubset$lat975 <= 50),]
}

table(hammersubset$shark)
# 
# hammersubset <- hammersubset %>%
#   dplyr::filter(
#     !meanMoveLocDist > 200000
#   )

## check which subset you are calculating
print(paste0("Now exploring dBBMMs for ", sp.f, " during ", tp.f,"."))

### ....................................................................................................
### [F] Calculate dBBMMs using movegroup package ----
### ....................................................................................................

# F1: Construct individual-level dBBMMs ----

movegroup::movegroup(
  data = hammersubset,
  ID = "shark",
  Datetime = "date", # DateTime
  Lat = "lat",
  Lon = "lon",
  # Group = NULL,
  dat.TZ = "UTC",
  # proj = sp::CRS("+proj=longlat +datum=WGS84"),
  # projectedCRS = "+init=epsg:32617", # https://epsg.io/32617 Bimini, Florida, ERROR "unused argument" commented out 20240827
  # sensor = "VR2W",
  moveLocError = hammersubset$meanMoveLocDist,
  timeDiffLong = 13, #original (TDL * 60), #adjust
  # Single numeric value. Threshold value in timeDiffUnits designating the length of long breaks in re-locations. Used for bursting a movement track into segments, thereby removing long breaks from the movement track. See ?move::bursted for details.
  timeDiffUnits = "hours",# original: "mins",
  # center = TRUE,
  buffpct = 10, #0.6, # Buffer extent for raster creation, proportion of 1.
  # rasterExtent = NULL,
  rasterCRS = sp::CRS("+proj=merc +datum=WGS84"),
  rasterResolution = 10000, # changed from 1000 to 10000 on 20240827
  # Single numeric value to set raster resolution - cell size in metres? 111000: 1 degree lat = 111km.
  # Tradeoff between small res = big file & processing time.
  # Should be a function of the spatial resolution of your receivers or positioning tags.
  # Higher resolution will lead to more precision in the volume areas calculations.
  # Try using 2*dbblocationerror.
  ##### Why did we choose 1000m?####
  dbbext = .3, # Ext param in the 'brownian.bridge.dyn' function in the 'move' package. Extends bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax). Default 0.3. - changed to 3 20240827
  # dbbwindowsize = 23,
  # writeRasterFormat = "ascii",
  # writeRasterExtension = ".asc",
  # writeRasterDatatype = "FLT4S",
  # absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
  # savedir = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/"), #SD
  savedir = file.path(saveloc, sp.f, tp.f), #VH
  alerts = TRUE
)

# F2: scale individual-level UDs/rasters to group-level UDs/rasters ----

## WARNING: MAKE SURE ALL INDIVIDUAL-LEVEL UDs OF INTEREST HAVE BEEN CREATED FIRST BEFORE RUNNING SCALERASTER()
movegroup::scaleraster(
  path = file.path(saveloc, sp.f, tp.f), 
  pathsubsets = file.path(saveloc), 
  # pattern = ".asc",
  # weighting = 1,
  # format = "ascii",
  # datatype = "FLT4S",
  # bylayer = TRUE,
  # overwrite = TRUE,
  # scalefolder = "Scaled",
  # weightedsummedname = "All_Rasters_Weighted_Summed",
  # scaledweightedname = "All_Rasters_Scaled_Weighted",
  crsloc = file.path(saveloc, sp.f, tp.f)
)

# F3: plot group-level UDs using movegroup for quick pattern exploration ----

movegroup::plotraster(
  x = file.path(saveloc, sp.f, tp.f, "Scaled","All_Rasters_Scaled_Weighted_UDScaled.asc"),
  crsloc = file.path(saveloc, sp.f, tp.f),
  trim = TRUE,
  myLocation = c(min(hammers$lon), min(hammers$lat), max(hammers$lon), max(hammers$lat)),
  # myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax) - default NULL, extents autocreated from data
  googlemap = TRUE,
  # gmapsAPI = NULL, # enter your Google maps API here, quoted character string.
  stadiaAPI = "6b11fbd4-2f7c-4f27-b2c6-2d0d4af6d849", # enter your Stadia maps API here, quoted character string.
  expandfactor = 1,
  mapzoom = 7,
  mapsource = "stadia",
  maptype = "stamen_toner_lite",
  contour1colour = "orange",
  contour2colour = "red",
  # plottitle = "Aggregated 95% and 50% UD contours",
  plotsubtitle = paste0("Scaled contours ", sp.f, " ", tp.f),
  # legendtitle = "Percent UD Contours",
  # plotcaption = paste0("movegroup, ", lubridate::today()),
  # axisxlabel = "Longitude",
  # axisylabel = "Latitude",
  # legendposition = c(0.9, 0.89),
  legendposition = c(1.12, .75),
  fontsize = 12,
  # fontfamily = "Times New Roman",
  # filesavename = paste0(lubridate::today(), "_dBBMM-contours.png"),
  savedir = file.path(saveloc, sp.f, tp.f,"Scaled", "Plot"),
  # receiverlats = NULL,
  # receiverlons = NULL,
  # receivernames = NULL,
  # receiverrange = NULL,
  # recpointscol = "black",
  # recpointsfill = "white",
  # recpointsalpha = 0.5,
  # recpointssize = 1,
  # recpointsshape = 21,
  # recbufcol = "grey75",
  # recbuffill = "grey",
  # recbufalpha = 0.5,
  # reclabcol = "black",
  # reclabfill = NA,
  # reclabnudgex = 0,
  # reclabnudgey = -200,
  # reclabpad = 0,
  # reclabrad = 0.15,
  # reclabbord = 0,
  surface = TRUE
)

# END OF CODE ----

# TODO:
# TODO1: automate entire movegroup-scaleraster-plotraster section so for all species and seasons











# F2: seasonal dBBMMs ----

# F2.1: Construct individual-level dBBMMs ----

##dBBMMhomeRange(
movegroup::movegroup(
  data = hammers.i, # data frame of data needs columns Lat Lon DateTime and optionally an ID and grouping columns.
  ID = "shark", # column name of IDs of individuals.
  Datetime = "date", # name of Datetime column. Must be in POSIXct format.
  Lat = "lat", # name of Lat & Lon columns in data.
  Lon = "lon",
  #Group = NULL, # name of grouping column in data. CURRENTLY UNUSED; MAKE USER DO THIS?
  #dat.TZ = "Etc/GMT+5", # timezone for as.POSIXct.
  #dat.TZ = "UTC", # timezone for as.POSIXct.
  #proj = sp::CRS("+proj=longlat +datum=WGS84"), # CRS for move function.
  projectedCRS = "+init=epsg:32617", # 32617 EPSG code for CRS for initial transform of latlon points; corresponds to rasterCRS zone
  #sensor = "unknown", # sensor for move function. Single character or vector with length of the number of coordinates. Optional.
  moveLocError = hammers.i$meandist, # location error in metres for move function. Numeric. Either single or a vector of lenth nrow data.
  #moveLocError = DET.i$meandist, # location error in metres for move function. Numeric. Either single or a vector of lenth nrow data.
  timeDiffLong = 13, #24*100, # threshold length of time in timeDiffUnits designating long breaks in relocations.
  timeDiffUnits = "hours", # units for time difference for move function.
  #center = TRUE, # center move object within extent? See spTransform.
  buffpct = .6, #3, # buffer extent for raster creation, proportion of 1.
  #rasterExtent = NULL, # if NULL, raster extent calculated from data, buffpct, rasterResolution. Else length 4 vector, c(xmn, xmx, ymn, ymx) decimal latlon degrees. Don't go to 90 for ymax
  # Doesn't prevent constraint to data limits (in plot anyway), but prevents raster clipping crash
  #rasterCRS = sp::CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"), # CRS for raster creation. 17
  rasterResolution = 1000, # 1 degree lat = 111km. Could maybe make a tad smaller. # numeric vector of length 1 or 2 to set raster resolution - cell size in metres? 111000: 1 degree lat = 111km
  #movemargin = 11,
  ##dbblocationerror = DET.i$meandist, # location.error param in brownian.bridge.dyn. Could use the same as moveLocError?
  dbbext = .3, # ext param in brownian.bridge.dyn. Extends bounding box around track. Numeric single (all edges), double (x & y), or 4 (xmin xmax ymin ymax). Default 0.3,
  dbbwindowsize = 23, # window.size param in brownian.bridge.dyn. The size of the moving window along the track. Larger windows provide more stable/accurate estimates of the brownian motion variance but are less well able to capture more frequent changes in behavior. This number has to be odd. A dBBMM is not run if total detections of individual < window size (default 31).
  writeRasterFormat = "ascii",
  writeRasterExtension = ".asc",
  writeRasterDatatype = "FLT4S",
  absVolumeAreaSaveName = "VolumeArea_AbsoluteScale.csv",
  savedir = paste0(saveloc, "seasonal/", seas.f, "/"),  # save outputs to a temporary directory (default) else.
  #savedir = paste0(saveloc, season.f,"/"),  # save outputs to a temporary directory (default) else.
  alerts = TRUE # audio warning for failures
)
beep(4)

# F2.2: Scale the rasters from 0 to 1, and weigh based on regions ----

## WARNING: MAKE SURE ALL INDIVIDUAL-LEVEL UDs OF INTEREST HAVE BEEN CREATED FIRST BEFORE RUNNING SCALERASTER()
# dBBMMhomeRange::scaleraster(
movegroup::scaleraster(
  path = paste0(saveloc, "seasonal/", seas.f, "/"), # same as savedir in dBBMM.build. Contains rasters to be scaled.
  pathsubsets = paste0(saveloc, "seasonal/"), # 'Parental folder': Location of ALL files created by dBBMM.build. No terminal slash.
  pattern = ".asc",
  #weighting = 1,
  # weighting to divide individual and summed-scaled rasters by, for unbalanced arrays
  format = "ascii",
  datatype = "FLT4S",
  bylayer = TRUE,
  overwrite = TRUE,
  scalefolder = "Scaled",
  #weightedsummedname = "All_Rasters_Weighted_Summed",
  scaledweightedname = "All_Rasters_Scaled_Weighted",
  crsloc = paste0(saveloc, "seasonal/", seas.f,"/"),
  # location of saved CRS Rds file from dBBMM.build.R. Should be same as path.
  returnObj = FALSE
)
beep(4)

### ....................................................................................................
### [G] Plot the contours ----
### ....................................................................................................

# G1: all year dBBMM contours ----

dir.create(paste0(saveloc, "allyear/Scaled/Plot"))

# ggmap::register_google()
movegroup::plotraster(
  # x = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/Scaled/All_Rasters_Scaled_Weighted_UDScaled.asc"), #SD
  # crsloc = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/"), #SD
  x = paste0(saveloc, "allyear/","Scaled/All_Rasters_Scaled_Weighted_UDScaled.asc"), #VH
  crsloc = paste0(saveloc, "allyear/"), #VH
  # trim = TRUE,
  # myLocation = NULL,
  googlemap = TRUE,
  # gmapsAPI = NULL,
  expandfactor = 1,
  mapzoom = 7,
  # mapsource = "google",
  # maptype = "satellite",
  # contour1colour = "red",
  # contour2colour = "orange",
  # plottitle = "Aggregated 95% and 50% UD contours",
  plotsubtitle = "Scaled contours. n = 10",
  # legendtitle = "Percent UD Contours",
  # plotcaption = paste0("movegroup, ", lubridate::today()),
  # axisxlabel = "Longitude",
  # axisylabel = "Latitude",
  legendposition = c(0.9, 0.89),
  # fontsize = 12,
  # fontfamily = "Times New Roman",
  # filesavename = paste0(lubridate::today(), "_dBBMM-contours.png"),
  # savedir = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/Scaled/Plot"), #SD
  savedir = paste0(saveloc, "allyear/Scaled/Plot/"), #VH
  # receiverlats = NULL,
  # receiverlons = NULL,
  # receivernames = NULL,
  # receiverrange = NULL,
  # recpointscol = "black",
  # recpointsfill = "white",
  # recpointsalpha = 0.5,
  # recpointssize = 1,
  # recpointsshape = 21,
  # recbufcol = "grey75",
  # recbuffill = "grey",
  # recbufalpha = 0.5,
  # reclabcol = "black",
  # reclabfill = NA,
  # reclabnudgex = 0,
  # reclabnudgey = -200,
  # reclabpad = 0,
  # reclabrad = 0.15,
  # reclabbord = 0,
  surface = TRUE
)

# G2: seasonal dBBMMs ----

dir.create(paste0(saveloc, "seasonal/summer/Scaled/Plot")); dir.create(paste0(saveloc, "seasonal/winter/Scaled/Plot"))

## group-level UDs scaled
for (thisseason in unique(season)) {
  ##### ISSUE####
  # movegroup said: "processing 7 of 7" i.e. not 8 of 8
  plotraster(
    # x = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/", thisshark, ".asc"), #SD
    # crsloc = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/"), #SD
    x = paste0(saveloc, "seasonal/", thisseason, "/Scaled/All_Rasters_Scaled_Weighted_UDScaled.asc"), #VH
    crsloc = paste0(saveloc, "seasonal/", thisseason, "/"), #VH
    googlemap = TRUE,
    expandfactor = 1,
    mapzoom = 7,
    plotsubtitle = "Scaled contours. n = 10",
    legendposition = c(0.9, 0.89),
    filesavename = paste0(lubridate::today(), "_", thisseason, "_dBBMM-contours.png"),
    # savedir = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/Scaled/Plot"))} #SD
    savedir = paste0(saveloc, "seasonal/", thisseason, "/Scaled/Plot"))
  } #VH

## individual-level UDs scaled - winter
for (thisshark in unique(hammer.i$shark)) {
  ##### ISSUE####
  # movegroup said: "processing 7 of 7" i.e. not 8 of 8
  plotraster(
    # x = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/", thisshark, ".asc"), #SD
    # crsloc = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/"), #SD
    x = paste0(saveloc, "seasonal/winter/Scaled/", thisshark, ".asc"), #VH
    crsloc = paste0(saveloc, "seasonal/winter/"), #VH
    googlemap = TRUE,
    expandfactor = 1,
    mapzoom = 7,
    plotsubtitle = "Scaled contours",
    legendposition = c(0.9, 0.89),
    filesavename = paste0(lubridate::today(), "_", thisshark, "_dBBMM-contours.png"),
    # savedir = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/Scaled/Plot"))} #SD
    savedir = paste0(saveloc, "seasonal/winter/Scaled/Plot"))
  } #VH

## inbdividual-level UDs scaled - winter
for (thisshark in unique(hammer.i$shark)) {
  ##### ISSUE####
  # movegroup said: "processing 7 of 7" i.e. not 8 of 8
  plotraster(
    # x = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/", thisshark, ".asc"), #SD
    # crsloc = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/"), #SD
    x = paste0(saveloc, "seasonal/summer/Scaled/", thisshark, ".asc"), #VH
    crsloc = paste0(saveloc, "seasonal/summer/"), #VH
    googlemap = TRUE,
    expandfactor = 1,
    mapzoom = 7,
    plotsubtitle = "Scaled contours",
    legendposition = c(0.9, 0.89),
    filesavename = paste0(lubridate::today(), "_", thisshark, "_dBBMM-contours.png"),
    # savedir = paste0(saveloc, "movegroup dBBMMs/", thissubset, "/", TDL, "h/Scaled/Plot"))} #SD
    savedir = paste0(saveloc, "seasonal/summer/Scaled/Plot"))
} #VH

## attempts at automation

# Automated Hammerhead Shark Movement Analysis
# This script runs through all species and seasons systematically

# Get unique species and seasons from the data
sp <- sort(unique(hammers$spp))
tp <- sort(unique(hammers$season))

print(paste("Species found:", paste(sp, collapse=", ")))
print(paste("Seasons found:", paste(tp, collapse=", ")))

# =============================================================================
# PHASE 1: Run movegroup() for ALL species and seasons
# =============================================================================
print("=== PHASE 1: Running movegroup() for all combinations ===")

for(s in 1:length(sp)) {
  sp.f <- sp[s]
  print(paste("Processing species:", sp.f, "(", s, "of", length(sp), ")"))
  
  for(t in 1:length(tp)) {
    tp.f <- tp[t]
    print(paste("  Processing season:", tp.f, "(", t, "of", length(tp), ")"))
    
    # Filter for species
    hammersubset <- filter(hammers, spp == sp.f)
    
    # Filter for season
    hammersubset <- filter(hammersubset, season == tp.f)
    
    # Check if we have data for this combination
    if(nrow(hammersubset) == 0) {
      print(paste("    No data found for", sp.f, "in", tp.f, "- skipping"))
      next
    }
    
    # Check relocations
    print(paste("    Relocations before filtering:"))
    print(table(hammersubset$shark))
    
    # Remove problematic latitude values
    check <- hammersubset[which(hammersubset$lat975 > 50),]
    if(nrow(check) == 0){
      hammersubset <- hammersubset
    } else {
      hammersubset <- hammersubset[which(hammersubset$lat975 <= 50),]
    }
    
    print(paste("    Relocations after filtering:"))
    print(table(hammersubset$shark))
    
    # Check if we still have data after filtering
    if(nrow(hammersubset) == 0) {
      print(paste("    No data remaining after filtering for", sp.f, "in", tp.f, "- skipping"))
      next
    }
    
    print(paste0("    Now calculating dBBMMs for ", sp.f, " during ", tp.f))
    
    # Run movegroup analysis
    tryCatch({
      movegroup::movegroup(
        data = hammersubset,
        ID = "shark",
        Datetime = "date",
        Lat = "lat",
        Lon = "lon",
        dat.TZ = "UTC",
        moveLocError = hammersubset$meanMoveLocDist,
        timeDiffLong = 13,
        timeDiffUnits = "hours",
        buffpct = 10,
        rasterCRS = sp::CRS("+proj=merc +datum=WGS84"),
        rasterResolution = 10000,
        dbbext = .3,
        savedir = file.path(saveloc, sp.f, tp.f),
        alerts = TRUE
      )
      print(paste("    ✓ Successfully completed movegroup for", sp.f, tp.f))
    }, error = function(e) {
      print(paste("    ✗ Error in movegroup for", sp.f, tp.f, ":", e$message))
    })
  }
}

print("=== PHASE 1 COMPLETE: All movegroup() analyses finished ===")

# =============================================================================
# PHASE 2: Run scaleraster() for ALL species and seasons
# =============================================================================
print("=== PHASE 2: Running scaleraster() for all combinations ===")

for(s in 1:length(sp)) {
  sp.f <- sp[s]
  print(paste("Scaling rasters for species:", sp.f, "(", s, "of", length(sp), ")"))
  
  for(t in 1:length(tp)) {
    tp.f <- tp[t]
    print(paste("  Scaling season:", tp.f, "(", t, "of", length(tp), ")"))
    
    # Check if the directory exists (i.e., movegroup was successful)
    raster_path <- file.path(saveloc, sp.f, tp.f)
    if(!dir.exists(raster_path)) {
      print(paste("    Directory not found for", sp.f, tp.f, "- skipping scaling"))
      next
    }
    
    tryCatch({
      movegroup::scaleraster(
        path = file.path(saveloc, sp.f, tp.f), 
        pathsubsets = file.path(saveloc), 
        crsloc = file.path(saveloc, sp.f, tp.f)
      )
      print(paste("    ✓ Successfully scaled rasters for", sp.f, tp.f))
    }, error = function(e) {
      print(paste("    ✗ Error in scaleraster for", sp.f, tp.f, ":", e$message))
    })
  }
}

print("=== PHASE 2 COMPLETE: All scaleraster() analyses finished ===")

# =============================================================================
# PHASE 3: Run plotraster() for ALL species and seasons
# =============================================================================
print("=== PHASE 3: Running plotraster() for all combinations ===")

for(s in 1:length(sp)) {
  sp.f <- sp[s]
  print(paste("Plotting rasters for species:", sp.f, "(", s, "of", length(sp), ")"))
  
  for(t in 1:length(tp)) {
    tp.f <- tp[t]
    print(paste("  Plotting season:", tp.f, "(", t, "of", length(tp), ")"))
    
    # Check if the directory exists (i.e., previous steps were successful)
    raster_path <- file.path(saveloc, sp.f, tp.f)
    if(!dir.exists(raster_path)) {
      print(paste("    Directory not found for", sp.f, tp.f, "- skipping plotting"))
      next
    }
    
    tryCatch({
      movegroup::plotraster(
        path = file.path(saveloc, sp.f, tp.f), 
        pathsubsets = file.path(saveloc), 
        crsloc = file.path(saveloc, sp.f, tp.f)
      )
      print(paste("    ✓ Successfully plotted rasters for", sp.f, tp.f))
    }, error = function(e) {
      print(paste("    ✗ Error in plotraster for", sp.f, tp.f, ":", e$message))
    })
  }
}

print("=== PHASE 3 COMPLETE: All plotraster() analyses finished ===")
print("=== ALL ANALYSES COMPLETE ===")

# Optional: Summary of what was processed
print("Summary:")
print(paste("Total species processed:", length(sp)))
print(paste("Total seasons processed:", length(tp)))
print(paste("Total combinations processed:", length(sp) * length(tp)))


## attempts at automation 2:
# Automated Hammerhead Shark Movement Analysis
# ONE-CLICK SOLUTION: Run this entire script once and it will:
# 1. Run movegroup() for all species/season combinations
# 2. Automatically run scaleraster() for all combinations  
# 3. Automatically run plotraster() for all combinations
# NO MANUAL INTERVENTION REQUIRED!

# Get unique species and seasons from the data
sp <- sort(unique(hammers$spp))
tp <- sort(unique(hammers$season))

print(paste("Species found:", paste(sp, collapse=", ")))
print(paste("Seasons found:", paste(tp, collapse=", ")))

# =============================================================================
# PHASE 1: Run movegroup() for ALL species and seasons
# =============================================================================
print("=== PHASE 1: Running movegroup() for all combinations ===")

for(s in 1:length(sp)) {
  sp.f <- sp[s]
  print(paste("Processing species:", sp.f, "(", s, "of", length(sp), ")"))
  
  for(t in 1:length(tp)) {
    tp.f <- tp[t]
    print(paste("  Processing season:", tp.f, "(", t, "of", length(tp), ")"))
    
    # Filter for species
    hammersubset <- filter(hammers, spp == sp.f)
    
    # Filter for season
    hammersubset <- filter(hammersubset, season == tp.f)
    
    # Check if we have data for this combination
    if(nrow(hammersubset) == 0) {
      print(paste("    No data found for", sp.f, "in", tp.f, "- skipping"))
      next
    }
    
    # Check relocations
    print(paste("    Relocations before filtering:"))
    print(table(hammersubset$shark))
    
    # Remove problematic latitude values
    check <- hammersubset[which(hammersubset$lat975 > 50),]
    if(nrow(check) == 0){
      hammersubset <- hammersubset
    } else {
      hammersubset <- hammersubset[which(hammersubset$lat975 <= 50),]
    }
    
    print(paste("    Relocations after filtering:"))
    print(table(hammersubset$shark))
    
    # Check if we still have data after filtering
    if(nrow(hammersubset) == 0) {
      print(paste("    No data remaining after filtering for", sp.f, "in", tp.f, "- skipping"))
      next
    }
    
    print(paste0("    Now calculating dBBMMs for ", sp.f, " during ", tp.f))
    
    # Run movegroup analysis
    tryCatch({
      movegroup::movegroup(
        data = hammersubset,
        ID = "shark",
        Datetime = "date",
        Lat = "lat",
        Lon = "lon",
        dat.TZ = "UTC",
        moveLocError = hammersubset$meanMoveLocDist,
        timeDiffLong = 13,
        timeDiffUnits = "hours",
        buffpct = 10,
        rasterCRS = sp::CRS("+proj=merc +datum=WGS84"),
        rasterResolution = 10000,
        dbbext = .3,
        savedir = file.path(saveloc, sp.f, tp.f),
        alerts = TRUE
      )
      print(paste("    ✓ Successfully completed movegroup for", sp.f, tp.f))
    }, error = function(e) {
      print(paste("    ✗ Error in movegroup for", sp.f, tp.f, ":", e$message))
    })
  }
}

print("=== PHASE 1 COMPLETE: All movegroup() analyses finished ===")

# =============================================================================
# PHASE 2: Run scaleraster() for ALL species and seasons
# =============================================================================
print("=== PHASE 2: Running scaleraster() for all combinations ===")

for(s in 1:length(sp)) {
  sp.f <- sp[s]
  print(paste("Scaling rasters for species:", sp.f, "(", s, "of", length(sp), ")"))
  
  for(t in 1:length(tp)) {
    tp.f <- tp[t]
    print(paste("  Scaling season:", tp.f, "(", t, "of", length(tp), ")"))
    
    # Check if the directory exists (i.e., movegroup was successful)
    raster_path <- file.path(saveloc, sp.f, tp.f)
    if(!dir.exists(raster_path)) {
      print(paste("    Directory not found for", sp.f, tp.f, "- skipping scaling"))
      next
    }
    
    tryCatch({
      movegroup::scaleraster(
        path = file.path(saveloc, sp.f, tp.f), 
        pathsubsets = file.path(saveloc), 
        crsloc = file.path(saveloc, sp.f, tp.f)
      )
      print(paste("    ✓ Successfully scaled rasters for", sp.f, tp.f))
    }, error = function(e) {
      print(paste("    ✗ Error in scaleraster for", sp.f, tp.f, ":", e$message))
    })
  }
}

print("=== PHASE 2 COMPLETE: All scaleraster() analyses finished ===")

# =============================================================================
# PHASE 3: Run plotraster() for ALL species and seasons
# =============================================================================
print("=== PHASE 3: Running plotraster() for all combinations ===")

for(s in 1:length(sp)) {
  sp.f <- sp[s]
  print(paste("Plotting rasters for species:", sp.f, "(", s, "of", length(sp), ")"))
  
  for(t in 1:length(tp)) {
    tp.f <- tp[t]
    print(paste("  Plotting season:", tp.f, "(", t, "of", length(tp), ")"))
    
    # Check if the directory exists (i.e., previous steps were successful)
    raster_path <- file.path(saveloc, sp.f, tp.f)
    if(!dir.exists(raster_path)) {
      print(paste("    Directory not found for", sp.f, tp.f, "- skipping plotting"))
      next
    }
    
    tryCatch({
      movegroup::plotraster(
        path = file.path(saveloc, sp.f, tp.f), 
        pathsubsets = file.path(saveloc), 
        crsloc = file.path(saveloc, sp.f, tp.f)
      )
      print(paste("    ✓ Successfully plotted rasters for", sp.f, tp.f))
    }, error = function(e) {
      print(paste("    ✗ Error in plotraster for", sp.f, tp.f, ":", e$message))
    })
  }
}

print("=== PHASE 3 COMPLETE: All plotraster() analyses finished ===")
print("=== ALL ANALYSES COMPLETE ===")

# Optional: Summary of what was processed
print("Summary:")
print(paste("Total species processed:", length(sp)))
print(paste("Total seasons processed:", length(tp)))
print(paste("Total combinations processed:", length(sp) * length(tp)))