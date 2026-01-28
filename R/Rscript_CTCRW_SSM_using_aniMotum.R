# TOP OF CODE ----
### ====================================================================================================
### Project:    Large bodied hammerhead complex in the western North Atlantic
### Analysis:   Processing and cleaning satellite telemetry data of fin-mounted SPOT tags for further steps
### Script:     ~2025_NWA_hammerhead_complex/R/Rscript2_CTCRW_SSM_SPOT_using_aniMotum.R
### Author:     Vital Heim
### Version:    1.0
### ====================================================================================================

### ....................................................................................................
### Content: this R script contains the code to filter, clean and process raw data collected by fin-
###          mounted Smart Position and Temperature (SPOT) transmitters.
###          This script contains the code fit a continuous time correlated random walk model (CTCRW
###          in a state space model (SSM) framework using sda-prefiltered data.
###          The code is based on the functions provided within the animotum (Jonsen et al. 2023) package:
###          > citation("aniMotum")
###          Ian Jonsen, W James Grecian, Lachlan Phillips, Gemma Carroll, Clive R. McMahon,
###          Robert G. Harcourt, Mark A. Hindell, and Toby A. Patterson (2023) aniMotum, an R
###          package for animal movement data: Rapid quality control, behavioural estimation
###          and simulation.  Methods in Ecology and Evolution DOI: 10.1111/2041-210X.14060
###          The script contains additional code for pre-fitting processing and cleaning of the data.
### ....................................................................................................

### ....................................................................................................
### [A] Setwd, paths and parameters ----
### ....................................................................................................

# A1: clear memory ----

rm(list = ls())

# A2: load necessary packages ----

# library("remotes")
## install if first time
# install.packages("tidyverse")
# install.packages("magrittr")
# install.packages("TMB", type = "source") # if package version inconsistency detected during loading of foieGras
# install.packages("aniMotum",
                    # repos = c("https://cloud.r-project.org",
                    # "https://ianjonsen.r-universe.dev"),
                    # dependencies = TRUE) #foiegras was removed from CRAN and replaced with Animotum on 12-12-2022
# install.packages("patchwork")
# install.packages("sf")
# install.packages("sp")
# install.packages("rnaturalearth")
# install.packages("rnaturalearthdata")
# remotes::install_github("ropensci/rnaturalearthhires")
# install.packages("ggplot2")
# install.packages("xts")
# install.packages("trip")
# install.packages("ggspatial")
# install.packages("virids")

## load
library(tidyverse)
library(magrittr)
library(aniMotum)
library(patchwork)
library(sf)
library(sp)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggplot2)
# library(xts)
# library(trip)
library(ggspatial)
library(viridis)

# A3: Specify needed functions

## NA

# A4: Specify data and saveloc ----

YOUR_IP <- "NA" # add your IP address, server name or similar if you connect via a shared drive

## Input data
dataloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2025_NWA_hammerhead_complex","Data_input", "Telemetry")

## Output data
saveloc <- file.path("/", YOUR_IP,"Science","Projects_current", "2025_NWA_hammerhead_complex","Data_output", "Telemetry", "CTCRW") # Adjust this

# A5: Define universal variables (e.g. for plotting) ----

# NA

# A6: Define universal options, variables, etc. (e.g. for plotting) ----

options(scipen=999) # so that R doesn't act up for pit numbers  

### ....................................................................................................
### [B] SSM with prefiltered data ----
### ....................................................................................................

# Remember: if you use prefiltered data, this means raw argos data was filtered to remove
# spurious detections based on  a speed-distance-angle algorithm. Based on the script/file either
# the sda filter from the package "argosfilter" or "trip" was used.
# In either case, the files were already adjusted for near duplicated observations so this does not
# need to be done again here in section B.

# Define first which sda filter you used

speedfilter <- "Argosfilter" # choices: "Argosfilter" or "TripSDA"

# B1: Import the filtered data ----

mydets_f <- list.files(path = dataloc,
                       pattern = paste0(speedfilter,"_.+.R"),
                       full.names = TRUE ) %>%
  purrr::map_dfr(readRDS) # make sure that there is only the most up to date speedfilter file in the directory

## check the data
# print(n = 1000, mydets_f %>% group_by(id) %>% dplyr::summarise(n = n()))
# A tibble: 69 × 2
# id         n
# <chr>  <int>
#   1 174515   129
# 2 175427   110
# 3 179468   798
# 4 179471   250
# 5 179472    38
# 6 180910   444
# 7 180912   876
# 8 180914   748
# 9 182746    36
# 10 183619  1048
# 11 183620   493
# 12 183621   537
# 13 183622   137
# 14 183624    72
# 15 198201   736
# 16 198202   481
# 17 198203   218
# 18 198204   217
# 19 198205   238
# 20 198206   299
# 21 210602   311
# 22 210603   789
# 23 210604    18
# 24 210605   264
# 25 210606   739
# 26 222134   258
# 27 222136    51
# 28 222137   877
# 29 222138   263
# 30 233536   724
# 31 233537   982
# 32 233538  1120
# 33 235278   466
# 34 235281   782
# 35 244597   687
# 36 244598   692
# 37 244600   616
# 38 244601  1197
# 39 244602   208
# 40 244603   103
# 41 244605  1451
# 42 244606  1991
# 43 244609   361
# 44 244610   388
# 45 244611   337
# 46 261302  2107
# 47 261303  1495
# 48 261304  1151
# 49 261305   619
# 50 261306   318
# 51 261309  1014
# 52 261742    82
# 53 264015  1609
# 54 264019   147
# 55 264021   356
# 56 264022   254
# 57 264023   261
# 58 264024   239
# 59 264025   468
# 60 264026   557
# 61 286887   269
# 62 286888   199
# 63 286889   280
# 64 286890   108
# 65 286895   104
# 66 286896   102
# 67 286897   193
# 68 286898    26
# 69 286899   105

# B2: housekeeping, tidy data and define column classes ----

## the fit_ssm and functions in the animotum package expect a data.frame or tibble in
## the following format if (!) error elipse information is available. If the data was collected
## using other methods than the CLS Argos' Kalman filter model other formats might be needed.
## Check: https://ianjonsen.github.io/aniMotum/articles/Overview.html

## We need:
## 'id'
## 'date'
## 'loc' -location class
## 'lon'
## 'lat'
## 'smaj' - error semi-major axis
## 'smin' - error Semi-minor axis
## 'eor' - error ellipse orientation

## Basic housekeeping
det_f <-mydets_f %>%
  dplyr::select( # select relevant columns, here: id, date, location class (lc), lon, lat,
    id,
    date,
    lc,
    lon,
    lat,
    smaj,
    smin,
    eor
  )

## check if there are Z locations left (should not!)
print(paste0("There are ", length(which(det_f$lc == "Z")), " Z-locations left."))

# B4: Known/tagging locations ----

## we already added tagging locations in a previous step. If you have no FastLoc data
## and no GPS location classes, then the tagging locations can be filtered
## by going for lc == G. If there is GPS data this section needs to be updated using the
## type == "user" argument.

## If this is not the case and the smaj, smin, eor data is already present from the previous
## script this can be commented out

# det_f <- det_f %>%
#   dplyr::mutate(
#     smaj = ifelse(lc == 'G', 50,
#                   smaj),
#     smin = ifelse(lc == 'G', 50,
#                   smin),
#     eor = ifelse(lc == 'G', 0,
#                  eor)
#   )

# B5: Deal with parametrically logical but biologically unreasonable locations ----

## Only execute this section if above speedfilters do not exclude all spurious
## detections satisfyingly.

## Plot filtered tracks individually. Plot the tracks of the desired speedfilter.

### Basemap
esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
                     'Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')
### Colors
nb.cols <- length(unique(det_f$id))
mycolors <- viridis(nb.cols)

### Tracks
sf_locs <- sf::st_as_sf(det_f, coords = c("lon","lat")) %>%
  sf::st_set_crs(4326)

sf_lines <- sf_locs %>%
  dplyr::arrange(id, date) %>%
  sf::st_geometry() %>%
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs$id))) %>%
  sf::st_cast("MULTILINESTRING") %>%
  sf::st_sf(id = as.factor(unique(det_f$id)))

sf_points <- sf_locs %>%
  dplyr::arrange(id, date) %>%
  sf::st_geometry() %>%
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs$id))) %>%
  sf::st_sf(id = as.factor(unique(det_f$id)))

for (i in 1:length(unique(det_f$id))){

  # plot the filtered locations by individual
  ggplot() +
    annotation_map_tile(type = esri_ocean,zoomin = 1,progress = "none") +
    layer_spatial(sf_points[i,], size = 0.5) +
    layer_spatial(sf_lines[i,], size = 0.75,aes(color = id)) +
    #scale_x_continuous(expand = expansion(mult = c(.6, .6))) +
    scale_colour_manual(values = mycolors[i]) +
    theme() +
    ggtitle("Speedfiltered (argosfilter::sdafilter()) Argos Location Paths",
            subtitle = paste("Sphyrna spp. from the U.S. Atlantic (n = ", length(unique(det_f$id)), ")"))

  # save it
  ggsave(file.path(saveloc,paste0("Speedfiltered_Argos_data_individual_", i, ".tiff")),
         width = 21, height = 15, units = "cm", device ="tiff", dpi=300)
}

## Notes code: if the code throws some errors/warning at you. However, there are not impacting
## the output, especially since the output only serves finding locations that are biologically
## unreasonable.

## Notes from comparing individual plots.
## Most individuals are fine. Some have locations on land in Andros, but given that
## the island is small in relation to available instrument accuray, these can be left in for now.
## Locations on land will later be dealt with by re-routing the tracks around land barriers

## 175427 has 1 spurious location at around -82.5°W and 28°N
## 2023-07-22 00:57:32 at lon. -36.06180 and lat: 34.09160 ###

## 179472 has two locations at around -85W
## one at 2019-11-16 00:03:57
## one at 2019-11-22 11:51:16

## 182746 has 1 location at -55 W
## at 2021-09-27 02:12:30

## 244611 has 1 location at -69.5 W

## 264024 has 1 location at -81.75W

## Filter out these segments manually
det_f <- det_f %>%
  filter(!(  
    #175427
    (id == "175427" & date == as.POSIXct("2025-02-28 15:18:10", tz = "UTC")) |
    #179472
    (id == "179472" & date == as.POSIXct("2019-11-16 00:03:57", tz = "UTC")) |
    (id == "179472" & date == as.POSIXct("2019-11-22 11:51:16", tz = "UTC")) |
    #182746
    (id == "182746" & date == as.POSIXct("2021-09-27 02:12:30", tz = "UTC")) |
    #235278
    (id == "235278" & date == as.POSIXct("2023-07-03 02:00:25", tz = "UTC")) |
    (id == "235278" & date == as.POSIXct("2023-07-02 01:11:49", tz = "UTC")) |
    (id == "235278" & date == as.POSIXct("2023-07-02 00:47:37", tz = "UTC")) |  
    #244611
    (id == "244611" & date == as.POSIXct("2024-02-17 15:55:20", tz = "UTC")) |
    #264024
    (id == "264024" & date == as.POSIXct("2025-02-14 22:40:38", tz = "UTC"))
  ))

## Visualise manually improved tracks again:
sf_locs_clean <- sf::st_as_sf(det_f, coords = c("lon","lat")) %>%
  sf::st_set_crs(4326)

sf_lines_clean <- sf_locs_clean %>%
  dplyr::arrange(id, date) %>%
  sf::st_geometry() %>%
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(det_f$id))) %>%
  sf::st_cast("MULTILINESTRING") %>%
  sf::st_sf(id = as.factor(unique(det_f$id)))

sf_points_clean <- sf_locs_clean %>%
  dplyr::arrange(id, date) %>%
  sf::st_geometry() %>%
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs_clean$id))) %>%
  sf::st_sf(id = as.factor(unique(det_f$id)))

#esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
#'Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')

### ALL
ggplot() +
  annotation_map_tile(type = esri_ocean,zoomin = 1,progress = "none") +
  layer_spatial(sf_points_clean, size = 0.5) +
  layer_spatial(sf_lines_clean, size = 0.75,aes(color = id)) +
  scale_x_continuous(expand = expansion(mult = c(.6, .6))) +
  scale_fill_manual(values = mycolors) +
  theme() +
  ggtitle("Speedfiltered (argosfilter::sdafilter()) and manually checked Argos Location Paths",
          subtitle = paste("Sphyrna spp. from the U.S. Atlantic (n = ", length(unique(det_f$id)), ")"))
ggsave(file.path(saveloc,"Speedfiltered_Argos_data_all_individuals_clean.tiff"),
       width = 21, height = 15, units = "cm", device ="tiff", dpi=150)

### IND
for (i in 1:length(unique(det_f$id))){

  # plot the filtered locations by individual
  ggplot() +
    annotation_map_tile(type = esri_ocean,zoomin = 1,progress = "none") +
    layer_spatial(sf_points_clean[i,], size = 0.5) +
    layer_spatial(sf_lines_clean[i,], size = 0.75,aes(color = id)) +
    #scale_x_continuous(expand = expansion(mult = c(.6, .6))) +
    scale_colour_manual(values = mycolors[i]) +
    theme() +
    ggtitle("Speedfiltered (argosfilter::sdafilter()) and manually checked Argos Location Paths",
            subtitle = paste("Sphyrna spp. from the U.S. Atlantic (n = ", length(unique(det_f$id)), ")"))

  # save it
  ggsave(file.path(saveloc,paste0("Speedfiltered_Argos_data_individual_", i, "_clean.tiff")),
         width = 21, height = 15, units = "cm", device ="tiff", dpi=150)
}


# B6: Calculate the time difference between detections in days and segment tracks ----

## Depending if all sharks survived and/or all tasg were deployed there might be NA values in the dataset
## check
check <- det_f[is.na(det_f$date),] # if these are NAs coming from post-release mortality sharks or undeployed tags, delete them
nrow(check) # if this is 0 you are good to continue

## Calculate the time difference between detections in days
det_f$tdiff.days <- unlist(tapply(det_f$date, INDEX = det_f$id,
                                  FUN = function(x) c(0, `units<-`(diff(x), "days"))))

## Find the maximum values
max(det_f$tdiff.days)
# [1] 149.5232

## Assign different segments if the time difference is larger than the cutoff
## we choose a cut off of 12 days
## Make sure to retain original ptt info and arrange df by timestamp
det_seg <- det_f %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(id = paste0(id,"_", 1+cumsum(tdiff.days >= 12))) %>%
  dplyr::mutate(shark = as.numeric(str_sub(id, start = 1, end = str_locate(id, "\\_")[,1] - 1))) %>%
  dplyr::arrange( # arrange by timestamp by individual so df can be used for fit_() functions
    shark,
    date
  )

## Let's check if there are sharks that have the tagging data as single segment
print(n = 1000,
      det_seg %>%
        dplyr::group_by(shark)%>%
        dplyr::slice(1:2)
)
## currently two sharks: 179472, 210603 have that issue
## we need to add the tagging event manually later on
## to do so we extract the tagging location of these sharks
# extra_tagging <- det_seg %>%
#   # dplyr::filter(
#   #   shark %in% c(235283, 261743)
#   # ) %>%
#   dplyr::group_by(
#     shark
#   ) %>%
#   dplyr::slice(1) # only get first row, i.e. tagging date

## Based on Logan et al. 2020, segments of <20 days should be removed
## Count number of rows per id
tracklengths <- det_seg %>%
  group_by(id) %>%
  summarise(
    num_locs = n(),
    start_date = min(as.Date(date)),
    end_date = max(as.Date(date)),
  ) %>%
  mutate(tracklength_in_days = as.numeric(end_date - start_date) + 1)

print(tracklengths, n = 999)

#View(tracklengths)

### Filter out short track segments
del_obs <- dplyr::filter(tracklengths, num_locs < 12 | tracklength_in_days < 20) # find suitable parameters
# del_obs <- dplyr::filter(tracklengths, tracklength_in_days < 15) # find suitable parameters

print(n = 999, del_obs); nrow(del_obs); sum(del_obs$num_locs)
# [1] 79 segments
# [1] 649 detections

det_seg_tl <- det_seg %>%
  inner_join(., tracklengths) %>%
  dplyr::filter( # remove tracks shorter than X days and/or less than Y total observations
    !(tracklength_in_days < 20 | num_locs < 12)
    # !(tracklength_in_days < 12)
  ) %>%
  dplyr::select( #remove unnecessary columns
    -start_date,
    -end_date
    #-tracklength_in_days
  )

## Find the distribution of time gaps between detections
#ddet_f_cleaner$tdiff.days[ddet_f_cleaner$tdiff.days < 0] <- 0
hist <- hist(det_seg_tl$tdiff.days
             , breaks = c(seq(from = 0, to = max(ceiling(det_seg_tl$tdiff.days)), by = .5))
             , plot = F)

hist$density <- hist$counts / sum(hist$counts)*100
hist$density

plot(hist, freq = F
     , ylim = c(0, 100)
     , col = "skyblue") # approximately 80 % of detections are less or equal to 12 hours apart. Use this info for track predictions.

## nearly 85% of detections are 12h or less apart, for SSM use 12 hour intervals

# B7: fit SSM using aniMotum ----

### make a map of your raw data before model fitting
world <- ne_countries(scale = 10, returnclass = "sf")

ggplot(data = world) +
  geom_sf() +
  geom_point(data = det_seg_tl, aes(x = lon, y = lat, colour = id), size = 2, shape = 20) +
  coord_sf(xlim = c(min(det_seg_tl$lon - 0.5), max(det_seg_tl$lon + 0.5)),
           ylim = c(min(det_seg_tl$lat - 0.5), max(det_seg_tl$lat + 0.5)), expand = F)

## Define parameters

## The fit_ssm() expects following parameters
# d = a data frame of observations including Argos KF error ellipse info (when present)
# vmax = max travel rate (m/s) passed to sda to identify outlier locations
# ang = angles (deg) of outlier location "spikes"
# distlim =lengths (m) of outlier location "spikes"
# spdf = (logical) turn trip::sda on (default; TRUE) or off
# min.dt = minimum allowable time difference between observations; dt <= min.dt will be ignored by the SSM
# pf = just pre-filter the data, do not fit the SSM (default is FALSE)
# model =  fit either a simple random walk ("rw") or correlated random walk ("crw") as a continuous-time process model
# time.step = options: 1) the regular time interval, in hours, to predict to; 2) a vector of prediction times, possibly not regular, must be specified as a data.frame with id and
# POSIXt dates; 3) NA - turns off prediction and locations are only estimated at observation times.
# scale = scale location data for more efficient optimization. This should rarely be needed (default = FALSE)
# emf = optionally supplied data.frame of error multiplication factors for Argos location quality classes. Default behaviour is to use the factors supplied in foieGras::emf()
# map =  a named list of parameters as factors that are to be fixed during estimation, e.g., list(psi = factor(NA))
# parameters = a list of initial values for all model parameters and unobserved states, default is to let sfilter specify these. Only play with this if you know what you are doing
# fit.to.subset = fit the SSM to the data subset determined by prefilter (default is TRUE)
# control = list of control settings for the outer optimizer (see ssm_control for details)
# inner.control = list of control settings for the inner optimizer (see MakeADFun for additional details)
# verbose = [Deprecated] use ssm_control(verbose = 1) instead, see ssm_control for details
# optim = [Deprecated] use ssm_control(optim = "optim") instead, see ssm_control for details
# optMeth = [Deprecated] use ssm_control(method = "L-BFGS-B") instead, see ssm_control for details
# lpsi [Deprecated] use ssm_control(lower = list(lpsi = -Inf)) instead, see ssm_control for details

## Data
## fit_ssm expects 'd' to be a dataframe or tiblle or sf-tibble (with projection info) with 5,7, or 8 columns
## The data should have  5 columns in the following order: "id", "date", "lc", "lon", "lat".
## Where "date" can be a POSIX object or text string in YYYY-MM-DD HH:MM:SS format
## Argos Kalman Filter (or Kalman Smoother) data should have 8 columns, including the above 5 plus
## "smaj", "smin", "eor" that contain Argos error ellipse variables (in m for "smaj", "smin" and deg for "eor").

## Do you want to calculate fitted observations across entire track?
entire_track <- "No" # change between "Yes" and "No", "No" means you are predicting observations at a given time interval using the segmented tracks

if (entire_track == "Yes"){ ## DO NOT change this
  dpf <- det_f
} else (dpf <- det_seg_tl)

## Remember: According to Lea et al. 2013

## "because each raw position has a different error field according to its Argos location class, we needed to
## decide the most probable location for each point within its error field. We achieved this by using a Bayesian
## state-space model (SSM) that adjusted the filtered tracks by producing regular positions based on the Argos location class,
## mean turning angle, and autocorrelation in speed and direction, producing the most probable track through the error fields"

## "Argos tracks only have locations for when the sharks were at the surface; consequently there is high
## variability in the number of locations in a given area, as a result of the shark’s varied surfacing behaviour
## rather than because of its actual location. This would introduce a bias into the analysis of time spent in different
## areas. To correct this bias, linear interpolation was used to normalise the transmission fre- quency by generating
## points at 12 hour intervals along track gaps of <20 days."

## Speedfilter
## For the prefiltering of the data we need a speedfilter. While other Sphyrna papers use the same speed
## as Vaudo et al. 2017 for mako sharks, this possibly is too high.
## Based on Ryan et al. 2015: https://link.springer.com/article/10.1007/s00227-015-2670-4
## we could use speed as function of FL, i.e. speed = 1xFL [m] * s^-1
## Or we could use modelled cruising speed from Payne et al. 2017, i.e. 2.1 m/s

## Prefilter
#vmax = 2.1 # based on Payne et al. 2016
#ang = c(15,25) # Vaudo et al. 2017, values are internal spikes, i.e. for values you need to 165 = 180 - 15, 155 = 180 - 25
#distlim = c(5000, 8000) # Vaudo et al. 2017
spdf = F
pf = F


## Model
model = "crw" # choose between rw, crw, mp. mpm and jpmp available in fit_mpm

## Fitted vs. predicted, i.e. choice of normalisation time-step
if(entire_track == "Yes"){
  time.step = NA # NA turns the time step off and estimates locations at observation times only
} else (time.step = 12)


## Optimizers
optim = "optim"
maxit = 2000
verbose = 2

## Specify SSM
mod.crw_pf <- aniMotum::fit_ssm(
  as.data.frame(dpf), # somehow getting errors if not ran with 'd' in data.frame class
  #vmax = vmax,
  #ang = ang,
  #distlim = distlim,
  spdf = spdf,
  min.dt = 0, # we have taken care of this during sda filtering steps
  pf = pf,
  model = model,
  time.step = time.step,
  #fit.to.subset = T,
  control = ssm_control(
    optim = optim,
    maxit = maxit,
    verbose = verbose), # shows parameter trace, 0: silent, 1: optimizer trace, 2: parametre trace (default)
  map = list(psi = factor(NA)) # this is needed for predicted model as of 2025-07-17
)
# saveRDS(mod.crw_pf, paste0(saveloc, "SSM_fitted.R"))
## we retain 17 segments

##Plot the model with normalised values
#tiff(paste0(saveloc, "SSM_MPM_output/SPOT_tracks_",spp.f,"_", model,"_",optim,"_",maxit, "iterations_fitted_with_argosfilter_output.tiff"),
#height = 30, width = 20, units = "cm", res = 150)
map(mod.crw_pf, what = if(entire_track == "Yes"){"f"} else ("p"), normalise = TRUE, silent = TRUE)
#dev.off()


# The `normalise` argument rescales the estimates to span the interval 0,1.
# Move persistence estimates from `fit_ssm()` tend to be smoothed more extremely compared
# those obtained from `fit_mpm()` and can lack contrast. Normalising the estimates provides a
# clearer view of changes in movement behaviour along tracks, highlighting regions where animals
# spend disproportionately more or less time. When fitting to a collection of individuals,
# the normalisation can be applied separately to each individual or collectively via
# `normalise = TRUE, group = TRUE`. In the latter case, the relative magnitudes of move persistence
# are preserved across individuals.

# B8: re-route path to account for land locations

## Some locations may still on land. Reroute them using route_path()
mod.crw_pf_rr <- route_path(mod.crw_pf,
                            what = if(entire_track == "Yes"){"fitted"} else ("predicted"),
                            map_scale = 10, # scale of rnaturalearth map to use for land mass - if you want 10, you need the package "rnaturalearthhires"
                            dist = 10000, # buffer distance (m) to add around track locations. The convex hull of these buffered locations defines the size of land polygon used to aid re-routing of points on land.
                            buffer = .5,  # buffer distance > 0 moves the track a bit further away from any land barriers, which can (sometimes) result in a better rerouted solution. Some trial and error with the buffer distance (units in km) is usually needed
                            centroids = TRUE, # setting centroids = TRUE essentially adds more resolution (& therefore more possible solutions) to the underlying mesh structure used to reroute the track.
                            append = T
                            )

## check if re-routing track changes nr. observations
check_1 <- grab(mod.crw_pf, what = if(entire_track == "Yes"){"fitted"} else ("predicted"))
check_1 <- check_1[,c(1:4)]
check_rr <- grab(mod.crw_pf_rr, what = "rerouted")
check_rr <- check_rr[,c(1:4)]
comparison <- anti_join(check_1, select(check_rr, -c(lon, lat))) # shows rows that are not present in check_rr , we lose 20 land locations
## plot deleted observations
ggplot(data = world) +
  geom_sf() +
  geom_point(data = comparison, aes(x = lon, y = lat, colour = id), size = 2, shape = 4) +
  coord_sf(xlim = c(min(comparison$lon - 0.5), max(comparison$lon + 0.5)),
           ylim = c(min(comparison$lat - 0.5), max(comparison$lat + 0.5)), expand = F) +
  labs(title = "Predicted locations without re-routed solutions",
       x = "Longitude",
       y = "Latitude"); ggsave(paste0(saveloc,"Predicted_locaations_no_rerouted_solutions.tiff"),
                               width = 21, height = 15, units = "cm", device ="tiff", dpi=150)


## regarding data loss issue: see response ianjonsen
## https://github.com/ianjonsen/aniMotum/discussions/67#discussioncomment-10418664

## map the locations to check if it improved
my.aes <- aes_lst(line = T,
                  conf = F,
                  mp = F,
                  obs = T)
m1 <- map(mod.crw_pf,
          what = if(entire_track == "Yes"){"fitted"} else ("predicted"),
          aes = my.aes,
          ext.rng = c(0.3, 0.1),
          silent = TRUE)

m2 <- map(mod.crw_pf_rr,
          what = "rerouted",
          aes = my.aes,
          ext.rng = c(0.3, 0.1),
          silent = TRUE)


### Comparison of locations between fitted vs. fitted & re-routed
# (m1 + labs(title = "SSM locations") | m2 + labs(title = "SSM re-routed locations")) &
#   theme(panel.grid = element_line(size=0.1, colour="grey60"))

### Detailview fitted & re-routed locations
m1 + labs(title = "Fitted/predicted locations") &
  theme(panel.grid= element_line (size=0.1, colour="grey60"))
ggsave(file.path(saveloc,paste0("Fitted_Predicted_locations_time_step_", time.step,"_h.tiff")), width = 15, height = 12, units = "cm", dpi = 150)

### Detailview fitted & re-routed locations
m2 + labs(title = "Re-routed locations") &
  theme(panel.grid= element_line (size=0.1, colour="grey60"))
ggsave(file.path(saveloc,paste0("Rerouted_locations_time_step_", time.step,"_h.tiff")), width = 15, height = 12, units = "cm", dpi = 150)

# B8: check goodness of fit of your ssm/mpm model ----

# Validating SSM fits to tracking data is an essential part of any analysis.
# SSM fits can be visualized quickly using the generic `plot` function on model fit objects.

# Both fitted and predicted locations (gold) can be plotted as 1-D time-series on top of the observations
# (blue) to visually detect any lack of fit. Observations that failed to pass the `prefilter` stage prior
# to SSM fitting (black x's) can be included (default) or omitted with the `outlier` argument.
# Uncertainty is displayed as a $\pm$ 2 SE envelope (pale gold) around estimates. A rug plot along the
#x-axis aids detection of data gaps in the time-series. Note, in second plot, the larger standard errors for
#predicted locations through small data gaps.

#plot(mod.crw_pf, what = "fitted", type = 2)
#plot(mod.crw_pf[1,], what = "predicted")

## As XY-plots (type 1) and 2-D tracks (type 2)
for (i in 1:2){

  plottype = i # change accordingly, save plottype 1 and 2

  IDs <- (mod.crw_pf_rr$id) ## change accordingly

  for (i in 1:length(IDs)){

    # open plot window - REMEMBER TO RENAME FITTED VS. PREDICTED
    png(file = file.path(saveloc,"Validation",paste0("Goodness_of_predicted_rerouted_", model, "_locations_", IDs[i],"_fitted_with_", optim, "_",maxit,"iterations_", speedfilter, "_filter_type_",plottype,"_plot.png")), res = 150, height = 20, width = 30, units = "cm")

    # plot the fitted locations by individual
    print(plot(mod.crw_pf_rr[i,], "rerouted", type = plottype, alpha = 0.1)) # "f" for fitted and type=2 for 2-D

    # save it
    dev.off()
  }
}

# Residual plots are important for validating models, but classical Pearson residuals,
# for example, are not appropriate for state-space models. Instead, a one-step-ahead prediction
# residual, provides a useful if computationally demanding alternative.
# In `aniMotum`, prediction residuals from state-space model fits are calculated using the
# `osar` function and can be visualized as time-series plots, Q-Q plots, or autocorrelation
# functions:

## calculate & plot residuals

for (i in 1:length(IDs)){

  # subset individual tracks
  res.s <- osar(mod.crw_pf_rr[i,]) ## change accordingly

  # open plot window
  png(file = file.path(saveloc,"Validation", paste0("Prediction_", model,"_fitted","_rerouted_residuals_for_validation_", IDs[i],"_fitted_with_", optim, "_",maxit,"iterations_", speedfilter, "_filter.png")), res = 150, height = 20, width = 30, units = "cm")

  #plot
  (plot(res.s, type = "ts") | plot(res.s, type = "qq")) /
          (plot(res.s, type = "acf") | plot_spacer())

  # save it
  dev.off()
}

# B9: export the model output ----

## what model did you calculate and want to export

m_type <- if(entire_track == "Yes"){"fitted"} else ("predicted") # choose between "Fitted", "Predicted", and "Rerouted" case sensitive

## fitted & predicted locations
### non-projected form
loc <- grab(mod.crw_pf, what = if(entire_track == "Yes"){"fitted"} else ("predicted"))

### projected form
loc.proj <- grab(mod.crw_pf, what = if(entire_track == "Yes"){"fitted"} else ("predicted"), as_sf = T)

## rerouted locations, either fitted or predicted
### non-projected form
locRR <- grab(mod.crw_pf_rr, what = "rerouted")

### projected form
locRR.proj <- grab(mod.crw_pf_rr, what = "rerouted", as_sf = T)

# B10: save the model output as .csv (non-projected form) and as R file (both) ----

# fitted & predicted
write.table(loc, file.path(saveloc, m_type, paste0("Data_aniMotum_CTCRW_output_",m_type,"_non-projected_with_",speedfilter, "_data.csv")),row.names=F,sep=",",dec=".")
saveRDS(loc, file.path(saveloc, m_type, paste0("Data_aniMotum_CTCRW_output_",m_type,"_non-projected_with_",speedfilter, "_data.rds")))
saveRDS(loc.proj, file.path(saveloc, m_type, paste0("Data_aniMotum_CTCRW_output_",m_type,"_projected_with_",speedfilter, "_data.rds")))

# rerouted - fitted - USE THIS IF YOU REROUTED *FITTED* LOCATIONS
write.table(locRR, file.path(saveloc, "rerouted_fitted", paste0("Data_aniMotum_CTCRW_output_",m_type, "_rerouted","_non-projected_with_",speedfilter, "_data.csv")),row.names=F,sep=",",dec=".")
saveRDS(locRR, file.path(saveloc, "rerouted_fitted", paste0("Data_aniMotum_CTCRW_output_",m_type,"_rerouted","_non-projected_with_",speedfilter, "_data.rds")))
saveRDS(locRR.proj, file.path(saveloc, "rerouted_fitted", paste0("Data_aniMotum_CTCRW_output_",m_type,"_rerouted","_projected_with_",speedfilter, "_data.rds")))

# re-routed - predicted - USE THIS IF YOU REROUTED *PREDICTED* LOCATIONS
write.table(locRR, file.path(saveloc, "rerouted_predicted", paste0("Data_aniMotum_CTCRW_output_",m_type, "_rerouted","_non-projected_with_",speedfilter, "_data.csv")),row.names=F,sep=",",dec=".")
saveRDS(locRR, file.path(saveloc, "rerouted_predicted", paste0("Data_aniMotum_CTCRW_output_",m_type,"_rerouted","_non-projected_with_",speedfilter, "_data.rds")))
saveRDS(locRR.proj, file.path(saveloc, "rerouted_predicted", paste0("Data_aniMotum_CTCRW_output_",m_type,"_rerouted","_projected_with_",speedfilter, "_data.rds")))

# for kmneans, dBBMMs, etc. - I.E. PREDICTED LOCS
saveRDS(locRR, file.path(dataloc, "kmeans", paste0("Data_aniMotum_CTCRW_output_",m_type,"_rerouted","_non-projected_with_",speedfilter, "_data.rds")))
saveRDS(locRR, file.path(dataloc, "dBBMMs", paste0("Data_aniMotum_CTCRW_output_",m_type,"_rerouted","_non-projected_with_",speedfilter, "_data.rds")))

### ....................................................................................................
### [C] Calculate CI for location estimates and convert to lat/lon ----
### ....................................................................................................

# The animotum output offers a "merc" location (x and y) as well as the standard error of predicted
# locations (in merc projection also). In case we want to calculate UDs or similar down the line,
# we need a CI for each loaction estimate that then can be used as a estimate for the location error.

# includes individual id, date, longitude, latitude, x and y (typically from the default Mercator
# projection) and their standard errors (`x.se`, `y.se` in km), `u`, `v` (and their standard errors,
# `u.se`, `v.se` in km/h) are estimates of signed velocity in the x and y directions. The `u`, `v`
# velocities should generally be ignored as their estimation uses time intervals between consecutive
# locations, whether they are observation times or prediction times. The columns `s` and `s.se`
# provide a more reliable 2-D velocity estimate, although standard error estimation is turned off
# by default as this generally increases computation time for the `crw` SSM. Standard error
# estimation for `s` can be turned on via the `control` argument to `fit_ssm`
# (i.e. `control = ssm_control(se = TRUE)`, see `?ssm_control` for futher details).


# To calculate a CI in decimal degress, we need to first calculate upper and lower (2.5% and 97.5%) limits
# of our locations. We can do this by adding/substracting 1.96xSE from the x and y coordinates.

# After that, we need to transform the coordinates from the output projection of foiegras to
# decimal degrees.

# D1: find the projection of your fitted ssm/mpm object ----

# chose your data
mod_data <- locRR
mod_data_proj <- locRR.proj

head(mod_data_proj)
# CRS:           +proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs

# D2: calculate upper and lower limits (2.5 and 97.5%) of location estimates ----

## Make new columns for longitude/x
mod_data$x.025 <- mod_data$x - (1.96*mod_data$x.se)
mod_data$x.975 <- mod_data$x + (1.96*mod_data$x.se)

## Make new columns for latitude/y
mod_data$y.025 <- mod_data$y - (1.96*mod_data$y.se)
mod_data$y.975 <- mod_data$y + (1.96*mod_data$y.se)

# D3: coordinates into a WGS84 lat/lon projection ----

## general
cord.merc <- SpatialPoints(cbind(mod_data$x, mod_data$y), proj4string=CRS("+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")) # define SpatialPoints object with output projection
cord.dec <- spTransform(cord.merc, CRS("+proj=longlat +datum=WGS84")) # transform the coordinates to a lat/lon projection
### physical check:
par(mfrow = c(1, 2))
plot(cord.merc, axes = TRUE, main = "Merc Coordinates", cex.axis = 0.95)
plot(cord.dec, axes = TRUE, main = "Lat-Lon Coordinates", col = "red", cex.axis = 0.95)
### create dataframe that can be added to mod.t object
coords.new <- as.data.frame(cord.dec)
colnames(coords.new) <- c("lon.new", "lat.new")

## Lower CI
cord.lowCI <- SpatialPoints(cbind(mod_data$x.025, mod_data$y.025), proj4string=CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs ")) # define SpatialPoints object with output projection
cord.025 <- spTransform(cord.lowCI, CRS("+proj=longlat +datum=WGS84")) # transform the coordinates to a lat/lon projection
### physical check:
par(mfrow = c(1, 2))
plot(cord.lowCI, axes = TRUE, main = "Merc Coordinates", cex.axis = 0.95)
plot(cord.025, axes = TRUE, main = "Lat-Lon Coordinates", col = "red", cex.axis = 0.95)
### create dataframe that can be added to mod.t object
CI025 <- as.data.frame(cord.025)
colnames(CI025) <- c("lon025", "lat025")

## upper CI
cord.upCI <- SpatialPoints(cbind(mod_data$x.975, mod_data$y.975), proj4string=CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs ")) # define SpatialPoints object with output projection
cord.975 <- spTransform(cord.upCI, CRS("+proj=longlat +datum=WGS84")) # transform the coordinates to a lat/lon projection
### physical check:
par(mfrow = c(1, 2))
plot(cord.upCI, axes = TRUE, main = "Merc Coordinates", cex.axis = 0.95)
plot(cord.975, axes = TRUE, main = "Lat-Lon Coordinates", col = "red", cex.axis = 0.95)
### create dataframe that can be added to mod.t object
CI975 <- as.data.frame(cord.975)
colnames(CI975) <- c("lon975", "lat975")

# D4: create df with lower and upper CI for lcoation estiamtes in [dd] ----

mod_CIs <- as.data.frame(cbind(mod_data, coords.new, CI025, CI975, mod_data$lon, mod_data$lat))

## double check that conversion was correct.
check <- mod_CIs
conversion.lon <- check$lon.new == check$mod.t.dd$lon
length(conversion.lon[conversion.lon == F]) # needs to be 0
conversion.lat <- check$lat.new == check$mod.t.dd$lat
length(conversion.lat[conversion.lat == F]) # needs to be 0
## if 0 twice continue

## create df and prepare for saving
## if SSM
mod.final <- dplyr::select(mod_CIs, id, date, lon.new, lon025, lon975, lat.new, lat025, lat975)
colnames(mod.final) <- c("id", "date", "lon", "lon025", "lon975", "lat", "lat025", "lat975")

## if MP
#mod.final <- dplyr::select(mod_CIs, id, date, lon.new, lon025, lon975, lat.new, lat025, lat975,g,g_normalized)
#colnames(mod.final) <- c("id", "date", "lon", "lon025", "lon975", "lat", "lat025", "lat975", "gamma_t", "gamma_t_normalized")

## this causes issues with kmeans in following scripts
# ## if you have predicted locations you need to add the tagging locations again
# tagging_selected <- extra_tagging %>%
#   dplyr::select(any_of(names(mod.final)))
# tagging_selected$date <- strptime(tagging_selected$date, format = "%Y-%m-%d %H", tz = "UTC")
#
# if (entire_track == "No")  {
#   mod.final <- mod.final %>%
#     dplyr::full_join(tagging_selected, by = c("id", "date"), suffix = c("", "_new")
#   ) %>%
#     dplyr::mutate(
#       lon = ifelse(is.na(lon), lon_new, lon),
#       lon025 = ifelse(is.na(lon025), lon, lon025), # if tagging date got removed, lon025 will be NA, make sure this is tagging loc
#       lon975 = ifelse(is.na(lon975), lon, lon975), # see above
#       lat = ifelse(is.na(lat), lat_new, lat),     # see above
#       lat025 = ifelse(is.na(lat025), lat, lat025), # see above
#       lat975 = ifelse(is.na(lat975), lat, lat975)  # see above
#     ) %>%
#     dplyr::select( # drop columns that are not needed
#       -lon_new,
#       -lat_new,
#       -shark
#     )
# } else (mod.final <- mod.final)
#
# ## make sure that the df is showing dates in ascending order
#
# mod.final %<>%
#   dplyr::mutate(shark = as.numeric(str_sub(id, start = 1, end = str_locate(id, "\\_")[,1] - 1))) %>%
#   dplyr::arrange(shark, date) %>%
#   dplyr::select(-shark)

### ....................................................................................................
### [E] Save movement data for further analyses ----
### ....................................................................................................

# if based on speedfilter data
if(entire_track == "Yes"){
  tracktype = "entire_track"
} else (tracktype = "segmented")

write.table(mod.final, file.path(saveloc, "rerouted_predicted", paste0("Data_aniMotum_CRW_output_",tracktype,"_rerouted","_proj_WGS84_converted_with_coord_CIs_with_",speedfilter, "_data.csv")),row.names=F,sep=",",dec=".")
saveRDS(mod.final, file.path(saveloc, "rerouted_predicted", paste0("Data_aniMotum_CRW_output_",tracktype,"_rerouted","_proj_WGS84_converted_with_coord_CIs_with_",speedfilter, "_data.rds")))

saveRDS(mod.final, file.path(dataloc, "kmeans", paste0("Data_aniMotum_CRW_output_",tracktype,"_rerouted","_proj_WGS84_converted_with_coord_CIs_with_",speedfilter, "_data.rds")))
saveRDS(mod.final, file.path(dataloc, "dBBMMs", paste0("Data_aniMotum_CRW_output_",tracktype,"_rerouted","_proj_WGS84_converted_with_coord_CIs_with_",speedfilter, "_data.rds")))

# END OF CODE ----