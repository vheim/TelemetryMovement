# TOP OF CODE ----
### ====================================================================================================
### Project:    Large bodied hammerhead complex in the western North Atlantic
### Analysis:   kmeans clustering for behavioral state estimation
### Script:     ~2025_NWA_hammerhead_complex/R/Rscript4_kmeans_clustering_behavioural_states.R
### Author:     Vital Heim (inspired by Simon Dedman)
### Version:    1.0
### ====================================================================================================

### ....................................................................................................
### Content: this code allows to define behavioural states along movement trajectories using the k-means 
### algorithm and cluster analysis as described in van Moorter et al. 2010:
### https://doi.org/10.2193/2009-155
### The code was largely sourced from Simon Dedman's GitHub repositorty:
### https://github.com/SimonDedman/SavingTheBlue
### ....................................................................................................

### ....................................................................................................
### [A] Ready environment, load packages ----
### ....................................................................................................

# A1: clear memory ----

rm(list = ls())

# A2: install and load necessary packages ----

## if first time
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("magrittr")
# install.packages("amt") #https://arxiv.org/pdf/1805.03227.pdf #sudo apt install libgsl-dev
# install.packages("RcppRoll")
# install.packages("lwgeom")
# install.packages("TropFishR") # VBGF 
# install.packages("data.table")
# install.packages("purrr")
# install.packages("beepr")
# install.packages("moments")
# install.packages("rgl")
# install.packages("clusterSim")
# install.packages("raster")
# install.packages("terra")
# install.packages("ncdf4")
# install.packages("stars")
# install.packages("rerddap")
# install.packages("rerddapXtracto")
# install.packages("maps")
# install.packages("mapdata")
# install.packages("RcolorBrewer")
# install.packages("maptools")
# install.packages("rasterVis")
# install.packages("ggplot2")
# install.packages("mapplots")
# remotes::install_github("SimonDedman/gbm.auto")# for cropmap
# install.packages("marmap")
# install.packages("scales")

## load packages and source needed functions
### general and general movement packages
library(tidyverse)
library(dplyr)
library(magrittr)
library(amt) #https://arxiv.org/pdf/1805.03227.pdf #sudo apt install libgsl-dev
library(RcppRoll)
library(lwgeom)
library(TropFishR) #VBGF
library(data.table)
library(lubridate)
library(purrr)
library(beepr)
# source("C:/Users/Vital Heim/switchdrive/Science/Rscripts/vanMoorter-et-al_2010/liRolling.R") # Simon's function for rolling Linearity Index values
### for data cleaning
library(moments)
### for vanmoorter
library(beepr)
library(rgl) # # sudo apt install libglu1-mesa-dev
library(clusterSim)
# source("C:/Users/Vital Heim/switchdrive/Science/Rscripts/vanMoorter-et-al_2010/p7_gap.statistic.r")
### for raster manipulation
library(raster)
library(terra)
library(ncdf4)
library(stars)
### for SST extraction
library(rerddap)
library(rerddapXtracto)
### optional, for nicer maps
# library(maps)
# library(mapdata)
# library(RcolorBrewer)
# library(maptools)
# library(rasterVis)
# library(ggplot2)
library(mapplots)
library(gbm.auto)
library(marmap)
library(rnaturalearth)
library(scales)

# A3: Specify needed functions ----

YOUR_IP <- "NA" # add your IP address, server name or similar if you connect via a shared drive

# *A3.1: range standardizations ----
## Function for range standardization
STrange <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

## Function to reverse-transform range standardization
rt_STrange <- function(x,y) {
  y*(max(x, na.rm = TRUE)-min(x, na.rm = TRUE)) + min(x, na.rm = TRUE)
}

# *A3.2: linearity index values function ----
source(file.path("/",YOUR_IP, "Science","Rscripts","vanMoorter-et-al_2010", "liRolling.R"))

# *A3.3: gap statistic function ----
source(file.path("/", YOUR_IP, "Science", "Rscripts","vanMoorter-et-al_2010","p7_gap.statistic.r" ))

# *A3.4: 2d barplot function ----
source(file.path("/",YOUR_IP, "Science","Rscripts","Barplot2dMap", "Barplot2dMap.R"))

# *A3.5: mapplot XYZ grid function ----
source(file.path("/",YOUR_IP, "Science","Rscripts","Mapplot","makeXYZ_HGerritsen.R"))

# A4: Specify data and saveloc ----

# YOUR_IP <- "NA" # add your IP address, server name or similar if you connect via a shared drive

## Input data
dataloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2025_NWA_hammerhead_complex","Data_input", "Telemetry")

## Output data
saveloc <- file.path("/", YOUR_IP,"Science","Projects_current", "2025_NWA_hammerhead_complex","Data_output", "Telemetry", "kmeans") # Adjust this

# A5: Define universal variables (e.g. for plotting) ----

## MAP APIs
## if you want to use stamen/stadia and/or google maps for plotting your UD rasters, you need to register your corresponding API key first
## for google maps this can be done in the movegroup::plotraster() function with the gmapsAPI = "YOUR_KEY_HERE" argument.
## For stamen/stadia maps you can register your key here (this will be updated in the movegroup package in the near future)

ggmap::register_stadiamaps("YOUR_API_KEY_HERE", write = F)

## Plotting colours

### define your plotting colors for kmeans behavioural tracks
behav_col <- c("#236684" , "#EDA904") # needs to be adjusted if more than 2 clusters in behavioural sates

# A6: Define universal options, variables, etc. (e.g. for plotting) ----

options(timeout = 3000) # manually increase time out threshold (needed when downloading basemap)
options(scipen=999) # so that R doesn't act up for pit numbers  
options(warn=1) #set this to two if you have warnings that need to be addressed as errors via traceback()
options(error = function() beep(9))  # give warning noise if it fails

### ....................................................................................................
### [B] Step length and turning angles ----
### ....................................................................................................

## remember - there are multiple species in your dataset. As per 20250718 we run the kmeans for each species
## separately.
### HOWEVER, this only becomes relevant starting with section [C] of this code.

# B1: Import needed data

## movement data
## use the regularised movement data output from the 12h CTCRW
hammers <- readRDS(file = file.path(dataloc, paste0("kmeans/Data_aniMotum_CRW_output_segmented_rerouted_proj_WGS84_converted_with_coord_CIs_with_Argosfilter_data.rds"))) |>
  mutate(shark = as.numeric(str_sub(id, start = 1, end = str_locate(id, "\\_")[,1] - 1)))

## metadata for length measurements
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

# B2: calculate step length and turning angle ----

## Step length will be calculated in km and relative to body length of the corresponding shark

# > li5day {5 day linearity index. 1: linear paths, 0: tortuous paths}
# > StepLengthKm {distance from previous day; Km}
# > StepLengthBL {distance from previous day; body lengths of that fish}
# > TurnAngleRelDeg {turning angle from previous day, relative to previous day at 0}
# > TurnAngleAzimDeg  {turning angle from previous day, azimuthal degrees}
# > FishLengthCm {fish length in cm}
hammers$li5day <- as.numeric(rep(NA, nrow(hammers)))
hammers$StepLengthKm <- as.numeric(rep(NA, nrow(hammers)))
hammers$StepLengthBL <- as.numeric(rep(NA, nrow(hammers)))
hammers$TurnAngleRelDeg <- as.numeric(rep(NA, nrow(hammers)))
hammers$TurnAngleAzimDeg <- as.numeric(rep(NA, nrow(hammers)))
hammers$Index <- 1:nrow(hammers)

if (!all(is.na(hammers$lat))) { # if not all lats are NA, i.e. there's something to be done
  # amt steps etc here
  #StepLength, TurnAngleRelDeg, TurnAngleAzimDeg
  #dfi remove NA dates and lats and lons rows
  # df_nona <- hammers[!is.na(hammers$Date),] 
  
  # *b2.1: omit rows with NA values for date, downsample to days only ----
  df_nona <- hammers
  df_nona[which(df_nona$lat > 90), "lat"] <- NA # lat over 90 breaks, impossible value, will fix upstream
  df_nona[which(df_nona$lat < -90), "lat"] <- NA # ditto
  df_nona[which(df_nona$lon > 180), "lon"] <- NA # lon over 180 breaks, impossible value, will fix upstream
  df_nona[which(df_nona$lon < -180), "lon"] <- NA # ditto
  df_nona <- df_nona[!is.na(df_nona$lat),] # omit rows with NA values for lat, downsample to days only
  df_nona <- df_nona[!is.na(df_nona$lon),] # omit rows with NA values for lon, downsample to days only
  fishlist <- unique(df_nona$shark)
  
  # "b2.2: loop id, calc li5day, make track ----
  for (i in fishlist) { #i <- fishlist[9]
    df_nonai <- df_nona[which(df_nona$shark == i),] #subset to each shark
    print(paste0(which(fishlist == i), " of ", length(fishlist), "; adding transit dive data to ", i))
    setDF(df_nonai) # else liRolling breaks
    if (nrow(df_nonai) > 5) df_nonai$li5day <- liRolling(x = df_nonai,
                                                         coords = c("lon", "lat"), # original, "lat", "lon"
                                                         roll = 5)
    ## 1: linear paths, 0: tortuous paths
    if (!nrow(df_nonai) > 1) next # track etc breaks with < 2 points
    track <- amt::mk_track(df_nonai,
                           .x = lon,
                           .y = lat,
                           .t = date, # was DateTimeUTCmin5
                           crs = 4326) %>%  #crs = sp::CRS("+init=epsg:4326")) %>%
      transform_coords(3395)# only needed here as we need projected crs in meters for SL but not for TA
      # transform_coords(6931)# only needed here as we need projected crs in meters for SL but not for TA
    stps <- steps(track)
    # "b2.3: step  length ----
    ## (sl;  in  CRSunits), 6931=m ## CHANGE THIS BASED ON YOUR transform_coords() epsg
    df_nonai$StepLengthKm <- c(NA, (stps$sl_/1000)) #add NA first since results start from second row
    
    # *b2.4: body lengths ----
    # use vonB in reverse to convert age to length per day in dfnona
    # VonBertalanffy parameters, Restrepo et al 2010:
    # Linf (cm): 314.9
    # k: 0.089
    # t0 (year): -1.13
    # variance Linf: 19.43 (not reqd?)
    # have age from length, want length from age
    # if age is missing, populate with NA instead of crashing
    
    # if (is.na(df_nonai$age[1])) {
    #   df_nonai$FishLengthCm <- rep(NA, nrow(df_nonai))
    #   df_nonai$StepLengthBL <- rep(NA, nrow(df_nonai))
    # } else {
    #   df_nonai$FishLengthCm <- TropFishR::VBGF(param = list(Linf = 314.9,
    #                                                         K = 0.089,
    #                                                         t0 = -1.13),
    #                                            t = df_nonai$age)
    
    
    # then divide StepLengthKm to BL
    # 1 km/day = 1000m/day = 100000cm/day
    df_nonai$StepLengthBL <- (df_nonai$StepLengthKm * 100000) / df_nonai$FishLengthCm
    
    # results is distance in body lengths per day.
    # Body lengths per second (mean value per day):
    # summary(df_nonai$StepLengthBL/24/60/60)
    
    # *b2.5: turning angles & track2 ----
    # turning angles (ta; in degrees;  notice that it cannot be calculated for steps that are not
    # preceded by a valid step), angles are between -pi and pi, *57.29578 to convert to +-180
    df_nonai$TurnAngleRelDeg <- c(NA, (stps$ta_*57.29578))
    
    # amt::movement_metrics
    # all could be added to metadata, one value per track
    # see behavr https://rethomics.github.io/behavr.html
    # straightness(track)
    # cum_dist(track)
    # tot_dist(track)
    # msd(track)
    # intensity_use(track)
    # sinuosity(track)
    # tac(track)
    ## Do I actually want 1 value per track though?? Up to 6 years, becomes meaningless.
    
    track2 <- amt::mk_track(df_nonai,
                            .x = lon,
                            .y = lat,
                            .t = date, # was DateTimeUTCmin5
                            crs = 4326) # sp::CRS("+init=epsg:4326")
    #TurnAngleAzimDeg: compass bearing of new direction, can compute for directionality, migration
    TurnAngleAzimDeg <- direction_abs(track2, full_circle = FALSE, zero_dir = "N", lonlat = TRUE, clockwise = TRUE) %>%
      as_degree
    TurnAngleAzimDeg <- TurnAngleAzimDeg[1:length(TurnAngleAzimDeg) - 1] #remove NA from back
    df_nonai$TurnAngleAzimDeg <- c(NA, TurnAngleAzimDeg) # add NA to front
    df_nonai %<>% dplyr::select(c(date, li5day, StepLengthKm, StepLengthBL, TurnAngleRelDeg, TurnAngleAzimDeg, FishLengthCm, shark))
    
    # *b2.6: save data to df ----
    setDT(hammers) # convert to data.table without copy
    setDT(df_nonai) # convert to data.table without copy
    # join and update "hammers" by reference, i.e. without copy
    # https://stackoverflow.com/questions/44930149/replace-a-subset-of-a-data-frame-with-dplyr-join-operations
    # If you want to return only df_nonai that have a matching hammers (i.e. rows where the key is in both tables), set the nomatch argument of data.table to 0.
    # nomatch isn't relevant together with :=, ignoring nomatch
    hammers[df_nonai, on = c("date", "shark"), li5day := i.li5day]
    hammers[df_nonai, on = c("date", "shark"), StepLengthKm := i.StepLengthKm]
    hammers[df_nonai, on = c("date", "shark"), StepLengthBL := i.StepLengthBL]
    hammers[df_nonai, on = c("date", "shark"), TurnAngleRelDeg := i.TurnAngleRelDeg]
    hammers[df_nonai, on = c("date", "shark"), TurnAngleAzimDeg := i.TurnAngleAzimDeg]
    hammers[df_nonai, on = c("date", "shark"), FishLengthCm := i.FishLengthCm]
    
    # kinda flat, low depth range. See weekly pdfs.
    # dives per Y hours < Z  (depends on definition of ‘dive’) = 0
    # distance covered per day > A? Do histogram of distances, then per
    # behaviour type to see if distance covered is sig higher when transitioning, e.g.
  } # close i
  hammers <- hammers[order(hammers[,"Index"]),] #reorder by index
  hammers %<>% dplyr::select(-Index) # remove unneeded Index col. Might need later?
  hammers <- as.data.frame(hammers)
} else {# close if (!all(is.na(hammers$lat)))
  print("all new days missing latitude data, can't get external data, nothing to do")
}

### ....................................................................................................
### [C] K-means Van Moorter et al. 2010 ----
### ....................................................................................................

## code based on van moorter et al. 2010 and Simon Dedman

## we need:
# p7_analysis_suppl_BL.R
# > kmeans2cluster {resident/transient	movement cluster based on body lengths}
# > kmeansBinary {0/1	movement cluster based on body lengths}
# >> clusterinfo.csv saved in saveloc

# C1: if multiple species, split data frame

spp_f <- "S.zygaena"
spp <- gsub("\\.", "", spp_f) # createa a character string to name files when saving, "." might cause issues downstream

hammers_i <- hammers %>%
  dplyr::filter(
    species %in% spp_f
  )

# C2: add columns for clusters ----

hammers_i$kmeans2cluster <- as.character(rep(NA, nrow(hammers_i)))
hammers_i$kmeansBinary <- as.integer(rep(NA, nrow(hammers_i)))

# C2: prepare steplength and turning angle data ----

## VanM paper: We log-transformed both activity measures and steplength to reduce positive skew. Furthermore, we standardized all variable values on their range (Steinley 2006a).
## Because we found no outliers (no data seemed outlying by visual inspection of the distribution, nor were any observations >3.3 SDs from the mean, which corresponds to a density of 0.001 in a normal distribution), we did not remove any data before data standardization.
## To demonstrate the importance of data preprocessing, we first performed a cluster analysis on the raw data, these data after log transformation and then after range standardization. In subsequent analyses, we used the transformed and standardized data, which is the recommended strategy.
## The raw data without any data preprocessing did not reveal any cluster structure; the gap statistic was maximum for one cluster (i.e., no cluster structure). After log transforming both activity measures and step lengths, we found a structure with 2 clusters.
## Steinley, D. 2006a. K-means clustering: a half-century synthesis. British Journal of Mathematical and Statistical Psychology 59:1–34.
## Standardisation on the range (not ((data - mean)/SD)) recommended but can't get that paper:
## Steinley, D. (2004a). Standardizing variables in K-means clustering. In D. Banks, L. House, F. R. McMorris, P. Arabie, & W. Gaul (Eds.), Classification, clustering, and data mining applications (pp. 53–60). New York: Springer.
## Suspect it's just data-mean

# *C2.1: test for positive skew, log transform if needed ----

## Steplengths
hist(hammers_i$StepLengthBL)
moments::skewness(hammers_i$StepLengthBL, na.rm = TRUE) # positive/negative skew, 0 = symmetrical. #13.20943
moments::kurtosis(hammers_i$StepLengthBL, na.rm = TRUE) # long tailed? Normal distribution = 3. #411.5077
hammers_i$StepLengthBLlog1p <- log1p(hammers_i$StepLengthBL) # log transform
hist(hammers_i$StepLengthBLlog1p)
moments::skewness(hammers_i$StepLengthBLlog1p, na.rm = TRUE) # -0.409 fine
moments::kurtosis(hammers_i$StepLengthBLlog1p, na.rm = TRUE) # 3.14 fine

## TurningAngles
hist(hammers_i$TurnAngleRelDeg)
moments::skewness(hammers_i$TurnAngleRelDeg, na.rm = TRUE) # positive/negative skew, 0 = symmetrical. #0.024
moments::kurtosis(hammers_i$TurnAngleRelDeg, na.rm = TRUE) # long tailed? Normal distribution = 3. #4.26
# hammers_i$TurnAngleRelDeglog1p <- hammers_i$TurnAngleRelDeg # only use this if no transformation is needed
hammers_i$TurnAngleRelDeglog1p <- log1p(abs(hammers_i$TurnAngleRelDeg))
hist(hammers_i$TurnAngleRelDeglog1p)
moments::skewness(hammers_i$TurnAngleRelDeglog1p, na.rm = TRUE) # -0.19 fine
moments::kurtosis(hammers_i$TurnAngleRelDeglog1p, na.rm = TRUE) # 3.14 fine

# *C2.2: Standardize variables ----

# As per van Morter et al. 2010: "Due to the role of the minimum and maximum value of a variable in range standardization,
# it is crucial to inspect the data for outliers (i.e., atypical data that are distant from the rest of the data) before
# performing this standardization."
# And: "outliers (no data seemed outlying by visual inspection of the distribution, nor were any observations >3.3 SDs
# from the mean, which corresponds to a density of 0.001 in a normal distribution).

## test outliers, remove >3.3 SDs from the mean
## Given that the definition of outliers depends on the mean, and we will later calculate the kmeans by individual,
## we need to test for outliers at the individual level as well.

## safety copy in case I mess up
df_i <- hammers_i

## Plot values with outliers
png(filename = file.path(saveloc, paste0("Kmeans-StepLength-TurnAngle-Scatter_", today(),"_all_data_",spp,".png")), width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
plot(df_i$StepLengthBLlog1, df_i$TurnAngleRelDeglog1p, xlab = "log(Step Length [BLs]) w/ outliers", ylab = "log(abs(Turn Angle Relative Degress w/ outliers))")
dev.off()
## StepLength Outliers
if (any(df_i$StepLengthBLlog1p > mean(df_i$StepLengthBLlog1p, na.rm = TRUE) + sd(df_i$StepLengthBLlog1p, na.rm = TRUE) * 3.3, na.rm = TRUE)) df_i$StepLengthBLlog1p[which(df_i$StepLengthBLlog1p > mean(df_i$StepLengthBLlog1p, na.rm = TRUE) + sd(df_i$StepLengthBLlog1p, na.rm = TRUE) * 3.3)] <- NA
if (any(df_i$StepLengthBLlog1p < mean(df_i$StepLengthBLlog1p, na.rm = TRUE) - sd(df_i$StepLengthBLlog1p, na.rm = TRUE) * 3.3, na.rm = TRUE)) df_i$StepLengthBLlog1p[which(df_i$StepLengthBLlog1p > mean(df_i$StepLengthBLlog1p, na.rm = TRUE) + sd(df_i$StepLengthBLlog1p, na.rm = TRUE) * 3.3)] <- NA
## TurningAngles Outliers
if (any(df_i$TurnAngleRelDeglog1p > mean(df_i$TurnAngleRelDeglog1p, na.rm = TRUE) + sd(df_i$TurnAngleRelDeglog1p, na.rm = TRUE) * 3.3, na.rm = TRUE)) df_i$TurnAngleRelDeg[which(df_i$TurnAngleRelDeg > mean(df_i$TurnAngleRelDeg, na.rm = TRUE) + sd(df_i$TurnAngleRelDeg, na.rm = TRUE) * 3.3)] <- NA
if (any(df_i$TurnAngleRelDeglog1p < mean(df_i$TurnAngleRelDeglog1p, na.rm = TRUE) - sd(df_i$TurnAngleRelDeglog1p, na.rm = TRUE) * 3.3, na.rm = TRUE)) df_i$TurnAngleRelDeg[which(df_i$TurnAngleRelDeg > mean(df_i$TurnAngleRelDeg, na.rm = TRUE) + sd(df_i$TurnAngleRelDeg, na.rm = TRUE) * 3.3)] <- NA
## Plot values without outliers
png(filename = file.path(saveloc, paste0("Kmeans-StepLength-TurnAngle-Scatter_", today(),"_no_outliers_all_data_", spp,".png")), width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
plot(df_i$StepLengthBLlog1, df_i$TurnAngleRelDeglog1p, xlab = "log(Step Length [BLs]) w/o outliers", ylab = "log(abs(Turn Angle Relative Degrees w/o outliers))")
dev.off()

## Now we can range-standardise the values for step lengths and turning angles
## We go with the range standardization, i.e. zi = (xi - min(x))/(max(x)-min(x)), suggested by Steinley (2004 and 2006)
df_i$StepLengthBLlog1pST <- STrange(df_i$StepLengthBLlog1p)
df_i$TurnAngleRelDeglog1pST <- STrange(df_i$TurnAngleRelDeglog1p) ## store cleaned and outlier free data for each individual in the list
## plot range standardised values
png(filename = file.path(saveloc, paste0("Standardised_Kmeans-StepLength-TurnAngle-Scatter_", today(),"_",spp,"_all.png")), width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
plot(df_i$StepLengthBLlog1pST, df_i$TurnAngleRelDeglog1pST, xlab = "Standardised log(Step Length [BLs])", ylab = "Standardised log(abs(Turn Angle Relative Degrees))")
dev.off()

# C3: calculate k-means ----

hammers_i <- df_i # convert back if above steps work correctly

setDT(hammers_i) # convert to data.table without copy

if (!all(is.na(hammers_i$lat))) { # if not all lats are NA, i.e. there's something to be done
  fishlist <- unique(hammers_i$shark) # already exists, is same
  for (i in fishlist) { # i <- fishlist[8]
    df_i <- hammers_i[which(hammers_i$shark == i),] #subset to each fish
    print(paste0(which(fishlist == i), " of ", length(fishlist), "; calculating KMeans clustering for ", i))
    if (nrow(df_i) < 5) next # skip else kmeans will break.
    setDF(df_i) # else things break

    # x <- data[,c("ACT1", "ACT2", "SL_METERS", "TURN_DEGRE")]
    # downsample df_i to days
    x <- df_i[!is.na(df_i$StepLengthBLlog1pST),] # omit NA rows
    x <- x[!is.na(x$TurnAngleRelDeglog1pST),] # ditto
    if (nrow(x) < 5) next # skip else kmeans will break.
    x <- x[,c("StepLengthBLlog1pST", "TurnAngleRelDeglog1pST")] # , "Date", "Index"
    
    ## *C3.1: determine GAP statistic for n clusters ----
    setDF(x) # else things break
    res <- data.frame(GAP = NA, s = NA, Wo = NA, We = NA)
    for (j in 1:5) { ##determine the GAP-statistic for 1 to 5 clusters. Changed from 10 to reduce compute time and we don't anticipate >4 XY XYT movement clusters
      if (j == 1) {  ##clall is the vector (in matrix format) of integers indicating the group to which each datum is assigned
        ones <- rep(1, nrow(x))
        clall <- matrix(ones)}
      if (j > 1) {
        cl1 <- kmeans(x, j, iter.max = 100)
        clall <- matrix(cl1$cluster)}
      g <- index.Gap.modif(x, clall, reference.distribution = "pc", B = 50, method = "k-means")
      res[j,] <- c(g$gap, g$s, g$Wo, g$We)
    } # close j
    
    ## *C3.2: plot GAP stat per clusters+SE ----
    par(mfrow = c(1,1))
    k <- seq(1:length(res$GAP))
    # png(filename = paste0(saveloc, "/KmeansClusters_", i, ".png"), width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
    png(filename = file.path(saveloc, "gap_plots", paste0("KmeansClusters_", i,"_",spp,"_",today(), ".png")), width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
    plot(k, res$GAP,
         xlab = expression(bold("Nr. clusters k")),
         ylab = expression(bold("GAP")),
         main = paste0("GAP statistic for shark id ",i, "of species ",spp_f), font.main = 2,
         type = "b")
    segments(k, c(res$GAP - res$s), k, c(res$GAP + res$s))
    kstar <- min(which(res$GAP[-length(res$GAP)] >= c(res$GAP - res$s)[-1])) # if none of the first set of values are individually smaller than their paired counterparts in the second set of values then which() produces all FALSEs and min() fails.
    kstar2 <- min(which(res$GAP[-length(res$GAP)] >= c(res$GAP - (1.96 * res$s))[-1])) # same
    if (!is.infinite(kstar2)) points(kstar2, res$GAP[kstar2], pch = 22, bg = "gray", cex = 1.25) #grey square box for tolerance2 # meaningless if kstar2 fails (now not run if so)
    if (!is.infinite(kstar)) points(kstar, res$GAP[kstar], col = "black", pch = 19, cex = 0.6) # black dot for tolerance1  # meaningless if kstar fails #change this to a different symbol?
    dev.off()
    
    ## *C3.3: plot var1 vs var2 clustering scatterplots ----
    # png(filename = paste0(saveloc, "/Kmeans-StepLength-TurnAngle-Scatter_", i, ".png"), width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
    png(filename = file.path(saveloc, "scatterplots", paste0("Kmeans-StepLength-TurnAngle-Scatter_", i, "_",spp,"_", today(),"_final.png")), width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
    plot(x$StepLengthBLlog1pST, x$TurnAngleRelDeglog1pST, xlab = "Standardised log(Step Length [BLs])", ylab = "Standardised log(abs(Turn Angle Relative Degrees))")
    dev.off()
    
    ## *C3.4: redo kmeans with selected number of clusters (kmeans output cl1 gets overwritten per i) ----
    ##### TODO make centers dynamic
    # See "Show how many clusters were chosen most commonly"
    # kstar & kstar2, might be different.
    kmeans2 <- kmeans(x, centers = 2, iter.max = 100) #run kmeans
    df_i[as.integer(names(kmeans2$cluster)), "kmeans2cluster"] <- kmeans2$cluster
    
    ## *C3.5: label transit/resident clusters algorithmically ----
    # kmeans2$centers
    # # StepLengthBLlog1p TurnAngleRelDeg
    # # 1    0.5816330       0.3355378
    # # 2    0.4993464       0.7817846
    # # these are standardised lognormalised values, not actual steplengths.
    # # Can get actual mean values since cluster bins are now added to extracted.
    cl1sl <- mean(df_i[which(df_i$kmeans2cluster == 1), "StepLengthBLlog1pST"], na.rm = T) # mean(extractedK1$StepLengthBLlog1p) #65.93499
    cl1ta <- mean(abs(df_i[which(df_i$kmeans2cluster == 1), "TurnAngleRelDeglog1pST"]), na.rm = T) # mean(extractedK1$TurnAngleRelDeg) #6.476207
    cl2sl <- mean(df_i[which(df_i$kmeans2cluster == 2), "StepLengthBLlog1pST"], na.rm = T) # mean(extractedK2$StepLengthBLlog1p) #39.47376
    cl2ta <- mean(abs(df_i[which(df_i$kmeans2cluster == 2), "TurnAngleRelDeglog1pST"]), na.rm = T) # mean(extractedK2$TurnAngleRelDeg) #74.15654

    df_i$kmeansBinary <- rep(NA, nrow(df_i))
    df_i$kmeansCharacter <- rep(NA, nrow(df_i))
    
    if (all(cl1sl > cl2sl, cl1ta < cl2ta, na.rm = TRUE)) { # cl1 longer & straighter
      df_i[which(df_i$kmeans2cluster == 1), "kmeansCharacter"] <- "transit" #replace 1&2 with named groups
      df_i[which(df_i$kmeans2cluster == 2), "kmeansCharacter"] <- "resident"
    } else if (all(cl1sl < cl2sl, cl1ta > cl2ta, na.rm = TRUE)) { #cl2 longer & straighter
      df_i[which(df_i$kmeans2cluster == 1), "kmeansCharacter"] <- "resident"
      df_i[which(df_i$kmeans2cluster == 2), "kmeansCharacter"] <- "transit"
    } else if (cl1sl > cl2sl) { # cl1 longer but also twistier
      if ((cl1sl / cl2sl) - (cl1ta / cl2ta) > 0) { #ratio of straightness > ratio of twistyness
        df_i[which(df_i$kmeans2cluster == 1), "kmeansCharacter"] <- "transit"
        df_i[which(df_i$kmeans2cluster == 2), "kmeansCharacter"] <- "resident"
      } else { #ratio of straightness < ratio of twistyness
        df_i[which(df_i$kmeans2cluster == 1), "kmeansCharacter"] <- "resident"
        df_i[which(df_i$kmeans2cluster == 2), "kmeansCharacter"] <- "transit"
      }
    } else if (cl1sl < cl2sl) { # cl1 shorter but also straighter
      if ((cl2sl / cl1sl) - (cl2ta / cl1ta) > 0) { #ratio of straightness > ratio of twistyness
        df_i[which(df_i$kmeans2cluster == 1), "kmeansCharacter"] <- "resident"
        df_i[which(df_i$kmeans2cluster == 2), "kmeansCharacter"] <- "transit"
      } else { #ratio of straightness < ratio of twistyness
        df_i[which(df_i$kmeans2cluster == 1), "kmeansCharacter"] <- "transit"
        df_i[which(df_i$kmeans2cluster == 2), "kmeansCharacter"] <- "resident"
      }
    } # close if elses block
    # see Blocklab/abft_diving/X_PlotsMisc/KMeans/clusterinfo.ods for worksheet looking at these.
    df_i[which(df_i$kmeansCharacter == "resident"), "kmeansBinary"] <- 0
    df_i[which(df_i$kmeansCharacter == "transit"), "kmeansBinary"] <- 1
    
    # kmeans2$withinss
    # [1] 4.705283 9.493627
    # kmeans2$size
    # [1] 137 291
    
    ## *C3.6: plot the scatter with color and centroids ----
    centroids2 <- kmeans2$centers # extract centroid locations for plotting
    png(filename = file.path(saveloc, "clusteridentity_plots", paste0("KmeansClusters_", i,"_",spp,"_",today(), "_cluster_identities.png")), width = 4*480, height = 4*480, units = "px", pointsize = 4*12, bg = "white", res = NA, family = "", type = "cairo-png")
    plot(df_i$StepLengthBLlog1pST, df_i$TurnAngleRelDeglog1pST,
         col = ifelse(df_i$kmeansCharacter == "resident", "skyblue4", "darkgoldenrod2"),
         pch = 16,
         xlab = expression(bold("Standardised log(step length in body lengths)")),
         ylab = expression(bold("Standardised log(abs(turning angle relative °)")))
    points(centroids2, col = "black", pch = 10, cex = 2, lwd = 4)
    dev.off()
    
    ## *C3.7: reverse transform -----
    #df_i$StepLengthBLlog1p <-  expm1(df_i$StepLengthBLlog1p + logmean) # old, original script
    df_i$StepLengthBL_rt <-  expm1(rt_STrange(hammers_i$StepLengthBLlog1p, df_i$StepLengthBLlog1pST))
    df_i$TurnAngleRelDeg_rt <-  expm1(rt_STrange(hammers_i$TurnAngleRelDeglog1p, df_i$TurnAngleRelDeglog1pST))
    
    ## *C3.8: save metadata clusterinfo.csv ----
    clusterinfoadd <- data.frame(nClustersTolerance1 = kstar,
                                 nClustersTolerance2 = kstar2,
                                 # actual data
                                 TransitClusterStepLengthBLlog1pSTMean = mean(df_i[which(df_i$kmeansCharacter == "transit"), "StepLengthBLlog1pST"], na.rm = T),
                                 TransitClusterTurnAngleAbsRelDeglog1pSTMean = mean(abs(df_i[which(df_i$kmeansCharacter == "transit"), "TurnAngleRelDeglog1pST"]), na.rm = T),
                                 ResidentClusterStepLengthBLlog1pSTMean = mean(df_i[which(df_i$kmeansCharacter == "resident"), "StepLengthBLlog1pST"], na.rm = T),
                                 ResidentClusterTurnAngleAbsRelDeglog1pSTMean = mean(abs(df_i[which(df_i$kmeansCharacter == "resident"), "TurnAngleRelDeglog1pST"]), na.rm = T),
                                 # reverse transformed data
                                 TransitClusterStepLengthBL_rt_Mean = mean(df_i[which(df_i$kmeansCharacter == "transit"), "StepLengthBL_rt"], na.rm = T),
                                 TransitClusterTurnAngleAbsRelDeg_rt_Mean = mean(abs(df_i[which(df_i$kmeansCharacter == "transit"), "TurnAngleRelDeg_rt"]), na.rm = T),
                                 ResidentClusterStepLengthBL_rt_Mean = mean(df_i[which(df_i$kmeansCharacter == "resident"), "StepLengthBL_rt"], na.rm = T),
                                 ResidentClusterTurnAngleAbsRelDeg_rt_Mean = mean(abs(df_i[which(df_i$kmeansCharacter == "resident"), "TurnAngleRelDeg_rt"]), na.rm = T),
                                 stringsAsFactors = FALSE)
    if (!exists("clusterinfo")) { #if this is the first i, create sensordates object
      clusterinfo <- clusterinfoadd
    } else { #else add to existing object
      clusterinfo <- rbind(clusterinfo,
                           clusterinfoadd,
                           stringsAsFactors = FALSE) # add to existing file
    } #close if else
    
    setDT(df_i) # convert to data.table without copy
    # join and update "df_i" by reference, i.e. without copy
    # https://stackoverflow.com/questions/44930149/replace-a-subset-of-a-data-frame-with-dplyr-join-operations
    hammers_i[df_i, on = c("date", "shark"), kmeans2cluster := i.kmeans2cluster]
    hammers_i[df_i, on = c("date", "shark"), kmeansCharacter := i.kmeansCharacter]
    hammers_i[df_i, on = c("date", "shark"), kmeansBinary := i.kmeansBinary]
  } #close i
  
  clusterinfo <- round(clusterinfo, 3)
  clusterinfo <- bind_cols(id = fishlist, clusterinfo)
  
  # write.csv(x = clusterinfo, file = paste0(saveloc, "/", today(), "_KmeansClusterinfo.csv"), row.names = FALSE) #SD
  write.csv(x = clusterinfo, file = file.path(saveloc, paste0(today(), "_individual_KmeansClusterinfo_",spp,".csv")), row.names = FALSE) #VH
} else {# close if (!all(is.na(alldaily$lat)))
  print("all new days missing latitude data, can't get external data, nothing to do")
}

# Show how many clusters were chosen most commonly
clustersvec <- c(clusterinfo$nClustersTolerance1, clusterinfo$nClustersTolerance2)
clustersvec <- clustersvec[!is.infinite(clustersvec)]
clustersvec %>%
  as_tibble() %>%
  group_by(value) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
# for all log transformed and all range standardies predicted movement data at 12h time interval
# SMOKARRAN
# # A tibble: 3 × 2
# value     n
# <dbl> <int>
#   1     3    32
# 2     2    20
# 3     1    12
# SLEWINI

#S SZYGAENA

## remove clusterinfo and clusterinfo add if you need to rerun kmeans for other species
rm(clusterinfo); rm(clusterinfoadd)

# could do this more systematically. Could also weight the Tolerance1/2 differently? Leave it for now, perfect enemy of good.
# See L325 TOT make centres dynamic

# C4: conclude kmeans calculations ----

setDF(hammers_i)
# hammers$StepLengthBLlog1p <-  expm1(hammers$StepLengthBLlog1p + logmean) # has been taken care of upstream
# saveRDS(object = hammers, file = paste0(saveloc, "/Hammers_KMeans.Rds")) #SD
## if with fitted data at original observation times
# dir.create(paste0(saveloc,"fitted"))
saveRDS(object = hammers_i, file = file.path(saveloc, paste0(spp,"_KMeans.Rds"))) #VH

## summarise steplength and TAs by cluster
mean_SlTa_trares <- hammers_i %>%
  dplyr::group_by(
    kmeansCharacter
  ) %>%
  dplyr::summarise(
    # SL log and ST
    meanSL_log1p_ST = mean(StepLengthBLlog1pST),
    sdSL_log1p_ST = sd(StepLengthBLlog1pST),
    # SL reverse T
    # meanSL_rt = mean(StepLengthBL_rt),
    # sdSL_rt = sd(StepLengthBL_rt),
    # SL BL raw
    meanSL_BL = mean(StepLengthBL),
    sdSL_BL = sd(StepLengthBL),
    # SL KM raw
    meanSL_KM = mean(StepLengthKm),
    sdSL_KM = sd(StepLengthKm),
    # TA log and ST
    meanTA_log_RST = mean(TurnAngleRelDeglog1pST),
    sdTA_log_RST = sd(TurnAngleRelDeglog1pST),
    # TA reverse T
    # meanTA_rt = mean(TurnAngleRelDeg_rt),
    # sdTA_rt = sd(TurnAngleRelDeg_rt),
    # TA rel degree raw
    meanTA = mean(abs(TurnAngleRelDeg)),
    sdTA = sd(abs(TurnAngleRelDeg))
  );mean_SlTa_trares

## save mean values
write.csv(mean_SlTa_trares, file.path(saveloc, paste0("Data_kmeans_cluster_summary_global_all_",spp,".csv")), row.names = F) #VH

## read clusterinfo for summaries if not yet loaded
clusters_info <- read.csv(file.path(saveloc, paste0(today(), "_individual_KmeansClusterinfo_",spp,".csv"))) #VH

## summarise
summary_all <- clusters_info %>% dplyr::summarise(
  ## log & range standardised
  meanSL_resident_log_RST = mean(ResidentClusterStepLengthBLlog1pSTMean),
  sdSL_resident_log_RST = sd(ResidentClusterStepLengthBLlog1pSTMean),
  meanTA_resident_log_RST = mean(ResidentClusterTurnAngleAbsRelDeglog1pSTMean),
  sdTA_resident_log_RST = sd(ResidentClusterTurnAngleAbsRelDeglog1pSTMean),
  meanSL_transit_log_RST = mean(TransitClusterStepLengthBLlog1pSTMean),
  sdSL_transit_log_RST = sd(TransitClusterStepLengthBLlog1pSTMean),
  meanTA_transit_log_RST = mean(TransitClusterTurnAngleAbsRelDeglog1pSTMean),
  sdTA_transit_log_RST = sd(TransitClusterTurnAngleAbsRelDeglog1pSTMean),
  ## backtransformed values
  meanSL_resident_rt = mean(ResidentClusterStepLengthBL_rt_Mean),
  sdSL_resident_rt = sd(ResidentClusterStepLengthBL_rt_Mean),
  meanTA_resident_rt = mean(ResidentClusterTurnAngleAbsRelDeg_rt_Mean),
  sdTA_resident_rt = sd(ResidentClusterTurnAngleAbsRelDeg_rt_Mean),
  meanSL_transit_rt = mean(TransitClusterStepLengthBL_rt_Mean),
  sdSL_transit_rt = sd(TransitClusterStepLengthBL_rt_Mean),
  meanTA_transit_rt = mean(TransitClusterTurnAngleAbsRelDeg_rt_Mean),
  sdTA_transit_rt = sd(TransitClusterTurnAngleAbsRelDeg_rt_Mean)
) %>%
  tidyr::pivot_longer(cols = everything(),  # Pivot all columns
                      names_to = "Variable",  # First column will store the column names
                      values_to = "Value" # Second column will store the corresponding values
  ); summary_all

write.csv(clusters_info, file.path(saveloc, paste0(today(), "_overall_Kmeans_cluster_SL_TA_summaries_",spp,".csv")), row.names = F) #VH

# rm(list = ls()) #remove all objects
# beep(8) #notify completion
# lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
# invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only = TRUE, unload = TRUE, force = TRUE))

### ....................................................................................................
### [D] 2D Barplot maps for kmeans ----
### ....................................................................................................

# D0: re-import data if needed ----
hammers_i <- readRDS(file = file.path(saveloc, paste0(spp,"_KMeans.Rds")))

# D1: create basemap and prepare plotting ----

## only run if you need new basemap - takes a long time!
cropmap <- gbm.auto::gbm.basemap(grids = hammers_i,
                                 gridslat = 5,
                                 gridslon = 2,
                                 savedir = saveloc)

## if basemap already downloaded
# cropmap <- sf::st_read(saveloc, "CroppedMap", "Crop_Map.shp")

## bathymetry map
bathysavepath <- file.path(saveloc,"getNOAAbathy")
# map2dbpSaveloc <- paste0(saveloc, "/2DbarplotMap/") #SD
map2dbpSaveloc <- file.path(saveloc, "2d_barplots")

# D3: create 2d barplot ----

for (i in c(0.25, 0.375, 0.5, 0.75, 0.875, 1)) {
  barplot2dMap(x = hammers_i |> tidyr::drop_na(kmeansCharacter),
               groupcol = "kmeansCharacter", # col name in x
               baseplot = cropmap,
               bathysavepath = bathysavepath,
               cellsize = c(i, i),
               # mycolours = c("#015475", "#F6D408"),
               mycolours = behav_col,
               legendloc = "topright",
               legendtitle = NULL,
               saveloc = paste0(map2dbpSaveloc,"/"),
               plotname = paste0(lubridate::today(), "_", spp, "_2DBarplot_Count_", i, "deg"))
}

# D4: plot kmeans clusters by individual by detection ----

# dir.create(file.path(saveloc, "individualPlots"))

# *D4.1: Prep background shapefile for individual 2d barplots

bg = ne_countries(scale = 10, continent = 'north america', returnclass = "sf") # needs to be adjusted depending where your study site is

### Prep EEZ shapefiles
# bah_eez <- read_sf("//Sharktank/Science/Data_raw/Shapefiles/Bahamas/Bahamas_EEZ_boundary.shp")
# st_crs(bah_eez)
# bah_eez <- st_transform(bah_eez, st_crs = proj4string(bathyR))
# # bah_eez_plot <- fortify(bah_eez)
# ## only if you want to use the eez shapefile and represent the entire area
# xlim_eez <- c(min(hammers$lon), -70.5105)
# ylim_eez <- c(20.3735, max(hammers$lat))

# *D4.2: create individual plots for each ptt id and save them ----

for (thisshark in unique(hammers_i$shark)){
  # filter your shark
  kplot_df <- hammers_i %>% dplyr::filter(shark == thisshark)
  kplot_df <- kplot_df[!is.na(kplot_df$kmeansCharacter),] # 2 clusters
  # kplot_df <- kplot_df[!is.na(kplot_df$kmeans2cluster),] # 3+ clusters
  
  #define factors
  kplot_df$kmeansCharacter <- as.factor(kplot_df$kmeansCharacter) # 2 clusters
  # kplot_df$kmeansCharacter <- as.factor(kplot_df$kmeans2cluster) # 3+ clusters
  ## create plot with dark themed background
  #p = basemap(dt, bathymetry = T, expand.factor = 1.2) + # for bathymetry with ggOceanMaps package
  p <- ggplot() +
    
    # lines and points
    geom_path(data = kplot_df,
              aes(x=lon,y=lat),
              alpha = 1, linewidth = 0.5)+
    geom_point(data = kplot_df,
               aes(x=lon,y=lat, group = kmeansCharacter, fill = kmeansCharacter, shape = kmeansCharacter), # 2 clusters
               # aes(x=lon,y=lat, group = kmeans2cluster, fill = kmeans2cluster, shape = kmeans2cluster), # 3+ clusters
               alpha = 0.8, size = 3.25, color = "black")+
    # if you create individually tailored plots, you will undoubtedly run into overlapping axis labels
    # depending on the range of your lon/lat values for a specific animal
    # here we can adjust the nr. breaks for each axis
    # breaks_x <- seq(min(kplot_df$lon), max(kplot_df$lon), length.out = 3)
    scale_x_continuous(
      breaks = seq(min(kplot_df$lon), max(kplot_df$lon), length.out = 3),
      labels = function(x) paste0(scales::number_format(accuracy = 0.1)(x), "°W")  # Append 'W' to each label as we lose the °W info when we scale the axis
    ) +
    
    # basemap
    geom_sf(data = bg, color = "black")+ # color is for border of geom object
    ## comment out the next three lines if you want to plot maps tailored to the movement extent of individual sharks
    ## and use L1187-1189 instead)
    # coord_sf(xlim = range(hammers$lon, na.rm = TRUE),
    #          ylim = range(hammers$lat , na.rm = TRUE),
    #          expand = T)+
    
    # bahamas eez shapefile
    # geom_sf(data = bah_eez, colour = "black", fill = NA, linewidth = .25) +
    ## for individuall tailored plots comment out L1184-1186 and use L1187-1189 instead
    # coord_sf(xlim = xlim_eez,
    #          ylim = ylim_eez+.25,
    #          expand = T)+
    coord_sf(xlim = c(min(kplot_df$lon)-0.2, max(kplot_df$lon)+0.2),
             ylim = c(min(kplot_df$lat)-0.2, max(kplot_df$lat)+0.2),
             expand = T)+
    
    # formatting
    # labs(x=NULL, y=NULL,
    #      fill = 'kmeans cluster',
    #      color = 'kmeans cluster')+
    scale_shape_manual(values = c(22,24))+ # 2 clusters
    # scale_shape_manual(values = c(22,24, 25))+ # 3+ clusters
    scale_color_manual(values = behav_col) +
    scale_fill_manual(values = behav_col) +
    
    # theme
    theme_light()+
    #theme(panel.background = element_rect(fill = "gray26", linewidth = 0.5, linetype = "solid", color = "black")) +
    theme(panel.grid = element_blank(), # remove grid lines
          legend.title = element_blank(), # remove legend title
          plot.title = element_text(face = "bold", size = 11.5),
          axis.text = element_text(size = 8)) +  # adjust size of lat/lon axis labels
    labs(x = "Longitude", y = "Latitude") +
    ggtitle(paste0("Kmeans derived movement states for shark id ", thisshark))
  p
  
  ggsave(file.path(saveloc, "individualPlots", paste0(today(), "_", spp, "_kmeans_",thisshark,".tiff")), width = 15, height = 10, units = "cm", dpi = 300)
}

# END OF CODE ----

# TODO: see below 
# TODO1: add environmental data to kmeans clusters and see if behavioural states are dependent on covariates
# TODO2: summarise behavioral states, swithces etc. by season or some sort of time frame
# TODO3: loop over multiple species rather than having to rerun individually

### ....................................................................................................
### [E] Add environmental data to kmeans output ----
### ....................................................................................................

### TBC THIS STILL NEEDS TO BE CODED!!!! NOT USABLE YET

# # We can add environmental data to our kmeans behavioural clusters for further analysis, if e.g.
# # different behaviours are correlated with different habitats, conditions, etc.
# 
# # Pre-step: import all dataframe with clusterinfo as we analysed separately by sex and species
# 
# ## Smokarran
# df_smok <- readRDS("C:/Users/Vital Heim/switchdrive/Science/Data/PhD_Chapter3_Output_files_(kmeans)/species_specific/Data_Kmeans_S.mokarran.R")
# write.csv(df_smok, paste0(saveloc, "Data_Kmeans_S.mokarran.csv"), row.names = F)
# ## Slewini
# df_slew <- readRDS("C:/Users/Vital Heim/switchdrive/Science/Data/PhD_Chapter3_Output_files_(kmeans)/species_specific/Data_Kmeans_S.lewini.R")
# write.csv(df_slew, paste0(saveloc, "Data_Kmeans_S.lewini.csv"), row.names = F)
# 
# ## combine
# df_all <- rbind(df_smok, df_slew)
# 
# ## all data
# # df_all <- list.files(path = "C:/Users/Vital Heim/switchdrive/Science/Data/PhD_Chapter3_Output_files_(kmeans)/species_specific/",pattern = "Data_K*.R") %>%
# #   purrr::map_dfr(readRDS)
# 
# # C1: import environmental data ----
# 
# # *C1.1: Bathymetry ----
# 
# bathy <- raster("C:/Users/Vital Heim/switchdrive/Science/Data/Bathymetry_maps_GMRT/GMRTv4_0_20221013topo_for_Sphyrna_SPOT_tracks.grd")
# 
# ## check crs
# proj4string(bathy)
# # > proj4string(bathy)
# # [1] "+proj=longlat +datum=WGS84 +no_defs"
# 
# ## find min and max values
# setMinMax(bathy)
# 
# ## we are working with ocean depth data, i.e. anyhting above 0 meter elevation we dont really need and can set to NA
# bathy <- clamp(bathy, upper = 0, value = F)
# 
# ## Extent of the raster
# ext(bathy)
# # > ext(bathy)
# # SpatExtent : -99.10546875, -65.513671875, 22.5693803462162, 42.9017222226332 (xmin, xmax, ymin, ymax)
# # > min(hammers$lat);max(hammers$lat);min(hammers$lon);max(hammers$lon)
# # [1] 23.42111
# # [1] 36.44029
# # [1] -94.75546
# # [1] -74.70941
# 
# # > bbox(bathy)
# # min       max
# # s1 -99.10547 -65.51367
# # s2  22.56938  42.90172 # make sure these values encompass the values of the hammers lon/lat
# 
# ## Visualise
# plot(bathy)
# # map("worldHires", col = "grey90", border = "grey50", fill = TRUE, add = TRUE)
# 
# # *C1.2: SST ----
# 
# # Available at url: https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.html
# 
# ## Extract SST - define parameters
# xpos <- hammers$lon # longitude across all tracks
# ypos <- hammers$lat # latitude across all tracks
# tpos <- hammers$date # date/time across all tracks/observations, you define if you extract daily, weekly, or monthly values
# 
# sstInfo <- rerddap::info("jplMURSST41") # this is your dataset ID from the website. Depending on daily, monthly, etc. datasets the id is different
# 
# ## Extract SST - extract analysed STT
# tic()
# sst <- rxtracto(sstInfo, parameter = 'analysed_sst',
#                 xcoord = xpos, ycoord = ypos, tcoord = tpos,
#                 #xlen = 0.01, ylen = 0.01, 
#                 progress_bar = TRUE)
# #sst2 <- griddap(murSST <- griddap(sstInfo, latitude = ypos, longitude = xpos, time = tpos, stride = c(1,0.2,0.2), fields = 'analysed_sst'))
# toc()  #takes 21 min to run - save and export next time
# 
# ## Extract SST - save the data
# saveRDS(sst, file = paste0(saveloc, "Extracted_SST_ERDDAP_jplMURSST41_daily_",spp,".RData"))
# 
# ## Extract SST - import previously extracted data
# #sst <- readDRS("C:/Users/Vital Heim/switchdrive/Science/Projects_and_Manuscripts/Bahamas_Hammerheads_2022/OutputData/Extracted_SST_ERDDAP_jplMURSST41_daily_BAH_Smok_2017-11-27_2021-11-21.RData")
# 
# ## Visualise SST
# plotTrack(sst,xpos,ypos,tpos,plotColor = "thermal", size = 2)
# 
# # *C1.3: Chlorophyll a ----
# 
# # # Available at url: https://coastwatch.pfeg.noaa.gov/erddap/griddap/nesdisVHNSQchlaDaily.html
# # 
# # ## Extract chlorophyll a - define parameters
# # xpos <- DET$lon # longitude across all tracks
# # ypos <- DET$lat # latitude across all tracks
# # tpos <- DET$date # date/time across all tracks/observations, you define if you extract daily, weekly, or monthly values
# # 
# # chlaInfo <- rerddap::info("erdMH1chla1day") # this is your dataset ID from the website. Depending on daily, monthly, etc. datasets the id is different
# # 
# # ## Extract chlorophyll a - extract analysed ChlA values
# # tic()
# # chla <- rxtracto(chlaInfo, parameter = 'chlorophyll',
# #                  xcoord = xpos, ycoord = ypos, tcoord = tpos,
# #                  #xlen = 0.2, ylen = 0.2, 
# #                  progress_bar = TRUE)
# # toc()  #takes 15 min to run - save and export next time
# # 
# # ## Extract chlorophyll a - save the data
# # saveRDS(chla, file = paste0(saveloc, "Extracted_ChlA_ERDDAP_erdMH1chla1day_daily_",spp,".RData"))
# # 
# # ## Extract chlorophyll a - import previously extracted data
# # #sst <- readDRS("C:/Users/Vital Heim/switchdrive/Science/Projects_and_Manuscripts/Bahamas_Hammerheads_2022/OutputData/Extracted_SST_ERDDAP_jplMURSST41_daily_BAH_Smok_2017-11-27_2021-11-21.RData")
# # 
# # ## Visualise chlorophyll a
# # plotTrack(chla,xpos,ypos,tpos,plotColor = colors$chlorophyll, size = 2)
# 
# # C2: append your environmental data to your movement data ----
# 
# # *C2.1: Bathymetry ---- 
# ## To append the environmental data such as bathymetry to our data, we need to create a shapefile of
# ## our locations
# 
# ## Define the crs of your spatial data
# mov_crs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# 
# ### with raster package
# ### Define your SpatialPointsDataFrame
# mov <- SpatialPointsDataFrame(hammers[,c(3,6)], hammers, proj4string = mov_crs)
# #plot(movR)
# mov_depths <- extract(bathy, mov)
# 
# ## add to your dataframe
# ### raster package
# hammers$depth <-mov_depths
# 
# ## classify depth as nearshore (>= - 30) and ofsshore (< -30m)
# 
# hammers <- hammers %>%
#   mutate(
#     shore = cut(depth,
#                 breaks = c(0,-30,-Inf),
#                 labels = c("off", "near"))
#   )
# 
# # # *C1.2: SST
# # 
# # DET$sst <- sst$`mean analysed_sst`
# # 
# # ## make a year-day and hour column and add it to the dataset so you can check if there is correlation between them and 
# # ## SST in case you want to use them as covariates
# # DET$yday <- yday(DET$date)
# # DET$hr <- hour(DET$date)
# # 
# # ## Is there a relationship between SST and yday?
# # plot(DET$yday, DET$sst)
# # #a non-linear pattern does appear to be present; let's stick w/ SST as driver and hr as a cyclical covariate
# # 
# # # *C1.3: Chlorophyll a
# # 
# # DET$chla <- chla$`mean chlorophyll`
# # 
# # 
# # 
# # rm(list = ls()) #remove all objects
# # 
# # beep(8) #notify completion
# # 
# # lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
# # invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only = TRUE, unload = TRUE, force = TRUE))
# 
# # C3: calculate proprtion of clusters within environmental conditions ----
# 
# ## depth/shore distance
# proportions <- hammers %>%
#   dplyr::group_by(
#     species,
#     sex,
#     kmeans2cluster,
#     shore
#   ) %>%
#   dplyr::summarise(
#     n = n()
#   ) %>%
#   dplyr::mutate(
#     frequency = n/sum(n)
#   )
# 
# ## save it
# write.csv(proportions, paste0(saveloc, "Proportions_of_detections_near_and_offshore_by_movement_states_", spp,".csv"), row.names = F)

### ....................................................................................................
### [F] Assign seasons and short summaries ----
### ....................................................................................................

### TBC - NEEDS TO BE CODED STILL!!! DONT USE YET

# F1: assign seasons ----

# ## Define seasons by month numbers (based on definitions in Kroetz et al. 2021)
# 
# swinter <- c(12,1,2) # december, january, february
# sspring <- c(3,4,5) # march, april, may
# ssummer <- c(6,7,8) # june, july, august
# sautumn <- c(9,10,11) # september, october, november
# 
# ## add seasons to kmeans df
# df_all <- df_all %>%
#   dplyr::mutate( #add a numeric month variable
#     month = as.numeric(format(date, format = "%m"))
#   ) %>%
#   dplyr::mutate( # create a new column for seasons based on month number
#     ., season = with(., case_when(
#       (month %in% swinter)  ~ "winter",
#       (month %in% sspring)   ~ "spring",
#       (month %in% ssummer)   ~ "summer",
#       (month %in% sautumn) ~ "autumn",
#       is.na(month) ~ "noseason",
#       TRUE ~ "noseason"
#     ))
#   ) %>%
#   dplyr::select( # cleanup
#     species,
#     shark,
#     sex,
#     date,
#     month,
#     season,
#     lon,
#     lat,
#     FishLengthCm,
#     li5day,
#     StepLengthKm,
#     StepLengthBL,
#     StepLengthBLlog1p,
#     StepLengthBLlog1pST,
#     TurnAngleAzimDeg,
#     TurnAngleRelDeg,
#     kmeans2cluster,
#     kmeansBinary,
#     kmeansCharacter
#   )
# 
# write.csv(df_all, paste0(saveloc, "Data_kmeans_all_tidy.csv"))
# saveRDS(df_all, paste0(saveloc, "Data_kmeans_all_tidy.R"))
# 
# # D2: create seasonal dataframes ----
# 
# ## S.mokarran
# ### Females
# k_spring_smokF <- df_all %>%
#   dplyr::filter(
#     season %in% c("spring"),
#     species %in% c("S.mokarran"),
#     sex %in% c("F")
#   ); write.csv(k_spring_smokF, paste0(saveloc, "Data_kmeans_Smok_spring_tidy_F.csv"))
# 
# k_summer_smokF <- df_all %>%
#   dplyr::filter(
#     season %in% c("summer"),
#     species %in% c("S.mokarran"),
#     sex %in% c("F")
#   ); write.csv(k_summer_smokF, paste0(saveloc, "Data_kmeans_Smok_summer_tidy_F.csv"))
# 
# k_autumn_smokF <- df_all %>%
#   dplyr::filter(
#     season %in% c("autumn"),
#     species %in% c("S.mokarran"),
#     sex %in% c("F")
#   ); write.csv(k_autumn_smokF, paste0(saveloc, "Data_kmeans_Smok_autumn_tidy_F.csv"))
# 
# k_winter_smokF <- df_all %>%
#   dplyr::filter(
#     season %in% c("winter"),
#     species %in% c("S.mokarran"),
#     sex %in% c("F")
#   ); write.csv(k_winter_smokF, paste0(saveloc, "Data_kmeans_Smok_winter_tidy_F.csv"))
# 
# ### Males
# k_spring_smokM <- df_all %>%
#   dplyr::filter(
#     season %in% c("spring"),
#     species %in% c("S.mokarran"),
#     sex %in% c("M")
#   ); write.csv(k_spring_smokM, paste0(saveloc, "Data_kmeans_Smok_spring_tidy_M.csv"))
# 
# k_summer_smokM <- df_all %>%
#   dplyr::filter(
#     season %in% c("summer"),
#     species %in% c("S.mokarran"),
#     sex %in% c("M")
#   ); write.csv(k_summer_smokM, paste0(saveloc, "Data_kmeans_Smok_summer_tidy_M.csv"))
# 
# k_autumn_smokM <- df_all %>%
#   dplyr::filter(
#     season %in% c("autumn"),
#     species %in% c("S.mokarran"),
#     sex %in% c("M")
#   ); write.csv(k_autumn_smokM, paste0(saveloc, "Data_kmeans_Smok_autumn_tidy_M.csv"))
# 
# k_winter_smokM <- df_all %>%
#   dplyr::filter(
#     season %in% c("winter"),
#     species %in% c("S.mokarran"),
#     sex %in% c("M")
#   ); write.csv(k_winter_smokM, paste0(saveloc, "Data_kmeans_Smok_winter_tidy_M.csv"))
# 
# ## S.lewini
# ### Females
# k_spring_slewF<- df_all %>%
#   dplyr::filter(
#     season %in% c("spring"),
#     species %in% c("S.lewini"),
#     sex %in% c("F")
#   ); write.csv(k_spring_slewF, paste0(saveloc, "Data_kmeans_Slew_spring_tidy_F.csv"))
# 
# k_summer_slewF <- df_all %>%
#   dplyr::filter(
#     season %in% c("summer"),
#     species %in% c("S.lewini"),
#     sex %in% c("F")
#   ); write.csv(k_summer_slewF, paste0(saveloc, "Data_kmeans_Slew_summer_tidy_F.csv"))
# 
# k_autumn_slewF <- df_all %>%
#   dplyr::filter(
#     season %in% c("autumn"),
#     species %in% c("S.lewini"),
#     sex %in% c("F")
#   ); write.csv(k_autumn_slewF, paste0(saveloc, "Data_kmeans_Slew_autumn_tidy_F.csv"))
# 
# k_winter_slewF <- df_all %>%
#   dplyr::filter(
#     season %in% c("winter"),
#     species %in% c("S.lewini"),
#     sex %in% c("F")
#   ); write.csv(k_winter_slewF, paste0(saveloc, "Data_kmeans_Slew_winter_tidy_F.csv"))
# 
# ### Males
# k_spring_slewM<- df_all %>%
#   dplyr::filter(
#     season %in% c("spring"),
#     species %in% c("S.lewini"),
#     sex %in% c("M")
#   ); write.csv(k_spring_slewM, paste0(saveloc, "Data_kmeans_Slew_spring_tidy_M.csv"))
# 
# k_summer_slewM <- df_all %>%
#   dplyr::filter(
#     season %in% c("summer"),
#     species %in% c("S.lewini"),
#     sex %in% c("M")
#   ); write.csv(k_summer_slewM, paste0(saveloc, "Data_kmeans_Slew_summer_tidy_M.csv"))
# 
# k_autumn_slewM <- df_all %>%
#   dplyr::filter(
#     season %in% c("autumn"),
#     species %in% c("S.lewini"),
#     sex %in% c("M")
#   ); write.csv(k_autumn_slewM, paste0(saveloc, "Data_kmeans_Slew_autumn_tidy_M.csv"))
# 
# k_winter_slewM <- df_all %>%
#   dplyr::filter(
#     season %in% c("winter"),
#     species %in% c("S.lewini"),
#     sex %in% c("M")
#   ); write.csv(k_winter_slewM, paste0(saveloc, "Data_kmeans_Slew_winter_tidy_M.csv"))
# 
# # D3: calculate proportions of movement in different states ----
# 
# ## by species, sex and season
# proportions <- df_all %>%
#   dplyr::group_by(
#     species,
#     season,
#     sex,
#     kmeansCharacter
#   ) %>%
#   dplyr::summarise(
#     n = n()
#   ) %>%
#   dplyr::mutate(
#     frequency = n/sum(n)
#   );proportions
# 
# ## now drop NA values coming from NA SL/TA
# proportions_nona <- df_all %>%
#   dplyr::filter(
#     !is.na(kmeansCharacter)
#   ) %>%
#   dplyr::group_by(
#     species,
#     season,
#     sex,
#     kmeansCharacter
#   ) %>%
#   dplyr::summarise(
#     n = n()
#   ) %>%
#   dplyr::mutate(
#     frequency = n/sum(n)
#   );proportions_nona
# 
# write.csv(proportions_nona, paste0(saveloc, "Data_proportions_movement_behaviours_by_species_sex_season.csv"))
# 
# # D4: summarise SL and TA by cluster
# 
# clusters_smok <- read.csv(paste0(saveloc, "2023-12-06_KmeansClusterinfo_S.mokarran.csv"))
# clusters_slew <- read.csv(paste0(saveloc, "2023-12-06_KmeansClusterinfo_S.lewini.csv"))
# 
# ## summarise
# ### Smokarran
# clusters_smok %>% dplyr::summarise(
#   meanSL_resident = mean(ResidentClusterStepLengthBLlog1pMean),
#   sdSL_resident = sd(ResidentClusterStepLengthBLlog1pMean),
#   meanTA_resident = mean(ResidentClusterTurnAngleRelDegAbsMean),
#   sdTA_resident = sd(ResidentClusterTurnAngleRelDegAbsMean),
#   meanSL_transit = mean(TransitClusterStepLengthBLlog1pMean),
#   sdSL_transit = sd(TransitClusterStepLengthBLlog1pMean),
#   meanTA_transit = mean(TransitClusterTurnAngleRelDegAbsMean),
#   sdTA_transit = sd(TransitClusterTurnAngleRelDegAbsMean)
# )
# # meanSL_resident sdSL_resident meanTA_resident sdTA_resident meanSL_transit sdSL_transit
# #       5336.738      2883.732        124.1423      18.33783        5347.05     2136.681
# # meanTA_transit sdTA_transit
# #       8.74231     6.743125
# 
# ### Slewini
# clusters_slew %>% dplyr::summarise(
#   meanSL_resident = mean(ResidentClusterStepLengthBLlog1pMean),
#   sdSL_resident = sd(ResidentClusterStepLengthBLlog1pMean),
#   meanTA_resident = mean(ResidentClusterTurnAngleRelDegAbsMean),
#   sdTA_resident = sd(ResidentClusterTurnAngleRelDegAbsMean),
#   meanSL_transit = mean(TransitClusterStepLengthBLlog1pMean),
#   sdSL_transit = sd(TransitClusterStepLengthBLlog1pMean),
#   meanTA_transit = mean(TransitClusterTurnAngleRelDegAbsMean),
#   sdTA_transit = sd(TransitClusterTurnAngleRelDegAbsMean)
# )
# # meanSL_resident sdSL_resident meanTA_resident sdTA_resident meanSL_transit sdSL_transit
# #       5617.829      2125.107        129.8071      13.58157       6512.136     2735.414
# # meanTA_transit sdTA_transit
# #          25.3     9.405318
