# TOP OF CODE ----
### ====================================================================================================
### Project:    Large bodied hammerhead complex in the western North Atlantic
### Analysis:   Animate movement tracks of tagged animals
### Script:     ~2025_NWA_hammerhead_complex/R/Rscript8_Animate_movement_data.R
### Author:     Vital Heim
### Version:    1.0
### ====================================================================================================

### ....................................................................................................
### Content: this script contains the code to visualse/animate time series location data
### such as animal tracking data.
### The script contains the following options to animate movement data
###
### Section D: Animate movement of a single animal in a fixed extent
### Section E: Animate movement of multiple animals with a dynamic extent that follow each track over
###            time.
###            
### Please see TODO list at end of script for open issues.
### ....................................................................................................

### ....................................................................................................
### [A] Ready environment, load packages ----
### ....................................................................................................

# A1: clear memory ----

rm(list = ls())

# A2: isntall and load necessary packages ----

## if first time
# install.packages("tidyverse")
# install.packages("magrittr")
# install.packages("rnaturalearth") # if first time
# install.packages("rnaturalearthdata") # if first time
# install.packages("gganimate") # if first time
# install.packages("gifski")
# install.packages("ncdf4")
# install.packages("stars")
# install.packages("raster")
# install.packages("terra")
# install.packages("tidyterra")
# install.packages("scales")
# install.packages("ggOceanMapsData")
# devtools::install_github("MikkoVihtakari/ggOceanMapsData") # required by ggOceanMaps
# install.packages("ggOceanMaps")
# install.packages("ggnewscale")
# install.packages("sf")
# install.packages("ggimage")
# remotes::install_github("SimonDedman/gbm.auto")
# install.packages("gbm.auto")
# install.packages("ggnewscale")
# install.packages("beepr")

## load packages
library(tidyverse)
library(magrittr)
library(rnaturalearth)
library(rnaturalearthdata)
library(gganimate)
library(gifski)
library(ncdf4)
library(stars)
library(raster)
library(terra)
library(tidyterra)
library(scales)
# library(ggOceanMapsData)
# library(ggOceanMaps")
library(ggnewscale)
library(sf)
# library(gbm.auto)
library(beepr)
# library("ggimage") # NOT READY YET

# A3: Specify needed functions ----

# *A3.1: function to synchronise starting years ----
shift_dates_custom_start <- function(x, id_col, date_col, start_year = 1980) {
  x %>%
    arrange({{ id_col }}, {{ date_col }}) %>%
    group_by({{ id_col }}) %>%
    mutate(
      # Get the first observation date for each ID
      first_date = first({{ date_col }}),
      
      # Calculate years, months, days, and time components
      year_diff = year({{ date_col }}) - year(first_date),
      
      # Create new date preserving month, day, and time but starting from specified year
      shifted_date = {{ date_col }} + years(start_year - year(first_date))
    ) %>%
    select(-first_date, -year_diff) %>%
    ungroup()
}

# A4: Specify data and saveloc ----

YOUR_IP <- "NA" # add your IP address, server name or similar if you connect via a shared drive

## Input data
dataloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2025_NWA_hammerhead_complex","Data_input", "Telemetry")
bathyloc <- file.path("/",YOUR_IP, "Science", "Data_raw", "Bathymetry_maps_GMRT", "GMRTv4_3_0_20250120topo_wider_NWA.grd")
misc_shapefileloc <- file.path("/",YOUR_IP, "Science", "Data_raw", "Shapefiles")

## Output data
saveloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2025_NWA_hammerhead_complex","Data_output", "Telemetry", "Animations")

# A5: Define universal variables (e.g. for plotting) ----

# Plotting variables

## bathymetry colors
shallow <- "#D3E5E8"
deep <- "#2B628B"
depth <- c(deep, shallow)

## speciescolours
slewcol <- "#EDA904"
smokcol <- "#70AB27"
szygcol <- "#E26306"

# A6: Define universal options, variables, etc. (e.g. for plotting) ----

options(timeout = 3000) # manually increase time out threshold (needed when downloading basemap)
options(scipen=999) # so that R doesn't act up for pit numbers  
options(warn=1) #set this to two if you have warnings that need to be addressed as errors via traceback()
options(error = function() beep(9))  # give warning noise if it fails

### ....................................................................................................
### [B] Data import ----
### ....................................................................................................

# B1: Movement data ----

## movement data
hammers <- readRDS(file = file.path(dataloc, "Animations", "Data_aniMotum_CRW_output_entire_track_rerouted_proj_WGS84_converted_with_coord_CIs_with_Argosfilter_data.rds")) |>
  dplyr::rename(
    shark = id
    ) %>%
  dplyr::mutate(
    shark = as.character(shark)
  )

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
  dplyr::mutate(
    shark = as.character(ptt_id)
    )

## combine movement data with metadata
hammers %<>% left_join(meta)  # , by = join_by(shark == id) # doesn't work naming columns, has gotten worse.

## VERY IMPORTANT: The movement data needs to have time stamps in ascending order,
## double check that this is true by ordering the dataframe
hammers %<>%
  dplyr::arrange(
    shark,
    date
  )

# B2: Basemap data and shapefiles ----

# *B2.1: Landmasses basemap ----
## Basemap data for landmasses- HQ
## Option 1 - required downloaded shapefile
### define resolution: 1=c, 2=l,3=i,4=h,5=f, 1:5 increasing quality
# res <- "f"
### read in worldmap
# world <- sf::st_read(dsn = paste0(shapefileloc,"/",res,"/GSHHS_", res, "_L1.shp"), layer = paste0("GSHHS_", res, "_L1"), quiet = TRUE) # read in worldmap

## World shapefile - LQ
world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")

# *B2.2: Bathymetry maps and other rasters ----

bathyR <- raster::raster(bathyloc)
## todo: use rast() rather than raster(), rast() is from the terra package, and terra has more options and will replace raster

#SpatExtent : -84.3068857607634, -74.5949697079866, 22.996942201359, 35.7020417853476 (xmin, xmax, ymin, ymax)
proj4string(bathyR) # from raster package
#> proj4string(bathy)
#[1] "+proj=longlat +datum=WGS84 +no_defs"

## find min and max values
raster::setMinMax(bathyR)
# minValue(bathyR); maxValue(bathyR)
### we are working with ocean depth data, i.e. anyhting above 0 meter elevation we dont really need and can set to NA
bathyR <- clamp(bathyR, upper = 50, useValues = F) # this is really just an aesthethic thing, needs to be adjsuted if you print a legend that goes between 0 to -8000 meters

## Extent of the raster
# bbox(bathyR)
# min       max
# s1 -99.10547 -65.51367
# s2  22.56938  42.90172

## Dimension
# dim(bathyR)
#> dim(bathy)
#[1] 2777 3822    1

## Resolution
# res(bathyR)
#[1] [1] 0.008789062 0.007321693

## Visualise
plot(bathyR)

## make a df needed for later plotting
raster.df <- as.data.frame(bathyR, xy = T)

## prep for TERRA PACKAGE.....
bathyterra <- terra::rast(bathyR)
plot(bathyterra)

# *B2.2: Other shapefiles, e.g. EEZ or closure boundaries ----

### Bahamian EEZ
# bah_eez <- read_sf("//Sharktank/Science/Data_raw/Shapefiles/Bahamas/Bahamas_EEZ_boundary.shp")
# sf::st_crs(bah_eez)
# bah_eez <- sf::st_transform(bah_eez, st_crs = proj4string(bathyR)) # bring to same CRS as bathymetry raster

### US State boundaries
usstates <- sf::st_read(dsn = file.path(misc_shapefileloc, "USA", "US_State_Boundaries", "US_State_Boundaries.shp"), quiet = TRUE) # read in worldmap

### Federal vs. state waters USA
# NEEDS UPDATE

### US EEZ
# NEEDS UPDATE

### ....................................................................................................
### [C] Data housekeeping and preparation ----
### ....................................................................................................

# C1: housekeeping DET dataframe ----

## filter the needed columns
hammers_f <- hammers %>%
  # dplyr::mutate( # define any additional variables that might be handy for plotting
  #   month = format(dateEST, format = "%b", tz = "US/Eastern"), # get rid of year for later season variable definitions
  #   season_bahamas = with(.,case_when(# summer/winter based on temperature by van Zinnicq Bergmann et al. 2022 ; summer =  1st Jun to 30th Nov, winter = 1st Dec to 31st May
  #     month %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct") ~ "wet",
  #     month %in% c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr") ~ "dry",
  #     TRUE ~ "seasonnonexistent")
  #   )
  # ) %>%
  dplyr::select(
    shark,
    date,
    lon,
    lat,
    species
  )

# C2: Process the data for animation ----

## Filter for one individual if you want
# ptt <- "200368"
# ind <- DET %>%
#   dplyr::filter(
#     id == ptt
#   )

# If you do not have standardised time stampes, e.g. 1 loc per day, and across multiple years,
# animating your df might be a bit of a pain as you will have long gaps, lags, etc.
# For now and until a better solution is found, we used indexing of the timestamps to create
# a df where each shark was "tagged" at the same time" and has the same time gap between each location estimate

## We add an index nr. first
det_anim <- hammers_f %>%
  group_by(
    shark
  ) %>%
  dplyr::mutate(
    Index = row_number()
  )

## we also add a shifted data so that all sharks were tagged in the same year
det_anim <- shift_dates_custom_start(det_anim, shark, date, 2020)

# C3: add image as datapoint for later plotting - NOT YET IMPLEMENTED ----
## add an image to the dataframe for later plotting - NOT READY YET
# df_all <- df_all %>%
  # dplyr::mutate(image = "C:/Users/Vital Heim/switchdrive/Science/Projects_and_Manuscripts/Andros_Hammerheads/InputData/hammerhead_shape.png")

### ....................................................................................................
### [D] Plot animated data with fixed plot window for individual animals or species ----
### ....................................................................................................

# D1: filter your dataframe by species if needed ----

## What species
sp.f <- "S.zygaena"

## filter
det_anim_spp <- det_anim %>%
  dplyr::filter(
    species %in% sp.f
  )

# D2: prepare parameters for static map ----

## plot extent
min(det_anim$lon);max(det_anim$lon)
min(det_anim$lat);max(det_anim$lat)

xmin <- ceiling(min(det_anim$lon)-1);xmin
xmax <- ceiling(max(det_anim$lon)+1);xmax
ymin <- floor(min(det_anim$lat)-1);ymin
ymax <- ceiling(max(det_anim$lat)+1);ymax

xlabs = seq(xmin, xmax, 5)
ylabs = seq(ymin, ymax, 5)

### ....................................................................................................
### [D] Plot animated data with fixed plot window for individual animals ----
### ....................................................................................................

## !! NEEDS UDPATE !!##

# D1: static map ----

# det_anim_spp %>% ggplot(aes(x = lon, y = lat)) + geom_path()

## create plot with bathymetry/raster background
#p = basemap(dt, bathymetry = T, expand.factor = 1.2) + # for bathymetry with ggOceanMaps package
st <- ggplot() +
  
  # bathymetry raster
  ## stars package
  # stars::geom_stars(data = stars::st_downsample(bathystar, 50) |> sf::st_transform(proj4string(bathyR)), inherit.aes = FALSE) + # choose your epsg code accordinlgy (here EPSG:3857 is for WGS 84 / Pseudo-Mercator -- Spherical Mercator, Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI) - personally don't like it takes ages even with downsample
  ## raster package
  # ggplot2::geom_raster(data = raster.df , aes(x = x, y = y, fill = layer)) +
  ## terra & tidyterra package - by far performs best
  tidyterra::geom_spatraster(data = bathyterra) +
  ## add a legend/gradient to your bathymetry raster
  ggplot2::scale_fill_gradientn(colors = depth, guide = "none") +   # guide = "none" prevents legend
  # # labs(x=NULL, y=NULL,
  #      fill = 'Depth [m]',
  #      color = 'Depth [m]')+
  # theme(panel.grid = element_blank(), legend.title = element_text(face = "bold")) +
  
  #coord_quickmap() +
  
  # start a new scale
  new_scale_fill() +
  
  # lines and points
  geom_path(data = det_anim_spp, # comment this out if you want to animate a wake effect later
            aes(x=lon,y=lat, group = shark),
            col = ifelse(det_anim_spp$species == "S.lewini", slewcol,
                         ifelse(det_anim_spp$species == "S.mokarran", smokcol,
                                ifelse(det_anim_spp$species == "S.zygaena", szygcol, "black"))),
            alpha = 1, linewidth = 1.5) +
  
  geom_point(data = det_anim_spp,
             aes(x=lon,y=lat, group = shark),
                 # group = seq_along(Index), # somehow needed if you want to keep data points, i.e. if you want keep_last = T in gganimate::animate()
                 fill=ifelse(det_anim_spp$species == "S.lewini", slewcol,
                             ifelse(det_anim_spp$species == "S.mokarran", smokcol,
                                    ifelse(det_anim_spp$species == "S.zygaena", szygcol, "black"))),
                 shape = 21,
             alpha = 1, size = 4) +
  
  # add a symbol rather than a simple data point
  # ggimage::geom_image(aes(image = image), size = 0.1)+ # NOT READY YET
  
  # basemap
  ggplot2::geom_sf(data = world, fill = "gray80", color = "black", size = 1, inherit.aes = F) +
  
  # raster
  #geom_polygon(data = r, aes(x = long, y = lat)) +
  
  # additional shapefiles
  ggplot2::geom_sf(data = usstates, colour = "black", fill = NA, size = .25) +
  
  # define plot limits
  ggplot2::coord_sf(xlim = c(xmin, xmax),
                    ylim = c(ymin+1, ymax),
                    expand = T) +
  
  # plot axis labels
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) +
  
  # formatting
  # scale_fill_viridis_c(option = "inferno")+
  # scale_color_viridis_c(option = "inferno")+
  
  ## continous variables
  # scale_fill_gradientn(colors = seasoncol) +
  # scale_color_gradientn(colors = seasoncol) +
  # scale_size_continuous(range = c(0.1,10))+
  
  ## categorical variables
  ### shape and legend for seasons
  # scale_shape_manual(values = seasonsym,
  #                    name= "Season", # sets legend name
  #                    drop = F) +
  ### fill colour for shark symbols
  # scale_fill_manual(name = "Shark-ID",
  #                   values = all_cols,
  #                   drop = F,
  #                   # guide = "none" # removes legend
  #                   guide = guide_legend(override.aes = list(fill=all_cols, shape = 22, color = "black")) # use this if you use a wake effect in the animation or if you want a dot legend for individuals
  # ) +
  ### colour and legend of shark track
  # scale_color_manual(name = "Shark-ID", # sets legend name
  #                    values = all_cols,
  #                    guide = "none", # use this if you don't want any tracks if you use a wake effect in the animation
  #                    drop = F) +
  
  # Define your theme aesthethics
  theme(panel.grid = element_blank(), #legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .75), # Add a black border around the plot
        # panel.background = element_rect(fill = "#C1D4E1"), # as alternative as GMRT background takes ages to animate
        # panel.background = element_rect(fill = "white"), # as alternative as GMRT background takes ages to animate
        text = element_text(family = "serif", face = "plain"), # all text to Times New Roman look-a-like
        plot.background = element_rect(fill = "white", color = "white"),    # Set overall plot background to white
        axis.text.x = element_text(size = 15), # change the font size of x.axis text
        axis.text.y = element_text(size = 15), # change the font size of y.axis text
        axis.title.x = element_blank(), # removes x axis title (i.e. "lon")
        axis.title.y = element_blank(), # removes y axis title (i.e. "lat")
        legend.title = element_text(size = 12, face = "bold"), # change the font size of the legend titles
        legend.text = element_text(size = 10),# change the font size of legend text
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(face = "bold", size = 14)
  ); st
  
  ## add a title
  # ggtitle("Great hammerheads of Andros Island, The Bahamas") +   # add title and subtitle
  
  ## animate your map
  # gganimate::transition_time(shifted_date) + # cannot be used with geom_path
  # gganimate::transition_reveal(along = shifted_date, keep_last = F) 
  # gganimate::shadow_mark(past = T, future = F)
  # gganimate::shadow_wake(wake_length = 0.2, size = 1, colour = slewcol)  # only has meaning if you do not have a geom_path() element - does not work nicely view_follow()
  # gganimate::view_follow() + # let the view follow the data in each step
  # gganimate::view_step(pause_length = 3, # NOT READY YET
  #                      step_length = 1,
  #                      nsteps = 5) + # This view is a bit like view_follow() but will not match the data in each frame. Instead it will switch between being static and zoom to the range of the data.
  # gganimate::view_zoom() #pause_length = 2,
  #                      step_length = 1,
  #                      nsteps = 7) + # in many ways equivalent to view_step() and view_step_manual() but instead of simply tweening the bounding box of each view it implement the smooth zoom and pan technique developed by Reach & North (2018).
  
  ## some more plot aesthethics
  # ease_aes('cubic-in-out')  # 'cubic-in-out' for a smoother appearance
  
# D2: animated map
  
setwd(file.path(saveloc))

## specify animation parameters
frames <- pull(det_anim_spp %>%
  count(shark) %>%
  summarise(max_rows = max(n)))
anim_dur <- 16 # tell me your desired animation duration
fps <- ceiling(floor(frames/4)/anim_dur)
  
## animate
anim = st +
 gganimate::transition_reveal(along = shifted_date, keep_last = F)+
    # view_follow() + # TODO: make coordinates and axes of plots dynamic so your animal can swim out of it and plot follow
    ease_aes('linear')
    # ggtitle("Storm - M - 315 cm sTL", subtitle = "Date: {frame_along}")
  
an <- gganimate::animate(anim, nframes = floor(frames/4), fps = fps, renderer = gifski_renderer(), width = 180, height = 200, units = "mm", res = 300)

## save it
anim_save(file.path(saveloc, paste0("Animated_tracks_",sp.f,"_with_bathymetry_",floor(frames/4),"_frames.gif")),
            animation = last_animation(), path = NULL, dpi = 300) #, width = 20, height = 15, units = "cm")

  
# END OF CODE ----
# FOR NOW - TO BE CONTINUED

# TODO LIST ----
# TODO 1: create animations with view_follow by individual within the same species, with diff. cols
# TODO 2: re-implement code with symbols and colours based on grouping variables

## create plot with dark themed background
#p = basemap(dt, bathymetry = T, expand.factor = 1.2) + # for bathymetry with ggOceanMaps package
# p <- ggplot() +
#
#   # basemap
#   geom_sf(data = bg)+
#   coord_sf(xlim = range(df_all$lon, na.rm = TRUE),
#            ylim = range(df_all$lat , na.rm = TRUE),
#            expand = T)+
#
#   # raster
#   #geom_polygon(data = r, aes(x = long, y = lat)) +
#
#   # lines and points
#   geom_path(data = df_all,
#             aes(x=lon,y=lat,group=id,color=Behavior),
#             alpha = 0.4, size = 1.5)+
#   geom_point(data = df_all,
#              aes(x=lon,y=lat,group=id,fill=Behavior),
#              alpha = 0.9, shape=21, size = 3)+
#
#   # formatting
#   #scale_fill_viridis_c(option = "inferno")+
#   #scale_color_viridis_c(option = "inferno")+
#   scale_fill_gradientn(colors = smok) +
#   scale_color_gradientn(colors = smok) +
#   scale_size_continuous(range = c(0.1,10))+
#   labs(x=NULL, y=NULL,
#        fill = 'Behaviour',
#        color = 'Behaviour')+
#   theme_dark()+
#   theme(panel.grid = element_blank(), legend.title = element_text(face = "bold"))
# p


# D2: animated map

setwd(paste0(saveloc,"/"))

## specify animation parameters
frames <- nrow(df_all)
anim_dur <- 12 # tell me your desired animation duration
fps <- ceiling(floor(frames/4)/anim_dur)

## animate
anim = p +
  gganimate::transition_reveal(along = date)+
  # view_follow() + # TODO: make coordinates and axes of plots dynamic so your animal can swim out of it and plot follow
  ease_aes('linear')+
  ggtitle("Storm - M - 315 cm sTL", subtitle = "Date: {frame_along}")

an <- gganimate::animate(anim, nframes = floor(frames/4), fps = fps, renderer = gifski_renderer(), width = 350, height = 200, units = "mm", res = 300)

## save it
anim_save(paste0(saveloc, "/Animated_track_Smok_",ptt,"_", min(df_all$date), "_to_", max(df_all$date), "_with_bathymetry_",floor(frames/4),"_frames.gif"),
          animation = last_animation(), path = NULL, dpi = 300) #, width = 20, height = 15, units = "cm")


### ....................................................................................................
### [E] Plot animated data with dynamic axes for multiple individuals ----
### ....................................................................................................

# E1: Plotting parameters ----

## Define your color palettes and shape symbols
### Individuals
# sort(unique(det_anim$id))
# [1] "183623" "200368" "200369" "209020" "222133" "235283" "244607" "244608", "261743"
all_cols <- c("#FFFFCC", #183623
              "#FFEDA0", #200368 same as 222133
              "#FED976", #200369 same as 244607
              "#FEB24C", #209020
              "#FD8D3C", #222133
              "#FC4E2A", #235283
              "#E31A1C", #244607
              "#BD0026", #244608
              "#800026" #261743
)

### Seasonal/grouping variables other than individuals
# summercol <- "#BDD149" # color for wet season
# wintercol <- "#579986" # color for dry season
seasonsym <- c("wet" = 21, "dry" = 25) # shapes for seasons

## Define colours for raster
shallow <- "#D3E5E8"
deep <- "#2B628B"
depth <- c(deep, shallow)

# to deal with odd labels and plot margins
# min(det_anim$lon);max(det_anim$lon)
# min(det_anim$lat);max(det_anim$lat)

## based on basemap
bbox <- st_bbox(basemap)
min_lon <- bbox["xmin"];min_lon
max_lon <- bbox["xmax"];max_lon
min_lat <- bbox["ymin"];min_lat
max_lat <- bbox["ymax"];max_lat

xlabs = seq(-90, -69, 2)
ylabs = seq(18, 35, 2)

## add a symbol of your study animal as datapoint - NOT IMPLEMENTED YET
# hammer_link <- "https://www.pngwing.com/en/free-png-zepta"
# df_all <- df_all %>%
#   mutate(
#     image = hammer_link
#   )

# E2: animated plot ----

dy <- ggplot() +

  # bathymetry raster
  ## stars package
  # stars::geom_stars(data = stars::st_downsample(bathystar, 50) |> sf::st_transform(proj4string(bathyR)), inherit.aes = FALSE) + # choose your epsg code accordinlgy (here EPSG:3857 is for WGS 84 / Pseudo-Mercator -- Spherical Mercator, Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI) - personally don't like it takes ages even with downsample
  ## raster package
  # ggplot2::geom_raster(data = raster.df , aes(x = x, y = y, fill = layer)) +
  ## terra & tidyterra package - by far performs best
  tidyterra::geom_spatraster(data = bathyterra) +
  ## add a legend/gradient to your bathymetry raster
  ggplot2::scale_fill_gradientn(colors = depth, guide = "none") +   # guide = "none" prevents legend
  # labs(x=NULL, y=NULL,
  #      fill = 'Depth [m]',
  #      color = 'Depth [m]')+
  # theme(panel.grid = element_blank(), legend.title = element_text(face = "bold")) +

   #coord_quickmap() +

  # start a new scale
  new_scale_fill() +

  # basemap
  ggplot2::geom_sf(data = basemap, fill = "gray80", color = "black", size = 1, inherit.aes = F) +

  # raster
  #geom_polygon(data = r, aes(x = long, y = lat)) +

  # additional shapefiles
  ## Bahamas EEZ
  ggplot2::geom_sf(data = bah_eez, colour = "black", fill = NA, size = 1) +
  ## US state boundaries
  # ggplot2::geom_sf(data = us_states, colour = "black", fill = NA, size = .25) +

  # define plot limits
  ggplot2::coord_sf(xlim = c(min(det_anim$lon)-0.5, max(det_anim$lon) + 0.5),
           ylim = c(min(det_anim$lat)-0.2, max(det_anim$lat)+0.2),
           expand = T) +

  # lines and points
  geom_path(data = det_anim, # comment this out if you want to animate a wake effect later
            aes(x=lon,y=lat, color = id),
            alpha = .75, linewidth = 1.1)+

  geom_point(data = det_anim,
             aes(x=lon,y=lat,
                 # group = seq_along(Index), # somehow needed if you want to keep data points, i.e. if you want keep_last = T in gganimate::animate()
                 fill=id,
                 shape = season),
             alpha = 1, size = 2.5) +

  # add a symbol rather than a simple data point
  # ggimage::geom_image(aes(image = image), size = 0.1)+ # NOT READY YET

  # plot axis labels
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) +

  # formatting
  # scale_fill_viridis_c(option = "inferno")+
  # scale_color_viridis_c(option = "inferno")+

  ## continous variables
  # scale_fill_gradientn(colors = seasoncol) +
  # scale_color_gradientn(colors = seasoncol) +
  # scale_size_continuous(range = c(0.1,10))+

  ## categorical variables
  ### shape and legend for seasons
  scale_shape_manual(values = seasonsym,
                     name= "Season", # sets legend name
                     drop = F) +
  ### fill colour for shark symbols
  scale_fill_manual(name = "Shark-ID",
                    values = all_cols,
                    drop = F,
                    # guide = "none" # removes legend
                    guide = guide_legend(override.aes = list(fill=all_cols, shape = 22, color = "black")) # use this if you use a wake effect in the animation or if you want a dot legend for individuals
                    ) +
  ### colour and legend of shark track
  scale_color_manual(name = "Shark-ID", # sets legend name
                     values = all_cols,
                     guide = "none", # use this if you don't want any tracks if you use a wake effect in the animation
                     drop = F) +

   # Define your theme aesthethics
  theme(panel.grid = element_blank(), #legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .75), # Add a black border around the plot
        # panel.background = element_rect(fill = "#C1D4E1"), # as alternative as GMRT background takes ages to animate
        # panel.background = element_rect(fill = "white"), # as alternative as GMRT background takes ages to animate
        text = element_text(family = "serif", face = "bold"), # all text to Times New Roman look-a-like
        plot.background = element_rect(fill = "white", color = "white"),    # Set overall plot background to white
        axis.text.x = element_text(size = 12), # change the font size of x.axis text
        axis.text.y = element_text(size = 12), # change the font size of y.axis text
        axis.title.x = element_blank(), # removes x axis title (i.e. "lon")
        axis.title.y = element_blank(), # removes y axis title (i.e. "lat")
        legend.title = element_text(size = 12, face = "bold"), # change the font size of the legend titles
        legend.text = element_text(size = 10),# change the font size of legend text
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(face = "bold", size = 14)
        ) +

  ## add a title
  ggtitle("Great hammerheads of Andros Island, The Bahamas") +   # add title and subtitle

  ## animate your map
  # gganimate::transition_time(Index) + # cannot be used with geom_path
  gganimate::transition_reveal(along = Index, keep_last = F) +
  # gganimate::shadow_wake(wake_length = 0.05) + # only has meaning if you do not have a geom_path() element - does not work nicely view_follow()
  gganimate::view_follow() + # let the view follow the data in each step
  # gganimate::view_step(pause_length = 3, # NOT READY YET
  #                      step_length = 1,
  #                      nsteps = 5) + # This view is a bit like view_follow() but will not match the data in each frame. Instead it will switch between being static and zoom to the range of the data.
  # gganimate::view_zoom() #pause_length = 2,
  #                      step_length = 1,
  #                      nsteps = 7) + # in many ways equivalent to view_step() and view_step_manual() but instead of simply tweening the bounding box of each view it implement the smooth zoom and pan technique developed by Reach & North (2018).

  ## some more plot aesthethics
  ease_aes('cubic-in-out')  # 'cubic-in-out' for a smoother appearance

# dy

# E3: animate your map ----

# animf <- dy +
#   gganimate::transition_reveal(along = date) +
#   gganimate::view_follow() + # let the view follow the data in each step
#   # gganimate::view_step(pause_length = 3, # NOT READY YET
#   #                      step_length = 1,
#   #                      nsteps = 5) + # This view is a bit like view_follow() but will not match the data in each frame. Instead it will switch between being static and zoom to the range of the data.
#   # # gganmiate::view_zoom(pause_length = 1,
#   #                      step_length = 1,
#   #                      nsteps = NULL) + # in many ways equivalent to view_step() and view_step_manual() but instead of simply tweening the bounding box of each view it implement the smooth zoom and pan technique developed by Reach & North (2018).
#
#   ## some more plot aesthethics
#   ease_aes('linear') +
#   ggtitle("Boo - F - 291 cm sTL", subtitle = "Date: {frame_along}") # add title and subtitle

# E4: save your animated map ----

## animation parameters
# frames <- nrow(df_all)
frames <- max(det_anim$Index)
anim_dur <- 25 # tell me your desired animation duration in seconds
fps <- ceiling(floor(frames/4)/anim_dur)

## animate with chosen parameters
gganimate::animate(dy, nframes = floor(frames/4), fps = fps, renderer = gifski_renderer(), width = 185, height = 145, units = "mm", res = 150)

##save it
anim_save(paste0(saveloc, "/Animated_tracks_hammerheads_of_Andros_dynamic_axis_follow_",anim_dur,"_seconds_bathy.gif"),
          animation = last_animation(), path = NULL) #, dpi = 200) #, width = 20, height = 15, units = "cm")

# END OF SCRIPT ----
#### TODO LIST ####
# TODO1: ----
# update animation section to retain original date/time info even when animating multiple individuals
# See this:
# C3: Align datetime info across usable time period for animation

## If you have data from animals tagged across multiple years, animating them
## might be a bit of a pain.
## For presentation pruposes it might be easier if you get them all within the same time period
## Since you might have animals tagged towards the end of the year that cross into the new year
## here we align them withing 24 months.

locs_aligned <- locs %>%
  mutate(month_aligned = cumsum(c(0,diff(month_num)<0))*12+month_num, .by = ptt)

## recreate a new date column based on aligned values

## month_date
locs_aligned$new_date <- as.numeric(paste0(locs_aligned$month_aligned, locs_aligned$day_num))

## create a dummy year variable
locs_aligned <- locs_aligned %>%
  mutate(# add a grouping variable for the shark's home island
    dummy_year = with(.,case_when(# summer/winter based on temperature by van Zinnicq Bergmann et al. 2022 ; summer =  1st Jun to 30th Nov, winter = 1st Dec to 31st May
      month_aligned %in% c(13:24) ~ "2024",
      month_aligned %in% c(1:12) ~ "2023",
      TRUE ~ "thereMightBeAnotherYear"))
  )

locs_aligned$pseudo_date <- paste0(locs_aligned$dummy_year, "-", locs_aligned$month_num, "-", locs_aligned$day_num, " ", locs_aligned$hours_time)
locs_aligned$pseudo_date <- as.Date(locs_aligned$pseudo_date, format = "%Y-%m-%d", tz = "UTC")
#locs_aligned$pseudo_date <- as.POSIXct(locs_aligned$pseudo_date, format = "%Y-%m-%d %H:%M", tz = "UTC")
#locs_aligned$pseudo_date = format(locs_aligned$pseudo_date, format = "%Y-%B-%d", tz = "UTC")
#locs_aligned$pseudo_date = as.Date(locs_aligned$pseudo_date, format = "%Y-%B-%d", tz = "UTC")

# C3: Process the data for animation

## Computed daily average positions and parameter values
locs_an = locs_aligned %>%
  group_by(ptt,pseudo_date) #%>%
#summarise(
#  lat = mean(lat, na.rm = TRUE),
#  lon = mean(lon, na.rm = TRUE)
#)

## create 'ideal' data with all combinations of data
#ideal = expand_grid(
# ptt = unique(locs_an$ptt),
#ideal_date = seq.Date(from = min(locs_an$pseudo_date), to = max(locs_an$pseudo_date), by = 1)
#)

## create complete dataset for plotting
#df_all = left_join(ideal,locs_an)
df_all = locs_an


## Computed daily average positions and parameter values
df_ave = ind %>%
  mutate(date=as.Date(time)) %>%
  group_by(id,date) #%>%
#summarise(
# lat = mean(lat, na.rm = TRUE),
#lon = mean(lon, na.rm = TRUE),
#behav = mean(Behavior, na.rm=TRUE)
#)

## create 'ideal' data with all combinations of data
#ideal = expand_grid(
# id = unique(df_ave$id),
#date = seq.Date(from = min(df_ave$date), to = max(df_ave$date), by = 1)
#)



###
# TODO2: implement bathyR code with terra package - COMPLETED 20250121 - conclusion: terra way faster use thsi ----
# terra package should potentially allow faster plotting, see here: https://github.com/r-spatial/stars/issues/503
# stars::geom_stars(data = stars::st_downsample(cropped_bathy, 2) |> sf::st_transform(3857), inherit.aes = FALSE) + # choose your epsg code accordinlgy (here EPSG:3857 is for WGS 84 / Pseudo-Mercator -- Spherical Mercator, Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI)
# ggplot2::geom_raster(data = bathyraster , aes(x = x, y = y, fill = layer))
# ggplot2::scale_fill_gradientn(colors = depth, guide = "none")   # guide = "none" prevents legend
### bathymetry raster - TODO implement using terra() or stars(), code with bathyloc or something
bathyR <- raster::raster("C:/Users/Vital Heim/switchdrive/Science/Data/Bathymetry_maps_GMRT/GMRTv4_0_20221013topo_for_Sphyrna_SPOT_tracks.grd")
### use rast() rather than raster(), rast() is from the terra package, and terra has more options and will replace raster
#SpatExtent : -84.3068857607634, -74.5949697079866, 22.996942201359, 35.7020417853476 (xmin, xmax, ymin, ymax)
proj4string(bathyR) # from raster package
#> proj4string(bathy)
#[1] "+proj=longlat +datum=WGS84 +no_defs"
### find min and max values
raster::setMinMax(bathyR)
minValue(bathyR); maxValue(bathyR)
### we are working with ocean depth data, i.e. anyhting above 0 meter elevation we dont really need and can set to NA
bathyR <- clamp(bathyR, upper = 0, useValues = F)
### make a df needed for later plotting
# bathy.df <- as.data.frame(bathyR, xy = T)
### make a stars object
bathystar <- stars::st_as_stars(bathyR)
# plot(bathystar, downsample = F)

# TODO3: add additional shapefile showing US State borders, EEZ and state/federal water borders ----
# TODO4: create a animation template with a wake and a fixed extent ----
