#Top of code----
### ====================================================================================================
### Project:    NA
### Analysis:   Visualising Utilisation Distributions (UDs), i.e. spatial probablity densities
### Script:     ~TelemetryMovement/Rscript_Visualise_SpaceUse_Rasters_UDs.R
### Author:     Vital Heim
### Version:    1.0
### ====================================================================================================

### ....................................................................................................
### Content: this R script contains the code to process and visualise space use rasters (i.e. rasters
### containing probability densities of an animal's whereabouts such as UDs) by plotting their values,
### specific quantiles or contours in spatial context.
### Examples using the R package ggplot2 and base R using a raster related packages are provided.
###
### This code has been heavily inspired by the movegroup::plotraster() function by Dedman & van Zinnicq
### Bergmann. For further information check:
### > citation("movegroup")
###   To cite package ‘movegroup’ in publications use:
###
###   Dedman, S, Bergmann vZ, MPM (2024). “movegroup:
###   Visualizing and Quantifying Space Use Data for Groups of Animals.” _In Prep._, *In Prep.*(In Prep.), In
###   Prep. doi:10.5281/zenodo.10805518
###   <https://doi.org/10.5281/zenodo.10805518>,
###   <https://github.com/SimonDedman/movegroup>.
### ....................................................................................................

### ....................................................................................................
### [A] Ready environment, load packages ----
### ....................................................................................................

# A1: clear memory ----

rm(list = ls())

# A2: install and load necessary packages

## if first time
# install.packages("raster")
# install.packages("sf")
# install.packages("sp")
# install.packages("rgdal")
# install.packages("rgeos")
# install.packages("stars")
# install.packages("ggplot2")
# install.packages("ggnewscale")
# install.packages("starsExtra")
# install.packages("terra")
# install.packages("magrittr")
# install.packages("beepr")

## load packages and source needed functions
### raster manipulation
library(raster)
library(sf)
library(sp)
library(rgdal)
library(rgeos)
library(stars)
library(starsExtra)
library(terra)
### plotting
library(ggplot2)
library(ggnewscale)
### general data manipulation
library(magrittr)
### error handling
library(beepr)

options(warn=1) #set this to two if you have warnings that need to be addressed as errors via traceback()

# A3: define file and savelocs ----

rasterloc <- file.path("C:","Users","Vital Heim", "switchdrive", "Science", "Data","MovementRasters","InputData","Andros_allyear")
crsloc <- file.path("C:","Users","Vital Heim", "switchdrive", "Science", "Data","MovementRasters","InputData","Andros_allyear")
shapefileloc <- file.path("C:","Users","Vital Heim", "switchdrive","Science","Data","Shapefiles","World", "gshhg-shp-2.3.7","GSHHS_shp")
saveloc <- file.path("C:","Users","Vital Heim", "switchdrive","Science","Data","MovementPlots")

# A4: define needed functions, universal variables etc.----

## NA

### ....................................................................................................
### [B] Prepare your data for plotting ----
### ....................................................................................................

# B1: import needed data ----

## movement raster
mov_raster <- stars::read_stars(file.path(rasterloc, "Andros_grouplevel_allyear_UDscaled.asc"))
### define projection of your raster
CRSdef <- readRDS(file.path(crsloc, "CRS.Rds"))
sf::st_crs(mov_raster) <- sp::proj4string(CRSdef)

## world shapefile - HQ
### define resolution: 1=c, 2=l,3=i,4=h,5=f, 1:5 increasing quality
# res <- "f"
### read in worldmap
# world <- sf::st_read(dsn = paste0(shapefileloc,"/",res,"/GSHHS_", res, "_L1.shp"), layer = paste0("GSHHS_", res, "_L1"), quiet = TRUE) # read in worldmap

## World shapefile - LQ
world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")

## bathymetry raster - TODO implement using terra() or stars(), code with bathyloc or something
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

## Additional shapefiles (e.g. EEZs, closures, etc.) - TODO code with EEZloc or similar
bah_eez <- sf::read_sf("C:/Users/Vital Heim/switchdrive/Science/Data/Shapefiles/Bahamas/Bahamas_EEZ_boundary.shp")

# B2: prepare raster data for plotting ----

## Trim your data for plotting the UD surface
ud_surface <- mov_raster
is.na(ud_surface[[1]]) <- ud_surface[[1]] == 0 # replace character pattern (0) with NA
# is.na(ud_surface[[1]]) <- ud_surface[[1]] < (max(ud_surface[[1]], na.rm = T)*0.05) # replace anything smaller than the 95 % prop. density
is.na(ud_surface[[1]]) <- ud_surface[[1]] < (max(ud_surface[[1]], na.rm = T)*0.01) # replace anything smaller than the 99 % prop. density

### remove all NA values
ud_surface <- starsExtra::trim2(ud_surface) # NA values were previously 0
### scale back to 100%
ud_surface[[1]] <- (ud_surface[[1]] / max(ud_surface[[1]], na.rm = TRUE)) * 100 # legend is 0:100%
plot(ud_surface)

## Create a raster df for later plotting - only needed if your raster is of type "raster", if you work with the "stars" package
## we deal with this directly in ggplot
# raster_df <- as.data.frame(rasterToPoints(ud_surface), xy = TRUE)
# colnames(raster_df) <- c("x", "y", "value")  # Rename the columns

# B3: extract contours for specific probability density levels ----

## contours for 95 and 50 % probability levels
### 95 %
contour95 <- stars::st_contour(x = mov_raster,
                               contour_lines = TRUE,
                               breaks = max(mov_raster[[1]],na.rm = TRUE) * 0.05
                               ) %>%
             sf::st_cast("POLYGON")
sf::st_crs(contour95) <- sp::proj4string(CRSdef)# set CRS

### 50 %
contour50 <- stars::st_contour(x = mov_raster,
                               contour_lines = TRUE,
                               breaks = max(mov_raster[[1]],na.rm = TRUE) * 0.5
                               ) %>%
             sf::st_cast("POLYGON")
sf::st_crs(contour50) <- sp::proj4string(CRSdef)# set CRS

### ....................................................................................................
### [C] Visualise your data using ggplot2 ----
### ....................................................................................................

# C1: define some plotting parameters ----

## contour colors
color95 <- "black"
color50 <- "black"

## contour linetypes
line95 <- "dashed"
line50 <- "solid"

## bathymetry colors
shallow <- "#D3E5E8"
deep <- "#2B628B"
depth <- c(deep, shallow)

## plot limits
## Option A: based on your ud raster
## find the bounding box, i.e. spatial extent of your ud_surface
bbox_ud_3857 <- st_bbox(ud_surface |> sf::st_transform(3857)) # Get the bounding box (spatial extent)
### Extract min/max longitude and latitude
min_lon <- bbox_ud_3857["xmin"]
max_lon <- bbox_ud_3857["xmax"]
min_lat  <- bbox_ud_3857["ymin"]
max_lat  <- bbox_ud_3857["ymax"]
### Print the results
min_lon
max_lon
min_lat
max_lat
## if manual or multiple plots with same extent
min_lon <- -9354962
max_lon <- -8318496
min_lat <- 2633249
max_lat <- 3495565

## re-convert bathymetry raster to rasterLayer to speed-up plotting
# ds_bathystar <- stars::st_downsample(bathystar, 2)
# cropped_bathy <- st_crop(ds_bathystar |> sf::st_transform(3857), bbox_ud_3857)
#### TODO: find efficient way to plot bathymetry raster
# C2: create a map

ggplot() +
  # plot the bathymetry raster
  # stars::geom_stars(data = stars::st_downsample(cropped_bathy, 2) |> sf::st_transform(3857), inherit.aes = FALSE) + # choose your epsg code accordinlgy (here EPSG:3857 is for WGS 84 / Pseudo-Mercator -- Spherical Mercator, Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI)
  # ggplot2::geom_raster(data = bathyraster , aes(x = x, y = y, fill = layer))
  # ggplot2::scale_fill_gradientn(colors = depth, guide = "none")   # guide = "none" prevents legend
  #
  # # start a new scale
  # ggnewscale::new_scale_fill() +

  # plot the UD surface
  stars::geom_stars(data = ud_surface |> sf::st_transform(3857), inherit.aes = FALSE) + # choose your epsg code accordinlgy (here EPSG:3857 is for WGS 84 / Pseudo-Mercator -- Spherical Mercator, Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI)
  ## adjust the legend and fill colours for the UD surface
  # viridis::scale_fill_viridis(
  #   option = "B", # A magma B inferno C plasma D viridis E cividis F rocket G mako H turbo
  #   alpha = 1, # 0:1
  #   begin = 0, # hue
  #   end = 1, # hue
  #   direction = 1, # colour order, 1 or -1
  #   discrete = FALSE, # false = continuous
  #   space = "Lab",
  #   na.value = "white",
  #   guide = "colourbar",
  #   aesthetics = "fill",
  #   name = "UD",
  #   #labels = ~ 100 - .mov_raster, # https://stackoverflow.com/questions/77609884/how-to-reverse-legend-labels-only-so-high-value-starts-at-bottom
  #   # values are 0-100 with 100=max in the centre but for proportion of time in UD we use % of max with 95% being 0.05 of max.
  #   # So we need to reverse the labels to convert usage density into proportion of time.
  #   position = "right"
  # ) +
  ## adjust the legend and fill colours for the UD surface
  ggplot2::scale_colour_gradient(
    low = "yellow",
    high = "red",
    space = "Lab",
    na.value = NA,
    guide = guide_colourbar(
      title = "UD %",
      barheight = unit(3, "cm"),
      barwidth = unit(.75, "cm"),
      ticks = T,
      frame.colour = "black", # make a box around gradient bar
      ticks.colour = "black",
      ticks.linewidth = .25,
      # reverse = T
    ),
    aesthetics = "fill",
    labels = ~ 100 - .x, # https://stackoverflow.com/questions/77609884/how-to-reverse-legend-labels-only-so-high-value-starts-at-bottom
    # values are 0-100 with 100=max in the centre but for proportion of time in UD we use % of max with 95% being 0.05 of max.
    # So we need to reverse the labels to convert usage density into proportion of time.
    # name = "UD %",
    position = "right"
  ) +

  # Add your 95 % UD contour
  ggplot2::geom_sf(data = stars::st_contour(x = mov_raster,
                                            contour_lines = TRUE,
                                            breaks = max(mov_raster[[1]],na.rm = TRUE) * 0.05) |> sf::st_transform(3857), # define 95 % contour and adjust crs accordingly
                  fill = NA,
                  inherit.aes = FALSE,
                  color = color95,
                  ggplot2::aes(linetype = "General use")) +

  # Add your 50 %  UD contour
  ggplot2::geom_sf(data = stars::st_contour(x = mov_raster,
                                            contour_lines = TRUE,
                                            breaks = max(mov_raster[[1]],na.rm = TRUE) * 0.5) |> sf::st_transform(3857), # define 95 % contour and adjust crs accordingly
                   fill = NA,
                   inherit.aes = FALSE,
                   color = color50,
                   ggplot2::aes(linetype = "Core use")) +

  # Specify linetype and legend for your contours
  ggplot2::scale_linetype_manual(name = "Space use",
                        values = c("General use" = "dashed", "Core use" = "solid"), # Custom linetypes
                        labels = c("Core use", "General use"), # custom line labels
                        guide = guide_legend(
                          keywidth = unit(1.0, "cm"),  # Width of the legend keys
                          keyheight = unit(0.75, "cm"), # Height of the legend keys
                          frame.colour = "black")
                        ) +

  # Add your 95 % and 50 % UD contours in one go
  # ggplot2::geom_sf(data = stars::st_contour(x = mov_raster,
  #                                           contour_lines = TRUE,
  #                                           breaks = c(max(mov_raster[[1]],na.rm = TRUE) * 0.05,
  #                                                      max(mov_raster[[1]],na.rm = TRUE) * 0.5))
  #                                           |> sf::st_transform(3857), # define 95 % contour and adjust crs accordingly
  #                  fill = NA,
  #                  inherit.aes = FALSE,
  #                  ggplot2::aes(colour = "")) +

  # Adjust line colours and line types to your needs


  # Plot the world shapefile
  ggplot2::geom_sf(data = world |> sf::st_transform(3857), fill = "gray90", color = "black", size = 0.5) +

  # Plot additional shapefiles if needed (e.g. for small islands)
  ## needs to be implemented

  # Define your plot limits
  coord_sf(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat)) +

  # Define the order of your legends
  # TBC

  # Clean up general non data related aesthetics
  # theme_minimal() +  # Start with a minimal theme
  theme(
  panel.grid.major = element_blank(),  # Remove major gridlines
  # panel.grid.minor = element_blank(), # Remove minor gridlines
  panel.background = element_rect(fill = "white", color = NA),  # Set background to white
  plot.background = element_rect(fill = "white", color = NA),    # Set overall plot background to white
  panel.border = element_rect(color = "black", fill = NA, linewidth = .75), # Add a black border around the plot
  axis.ticks.length = unit(0.1, "cm"),  # Set length of axis ticks
  text = element_text(family = "serif"), # all text to Times New Roman look-a-like
  axis.text.x = element_text(size = 12), # change the font size of x.axis text
  axis.text.y = element_text(size = 12), # change the font size of y.axis text
  legend.title = element_text(size = 12, face = "bold"), # change the font size of the legend titles
  legend.text = element_text(size = 10), # change the font size of legend text
  # legend.position = "bottom"
  ) +
  labs(
  # x = "Latitude", y = "Longitude"
  )

ggsave(paste0(saveloc, "/UDmap_Andros_grouplevel_allyear.tiff"), width = 21, height = 21, units = "cm", dpi = 300)

#### TO DO List ggplot2 ####

# clean up axes (for latlon aesthetics I prefer degree-minutes-seconds notations, and tick marks on all four side)
# change position, names, size, orientation etc. of legends
# check if ggplot contour legends can be boxes as contours are technically showing areas
# add distance bar and add north arrow
# add option to include additional rasters such as bathymetry raster
# add option to plot 95 and 50 % space use areas as filled out polygons without the continuous UD value scale
# find/add higher res shapefile for land masses

### ....................................................................................................
### [D] Visualise your data using base R,s stars and terra  ----
### ....................................................................................................

# D0: transform your stars raster to a terra raster and prepare for plotting

ud_terra <- rast(ud_surface) |> terra::project("EPSG:4326")

# D1: define some plotting parameters ----

## ud_surface colours
colfuncUD <- colorRampPalette(c("yellow", "red"))

## contour colors
color95 <- "black"
color50 <- "black"

## contour linetypes
line95 <- "dashed"
line50 <- "solid"

## plot limits
## Option A: based on your ud raster
## find the bounding box, i.e. spatial extent of your ud_surface
bbox_ud_4326 <- st_bbox(terra::rast(ud_surface) |> terra::project("EPSG:4326")) # Get the bounding box (spatial extent)
### Extract min/max longitude and latitude
min_lon <- bbox_ud_4326["xmin"]
max_lon <- bbox_ud_4326["xmax"]
min_lat  <- bbox_ud_4326["ymin"]
max_lat  <- bbox_ud_4326["ymax"]
### Print the results
print(c(min_lon,max_lon,min_lat,max_lat))


# D2: create a map

# *D2.1: ud surface ----

## define plotting parameters and legend parameters for your rasters
plgUD = list(#ext=c(max_lon + .1, max_lon + .15, min_lat, max_lat),
             title = "UD %",
             title.cex = 0.9, cex = 0.7, shrink=0)
## define general plotting and aesthetics parameters
pax <- list(side=1:4, lab=c(2,3), tick=1:4, retro=TRUE)

## Plot ud raster
terra::plot(ud_terra, plg=plgUD, pax=pax, las=1, col=colfuncUD(20))

# *D2.2: contours ----

## add 95 % contour
terra::plot(st_geometry(stars::st_contour(x = mov_raster,
                                          contour_lines = TRUE,
                                          breaks = max(mov_raster[[1]],na.rm = TRUE) * 0.05) |> sf::st_transform(4326)),
            add = T, col = "black", lwd = 2, lty = 3, key.pos = 3) # lty=1:solid, lty=2:dashed, lty=3:dotted, lty=4:dot-dash, lty=5:long dash, lty=6:two-dash

## add 50 % contour
terra::plot(st_geometry(stars::st_contour(x = mov_raster,
                                          contour_lines = TRUE,
                                          breaks = max(mov_raster[[1]],na.rm = TRUE) * 0.5) |> sf::st_transform(4326)),
     add = T, col = "black", lwd = 2, lty = 1)

# *D2.3: shapefiles of landmasses

terra::plot(terra::vect(world |> sf::st_transform(4326)), # terra natively uses SpatVector objects for shapefiles
            add = T, col= "gray90", border = "black", lwd = 1.5)

# *D2.4: general plot aesthetics

## Reintroduce plot borders that got overlaid with shapefile
### define extent of your box
ext <- ext(terra::rast(ud_surface) |> terra::project("EPSG:4326"))
### Add a box around the extent
rect(xleft = ext[1], xright = ext[2], ybottom = ext[3], ytop = ext[4],
     border = "black", lwd = 2)

## scalebar
sbar(400, type = "bar", divs = 4)


######## testing area -----

# > class(contour95)
# [1] "sf"         "data.frame"

terra::plot(st_geometry(contour95 |> sf::st_transform(4326)),
            add = T, col = "black", lwd = 2, lty = 3, key.pos = 3) # lty=1:solid, lty=2:dashed, lty=3:dotted, lty=4:dot-dash, lty=5:long dash, lty=6:two-dash

test <- contour95 |> sf::st_transform(4326)
class(test)

testvect <- terra::vect(test)
testvect
names(testvect) <- "general"
testvect

test2 <- contour50 |> sf::st_transform(4326)
class(test)

testvect2 <- terra::vect(test2)
names(testvect2) <- "core"

spat_vector1 <- testvect
spat_vector2 <- testvect2

# Check the attribute names of each SpatVector
names(spat_vector1)
names(spat_vector2)

# Find missing columns in spat_vector1 and add them
missing_cols1 <- setdiff(names(spat_vector2), names(spat_vector1))
if(length(missing_cols1) > 0){
  for(col in missing_cols1){
    spat_vector1[[col]] <- NA  # Add missing columns with NA values
  }
}

# Find missing columns in spat_vector2 and add them
missing_cols2 <- setdiff(names(spat_vector1), names(spat_vector2))
if(length(missing_cols2) > 0){
  for(col in missing_cols2){
    spat_vector2[[col]] <- NA  # Add missing columns with NA values
  }
}

contours <- rbind(spat_vector1, spat_)

plot(contours, "names")

v <- vect(system.file("ex/lux.shp", package="terra"))

## TO DO base R and terra ####
# legend position ud surface and aesthethics
# legend of ud contours
# specific retro labels if you are degree minutes type-of-person
# 95 and 50 % UDs as areas
# legend for UD contours
# north arrow and distance bar

