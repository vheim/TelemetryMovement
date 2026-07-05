#TOP OF CODE----
### ====================================================================================================
### Project:    NA
### Analysis:   Visualising Utilisation Distributions (UDs), i.e. spatial probablity densities
### Script:     ~/Project_code/TelemetryMovement/R/Rscript_UD_space_use_raster_visualisation.R
### Author:     Vital Heim
### Version:    1.0
### ====================================================================================================

### ....................................................................................................
### Content: this R script contains the code to process and visualise space use rasters (i.e. rasters
### containing probability densities of an animal's whereabouts such as UDs) by plotting their values,
### specific quantiles or contours in spatial context.
###
### Examples using the R package ggplot2 and base R using a raster related packages are provided.
### ....................................................................................................

### ....................................................................................................
### [A] Ready environment, load packages ----
### ....................................................................................................

# A1: clear memory ----

## clear memory
rm(list = ls())

## if you are working on a network server
# YOUR_IP <- "NA" # add your IP address, server name or similar if you connect via a shared drive
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  use_server <- rstudioapi::showQuestion(
    title = "Network Server",
    message = "Are you connecting via a network server/NAS?",
    ok = "Yes", cancel = "No"
  )

  if (use_server) {
    YOUR_IP <- rstudioapi::showPrompt(
      title = "Server Address",
      message = "Enter your IP address or server name:"
    )
    if (is.null(YOUR_IP) || nchar(trimws(YOUR_IP)) == 0)
      stop("No IP address entered. Please restart and try again.")
  } else {
    YOUR_IP <- ""
  }
} else {
  YOUR_IP <- ""
}

# A2: install and load necessary packages

## if first time
# install.packages("move")
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
### ud calculations
library(move)
### raster manipulation
library(raster)
library(sf)
library(sp)
# library(rgdal)
# library(rgeos)
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

# A3: Specify needed functions ----

## FUNCTION1: enter and format species name for multiple species projects ----
format_species_name <- function() {
  # Request species name from user with colored prompt
  cat("\033[33mEnter species name (e.g., 'S. mokarran' or 'Sphyrna mokarran'): \033[0m")
  species_input <- readline()

  # Remove dots first
  species_input <- gsub("\\.", "", species_input)

  # Split by spaces
  parts <- trimws(unlist(strsplit(species_input, "\\s+")))

  # Extract first letter of genus and full species name
  if (length(parts) >= 2) {
    genus_initial <- toupper(substr(parts[1], 1, 1))
    species_name <- tolower(parts[2])
    spp <- paste0(genus_initial, species_name)
  } else {
    # If no space (already formatted), just clean it
    spp <- gsub("\\s+", "", species_input)
    spp <- paste0(toupper(substr(spp, 1, 1)), tolower(substr(spp, 2, nchar(spp))))
  }

  # Display formatted species name and confirmation in orange/yellow
  # cat("\033[33m==========================================\033[0m\n")
  cat("\033[33mFormatted species name: ", spp, "\033[0m\n", sep = "")
  cat("\033[33mUD maps will be generated for: ", spp, "\033[0m\n", sep = "")
  # cat("\033[33m==========================================\033[0m\n")

  # Return the formatted name
  return(spp)
}

# A4: Specify data and saveloc ----

## Project folder
projloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA")

# Input data
rasterloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA","Data_input", "UD_visualisation")
crsloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA","Data_input", "UD_visualisation")
world_shapefileloc <- file.path("/",YOUR_IP, "Science", "Data_raw", "Shapefiles", "World", "gshhg-shp-2.3.7", "GSHHS_shp")
misc_shapefileloc <- file.path("/",YOUR_IP, "Science", "Data_raw", "Shapefiles")
bathyloc <- file.path("/",YOUR_IP, "Science", "Data_raw", "Bathymetry_maps_GMRT", "GMRTv4_3_0_20250120topo_wider_NWA.grd")

# Output data
saveloc <- file.path("/",YOUR_IP, "Science","Projects_current", "2026_IUCN_ISRA_NWA","Data_output", "UD_visualisation")
### check if saveloc already present
if (!dir.exists(saveloc)) {
  dir.create(saveloc, recursive = TRUE)
  cat("New saveloc was created at:", saveloc, "\n")
} else {
  cat("Saveloc is already present","\n")
}

# A5: Define universal variables (e.g. for plotting) ----

# Plotting variables

## contour colors
color95ct <- "black"
color50ct <- "black"

## contour linetypes
line95ct <- "dashed"
line50ct <- "solid"

## UD surface
lowUD <- "#EACA08"
highUD <-  "#E40200"

## bathymetry colors
shallow <- "#D3E5E8"
deep <- "#2B628B"
depth <- c(deep, shallow)

# A6: Define universal options, variables, etc. (e.g. for plotting) ----

options(timeout = 3000) # manually increase time out threshold (needed when downloading basemap)
options(scipen=999) # so that R doesn't act up for pit numbers
options(warn=1) #set this to two if you have warnings that need to be addressed as errors via traceback()
options(error = function() beep(9))  # give warning noise if it fails

### ....................................................................................................
### [B] Prepare your data for plotting GROUP-LEVEL ----
### ....................................................................................................

# B0: tell the script what species (or other groupings such season) you want to plot

spp <- format_species_name()

# B1: import needed data ----

## UD raster w/stars
# mov_raster <- stars::read_stars(file.path(rasterloc, spp, "All_Rasters_Scaled_Weighted_UDScaled.asc")) # should not have to change this
### define projection of your raster
# CRSdef <- readRDS(file.path(crsloc, spp, "CRS.Rds")) # should not have to change this
# sf::st_crs(mov_raster) <- sp::proj4string(CRSdef)

## UD raster w/ raster
mov_raster <- raster::raster(file.path(rasterloc, spp, "All_Rasters_Scaled_Weighted_UDScaled.asc")) # should not have to change this
### define projection of your raster
CRSdef <- readRDS(file.path(crsloc, spp, "CRS.Rds")) # should not have to change this
raster::crs(mov_raster) <- CRSdef

## world shapefile - HQ
### define resolution: 1=c, 2=l,3=i,4=h,5=f, 1:5 increasing quality
# res <- "f"
### read in worldmap
# world <- sf::st_read(dsn = paste0(shapefileloc,"/",res,"/GSHHS_", res, "_L1.shp"), layer = paste0("GSHHS_", res, "_L1"), quiet = TRUE) # read in worldmap

## World shapefile - LQ
world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")

## bathymetry raster
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
bathyR <- clamp(bathyR, upper = 10, useValues = F) # this is really just an aesthethic thing, needs to be adjsuted if you print a legend that goes between 0 to -8000 meters

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

## Additional shapefiles (e.g. EEZs, closures, etc.) - TODO code with EEZloc or similar
# bah_eez <- sf::read_sf("C:/Users/Vital Heim/switchdrive/Science/Data/Shapefiles/Bahamas/Bahamas_EEZ_boundary.shp")

### US State boundaries
usstates <- sf::st_read(dsn = file.path(misc_shapefileloc, "USA", "US_State_Boundaries", "US_State_Boundaries.shp"), quiet = TRUE) # read in worldmap

### US EEZ
eez_usa <- sf::st_read(dsn = file.path(misc_shapefileloc, "USA", "US_EEZ","eez.shp"), quiet = TRUE) # read in worldmap

# B2: prepare raster data for plotting ----

## first make a new UD object to plot volumetric UD later
newUD <- new(".UD", mov_raster)

## calculate volumetric UD
volUD <- move::getVolumeUD(newUD)
raster::crs(volUD) <- CRSdef # reset CRS - safety

# convert your volumetric UD object into a stars raster for easier contour plotting later
ud_stars <- stars::st_as_stars(volUD)
# add your crs info, function from sp
sf::st_crs(ud_stars) <- sp::proj4string(CRSdef) # this should be used for volumetric UD calcs, i.e. ud_stars

## Trim your data for plotting the UD surface
### tell the script if you want to trim (prep for movegroup inclusion)
trim <- T
### define the % UD you want to plot (e.g. 99%) as surface on continuous scale
surfaceUDpct <- .99
### Trim your UD raster surface accordingly
y <- ud_stars # make dupe object else removing all data < 0.05 means the 0.05 contour doesn't work in ggplot
if (trim) { # trim raster extent to data, remove 0 values and crop to 99% UD
  is.na(y[[1]]) <- y[[1]] == 0 # replace char pattern (0) in whole df/tbl with NA
  # is.na(y[[1]]) <- y[[1]] < (max(y[[1]], na.rm = TRUE) * 0.05) # replace anything < 95% contour with NA since it won't be drawn - DEPRECATED 20260626 as based on raw probabilities
  is.na(y[[1]]) <- as.vector(y[[1]]) > surfaceUDpct # mask all cells outside the 95% UD
  # does this mean subsequent 50 & 95%'s are 95% of 95%?
  # No: Y is only used for the UD surface, NOT the contours, that's based on entire UD surface
}
y <- starsExtra::trim2(y) # remove NA columns, which were all zero columns. This changes the bbox accordingly
y[[1]] <- (y[[1]] / max(y[[1]], na.rm = TRUE)) * 100 # convert from raw values to 0:100 scale for later colour display, legend breaks will be corrected later so they correspond to absolute UD values
# plot(y)

### Since we trim the raster we need to prepare breaks and labels for later plotting
legendbreaks <- c(0, 25/surfaceUDpct, 50/surfaceUDpct, 75/surfaceUDpct, 100)
legendlabels <- c("0", "25", "50","75", paste0(surfaceUDpct*100))

## Create a raster df for later plotting - only needed if your raster is of type "raster", if you work with the "stars" package
## we deal with this directly in ggplot
# raster_df <- as.data.frame(rasterToPoints(ud_surface), xy = TRUE)
# colnames(raster_df) <- c("x", "y", "value")  # Rename the columns

# B3: extract contours for specific probability density levels ----

# ## contours for raw 95 and 50 % probability levels - DEPRECATED 20260705
# ### 95 %
# contour95 <- stars::st_contour(x = mov_raster,
#                                contour_lines = TRUE,
#                                breaks = max(mov_raster[[1]],na.rm = TRUE) * 0.05
#                                ) %>%
#              sf::st_cast("POLYGON")
# sf::st_crs(contour95) <- sp::proj4string(CRSdef)# set CRS
#
# ### 50 %
# contour50 <- stars::st_contour(x = mov_raster,
#                                contour_lines = TRUE,
#                                breaks = max(mov_raster[[1]],na.rm = TRUE) * 0.5
#                                ) %>%
#              sf::st_cast("POLYGON")
# sf::st_crs(contour50) <- sp::proj4string(CRSdef)# set CRS



### ....................................................................................................
### [C] Visualise your data using ggplot2 ----
### ....................................................................................................

# C1: define some plotting parameters ----

## plot limits
## Option A: based on your ud raster
## find the bounding box, i.e. spatial extent of your ud_surface
bbox_ud_3857 <- st_bbox(y |> sf::st_transform(3857)) # Get the bounding box (spatial extent) #
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
## SLEWINI
# min_lon <- -11175792
# max_lon <- -7624101
# min_lat <-2552507
# max_lat <- 5017212
## SMOKARRAN
# min_lon <- -10365241
# max_lon <- -8344251
# min_lat <- 2733758
# max_lat <- 4385954
## SZYGAENA
# min_lon <- -9830292
# max_lon <- -7845821
# min_lat <- 2665792
# max_lat <- 5138885

## re-convert bathymetry raster to rasterLayer to speed-up plotting
# ds_bathystar <- stars::st_downsample(bathystar, 2)
# cropped_bathy <- st_crop(ds_bathystar |> sf::st_transform(3857), bbox_ud_3857)
#### TODO: find efficient way to plot bathymetry raster in high res

## contour levels - define them (prep for function release movegroup)
generalUDpct <- .95
coreUDpct <- .5
## contour legend labels, driven by the users chosen ud pcts (prep for function movegroup)
## NOTE to Vital: if the user specifies the contours to be 95 and 50%, this would match
## generally used definitions of general and core space use. So it seems nicer, if these expressions
## are used as labels for easier udnerstanding by non-ecologists. However, the user might select
## some non-conventially used percentages. So we want to label accordingly. We can do this by running this function first:
## Function to name text at the conventional cutoffs, numeric otherwise
udlab <- function(pct) {
  if (isTRUE(all.equal(pct, .95))) return("General use")
  if (isTRUE(all.equal(pct, .50))) return("Core use")
  paste0(round(pct * 100, 2), "% UD")
}
## and now update labels
# genlab  <- paste0(generalUDpct * 100, "% UD")   # "95% UD"
# corelab <- paste0(coreUDpct    * 100, "% UD")   # "50% UD"
genlab  <- udlab(generalUDpct)   # .95 -> "General use", .91 -> "91% UD"
corelab <- udlab(coreUDpct)      # .50 -> "Core use",    .45 -> "45% UD"

# C2: create a map

g <- ggplot() +
  # plot the bathymetry raster
  ## terra & tidyterra package - by far performs best
  tidyterra::geom_spatraster(data = bathyterra |> terra::project("EPSG:3857")) +
  ## add a bathymetry legend
  ggplot2::scale_fill_gradientn(
    colors   = depth,
    name     = "Depth [m]",
    na.value = NA,
    guide = ggplot2::guide_colourbar(
      barheight       = ggplot2::unit(3, "cm"),
      barwidth        = ggplot2::unit(0.75, "cm"),
      frame.colour    = "black",
      ticks.colour    = "black",
      ticks.linewidth = 0.25,
      order           = 3          # forces this beneath the UD legend
    )
  ) +
  ## or suppress a legend/gradient to your bathymetry raster
  # ggplot2::scale_fill_gradientn(colors = depth, guide = "none")    # guide = "none" prevents legend
  # labs(x=NULL, y=NULL,
  #      fill = 'Depth [m]',
  #      color = 'Depth [m]')+
  # theme(panel.grid = element_blank(), legend.title = element_text(face = "bold")) +

  # start a new scale
  new_scale_fill() +
  # # start a new scale
  # ggnewscale::new_scale_fill() +

  # plot the UD surface
  stars::geom_stars(data = y |> sf::st_transform(3857), inherit.aes = FALSE, # choose your epsg code accordinlgy (here EPSG:3857 is for WGS 84 / Pseudo-Mercator -- Spherical Mercator, Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI)
                    color = NA) + # set color to NA so that you dont have border lines for each raster cell
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
    low = highUD,
    high = lowUD,
    space = "Lab",
    # na.value = NA,
    guide = guide_colourbar(
      title = "UD %",
      barheight = unit(3, "cm"),
      barwidth = unit(.75, "cm"),
      ticks = T,
      frame.colour = "black", # make a box around gradient bar
      ticks.colour = "black",
      ticks.linewidth = .25,
      reverse = T,
      order = 2
    ),
    aesthetics = "fill",
    # labels = ~ ifelse(100 - .x == 0, "", 100 - .x), # https://stackoverflow.com/questions/77609884/how-to-reverse-legend-labels-only-so-high-value-starts-at-bottom
    # # values are 0-100 with 100=max in the centre but for proportion of time in UD we use % of max with 95% being 0.05 of max.
    # # So we need to reverse the labels to convert usage density into proportion of time. And we remove the 0 values as a 0% UD seems counterintutitve to inexperienced readers
    # # name = "UD %",
    labels = legendlabels,
    breaks = legendbreaks,
    position = "right"
  ) +

  # Add your 95 % UD contour
  ## raw probability - DEPRECATED 202607
  # ggplot2::geom_sf(data = stars::st_contour(x = mov_raster,
  #                                           contour_lines = TRUE,
  #                                           breaks = max(mov_raster[[1]],na.rm = TRUE) * 0.05) |> sf::st_transform(3857), # define 95 % contour and adjust crs accordingly
  #                 fill = NA,
  #                 inherit.aes = FALSE,
  #                 color = color95,
  #                 linewidth = .6, # default 0.5
  #                 ggplot2::aes(linetype = "General use")) +
  # Add your 95 % and 50 % UD contours in one go - DEPRECATED 202607
  # ggplot2::geom_sf(data = stars::st_contour(x = mov_raster,
  #                                           contour_lines = TRUE,
  #                                           breaks = c(max(mov_raster[[1]],na.rm = TRUE) * 0.05,
  #                                                      max(mov_raster[[1]],na.rm = TRUE) * 0.5))
  #                                           |> sf::st_transform(3857), # define 95 % contour and adjust crs accordingly
  #                  fill = NA,
  #                  inherit.aes = FALSE,
  #                  ggplot2::aes(colour = "")) +
  ## 95% volumetric UD
  ggplot2::geom_sf(data = stars::st_contour(x = ud_stars,
                                            contour_lines = TRUE,
                                            breaks = generalUDpct) |>
                     sf::st_transform(3857),
                   fill = NA,
                   inherit.aes = FALSE,
                   ggplot2::aes(colour = genlab, linetype = genlab)
  ) +
  # Add your 50 %  UD contour
  ## ## raw probability - DEPRECATED 202607
  # ggplot2::geom_sf(data = stars::st_contour(x = mov_raster,
  #                                           contour_lines = TRUE,
  #                                           breaks = max(mov_raster[[1]],na.rm = TRUE) * 0.5) |> sf::st_transform(3857), # define 95 % contour and adjust crs accordingly
  #                  fill = NA,
  #                  inherit.aes = FALSE,
  #                  color = color50,
  #                  linewidth = .6, # default 0.5
  #                  ggplot2::aes(linetype = "Core use")) +
  ## 50% volumetric UD
  ggplot2::geom_sf(data = stars::st_contour(x = ud_stars,
                                            contour_lines = TRUE,
                                            breaks = coreUDpct) |>
                     sf::st_transform(3857),
                   fill = NA,
                   inherit.aes = FALSE,
                   ggplot2::aes(colour = corelab, linetype = corelab)
  ) +

  # Specify linetype and legend for your contours
  # ggplot2::scale_linetype_manual(name = "Space use",
  #                       values = c("General use" = "21", "Core use" = "solid"), # Custom linetypes: dashed = 44, dotted = 13, tight dots = 31, tight dashes = 42, very short dashes, very short gaps = 21
  #                       labels = c("Core use", "General use"), # custom line labels
  #                       guide = guide_legend(
  #                         keywidth = unit(1.0, "cm"),  # Width of the legend keys
  #                         keyheight = unit(0.75, "cm"), # Height of the legend keys
  #                         frame.colour = "black")
  #) +
  ggplot2::scale_colour_manual(
    name   = "Space use",
    breaks = c(genlab, corelab),
    values = setNames(c(color95ct, color50ct), c(genlab, corelab)),
    guide  = ggplot2::guide_legend(order = 1) # puts this legend to the top
  ) +
  ggplot2::scale_linetype_manual(
    name   = "Space use",
    breaks = c(genlab, corelab),
    values = setNames(c(line95ct, line50ct), c(genlab, corelab)),
    guide  = ggplot2::guide_legend(
      order        = 1, # puts this at the top of the legends
      keywidth     = ggplot2::unit(1.0, "cm"),
      keyheight    = ggplot2::unit(0.75, "cm"),
      override.aes = list(linewidth = 0.6)
    )
  ) +

  # Plot the world shapefile
  ggplot2::geom_sf(data = world |> sf::st_transform(3857), fill = "gray90", color = "black", size = 0.5) +

  # Plot additional shapefiles if needed (e.g. for small islands)
  ## US State boundaries
  ggplot2::geom_sf(data = usstates |> sf::st_transform(3857), fill = NA, color = "black", size = 0.5) +
  ## US EEZ
  ggplot2::geom_sf(data = eez_usa, colour = "black", fill = NA, size = .25) +

  # Define your plot limits
  coord_sf(xlim = c(min_lon+00000, max_lon-00000), ylim = c(min_lat+00000, max_lat-00000)) +

  # Define the order of your legends
  # done above specifically for each legend, not ideal, but seems to be the only way
  # if we have mutliple 'fill' legends. Something to look into with package updates etc. in future.

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
  ); print (g)

ggsave(file.path(saveloc, paste0("UDmap_",spp,"_group_level.tiff")), width = 21, height = 21, units = "cm", dpi = 300)

#### TO DO List ggplot2 ####

# clean up axes (for latlon aesthetics I prefer degree-minutes-seconds notations, and tick marks on all four side)
# change position, names, size, orientation etc. of legends
# check if ggplot contour legends can be boxes as contours are technically showing areas
# add distance bar and add north arrow
# add option to include additional rasters such as bathymetry raster
# add option to plot 95 and 50 % space use areas as filled out polygons without the continuous UD value scale
# find/add higher res shapefile for land masses

### ....................................................................................................
### [D] Prepare your data for plotting at INDIVIDUAL-LEVEL ----
### ....................................................................................................

# D1: import individual-level UD rasters

## Remember: based on your dBBMM script, all scaled, individual-level rasters should have a filename that
## starts with a capital X, i.e. "X123456". We can use this to import all files at once without importing
## group-lelve rasters within the same directory

ind_rasters <- list.files(path = file.path(rasterloc, spp),
                          pattern = "^X.*\\.asc$",
                          full.names = TRUE)

## check file import
if (length(ind_rasters) == 0) {
  stop("No raster files starting with 'X' found in the directory")
}; print(paste("Found", length(ind_rasters), "raster files to process"))

# D2: import all files that are needed for the plot but only once (e.g. bathyraster, CRS, etc.)

## If you ran the group-level code, all files have been previously imported. Check above.

### ....................................................................................................
### [E] Visualise your data using ggplot2 at INDIVIDUAL-LEVEL absolute scale ----
### ....................................................................................................

# E1: Prepare plot limits ----

## Depending on your choice, the plots can have a fixed extent across all individuals, or we can plat
## the rasters so that the extent of the plot changes depending on the space use range of your animal

## IMPORTANT!! if your data covers a very large extents (e.g. several UTM zones), dynamic extents can be
## be problematic. In this case choose fixed extent below.
## Options to deal with this issue will be added in the future.

## User specificed plot limits choice and calculation of bounding box
valid_input <- FALSE
while (!valid_input) {
  cat("\033[33mEnter your choice (1 - DYNAMIC or 2 - FIXED extent): \033[0m")
  user_choice <- readline()

  if (user_choice == "1") {
    bbox_option <- "dynamic"
    valid_input <- TRUE
    cat(paste("\nYou selected:", toupper(bbox_option), "map extent option\n"))
  } else if (user_choice == "2") {
    bbox_option <- "fixed"
    valid_input <- TRUE
    cat(paste("\nYou selected:", toupper(bbox_option), "map extent option\n"))

    # Set up fixed bbox using previously defined values
    cat("\nUsing previously defined group-level fixed bbox values\n")

    # Check if the variables exist from previous code section
    if (!exists("min_lon") || !exists("max_lon") || !exists("min_lat") || !exists("max_lat")) {
      stop("Error: Fixed bbox values (min_lon, max_lon, min_lat, max_lat) not found.\nPlease ensure these are defined in previous code for group-level UD maps.")
    }

    # Store the previously defined values as fixed bbox
    fixed_min_lon <- min_lon
    fixed_max_lon <- max_lon
    fixed_min_lat <- min_lat
    fixed_max_lat <- max_lat

    cat(paste("Fixed bbox values:\n"))
    cat(paste("  X range:", round(fixed_min_lon, 2), "to", round(fixed_max_lon, 2), "\n"))
    cat(paste("  Y range:", round(fixed_min_lat, 2), "to", round(fixed_max_lat, 2), "\n"))
  } else {
    cat("Invalid choice. Please enter 1 or 2.\n")
  }
}

# E2: Plot individual-level rasters ----

## create new directory for individual maps first
individual_dir <- file.path(saveloc, paste0(spp, "_individuals"))
if (!dir.exists(individual_dir)) {
  dir.create(individual_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\033[32mCreated new directory for individual maps:\033[0m\n")
  cat("\033[32m", individual_dir, "\033[0m\n")
} else {
  cat("\033[32mUsing existing directory for individual maps:\033[0m\n")
  cat("\033[32m", individual_dir, "\033[0m\n")
}

for (i in seq_along(ind_rasters)) {

  cat(paste("Processing file", i, "of", length(ind_rasters), "\n"))
  cat(paste("File:", basename(ind_rasters[i]), "\n"))

  # Extract shark ID from filename (remove X prefix and .asc extension)
  shark_id <- tools::file_path_sans_ext(basename(ind_rasters[i]))
  shark_id <- sub("^X", "", shark_id)  # Remove leading X

  # Read the current raster
  # mov_raster_i <- stars::read_stars(ind_rasters[i])
  # sf::st_crs(mov_raster_i) <- sp::proj4string(CRSdef)
  mov_raster_i <- raster::raster(ind_rasters[i])   # RasterLayer, needed for move::getVolumeUD()
  raster::crs(mov_raster_i) <- CRSdef

  # NOTE: When writing rasters to disk with dataType FLTS4 there is a risk of accuracy loss.
  # This can cause an error 'The used raster is not a UD (sum unequal to 1)'
  # This can be resolved by re-standardizing individual rasters
  mov_raster_i <- mov_raster_i / raster::cellStats(mov_raster_i, "sum", na.rm = TRUE)

  # Calculate the volumetric UD
  newUD_i <- new(".UD", mov_raster_i)
  volUD_i    <- move::getVolumeUD(newUD_i)
  raster::crs(volUD_i) <- CRSdef

  # Convert volumetric UD into a stars raster for easier contour plotting
  ud_stars_i <- stars::st_as_stars(volUD_i)
  sf::st_crs(ud_stars_i) <- sp::proj4string(CRSdef)

  # Trim your data for plotting the UD surface
  y_i <- ud_stars_i
  if (trim) {
    is.na(y_i[[1]]) <- y_i[[1]] == 0
    is.na(y_i[[1]]) <- as.vector(y_i[[1]]) > surfaceUDpct
  }
  y_i <- starsExtra::trim2(y_i)
  y_i[[1]] <- (y_i[[1]] / max(y_i[[1]], na.rm = TRUE)) * 100

  # Define plot limits based on chosen option
  if (bbox_option == "fixed") {
    # Use fixed bbox for all maps
    min_lon <- fixed_min_lon
    max_lon <- fixed_max_lon
    min_lat <- fixed_min_lat
    max_lat <- fixed_max_lat
    cat("Using fixed bbox for all maps\n")
  } else {
    # Use individual bbox for each map
    bbox_ud_3857 <- st_bbox(y_i |> sf::st_transform(3857))
    min_lon <- bbox_ud_3857["xmin"]
    max_lon <- bbox_ud_3857["xmax"]
    min_lat <- bbox_ud_3857["ymin"]
    max_lat <- bbox_ud_3857["ymax"]
    cat("Using individual bbox for this map.\n")
  }

  cat(paste("  Plot extent - X:", round(min_lon, 2), "to", round(max_lon, 2), "\n"))
  cat(paste("  Plot extent - Y:", round(min_lat, 2), "to", round(max_lat, 2), "\n"))

  # Create the plot
  p <- ggplot() +
    # Bathymetry raster
    tidyterra::geom_spatraster(data = bathyterra |> terra::project("EPSG:3857"),
                               maxcell = 1e7) + # Increase max cells for higher resolution
    ## add a bathymetry legend
    ggplot2::scale_fill_gradientn(
      colors   = depth,
      name     = "Depth [m]",
      na.value = NA,
      guide = ggplot2::guide_colourbar(
        barheight       = ggplot2::unit(3, "cm"),
        barwidth        = ggplot2::unit(0.75, "cm"),
        frame.colour    = "black",
        ticks.colour    = "black",
        ticks.linewidth = 0.25,
        order           = 3          # forces this beneath the UD legend
      )
    ) +
    ## or supress the legend
    # ggplot2::scale_fill_gradientn(colors = depth, guide = "none") +

    # Start new scale
    new_scale_fill() +

    # UD surface
    stars::geom_stars(data = y_i |> sf::st_transform(3857),
                      inherit.aes = FALSE,
                      color = NA) +

    # UD surface fill scale
    ggplot2::scale_colour_gradient(
      low = highUD,
      high = lowUD,
      space = "Lab",
      na.value = NA,
      guide = guide_colourbar(
        title = "UD %",
        barheight = unit(3, "cm"),
        barwidth = unit(.75, "cm"),
        ticks = T,
        frame.colour = "black",
        ticks.colour = "black",
        ticks.linewidth = .25,
        reverse = T,
        order = 2
      ),
      aesthetics = "fill",
      # labels = ~ ifelse(100 - .x == 0, "", 100 - .x),
      labels = legendlabels,
      breaks = legendbreaks,
      position = "right"
    ) +

    # 95% UD contour
    ggplot2::geom_sf(data = stars::st_contour(x = ud_stars_i,
                                              contour_lines = TRUE,
                                              # breaks = max(mov_raster_i[[1]], na.rm = TRUE) * 0.05
                                              breaks = generalUDpct) |>
                       sf::st_transform(3857),
                     fill = NA,
                     inherit.aes = FALSE,
                     # color = color95,
                     # linewidth = .6, # default = .5
                     # ggplot2::aes(linetype = "General use"))
                     ggplot2::aes(colour = genlab, linetype = genlab)
    ) +

    # 50% UD contour
    ggplot2::geom_sf(data = stars::st_contour(x = ud_stars_i,
                                              contour_lines = TRUE,
                                              # breaks = max(mov_raster_i[[1]], na.rm = TRUE) * 0.5
                                              breaks = coreUDpct) |>
                       sf::st_transform(3857),
                     fill = NA,
                     inherit.aes = FALSE,
                     # color = color50,
                     # linewidth = .6, # default = .5
                     # ggplot2::aes(linetype = "Core use"))
                     ggplot2::aes(colour = corelab, linetype = corelab)
    ) +

    # Linetype scale
    # ggplot2::scale_linetype_manual(
    #   name = "Space use",
    #   values = c("General use" = "21", "Core use" = "solid"), # Custom linetypes: dashed = 44, dotted = 13, tight dots = 31, tight dashes = 42, very short dashes, very short gaps = 21
    #   labels = c("Core use", "General use"),
    #   guide = guide_legend(
    #     keywidth = unit(1.0, "cm"),
    #     keyheight = unit(0.75, "cm"),
    #     frame.colour = "black"
    #   )
    # ) +
    ggplot2::scale_colour_manual(
      name   = "Space use",
      breaks = c(genlab, corelab),
      values = setNames(c(color95ct, color50ct), c(genlab, corelab)),
      guide  = ggplot2::guide_legend(order = 1)
    ) +
    ggplot2::scale_linetype_manual(
      name   = "Space use",
      breaks = c(genlab, corelab),
      values = setNames(c(line95ct, line50ct), c(genlab, corelab)),
      guide  = ggplot2::guide_legend(
        order = 1,
        keywidth = ggplot2::unit(1.0,"cm"),
        keyheight = ggplot2::unit(0.75,"cm"),
        override.aes = list(linewidth = 0.6)
      )
    ) +

    # World shapefile
    ggplot2::geom_sf(data = world |> sf::st_transform(3857),
                     fill = "gray90", color = "black", size = 0.5) +

    # Additional shapefiles
    ## US State boundaries
    ggplot2::geom_sf(data = usstates |> sf::st_transform(3857), fill = NA, color = "black", size = 0.5) +
    ## US EEZ
    ggplot2::geom_sf(data = eez_usa, colour = "black", fill = NA, size = .25) +

    # Plot limits
    coord_sf(xlim = c(min_lon + 00100, max_lon - 00100),
             ylim = c(min_lat + 00100, max_lat - 00100)) +

    # Theme
    theme(
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = .75),
      axis.ticks.length = unit(0.1, "cm"),
      text = element_text(family = "serif"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    labs() +

    # Add title with shark ID
    ggtitle(paste("Shark ID:", shark_id))

  # Create output filename (remove .asc extension and add .png)
  output_filename <- file.path(individual_dir,
                               paste0("UD_map_", tools::file_path_sans_ext(basename(ind_rasters[i])),
                                      "_individual-level.tiff"))

  # Save your individual map
  ggsave(filename = output_filename, plot = p, width = 21, height = 21, units = "cm", dpi = 300)

  cat(paste("Plot saved to:", output_filename, "\n"))

  # Optional: display the plot in RStudio
  # print(p)

  # Clean up large objects for this iteration
  rm(mov_raster_i, ud_stars_i, newUD_i, volUD_i, p)
  gc()
}; cat("\033[32mAll rasters processed successfully!\033[0m\n")

# END OF CODE ----

## TODO:
## TODO1 indivudal-level: clean up loop and comment out -----
## TODO2: add a title/label for the ID of the animal to the individual plot ----
## TODO3: based on your extent choice, create subdirectory and save plots accordingly in dedicated directory ----
## TODO4: current problem with extents if plotting both: ----
##        if you plot individual level plots with dynamic extent and then go back and plot them with fixed as well, you end up with extent from last individual raster. makes no sense. need to extent the loop that first all inbdividual maps are plotted with fixed and then dynamic extent and also combine it with the targeted saving in specific directories from to do 3.
## TODO5: deal with finding EPSG for data with large xtent
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
# END of script ----

#TODO-List
## TODO1: ----
# legend position ud surface and aesthethics
# legend of ud contours
# specific retro labels if you are degree minutes type-of-person
# 95 and 50 % UDs as areas
# legend for UD contours
# north arrow and distance bar
# remove 0% from UD scale bar
# add a background raster using terra similar to animate movement script
# TODO2: ----
# automate plotting/plot selection at the beginning
# ADD US EEZ shapefile

