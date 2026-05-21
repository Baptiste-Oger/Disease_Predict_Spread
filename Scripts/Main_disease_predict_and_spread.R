
## Settings ----------------------------------------------------------------
Sys.setlocale("LC_TIME", "C")

# Loading packages
library(sf)
library(dplyr)
library(lubridate)
library(giscoR)
library(here)


# Set paths
source(file.path(here(), "Config.R"))

# Source code
source(file.path(Path_script_source, "Source_code_sat_RF_inflated.R")) 
source(file.path(Path_script_source, "Source_code_DPC.R")) 
source(file.path(Path_script_source, "Source_code_Smoothing_kernel.R")) 


## Parameters --------------------------------------------------------------

# Target year used for validation
validation_year <- 2020  # 2020 or 2023

# Proportion of observations from the validation year included in training
train_prop_in_validation_year <- 0.2 

# Time step for dpc prediction
time_step_dpc <- 7

# Minimum number of values to fit a DPC 
min_obs_dpc <- 3


## Data prep ---------------------------------------------------------------

DD_raw <- st_read(file.path(Path_data, "jaunisse_satellite_2023-10-12_anonym.shp"))
names(DD_raw)[names(DD_raw) != "geometry"] <- as.character(read.table(file.path(Path_data, "jaunisse_satellite_2023-10-12.names.txt"))[, 2])

DD_L93 <- st_transform(DD_raw, 2154) # Lambert 93 projection (for direct distances)

# Additional time related covariates
DD_L93$date <- ymd(DD_L93$date)
DD_L93$month <- month(DD_L93$date)
DD_L93$day_number <- yday(DD_L93$date)

# Variable to be predicted
Y_total <- DD_L93$gravite

# Covariates used for predictions
pred_covariates <- c(2,3,   # "variete" & campagne (year)
                     14:15, # "SAFRAN_X" & "SAFRAN_Y"
                     17,    # "date_satellite"
                     18:46, # All indices
                     47:56, # All S2 bands
                     87:88, # 60:81, # S2 metadata (angles)
                     82:83, # "timelag" & "j_notation" 
                     89,90) # day_number month
X_total <- st_drop_geometry(DD_L93)[, pred_covariates]

# Information about each individual 
vect_id <- DD_L93$uuid # Location id
vect_date <- DD_L93$date # date
vect_year <- lubridate::year(DD_L93$date) # year
vect_geometry <- DD_L93$geometry # coordinates

# Creating a vector with corresponding departments (sampling unit)
france_sf <- gisco_get_nuts(country = "FR", 
                            resolution = "1",  # 1 = high resolution
                            nuts_level = 3,  # NUTS-3 = departments
                            cache = T, update_cache = T)

france_sf <- st_transform(france_sf, st_crs(DD_L93))
vect_dpt <- st_join(DD_L93, france_sf, join = st_within)$NUTS_NAME


## Stacked hurdle model ----------------------------------------------------

# TRUE / FALSE vector to select training sample
train_select <- Select_training_data(vect_id = vect_id,
                                     vect_year = vect_year,
                                     sampling_group = vect_dpt,
                                     predict_year = validation_year,
                                     prop = train_prop_in_validation_year,
                                     randomness_seed = 10)

# Stacked hurdle model (training and prediction)
result_stacked <- stacked_hurdle_model(train_select = train_select,
                                       X_total = X_total,
                                       Y_total = Y_total,
                                       vect_year = vect_year,
                                       geometry = vect_geometry)



## Reshape data  -----------------------------------------------------------

# Reshape data 
df <- st_as_sf(data.frame(uuid = vect_id,
                          date = vect_date,
                          year = vect_year,
                          y_obs = Y_total, # observation
                          y_sat = result_stacked$pred_total_F$pred, # prediction
                          geometry = vect_geometry,
                          group = vect_dpt))

# For ids used in training, observation is used instead of prediction
df$y_sat[train_select] <- df$y_obs[train_select]

# Subset data for the validation year
df_year <- df[df$year == validation_year,]

grid <- read.csv(file.path(Path_data, paste0("grid_", validation_year, ".csv")))


## Spatio-temporal semi parametric spread model ----------------------------

# Applied to observation (Benchmark)
resultat_obs <- Spatio_temporal_spread(y_col = df_year$y_obs,
                                       id_col = df_year$uuid,
                                       date_col = df_year$date,
                                       geom_col = df_year$geometry,
                                       grp_col = df_year$group,
                                       time_step_dpc = time_step_dpc,
                                       min_obs_dpc = min_obs_dpc,
                                       distance_interval = c(1, 100) * 10^3,
                                       grid = grid,
                                       nrepeat = NULL)

# Applied to local-scale prediction
resultat_sat <- Spatio_temporal_spread(y_col = df_year$y_sat,
                                       id_col = df_year$uuid,
                                       date_col = df_year$date,
                                       geom_col = df_year$geometry,
                                       grp_col = df_year$group,
                                       time_step_dpc = time_step_dpc,
                                       min_obs_dpc = min_obs_dpc,
                                       distance_interval = c(1, 100) * 10^3,
                                       grid = grid,
                                       nrepeat = NULL)

# Store results

results <- list(res_obs = resultat_obs$res,
                res_sat = resultat_sat$res,
                time_seq = resultat_obs$time_seq,
                grid = grid)


## Plot disease dynamics ---------------------------------------------------

# Define a common color scale
max_scale <- ceiling(max(c(max(results$res_obs), max(results$res_sat)))*10)/10

# Define spatial bounding box 
bbox_dd <- st_bbox(c(xmin = 450000, xmax = 880000, ymin = 6690000, ymax = 7100000), 
                   crs = st_crs(2154))  

# Get French regions 
fr_regions <- gisco_get_nuts(country = "FR", resolution = "1", nuts_level = 1,
                             cache = TRUE, update_cache = TRUE)
fr_regions <- st_transform(fr_regions, 2154)

# Get neighboring countries for map context
neigh_countries <- c("BE", "LU", "DE", "CH", "IT", "ES","NL")
neigh_regions <- giscoR::gisco_get_countries(country = neigh_countries, resolution = "1")
neigh_regions <- st_transform(neigh_regions, 2154)

# Create a polygon representing the sea
sea_poly <- st_difference(st_make_valid(st_buffer(st_as_sfc(st_bbox(bbox_dd)), 1e5)),
                          st_union(st_geometry(fr_regions),
                                   st_geometry(neigh_regions)))
sea_poly <- st_as_sf(sea_poly)

# Prepare region labels
fr_regions$region_label <- fr_regions$NUTS_NAME

centroids <- st_centroid(fr_regions)
fr_regions$label_x <- st_coordinates(centroids)[,1]
fr_regions$label_y <- st_coordinates(centroids)[,2]

# Adjust to make it visible
fr_regions$label_x[fr_regions$region_label == "Normandie"] <- fr_regions$label_x[fr_regions$region_label == "Normandie"] + 35000
fr_regions$label_x[fr_regions$region_label == "Grand Est"] <- fr_regions$label_x[fr_regions$region_label == "Grand Est"] - 90000


# Selected weeks to plot (first, intermediate and last week of the season)
week1 <- 1
week2 <- floor(length(results$time_seq)/2)
week3 <- length(results$time_seq)


## Plot --------------------------------------------------------------------

# Function to generate a single map
plot_map <- function(temp, week, color_high, time_seq, grid) {
  
  ggplot() +
    geom_sf(data = sea_poly, fill = "#d6ecff", color = NA) +
    geom_sf(data = neigh_regions, fill = "grey90", color = "grey70", size = 0.3) +
    geom_sf(data = fr_regions, fill = "white", color = "black", size = 0.4) +
    
    geom_raster(data = grid, aes(x = x, y = y, fill = temp[, week]), alpha = 0.7) +
    
    geom_text(data = fr_regions,
              aes(x = label_x, y = label_y, label = region_label),
              size = 4, fontface = "bold") +
    
    scale_fill_gradient(low = "white", high = color_high, limits = c(0, max_scale)) +
    
    geom_sf(data = fr_regions, fill = NA, color = "black", size = 0.4) +
    
    coord_sf(crs = 2154,
             xlim = c(bbox_dd["xmin"], bbox_dd["xmax"]),
             ylim = c(bbox_dd["ymin"], bbox_dd["ymax"])) +
    
    annotation_north_arrow(location = "tl",
                           which_north = "true",
                           style = north_arrow_fancy_orienteering,
                           height = unit(1.2, "cm"),
                           width = unit(1.2, "cm"),
                           pad_x = unit(0.5, "cm"),
                           pad_y = unit(0.5, "cm")) +
    
    annotation_scale(location = "bl",
                     bar_cols = c("grey80", "white"),
                     text_cex = 0.8,
                     width_hint = 0.25) +
    
    theme_minimal() +
    labs(fill = "Severity",
         subtitle = paste("Week of", time_seq[week])) +
    
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_blank(),
          plot.subtitle = element_text(hjust = 0.5, size = 12))
}

# Generate maps
rst_obs1 <- plot_map(temp = results$res_obs, week = length(results$time_seq), "darkred", 
                     time_seq = results$time_seq, grid = results$grid)





