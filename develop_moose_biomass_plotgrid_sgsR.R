# ------------------------------------------------------------------------------
# develop sample locations for Nazko 
#
# author: Brianne Boufford
# date created: June 20, 2026
# last edit: August 3, 2026
#
# -----------------------------------------------------------------------------

# load libraries 
library(terra) 
library(sf)
library(dplyr)
library(sgsR) 
library(exactextractr)

# set paths 
if (TRUE){
  project_path <- file.path("F:", "Nazko")
  
  # ntems data
  ntems_path <- file.path(project_path, "data", "src", "ntems", "10S")
  
  change_annual_path <- list.files(file.path(ntems_path, "change_annual"), 
                                   pattern = "2023_change_annual.dat$",
                                   full.names = TRUE)
  change_attribution_path <- list.files(file.path(ntems_path, "change_attribution"), 
                                        pattern = ".dat$",
                                        full.names = TRUE)
  fr_path <- list.files(file.path(ntems_path, "change_metrics"), 
                        pattern = "First_Change_Year.dat$",
                        full.names = TRUE)
  lst_path <- list.files(file.path(ntems_path, "change_metrics"), 
                         pattern = "Last_Change_Year.dat$",
                         full.names = TRUE)
  age_path <- list.files(file.path(ntems_path, "age"), 
                         pattern = ".dat$",
                         full.names = TRUE)
  lc_path <- list.files(file.path(ntems_path, "VLCE2.0"), 
                        pattern = "2023_v20_v20.dat$",
                        full.names = TRUE)
  fc_path <- list.files(file.path(ntems_path, "structure", "percentage_first_returns_above_2m"), 
                        pattern = "2023.dat$",
                        full.names = TRUE)
  h_path <- list.files(file.path(ntems_path, "structure", "elev_p95"), 
                       pattern = "2023.dat$",
                       full.names = TRUE)
  
  large_roads_path <- file.path(project_path, "data", "src", "roads", "large_roads_buffer_1km.shp")
  small_roads_path <- file.path(project_path, "data", "src", "roads", "small_roads_buffer_1km.shp")
  
  dpar_roads_path <- file.path(project_path, "data", "src", "roads", "DPAR_clip.shp")
  
  sa_path <- file.path(project_path, "data", "src", "ntems_boundary", "nz_buffer_51km.shp")
  als_grid_path <- file.path(project_path, "data", "bclidar", "aoi_las_grid.shp")
  
  sbps_path <- file.path(project_path, "data", "src", "Nazko_boundaries", "SBPSDC_nz_buffer_51km.shp")
  
}

# read data 
if (TRUE) {
  change_annual <- rast(change_annual_path)
  change_attribution <- rast(change_attribution_path)
  fr <- rast(fr_path)
  lst <- rast(lst_path)
  age <- rast(age_path)
  lc <- rast(lc_path)
  sa <- st_read(sa_path)
  fc <- rast(fc_path)
  h <- rast(h_path)
  als_grid <- st_read(als_grid_path) %>% 
    select(c("year", "geometry")) %>% 
    st_transform(., crs(fc))
  lg_roads <- st_read(large_roads_path) %>%
    st_transform(., crs(fc))
  sm_roads <- st_read(small_roads_path) %>% 
    st_transform(., crs(fc))
  dpar_roads <- st_read(dpar_roads_path)
  sbps <- st_read(sbps_path)
}

# ------------------------------------------------------------------------------
# resample all layers to 60m 
# plot sample design includes satellite plots 30m from plot centre so 
# overall plot radius must be 60m 
# ALS data was collected from 2023-2025, newest disturbance data is from 2023 
# so any distrubance after 2023 is not reflected in the data stratification
# 
# resample VLCE data to 60m using modal bc moose plots ~ 60m pixel size (r = 30m)

# ------------------------------------------------------------------------------
# LANDCOVER 
# ------------------------------------------------------------------------------
lc_60 <- terra::aggregate(lc, 
                          fact = 2, 
                          fun = "modal")
### WRITE
writeRaster(lc_60, 
            file.path(project_path, "data", "src", "vlce_60m.tif"),
            overwrite = TRUE)

# ------------------------------------------------------------------------------
# ROADS
# ------------------------------------------------------------------------------
# filter out undesirable road classes 
dpar_sub <- dpar_roads[!dpar_roads$RD_CLASS %in% c("service", "recreation", "restricted", "trail", "unclassified", "driveway", "lane",
                                                  "skid"), ]
# filter otu undesirable road surface types 
dpar_sub <- dpar_sub[!dpar_sub$RD_SURFACE %in% c("overgrown", "decommissioned"), ]

# remove roads without names (should be none in DPAR dataset)
dpar_sub <- dpar_sub[!is.na(dpar_sub$NAME_FULL), ]

# reproject to CRS of landcover data
dpar_sub <- st_transform(dpar_sub, 
                         crs = crs(lc_60))

### WRITE 
st_write(dpar_sub, 
         file.path(project_path, "data", "src", "roads", "DPAR_subset_aug4.shp"), 
         append = FALSE)

# ------------------------------------------------------------------------------
# ALS GRID
# ------------------------------------------------------------------------------
# rasterize als year of collection 
als_grid_rast <- rasterize(als_grid, 
                      fc, 
                      field = "year") %>% 
  terra::project(., lc_60)

# ------------------------------------------------------------------------------
# MATURE AGE
# ------------------------------------------------------------------------------
# any pixel where last change year was before 1986 is >40 years old and therefore 
# classified as mature
# note: cant use "mode" as function in reprojection so calculating the mean 
# then using a filter for the majority 
mat_mask <- lst < 1986 

# take the mean 
mat_mask <- mat_mask %>% 
  terra::project(., lc_60, mean)

# 60m pixels with 2/4 30m pixels mature will be considered mature 
mat_mask[mat_mask > 0.5] <- 1
mat_mask[mat_mask < 1] <- 0

### WRITE
writeRaster(mat_mask, 
            file.path(project_path, "data", "src", "mature_mask_gt40_aug3.tif"),
            overwrite = TRUE)

# ------------------------------------------------------------------------------
# BEC ZONE (SBPSdc)
# ------------------------------------------------------------------------------

# sub boreal pine spruce region only 
sbps <- st_transform(sbps , crs = crs(lc_60))

### WRITE 
st_write(sbps, 
            file.path(project_path, "data", "src", "SBPS.shp"))

# ------------------------------------------------------------------------------
# WETLANDS 
# - undisturbed / not disturbed in last 40 years 
# - within SBPSdc
# - within ALs coverage 
# - Wetland or Wetland-Treed landcover 
# ------------------------------------------------------------------------------
# 
wetland_mask <- lc_60 == "Wetland" | lc_60 == "Wetland-Treed"
wl <- mask(lc_60, wetland_mask, maskvalues = FALSE) %>% 
  mask(., mat_mask, maskvalues = FALSE) %>%  
  mask(., mask = sbps) %>% 
  mask(., mask = als_grid_rast, maskvalues = NA)

### WRITE 
writeRaster(wl, 
            file.path(project_path, "data", "interm", "wetland_pixels_aug3.tif"),
            overwrite = TRUE)

# ------------------------------------------------------------------------------
# BROADLEAF 
# - undisturbed / not disturbed in last 40 years 
# - within SBPSdc
# - within ALs coverage 
# - broadleaf landcover 
# ------------------------------------------------------------------------------
bl_mask <- lc_60 == "Broadleaf"
bl <- mask(lc_60, bl_mask, maskvalues = FALSE) %>%
  mask(., mat_mask, maskvalues = FALSE) %>%
  mask(., mask = sbps)  %>% 
  mask(., mask = als_grid_rast, maskvalues = NA)

writeRaster(bl, 
            file.path(project_path, "data", "interm", "broadleaf_pixels_aug3.tif"),
            overwrite = TRUE)

# ------------------------------------------------------------------------------
# CONIFEROUS FOREST - mature and disturbance
# ------------------------------------------------------------------------------
# coniferous forest pixels 
conif_mask <- lc_60 == "Coniferous"

# binary layer T/F for areas that experienced fire 
fire_mask <- change_attribution == "Fire"

# binary gets converted to 0-1 scale with reprojection 
fire_mask <- fire_mask %>% 
  terra::project(., lc_60)

# apply modal filter 
fire_mask[fire_mask > 0.5] <- 1
fire_mask[fire_mask < 1] <- 0

fire_mask[fire_mask == 0] <- NA

# do the same for harvesting
harvest_mask <- change_attribution == "Harvesting"

harvest_mask <- harvest_mask %>% 
  terra::project(., lc_60)

harvest_mask[harvest_mask > 0.5] <- 1
harvest_mask[harvest_mask < 1] <- 0

harvest_mask[harvest_mask == 0] <- NA

# make mature coniferous layer 
conif_mat <- mask(lc_60, conif_mask, maskvalues = FALSE) %>%
  mask(., mat_mask, maskvalues = FALSE) %>%
  mask(., mask = sbps)  %>% 
  mask(., mask = als_grid_rast, maskvalues = NA)

### WRITE 
writeRaster(conif_mat, 
            file.path(project_path, "data", "interm", "mature_coinf_pixels_aug3.tif"),
            overwrite = TRUE)

# ------------------------------------------------------------------------------
# recontruct age 
# ------------------------------------------------------------------------------

# for leaflet 
lst_60 <- terra::project(lst, lc_60, method = max)
writeRaster(lst_60, 
            file.path(project_path, "data", "interm", "last_change_year_aug3.tif"),
            overwrite = TRUE)

# also for leaflet 
ch_att_60 <- terra::project(change_attribution, 
                            lc_60)
ch_att_60[!(ch_att_60 == "Fire" | ch_att_60 == "Harvesting")] <- NA
writeRaster(ch_att_60, 
            file.path(project_path, "data", "interm", "change_attribution_fire_harvest_aug3.tif"),
            overwrite = TRUE)

# subtract last change year from
time_diff <- 2026 - lst
time_diff[time_diff == 2026] <- 100
time_diff <- time_diff %>%
  terra::project(., lc_60, min)

r1 <- time_diff <= 10
r2 <- time_diff > 10 & time_diff <= 20
r3 <- time_diff > 20 & time_diff <= 30
r4 <- time_diff > 30 & time_diff <= 40

# recovering forest mask 
# filtering out unsuitable landcover types for post disturbance regions 
recovery_mask <- !(lc_60 == "Water") & !(lc_60 == "Rock/Rubble") & !(lc_60 == "Wetland") & !(lc_60 == "Wetland-Treed") & !(lc_60 == "Broadleaf")

writeRaster(recovery_mask, 
            file.path(project_path, "data", "interm", "recovery_conif_pixels.tif"),
            overwrite = TRUE)

# fire
# 0-10 years
fire_r1 <- mask(lc_60, recovery_mask, maskvalues = FALSE) %>%
  mask(., fire_mask, maskvalues = NA) %>%
  mask(., r1, maskvalues = FALSE) %>%
  mask(., sbps) %>%
  mask(., als_grid_rast, maskvalues = NA)

# 10-20 years
fire_r2 <- mask(lc_60, recovery_mask, maskvalues = FALSE) %>%
  mask(., fire_mask, maskvalues = NA) %>%
  mask(., r2, maskvalues = FALSE) %>%
  mask(., mask = sbps) %>%
  mask(., als_grid_rast, maskvalues = NA)

# 20-30 years
fire_r3 <- mask(lc_60, conif_mask, maskvalues = FALSE) %>%
  mask(., fire_mask, maskvalues = NA) %>%
  mask(., r3, maskvalues = FALSE) %>%
  mask(., mask = sbps) %>%
  mask(., als_grid_rast, maskvalues = NA)

# 30-40 years
fire_r4 <- mask(lc_60, conif_mask, maskvalues = FALSE) %>%
  mask(., fire_mask, maskvalues = NA) %>%
  mask(.,r4, maskvalues = FALSE) %>%
  mask(., mask = sbps) %>%
  mask(., als_grid_rast, maskvalues = NA)


# harvest
# 0-10 years 
harvest_r1 <- mask(lc_60, recovery_mask, maskvalues = FALSE) %>%
  mask(., harvest_mask, maskvalues = NA) %>%
  mask(., r1, maskvalues = FALSE) %>%
  mask(., mask = sbps) %>%
  mask(., als_grid_rast, maskvalues = NA)

# 10-20 
harvest_r2 <- mask(lc_60, recovery_mask, maskvalues = FALSE) %>%
  mask(., harvest_mask, maskvalues = NA) %>%
  mask(.,r2, maskvalues = FALSE) %>%
  mask(., mask = sbps) %>%
  mask(., als_grid_rast, maskvalues = NA)

# 20-30 
harvest_r3 <- mask(lc_60, conif_mask, maskvalues = FALSE) %>%
  mask(., harvest_mask, maskvalues = NA) %>%
  mask(.,r3, maskvalues = FALSE) %>%
  mask(., mask = sbps) %>%
  mask(., als_grid_rast, maskvalues = NA)

# 30-40
harvest_r4 <- mask(lc_60, conif_mask, maskvalues = FALSE) %>%
  mask(., harvest_mask, maskvalues = NA) %>%
  mask(.,r4, maskvalues = FALSE) %>%
  mask(., mask = sbps) %>%
  mask(., als_grid_rast, maskvalues = NA)

# ------------------------------------------------------------------------------
# sample with sgsR
# ------------------------------------------------------------------------------
set.seed(1113)

bl_samp <- sample_srs(bl,
                      nSamp = 10, # number of samples
                      mindist = 90, # bwteen samples 
                      access = dpar_sub, # subset roads layer
                      buff_inner = 60, # at least 1 pixel from road
                      buff_outer = 300, # max 200m from road
                      plot = TRUE)
bl_samp$TYPE <- "broadleaf"

wl_samp <- sample_srs(wl,
                      nSamp = 10, 
                      mindist = 90, 
                      access = dpar_sub, 
                      buff_inner = 60, 
                      buff_outer = 300, 
                      plot = TRUE)
wl_samp$TYPE <- "wetland"

conif_samp <- sample_srs(conif_mat,
                         nSamp = 10, 
                         mindist = 90, 
                         access = dpar_sub, 
                         buff_inner = 60, 
                         buff_outer = 300, 
                         plot = TRUE)
conif_samp$TYPE <- "Mature Coniferous"

fire_r1_samp <- sample_srs(fire_r1,
                         nSamp = 8, 
                         mindist = 90, 
                         access = dpar_sub, 
                         buff_inner = 60, 
                         buff_outer = 300, 
                         plot = TRUE)
fire_r1_samp$TYPE <- "Fire 0-10"

fire_r2_samp <- sample_srs(fire_r2,
                           nSamp = 8, 
                           mindist = 90, 
                           access = dpar_sub, 
                           buff_inner = 60, 
                           buff_outer = 300, 
                           plot = TRUE)
fire_r2_samp$TYPE <- "Fire 11-20"

fire_r3_samp <- sample_srs(fire_r3,
                           nSamp = 8, 
                           mindist = 90, 
                           access = dpar_sub, 
                           buff_inner = 60, 
                           buff_outer = 300, 
                           plot = TRUE)
fire_r3_samp$TYPE <- "Fire 21-30"

fire_r4_samp <- sample_srs(fire_r4,
                           nSamp = 8, 
                           mindist = 90, 
                           access = dpar_sub, 
                           buff_inner = 60, 
                           buff_outer = 300, 
                           plot = TRUE)
fire_r4_samp$TYPE <- "Fire 31-40"

harvest_r1_samp <- sample_srs(harvest_r1,
                           nSamp = 8, 
                           mindist = 90, 
                           access = dpar_sub, 
                           buff_inner = 60, 
                           buff_outer = 300, 
                           plot = TRUE)
harvest_r1_samp$TYPE <- "Harvest 0-10"

harvest_r2_samp <- sample_srs(harvest_r2,
                              nSamp = 8, 
                              mindist = 90, 
                              access = dpar_sub, 
                              buff_inner = 60, 
                              buff_outer = 300, 
                              plot = TRUE)
harvest_r2_samp$TYPE <- "Harvest 11-20"

harvest_r3_samp <- sample_srs(harvest_r3,
                              nSamp = 8, 
                              mindist = 90, 
                              access = dpar_sub, 
                              buff_inner = 60, 
                              buff_outer = 300, 
                              plot = TRUE)
harvest_r3_samp$TYPE <- "Harvest 21-30"

harvest_r4_samp <- sample_srs(harvest_r4,
                              nSamp = 8, 
                              mindist = 90, 
                              access = dpar_sub, 
                              buff_inner = 60, 
                              buff_outer = 300, 
                              plot = TRUE)
harvest_r4_samp$TYPE <- "Harvest 31-40"

all_samp <- rbind(wl_samp, 
                  bl_samp) %>%
  rbind(fire_r1_samp) %>%
  rbind(fire_r2_samp) %>%
  rbind(fire_r3_samp) %>%
  rbind(fire_r4_samp) %>%
  rbind(harvest_r1_samp) %>%
  rbind(harvest_r2_samp) %>%
  rbind(harvest_r3_samp) %>%
  rbind(harvest_r4_samp)%>%
  rbind(conif_samp)

# need to stack LANDCOVER, AGE, CHANGE ATTRIBUTION
ch_att <- change_attribution %>% 
  terra::project(., lc_60)
  
r_stack <- c(lc_60, time_diff, ch_att_60)
names(r_stack) <- c("landcover", "age", "change_attribution")
all_samp_w_attributes <- terra::extract(r_stack, vect(all_samp))

all_samp <- cbind(all_samp, all_samp_w_attributes)

# WRITE OUT SAMPLE LOCATIONS
st_write(all_samp, 
         file.path(project_path, "data", "sgsR", "plot_locations_Aug4_2026.shp"),
         append = FALSE)

# ------------------------------------------------------------------------------
# after making changes, read in and reassign attributes
# ------------------------------------------------------------------------------
plots <- st_read(file.path(project_path, "data", "sgsR", "plot_locations_Aug4_2026.shp"))
plots <- plots[, c("TYPE", "Notes", "geometry")]
plot_attributes <- terra::extract(r_stack, vect(plots))
plots_all <- cbind(plots, plot_attributes)

st_write(plots_all, 
         file.path(project_path, "data", "sgsR", "updated_plot_locations_Aug4_2026.shp"),
         append = FALSE)


#---------- end here 
# check out places with more than one change year (actually not needed but keeping code bite)
stand_replacing_mask <- fr
stand_replacing_mask[change_attribution == "Fire" | change_attribution == "Harvesting"] <- 1
stand_replacing_mask[stand_replacing_mask != 1] <- NA

more_than_one_cy <- fr 
more_than_one_cy[more_than_one_cy == lst] <- NA

more_than_one_cy_masked <- more_than_one_cy*stand_replacing_mask

writeRaster(more_than_one_cy_masked, 
            file.path(project_path, "data", "sgsR", "more_than_one_cy_SR.tif"),
            overwrite = TRUE)
