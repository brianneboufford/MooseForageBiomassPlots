# subsetting BEC zone to AOI
# June 2nd, 2026

library(sf)
library(dplyr)

# set wd to Nazko data dir
setwd(file.path("F:", "Nazko"))

# path to nazko boundary files
path_to_nz_data <- file.path("F:", "Nazko", "data", "src", "Nazko_boundaries")

# path to community boundary area
path_to_bound <- file.path(path_to_nz_data, "Nazko_Community_Area.shp")

# path to BEC zones 
path_to_bec <- file.path("F:", "_BCLayers", "BEC_zones_2024", "BEC_BIOGEOCLIMATIC_POLY", "BEC_POLY_polygon.shp")

# NASKO COMMUNITY AREA
nz <- st_read(path_to_bound)

# buffer by 50km 
nz_buffer <- st_buffer(nz, dist = 51000)
st_write(nz_buffer, 
         file.path(path_to_nz_data, "nz_buffer_51km.shp"), 
         append = FALSE)

# ALL BEC
bec_bc <- st_read(path_to_bec)

# JUST BECS of interests 
sbpsdc <- subset(bec_bc, MAP_LABEL %in% c("SBPSdc"))
st_write(sbpsdc, 
      file.path(path_to_nz_data, "SBPSDC_BEC_SUBZONE.shp"), 
      append = FALSE)

# get intersection of SBPSdc and nazko 50km buffer
nz_sbpsdc <- st_intersection(sbpsdc, 
                     nz_buffer)
st_write(nz_sbpsdc, 
         file.path(path_to_nz_data, "SBPSDC_nz_buffer_51km.shp"), 
         append = FALSE)

                     
# look back in my notes but I think it is 
# - 50 km from Nazko community
# - < 1km from a road
# - 30 plots (0-40 yrs post harvest) (NTEMS)
# - 30 plots (0-40 yrs post wildfire) (NTEMS)
# - 5-10 plots mature stands and maybe some old seral (NTEMS)
# 5-10 deciduous (NTEMS)
# 5-10 wetlands (NTEMS)

# for model building 
# HLS spectral
# ALS structure
# NTEMS 


