##########################################################################################################
# Create sub-site ROIs from center points.
# Outputs: ROI polygon shapefiles and base rasters.
##########################################################################################################
library(rjson) 
library(terra) 
library(sf)
library(raster) 
library(stringr)
library(lubridate)
library(dplyr)


#########################################################
params <- rjson::fromJSON(file='/../Flowering/input/PFP_Parameters.json')
print(params$setup$rFunctions)
source(params$setup$rFunctions)


#########################################################
#### setup

## directory
curdir <- params$setup$workDir
rawdir <- params$setup$dataDir
spatdatdir <- file.path(curdir, "spatdat")
jsondir <- file.path(spatdatdir, "geojson")

## sites, sub-sites, ids
fl <- list.files(path=spatdatdir, pattern="cp.*\\.shp$", recursive=T, full.names=T)  
ssites <- sapply(strsplit(basename(fl), "_"), "[[",1)
sps <- sapply(strsplit(basename(fl), "_"), "[[",2)
ids <- paste(ssites, sps, sep="_")

## square side length
sl <- 400    # unit: m


#########################################################
#### create sub-site ROIs

for (s in seq(ssites)) {
  
  # sub-site name
  ssite <- ssites[s]
  
  # id
  id <- ids[s]
  
  # site name
  site <- gsub("[0-9]+", "", ssite)
  
  # base tif
  tgdir <- list.dirs(file.path(rawdir,site), full.names=T, recursive=F)
  tiflist <- list.files(path=tgdir, pattern = "SR_harmonized.*\\.tif$", recursive=T, full.names=T) 
  basetif <- tiflist[which.max(file.info(tiflist)$size)[1]]
  basetif <- rast(basetif)
  r <- basetif[[1]]

  # center point of sub-site
  cpfl <- list.files(spatdatdir, pattern= paste0(id, ".*cp.*\\.shp$"), recursive=T, full.names=T)
  cp <- st_read(cpfl)
  cp_utm <- st_transform(cp, crs(r))
  v <- vect(cp_utm)
  cell_id <- terra::cellFromXY(r, terra::crds(v))
  cell_xy <- terra::xyFromCell(r, cell_id)
  out <- data.frame(
    id = ids[s],
    cell = cell_id,
    cell_x = cell_xy[, 1],
    cell_y = cell_xy[, 2]
  )

  # sub-site ROI polygon
  res_m <- res(r)[1]
  half <- sl / 2
  n <- half %/% res_m   
  xmin <- out$cell_x - (n + 0.5) * res_m
  xmax <- out$cell_x + (n + 0.5) * res_m
  ymin <- out$cell_y - (n + 0.5) * res_m
  ymax <- out$cell_y + (n + 0.5) * res_m
  samp_plg <- sf::st_as_sfc(
    mapply(
      function(a, b, c, d) sf::st_polygon(list(rbind(
        c(a, c), c(b, c), c(b, d), c(a, d), c(a, c)
      ))),
      xmin, xmax, ymin, ymax,
      SIMPLIFY = FALSE
    ),
    crs = sf::st_crs(terra::crs(r))
  )
  samp_plg_sf <- sf::st_sf(
    id = out$id,
    cell = out$cell,
    geometry = samp_plg
  )
  st_write(
    samp_plg_sf,
    file.path(spatdatdir, "sample", "plg_shp", paste0(ids[s], "_plg.shp")),
    delete_layer = TRUE
    )
  
  # sub-site ROI raster
  samp_plg_v <- terra::vect(samp_plg_sf)
  samp_r <- terra::crop(r, samp_plg_v)
  samp_r <- terra::mask(samp_r, samp_plg_v)
  values(samp_r) <- 1
  names(samp_r) <- "roi"
  writeRaster(
    samp_r,
    file.path(spatdatdir, "sample", "base_img_tif", paste0(ids[s],".tif")),
    overwrite = TRUE
  )

}

