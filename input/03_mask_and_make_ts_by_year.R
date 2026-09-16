##########################################################################################################
# Apply quality filtering and extract sub-daily PlanetScope time series for each sub-site.
# Outputs: red, green, blue, NIR, EVI2, and EBI time-series matrices.
##########################################################################################################
library(rjson) 
library(terra)
library(sf)
library(raster) 
library(stringr)
library(lubridate)
library(dplyr)
library(tidyr)


#########################################################
args <- commandArgs()
print(args)

numSite <- as.numeric(substr(args[3],1,3))
yy      <- as.numeric(substr(args[3],4,7)) 
# numSite <- 1; yy <- 2025


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
outputdir <- file.path(curdir, "output")

## temporary file settings
tmp <- "/../Flowering/_tmp"
if (!dir.exists(tmp)) dir.create(tmp, recursive = TRUE)
prjname <- basename(normalizePath(params$setup$workDir))

## site, sub-site, id
fl <- list.files(path=file.path(spatdatdir, "sample"), 
                 pattern="cp.*\\.shp$", recursive=T, full.names=T)  
ssites <- sapply(strsplit(basename(fl), "_"), "[[",1)
sps <- sapply(strsplit(basename(fl), "_"), "[[",2)
ids <- paste(ssites, sps, sep="_")
ssite <- ssites[numSite]
sp <- sps[numSite]
id <- ids[numSite]
print(id)


#########################################################
#### extract sample pixels within the sub-site ROI

## extract points from the base raster
tgtiflist <- list.files(path = file.path(spatdatdir,"sample"),
                        pattern = paste0(id, ".*\\.tif$"),
                        full.names = TRUE, recursive=T)
obj <- rast(tgtiflist)
xy <- xyFromCell(obj, which(!is.na(values(obj))))
pts_vect <- vect(xy, crs = crs(obj))
pts_vect$subsite <- ssite
pts_vect$species <- sp
pts_vect$id <- id
pts_vect$pid <- seq_len(nrow(pts_vect))

# # save as shp file (optional)
# out.file <- file.path(spatdatdir, "sample", "pts_shp", paste0(id,"_sampled_pts", ".shp"))
# if (!dir.exists(dirname(out.file))) {
#   dir.create(dirname(out.file), recursive = TRUE)
# }
# terra::writeVector(pts_vect, out.file, overwrite = TRUE)


##############################################################
#### calculate vegetation and flowering indices

## EBI
f_EBI <- function(red, green, blue){
  tryCatch(
    {
      ebi <- (red+green+blue)/((green/blue)*(red-blue+1))
      return(ebi)
    }, error = function(e) NULL
  )
}

## EVI2
f_EVI2 <- function(nir, red){
  tryCatch(
    {
      evi2 <- 2.5*(nir-red)/(nir+(2.4*red)+1)
      return(evi2)
    }, error = function(e) NULL
  )
}

##############################################################
#### mask pixels and make time series (by year) 

print(ssite)
print(id)
site <- gsub("[0-9]+", "", ssite)
print(site)

## roi
tgplg <- list.files(file.path(spatdatdir,"sample","plg_shp"),
                    pattern=paste0(id, ".*\\.shp$"), full.names=T)
roi <- terra::vect(tgplg)

## id data frame
idfl <- list.files(path=rawdir, pattern=paste0("rawdata_list.*", site, ".*\\.csv$"), recursive=T, full.names=T)
iddf <- read.csv(idfl)


terraOptions(tempdir = tmp)
rasterOptions(tmpdir = tmp)

subiddf <- iddf[iddf$year==yy,]
subiddf <- subiddf[!is.na(subiddf$item_type),]

tiflist <- paste(subiddf$dir, subiddf$tif, sep="/")
udmlist <- paste(subiddf$dir, subiddf$udm, sep="/")

times <- format(
  as.POSIXct(subiddf$kst_time, tz = "Asia/Seoul"),
  "%Y%m%d_%H%M%S"
)

## output matrix setup
nr <- dim(pts_vect)[1]
nc <- length(times)
b1 <- matrix(NA_real_, nrow = nr, ncol = nc) # red
b2 <- matrix(NA_real_, nrow = nr, ncol = nc) # green
b3 <- matrix(NA_real_, nrow = nr, ncol = nc) # blue
b4 <- matrix(NA_real_, nrow = nr, ncol = nc) # nir
v1 <- matrix(NA_real_, nrow = nr, ncol = nc) # EVI2
v2 <- matrix(NA_real_, nrow = nr, ncol = nc) # EBI

for (i in seq(tiflist)) {
  
  # i=1
  tiffn <- tiflist[i]
  udmfn <- udmlist[i]
  
  cat("processing: i=",i,", ", tiffn, "\n")
  
  bandall <- tryCatch(rast(tiffn), error = function(e) NULL)
  udmall <- tryCatch(rast(udmfn), error = function(e) NULL)
  
  if (is.null(bandall) || is.null(udmall)) {
    next
  }
  
  udmb1 <- udmall$clear # 0, 1
  udmb8 <- udmall$udm1 # 0-255
  
  clbandall <- tryCatch({
    terra::mask(bandall, udmb1, maskvalues=0,
                filename=tempfile(
                  pattern = sprintf("%s_%s_%s_", prjname, site, "clear"),
                  tmpdir  = terraOptions()$tempdir,
                  fileext = ".tif"
                ),
                overwrite = TRUE)
  }, error = function(e) {
    stop(e)
  })
  
  clbandall <- tryCatch({
    badmask <- udmb8 != 0
    terra::mask(clbandall, badmask, maskvalues=TRUE,
                filename=tempfile(
                  pattern = sprintf("%s_%s_%s_", prjname, site, "masked_bad"),
                  tmpdir  = terraOptions()$tempdir,
                  fileext = ".tif"
                ),
                overwrite = TRUE)
  }, error = function(e) {
    stop(e)
  })
  
  clbandall <- clbandall/10000
  
  if (!identical(crs(roi), crs(clbandall))) {
    roi <- project(roi, terra::crs(clbandall))
  }

  interarea <- terra::intersect(terra::ext(clbandall), terra::ext(roi))
  
  if (!is.null(interarea)) {
    clbandall <- crop(clbandall, roi)
    clbandall <- mask(clbandall, roi)
    
    blue <- clbandall$blue
    green <- clbandall$green
    red <- clbandall$red
    nir <- clbandall$nir
    EVI2 <- f_EVI2(nir, red)
    EBI <- f_EBI(red, green, blue)
    
  } else {
    next
  }
  
  bstack <- c(red, green, blue, nir, EVI2, EBI)
  names(bstack) <- c("red", "green", "blue", "nir", "EVI2", "EBI")
  
  ex_df <- terra::extract(bstack, pts_vect, df = TRUE)
  
  b1[,i] <- ex_df$red
  b2[,i] <- ex_df$green
  b3[,i] <- ex_df$blue
  b4[,i] <- ex_df$nir
  v1[,i] <- ex_df$EVI2
  v2[,i] <- ex_df$EBI
  
}

red <- b1    # red
green <- b2  # green
blue <- b3   # blue
nir <- b4    # nir
evi2 <- v1   # EVI2
ebi <- v2    # EBI

out.file <- file.path(outputdir, id, "tsmat", paste0(id, "_sub_daily_ts_", yy, ".rda"))
if (!dir.exists(dirname(out.file))) {
  dir.create(dirname(out.file), recursive = TRUE)
}

save(times, red, green, blue, nir, evi2, ebi, file = out.file)  


