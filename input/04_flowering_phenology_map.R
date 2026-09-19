##########################################################################################################
# Detect pixel-level green-up and flowering phenology from PlanetScope time series.
# Outputs: green-up and flowering phenology maps for multiple amplitude thresholds.
##########################################################################################################

library(sf) 
library(terra) 
library(raster) 
library(stringr)
library(lubridate)
library(dplyr)
library(tidyr)
library(purrr)
library(future)
library(furrr)
library(zoo)


###############################
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

## site, sub-site, id
fl <- list.files(path=spatdatdir, pattern="cp.*\\.shp$", recursive=T, full.names=T)  
ssites <- sapply(strsplit(basename(fl), "_"), "[[",1)
sps <- sapply(strsplit(basename(fl), "_"), "[[",2)
ids <- paste(ssites, sps, sep="_")
ssite <- ssites[numSite]
sp <- sps[numSite]; sg <- sp
print(sp)
id <- ids[numSite]
print(id)

## year
yr <- yy

## id data frame
print(ssite)
site <- gsub("[0-9]+", "", ssite)
print(site)
idfl <- list.files(path=rawdir, pattern=paste0("rawdata_list.*", site, ".*\\.csv$"), recursive=T, full.names=T)
iddf <- read.csv(idfl)

## index
VI <- "EVI2"  
FI <- "EBI"  

## detection parameter
phenopar    <- params$phenology_parameters

## haze threshold
hthr <- phenopar$hazeThresh

fampthr <- seq(0, 0.03, 0.005)


##################################################################
#### input spatial data

## base tif
tgtiflist <- list.files(path = file.path(spatdatdir, "sample", "base_img_tif"),
                        pattern = paste0(id, ".*\\.tif$"),
                        full.names = TRUE, recursive=T)
r <- rast(tgtiflist)
nc <- ncell(r)


##################################################################
#### aggregate to daily data

print(paste0("year: ", yr))

## input sub-daily data
tsfl <- list.files(path = file.path(outputdir, id, "tsmat"),
                   pattern =  "daily_ts",
                   recursive= T,
                   full.names = T)
tsfl <- tsfl[grepl("sub", tsfl)] # sub-daily data
tsflyr <- tsfl[grepl(yr, basename(tsfl))]
print(tsflyr)
load(tsflyr)     # input: times, red, green, blue, nir, evi2, ebi


## filter hazy images
hazetimes <- iddf %>% filter(year==yr) %>% 
  filter(haze_percent>hthr) %>% 
  pull(kst_time) %>% as.POSIXct() %>% format("%Y%m%d_%H%M%S")
hazeidx <- match(hazetimes, times)

red[,hazeidx]<-NA; green[,hazeidx]<-NA; blue[,hazeidx]<-NA; nir[,hazeidx]<-NA
evi2[,hazeidx]<-NA; ebi[,hazeidx]<-NA
cleartimes <- times
cleartimes[c(hazeidx)] <- NA
print(paste0("Num of haze-filtered images: ", length(hazeidx)))


## additional filter: image-level high EBI images
img_qt  <- apply(red+green+blue, 2, quantile, probs = 0.01, na.rm = TRUE)
thr_qt <- median(img_qt, na.rm = TRUE) + 2 * mad(img_qt, na.rm = TRUE)
qtidx <- which(img_qt > thr_qt)
cleartimes[qtidx]
red[, qtidx]   <- NA
green[, qtidx] <- NA
blue[, qtidx]  <- NA
nir[, qtidx]   <- NA
evi2[, qtidx]  <- NA
ebi[, qtidx]   <- NA
print(paste0("Num of q01-filtered images: ", length(qtidx)))
cleartimes[c(qtidx)] <- NA


## update clear times after both filters
dates <- unique(substr(cleartimes[!is.na(cleartimes)],1,8))
print(paste0("Num of dates: ", length(dates)))


## aggregate
yyyymmdd <- substr(times, 1, 8)   
to_daily <- function(mat, yyyymmdd, dates) {
  daily_list <- lapply(dates, function(d) {
    # d <- dates[1]
    # mat <- evi2
    cols <- which(yyyymmdd == d)
    if (length(cols) == 1) {
      mat[, cols, drop = FALSE]
    } else {
      matrix(rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE), ncol = 1)
    }
  })
  out <- do.call(cbind, daily_list)
  out
}
RED   <- to_daily(red,  yyyymmdd, dates)
GREEN <- to_daily(green,yyyymmdd, dates)
BLUE  <- to_daily(blue, yyyymmdd, dates)
NIR   <- to_daily(nir,  yyyymmdd, dates)
EVI2  <- to_daily(evi2, yyyymmdd, dates)
EBI   <- to_daily(ebi,  yyyymmdd, dates)

RED[is.nan(RED)] <- NA
GREEN[is.nan(GREEN)] <- NA
BLUE[is.nan(BLUE)] <- NA
NIR[is.nan(NIR)] <- NA
EVI2[is.nan(EVI2)] <- NA
EBI[is.nan(EBI)]  <- NA


###################################################################
#### detect flowering phenology

dates <- as.Date(dates, format = "%Y%m%d")
vi <- get(VI)
fi <- get(FI)

vgonset  <- rep(NA_real_, nc)
vgmid    <- rep(NA_real_, nc)
vgend    <- rep(NA_real_, nc)
vfstart  <- rep(NA_real_, nc)
vfpeak   <- rep(NA_real_, nc)
vfend    <- rep(NA_real_, nc)
vfamp    <- rep(NA_real_, nc)
vvamp    <- rep(NA_real_, nc)
vvmax    <- rep(NA_real_, nc)
vvmin    <- rep(NA_real_, nc)

for (cid in 1:nc) {
  
  vi_sub <- vi[cid,]
  fi_sub <- fi[cid,]
  vi_sub <- data.frame(date=dates, value=vi_sub)
  fi_sub <- data.frame(date=dates, value=fi_sub)

  ## extract green-up phenophase
  gdetecs <- detect_greenup(dfvi=vi_sub, dffi=fi_sub, phenopar)
  green_date <- gdetecs$phase
  vamp <- gdetecs$amp
  vmax <- gdetecs$max
  vmin <- gdetecs$min
  gonset <- green_date$date[1]
  gmid <- green_date$date[2]
  gend <- green_date$date[3]
  
  ## smoothing and find flowering peak
  fitresult <- detect_flower(dfvi=vi_sub, dffi=fi_sub, sg, phenopar)
  fpeak <- fitresult$peak$peak_date
  fstart <- fitresult$peak$start_date
  fend <- fitresult$peak$end_date
  famp <-  fitresult$peak$peak_amp

  vgonset[cid]  <- yday(gonset)
  vgmid[cid]    <- yday(gmid)
  vgend[cid]    <- yday(gend)
  vfstart[cid]  <- yday(fstart)
  vfpeak[cid]   <- yday(fpeak)
  vfend[cid]    <- yday(fend)
  vfamp[cid]    <- famp
  vvamp[cid]    <- vamp
  vvmax[cid]    <- vmax
  vvmin[cid]    <- vmin

}

lyrnames <- c("Gonset", "Gmid", "Gend",
              "Fstart", "Fpeak", "Fend",
              "Famp", "Vamp",
              "Vmax", "Vmin")
out.files <- file.path(outputdir, id, "map", "phe", 
                       paste0("haze_", phenopar$hazeThresh,
                              "_spar_",phenopar$WsplineSpar,
                              "_famp_",phenopar$FIampThresh),
                       paste0(id,"_",yr,"_",lyrnames,"_",paste0("famp","_",phenopar$FIampThresh),
                              ".tif"))
dir.create(unique(dirname(out.files)), recursive = TRUE, showWarnings = FALSE)

writeRaster(setValues(rast(r), vgonset), out.files[1], overwrite = TRUE)
writeRaster(setValues(rast(r), vgmid),   out.files[2], overwrite = TRUE)
writeRaster(setValues(rast(r), vgend),   out.files[3], overwrite = TRUE)
writeRaster(setValues(rast(r), vfstart), out.files[4], overwrite = TRUE)
writeRaster(setValues(rast(r), vfpeak),  out.files[5], overwrite = TRUE)
writeRaster(setValues(rast(r), vfend),   out.files[6], overwrite = TRUE)
writeRaster(setValues(rast(r), vfamp),   out.files[7], overwrite = TRUE)
writeRaster(setValues(rast(r), vvamp),   out.files[8], overwrite = TRUE)
writeRaster(setValues(rast(r), vvmax),   out.files[9], overwrite = TRUE)
writeRaster(setValues(rast(r), vvmin),   out.files[10], overwrite = TRUE)


#### apply flowering amplitude thresholds
for (thr in fampthr) {
  
  if (thr==phenopar$FIampThresh) {
    next
  }
  
  out.files <- file.path(outputdir, id, "map", "phe", 
                         paste0("haze_", phenopar$hazeThresh,
                                "_spar_",phenopar$WsplineSpar,
                                "_famp_",thr),
                         paste0(id,"_",yr,"_",lyrnames,"_",paste0("famp","_",thr),
                                ".tif"))
  dir.create(unique(dirname(out.files)), recursive = TRUE, showWarnings = FALSE)

  vgonsettp <- vgonset
  vgmidtp   <- vgmid
  vgendtp   <- vgend
  vfstarttp <- vfstart
  vfpeaktp  <- vfpeak
  vfendtp   <- vfend
  vfamptp   <- vfamp
  vvamptp   <- vvamp
  vvmaxtp   <- vvmax
  vvmintp   <- vvmin

  idx <- vfamp < thr

  vgonsettp[idx] <- NA
  vgmidtp[idx]   <- NA
  vgendtp[idx]   <- NA
  vfstarttp[idx] <- NA
  vfpeaktp[idx]  <- NA
  vfendtp[idx]   <- NA
  vfamptp[idx]   <- NA
  vvamptp[idx]   <- NA
  vvmaxtp[idx]   <- NA
  vvmintp[idx]   <- NA
  
  writeRaster(setValues(rast(r), vgonsettp), out.files[1],  overwrite = TRUE)
  writeRaster(setValues(rast(r), vgmidtp),   out.files[2],  overwrite = TRUE)
  writeRaster(setValues(rast(r), vgendtp),   out.files[3],  overwrite = TRUE)
  writeRaster(setValues(rast(r), vfstarttp), out.files[4],  overwrite = TRUE)
  writeRaster(setValues(rast(r), vfpeaktp),  out.files[5],  overwrite = TRUE)
  writeRaster(setValues(rast(r), vfendtp),   out.files[6],  overwrite = TRUE)
  writeRaster(setValues(rast(r), vfamptp),   out.files[7],  overwrite = TRUE)
  writeRaster(setValues(rast(r), vvamptp),   out.files[8],  overwrite = TRUE)
  writeRaster(setValues(rast(r), vvmaxtp),   out.files[9],  overwrite = TRUE)
  writeRaster(setValues(rast(r), vvmintp),   out.files[10], overwrite = TRUE)

}

