##########################################################################################################
# Identify PlanetScope data for each site and compile scene-level metadata.
# Output: rawdata_list_<site>.csv
##########################################################################################################
library(rjson) 
library(jsonlite) 
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

## sites
jsonlist <- list.files(path=jsondir)  
sites <- sub("\\.geojson$", "", jsonlist)


#########################################################
#### identify PlanetScope data for each site

## check raw data

for (s in 1:length(sites)) {

  site <- sites[s]
  
  tgdir <- list.dirs(path=file.path(rawdir, gsub("[0-9]+", "", site)), 
                     full.names=T,
                     recursive=F)
  
  tiflist <- list.files(path=tgdir, 
                        recursive=T, 
                        full.names=T, 
                        pattern = "SR_harmonized.*\\.tif$") 

  udmlist <- list.files(path=tgdir,
                        recursive=T,
                        full.names=T,
                        pattern = glob2rx("*udm*.tif$*"))
  udmlistl <- lapply(tiflist, function(x) {
    str <- sub("_3B_AnalyticMS.*$", "", x)
    f <- grep(str, udmlist, value = TRUE)
    if (length(f) == 0) NA_character_ else f
  })
  udmlist <- unlist(udmlistl)
  
  jsonlist <- list.files(path=tgdir,
                         pattern = glob2rx("*meta*.json$*"),
                         full.names = TRUE,
                         recursive=TRUE)
  jsonlistl <- lapply(tiflist, function(x) {
    str <- sub("_3B_AnalyticMS.*$", "", x)
    f <- grep(str, jsonlist, value = TRUE)
    if (length(f) == 0) NA_character_ else f
  })
  jsonlist <- unlist(jsonlistl)
  
  keep <- !(is.na(udmlist) | is.na(jsonlist))
  tiflist   <- tiflist[keep]
  udmlist   <- udmlist[keep]
  jsonlist  <- jsonlist[keep]
  
  keep      <- !duplicated(basename(tiflist))
  tiflist   <- tiflist[keep]
  udmlist   <- udmlist[keep]
  jsonlist  <- jsonlist[keep]

  length(jsonlist)
  acdates <- sapply(strsplit(basename(tiflist),"_"), "[[", 1)
  actimes <- sapply(strsplit(basename(tiflist),"_"), "[[", 2)

  idlist <- lapply(jsonlist, function(jsonfn) {
    tryCatch(
      {
     
      json <- jsonlite::fromJSON(jsonfn)
      
      it <- json$properties$item_type
      inst <- json$properties$instrument
      pr <- json$properties$pixel_resolution
      gsd <- json$properties$gsd
      light_haze <- json$properties$light_haze_percent
      cloud <- json$properties$cloud_percent
      
      acdatetime <- json$properties$acquired
      utctime <- lubridate::ymd_hms(acdatetime, tz = "UTC")
      ksttime <- lubridate::with_tz(utctime, tzone = "Asia/Seoul")
      
      list(
        item_type = it,
        instrument = inst,
        pixel_resolution = pr,
        gsd = gsd,
        kst_time = ksttime,
        haze_percent = light_haze,
        cloud_percent = cloud
      )
      
    }, error = function(e) {
      message("fail to open json file: ", jsonfn, " → ", conditionMessage(e))
      list(
        item_type = NA,
        instrument = NA,
        pixel_resolution = NA,
        gsd = NA,
        kst_time = NA,
        haze_percent = NA,
        cloud_percent = NA
      )
    }
    )
  }
  )  

idlist <- lapply(idlist, function(x) {
  x <- lapply(x, function(v) {
    if (is.null(v)) return(NA)
    if (length(v) > 1) return(v[1])
    return(v)
  })
  return(x)
})
iddf <- do.call(rbind, lapply(idlist, as.data.frame))
iddf$dir <- dirname(tiflist) 
iddf$tif <- basename(tiflist)
iddf$udm <- basename(udmlist)
iddf$meta <- basename(jsonlist)
iddf$utc_acdate <- acdates
iddf$utc_actime <- actimes
iddf$year <- as.numeric(substr(iddf$utc_acdate,1,4))
iddf$site <- site
iddf$kst_time <- format(iddf$kst_time, "%Y-%m-%d %H:%M:%S")

write.csv(iddf, paste0(rawdir,"/rawdata_list_", site, ".csv"), row.names=F)

}

## summarize data availability
DF <- data.frame()
idfllist <- list.files(path=rawdir, pattern="\\.csv$", full.names=T)
for (d in idfllist) {
 df <- read.csv(d)
 df <- df %>% dplyr::group_by(site) %>% summarise(start=min(kst_time, na.rm=T), end=max(kst_time, na.rm=T))
 DF <- rbind(DF, df)
 }



