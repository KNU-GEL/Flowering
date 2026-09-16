###############################
library(rjson)
params <- rjson::fromJSON(file='/../Flowering/input/PFP_Parameters.json')

geojsonList <- list.files(path=params$setup$geojsonDir)
sites <- 1:length(geojsonList)

cpl <- list.files(path=file.path(params$setup$workDir,"spatdat", "sample"),
                  pattern="cp.*\\.shp$", recursive=T, full.names=T)
ssites <- 1:length(cpl)


###############################
### 00 Order and download PlanetScope imagery for study areas defined by GeoJSON files


###############################
### 01 Identify PlanetScope imagery

system(paste('qsub -V -m n -pe omp 2 -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_01.sh ', sep=''))


###############################
### 02 Define sub-site areas from center points and create polygon/raster files

system(paste('qsub -V -m n -pe omp 2 -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_02.sh ', sep=''))


###############################
### 03 Quality filtering and sub-daily time-series construction

setwd(paste0(params$setup$logDir,'03'))
for(numSite in ssites){
  nn <- sprintf('%03d',numSite)
  for(yy in params$setup$phenStartYr:params$setup$phenEndYr){
    system(paste('qsub -V -m n -pe omp 2 -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_03.sh ',nn,yy,sep=''))
  }
}


###############################
### 04 Generate flowering phenology maps

setwd(paste0(params$setup$logDir,'04'))
for(numSite in ssites){
  nn <- sprintf('%03d',numSite)
  for(yy in params$setup$phenStartYr:params$setup$phenEndYr){
    system(paste('qsub -V -m n -pe omp 2 -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_04.sh ',nn,yy,sep=''))
  }
}

