# PlanetScope Flowering Phenology Detection

This repository contains the code used to detect tree flowering phenology from PlanetScope time series.

The workflow combines EVI2-derived canopy development stages with flowering signals derived from the Enhanced Bloom Index (EBI). Phenology-based detection windows are used to identify flowering peaks for tree species that flower before, during, or after canopy development.

## Workflow

1. Identify PlanetScope imagery and organize scene-level metadata.
2. Define sub-site regions of interest (ROIs).
3. Apply quality filtering and construct sub-daily PlanetScope time series.
4. Aggregate observations to daily time series.
5. Detect green-up phenophases from EVI2.
6. Define flowering detection windows according to flowering strategy.
7. Detect flowering peaks from EBI.
8. Generate pixel-level green-up and flowering phenology maps.

## Repository structure

```text
Flowering/
├── input/
│   ├── 01_identify_planetscope_data.R
│   ├── 02_make_sample.R
│   ├── 03_mask_and_make_ts_by_year.R
│   ├── 04_flowering_phenology_map.R
│   ├── PFP_Functions.R
│   ├── PFP_Parameters.json
│   └── run_script_*.sh
├── raw/
│   └── Example metadata file
└── spatdat/
    ├── geojson/
    └── sample/
```

## Main indices

**EVI2** (Jiang et al., 2008)

Used to characterize annual canopy development and derive green-up phenophases.

**EBI** (Chen et al., 2019)

Used to characterize flowering-related spectral signals and identify flowering peaks within phenology-based detection windows.

## Input data

The workflow was developed for PlanetScope Surface Reflectance imagery with associated quality-mask information.

PlanetScope imagery is not distributed with this repository because of data licensing restrictions. Users must obtain PlanetScope imagery separately.

Example spatial input files and metadata structures are included to illustrate the required input format.

## Parameters

Processing and flowering-detection parameters are specified in:

```text
input/PFP_Parameters.json
```

The main functions for green-up and flowering detection are provided in:

```text
input/PFP_Functions.R
```

## Usage

The processing workflow is organized sequentially from `01` to `04`.

```text
01_identify_planetscope_data.R
02_make_sample.R
03_mask_and_make_ts_by_year.R
04_flowering_phenlogy_map.R
```

The shell scripts in the `input/` directory are used for batch processing on a computing cluster.

File paths and processing parameters should be updated in `PFP_Parameters.json` before running the workflow.

## Citation

Kim, Sukyung and Kim, Kyu Rang and Han, Young Jong and Kim, Hyun Seok and Moon, Minkyu, Detecting tree flowering across contrasting canopy development stages using PlanetScope imagery. Available at SSRN: https://ssrn.com/abstract=7484765 or http://dx.doi.org/10.2139/ssrn.7484765


