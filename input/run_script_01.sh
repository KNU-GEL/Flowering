#!/bin/bash

# module load python3/3.7.7
# module load gdal/3.1.2
# module load R 

echo Submitting $1
R --vanilla < /../Flowering/input/01_identify_planetscope_data.R


