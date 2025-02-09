# LAWS-HiC

LAWS-HiC is a post-processing tool designed to enhance the results of Hi-C peak callers. Using a *locally adaptive weighting and screening* (LAWS) approach, LAWS-HiC adjusts p-values by accounting for the local structure of data, addressing the common assumption of statistical independency among interactions made by most peak callers. This method enables more accurate identification of significant chromatin interactions.

## Installation
LAWS-HiC now is avaliable as an R package. This package can be installed from GitHub:
```r
install.packages("devtools")
devtools::install_github("lb-zhou/LAWSHiC")
```

## Usage
Please check `example.Rmd` or `example.html` for detailed instruction.