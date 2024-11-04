# LAWS-HiC

LAWS-HiC is a post-processing tool designed to enhance the results of Hi-C peak callers. Using a *locally adaptive weighting and screening* (LAWS) approach, LAWS-HiC adjusts p-values by accounting for the local structure of data, addressing the common assumption of statistical independence among interactions made by most peak callers. This method enables more accurate identification of significant chromatin interactions.

For a detailed example of how to use LAWS-HiC, please see `example.html` or `example.Rmd'.