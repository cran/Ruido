# Ruido: An R Package for Profiling Background Noise and Calculating Soundscape Saturation

<img src="man/figures/ruidoIcon.png" alt="Icon of Ruido" align="right" height="300"/>

`Ruido` is an `R` package that aims to provide a simple and accessible framework for calculating less common soundscape metrics that describes noise dynamics. It provides accessible tools for calculating less common, but ecologically meaningful soundscape metrics, helping researchers move beyond standard and classic indices.

### The package implements methods to estimate:

- **Background Noise (BGN)** and **Soundscape Power (POW)**, following Towsey et al. (2014)
- **Soundscape Saturation (SAT)**, following Burivalova et al. (2021)

These metrics can be used to explore acoustic complexity, biotic activity, and environmental disturbance, making `Ruido` useful for ecological monitoring, bioacoustic surveys, or experimental soundscape studies.

## Installation

### CRAN Download:

``` r
install.packages("Ruido")
library(Ruido)
```

### Github Download:

``` r
devtools::install_github("Arthurigorr/Ruido")
library(Ruido)
```
