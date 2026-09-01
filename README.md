# Ruido: An R Package for Calculating Spectral Metrics for Ecoacoustics Research

<img src="man/figures/ruidoIcon.png" alt="Icon of Ruido" align="right" height="300"/>

[![CRANStatusBadge](https://www.r-pkg.org/badges/version-ago/Ruido)](https://cran.r-project.org/package=Ruido) ![packageDownloads](https://cranlogs.r-pkg.org/badges/grand-total/Ruido?color=blue) [![R-CMD-check](https://github.com/Arthurigorr/Ruido/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Arthurigorr/Ruido/actions/workflows/R-CMD-check.yaml) ![lastGitCommit](https://img.shields.io/github/last-commit/Arthurigorr/Ruido) [![Codecov test coverage](https://codecov.io/gh/Arthurigorr/Ruido/graph/badge.svg)](https://app.codecov.io/gh/Arthurigorr/Ruido)

`Ruido` is an `R` package that aims to provide a simple and accessible framework for calculating spectral soundscape metrics that describes noise dynamics. It provides accessible tools for calculating less common, but ecologically meaningful soundscape metrics, helping researchers move beyond standard and classic indices.

### Currently, the package implements methods to estimate:

-   **Background Noise (BGN)** and **Soundscape Power (POW)**, following Towsey et al. 2017
-   **Soundscape Saturation (SAT)**, following Burivalova et al. 2018
-   **Spectral Acoustic Complexity Index (ACI)**, following Pieretti, et al. 2011
-   **Spectral Temporal Entropy Index (ENT)**, following Towsei et al. 2017

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

## Examples:

To illustrate the package's use, we are going to use the recordings available at: <https://zenodo.org/records/17243660>. <br>Use <https://zenodo.org/records/17575795> instead if you want to download lighter files.

If you wish to temporary download the files using R to follow the examples, run:

``` r
options(timeout = 500)

dir <- tempdir()
recName <- paste0("GAL24576_20250401_", sprintf("%06d", seq(0, 230000, by = 10000)), ".wav")
recDir <- paste(dir, recName, sep = "/")

for (rec in recName) {
  print(rec)
  url <- paste0("https://zenodo.org/records/17243660/files/",
                rec,
                "?download=1")
  download.file(url, destfile = paste(dir, rec, sep = "/"), mode = "wb")
}
```

These examples use ggplot2 and patchwork to plot their results. Before running them, first run:

``` r
library(ggplot2)
library(patchwork)
```

### Background Noise (BGN) and Soundscape Power (POW)

``` r
BGN_POW <- lapply(recDir, bgNoise)

time <- sapply(strsplit(recName, "_"), function(x)
  paste(substr(x[3], 1, 2), substr(x[3], 3, 4), substr(x[3], 5, 6), sep = ":"))
date <- sapply(strsplit(recName, "_"), function(x)
  paste(substr(x[2], 1, 4), substr(x[2], 5, 6), substr(x[2], 7, 8), sep = "-"))

dateTime <- as.POSIXct(paste(date, time))
sampRate <- BGN_POW[[1]]@sampRate
kHz <- cumsum(c(0, rep(sampRate / 6, 6))) / 1000
breaks <- round(c(1, cumsum(rep(256 / 6, 6))))

timeLabels <- time[c(1, 7, 13, 19, 24)]
timeBreaks <- as.character(dateTime[c(1, 7, 13, 19, 24)])

plotList <- list()
plotN <- 1

for (ind in c("BGN", "POW")) {
  for (cha in c("left", "right")) {
    core <- do.call(cbind, lapply(BGN_POW, function(x) {
      x@values[[cha]][[ind]]
    }))
    
    sDim <- dim(core)
    
    coreDf <- data.frame(
      TIME = as.character(rep(dateTime, each = sDim[1] * 3) + rep(rep(c(0, 60, 120), each = sDim[1]), sDim[2] / 3)),
      SPEC = rep(seq(sDim[1]), sDim[2]), VAL = c(unlist(core))
    )
    
    plotList[[plotN]] <- ggplot(coreDf, aes(x = TIME, y = SPEC, fill = VAL)) +
      geom_tile() +
      theme_classic() +
      scale_y_continuous(expand = c(0, 0), labels = kHz, breaks = breaks) +
      scale_x_discrete(expand = c(0, 0), labels = timeLabels, breaks = timeBreaks) +
      scale_fill_viridis_c(option = "magma", name = ind) +
      labs( x = "Time of Day", y = "Frequency (kHz)", title = paste(ind, "in the", cha, "channel")
      )
    
    plotN <- plotN + 1
    
  }
}

plotList[[1]] + plotList[[2]] + plotList[[3]] + plotList[[4]]
```

<img src="man/figures/BGN&amp;POW.png" align="center"/>

### Soundscape Saturation (SAT)

``` r
sat <- soundSat(dir)
SAT <- sat$values

satForPlot <- cbind(
  aggregate(SAT ~ AUDIO + CHANNEL, data = SAT, sd),
  aggregate(SAT ~ AUDIO + CHANNEL, data = SAT, mean)$SAT,
  TIME = rep(substr(time, 1, 5), 2)
)
colnames(satForPlot)[c(3, 4)] <- c("sdSAT", "meanSAT")

ggplot(
  satForPlot,
  aes(x = TIME, y = meanSAT * 100, group = CHANNEL, fill = CHANNEL,
      ymin = pmax(meanSAT - sdSAT, 0) * 100, ymax = pmin(meanSAT + sdSAT, 100) * 100
  )
) +
  geom_ribbon(alpha = 0.5) +
  geom_line() +
  geom_point() +
  theme_classic() +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_x_discrete( expand = c(0, 0), breaks = c("00:00", "06:00", "12:00", "18:00", "23:00")
  ) +
  labs(y = "Soundscape Saturation (%)") +
  theme(
    axis.title.x = element_blank(), axis.text = element_text(size = 15),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15)
  ) +
  guides(fill = guide_legend(title = "Side"))
```

<img src="man/figures/SAT.png" align="center"/>

### Acoustic Activity

``` r
act <- multActivity(dir, powthr = sat$powthresh, bgnthr = sat$bgnthresh / 100)

time <- sapply(strsplit(recName, "_"), function(x)
  paste(substr(x[3], 1, 2), substr(x[3], 3, 4), substr(x[3], 5, 6), sep = ":"))
date <- sapply(strsplit(recName, "_"), function(x)
  paste(substr(x[2], 1, 4), substr(x[2], 5, 6), substr(x[2], 7, 8), sep = "-"))

dateTime <- as.POSIXct(paste(date, time))
sampRate <- act$info$SAMPRATE[[1]]
kHz <- cumsum(c(0, rep(sampRate / 6, 6))) / 1000
breaks <- round(c(1, cumsum(rep(256 / 6, 6))))

timeLabels <- time[c(1, 7, 13, 19, 24)]
timeBreaks <- as.character(dateTime[c(1, 7, 13, 19, 24)])

plotList <- list()
plotN <- 1

for (cha in c("left", "right")) {
  actCurrent <- act$values[, act$info$CHANNEL == cha]
  actCurrentDF <- data.frame(
    TIME = as.character(rep(dateTime, each = sDim[1] * 3) + rep(rep(c(0, 60, 120), each = sDim[1]), sDim[2] / 3)),
    SPEC = rep(seq(sDim[1]), sDim[2]),
    VAL = factor(c(unlist(actCurrent)), levels = c(0, 1))
  )
  
  plotList[[plotN]] <- ggplot(actCurrentDF, aes(x = TIME, y = SPEC, fill = VAL)) +
    geom_tile() +
    theme_classic() +
    scale_y_continuous(expand = c(NA, NA), labels = kHz, breaks = breaks) +
    scale_x_discrete(expand = c(0, 0), labels = timeLabels, breaks = timeBreaks) +
    scale_fill_manual(values = c("white", "black"), labels = c("Inactive", "Active")) +
    guides(fill = guide_legend(title = "Acoustic Activity")) +
    labs(
      x = "Time of Day",
      y = "Frequency (kHz)",
      title = paste("Acoustic Activity in the", cha, "channel")
    )
  
  plotN <- plotN + 1
  
}

plotList[[1]] + plotList[[2]] +
  plot_layout(guides = "collect")
```

<img src="man/figures/ACTIVITY.png" align="center"/>

## References

-   Burivalova, Z., Towsey, M., Boucher, T., Truskinger, A., Apelis, C., Roe, P., & Game, E. T. (2018). Using soundscapes to detect variable degrees of human influence on tropical forests in Papua New Guinea. Conservation Biology, 32(1), 205-215. <https://doi.org/10.1111/cobi.12968>
-   Pieretti, N., Farina, A., & Morri, D. (2011). A new methodology to infer the singing activity of an avian community: The Acoustic Complexity Index (ACI). Ecological Indicators, 11(3), 868–873. <https://doi.org/10.1016/j.ecolind.2010.11.005>
-   Towsey, M. W. (2017). The calculation of acoustic indices derived from long-duration recordings of the natural environment. In eprints.qut.edu.au. <https://eprints.qut.edu.au/110634/>
