getSampleBins <- function(samples, samp.rate, binSize) {
  b <- seq(1, samples, by = samp.rate * binSize)
  e <- pmin(b + samp.rate * binSize - 1, samples)

  keepThese <- ((samp.rate * binSize) * 0.1) < e - b

  if (length(b) == length(e)) {
    data.frame(b, e)[keepThese, ]
  } else {
    data.frame(b, e = c(e, samples))
  }

}
