#' @title Single Soundscape Saturation Index
#'
#' @param soundfile tuneR Wave object or path to a valid audio
#' @param channel channel where the background noise values will be extract from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`.
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Defaults to `-90`, as set by Towsey 2017.
#' @param targetSampRate sample rate of the audios. Defaults to `NULL` to not change the sample rate. This argument is only used to down sample the audio.
#' @param wl window length of the spectrogram. Defaults to `512`.
#' @param window window used to smooth the spectrogram. Defaults to `signal::hammning(wl)`. Switch to `signal::hanning(wl)` if to use hanning instead.
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges`", `"scott"` and `100`. Defaults to `"FD"`.
#' <br>Can also be set to any number to limit or increase the amount of breaks.
#' @param powthr a single value to evaluate the activity matrix for Soundscape Power (in %dB). Defaults to `10`.
#' @param bgnthr a single value to evaluate the activity matrix for Background Noise (in %). Defaults to `0.8`
#' @param beta how BGN thresholds are calculated. If TRUE, BGN thresholds are computed using all recordings combined.
#'
#' @export
#' @returns A vector containing the saturation values for all time bins of the inputted file
#' @details  Soundscape Saturation (`SAT`) is a measure of the proportion of frequency bins that are acoustically active in a determined window of time. It was developed by Burivalova et al. 2018 as an index to test the acoustic niche hypothesis.
#' To calculate this function, first we need to generate an activity matrix for each time bin of your recording with the following formula:
#'
#'\deqn{a_{mf} = 1\  if (BGN_{mf} > \theta_{1})\  or\  (POW_{mf} > \theta_{2});\  otherwise,\  a_{mf} = 0,}
#'
#'Where \eqn{\theta_{1}} is the threshold of BGN values and \eqn{\theta_{2}} is a threshold of dB values.
#'Since we define a single threshold for both in this function, we don't have to worry about generating a saturation value for many different combinations.
#'For the selected threshold a soundscape saturation measure will be taken with the following formula:
#'
#'\deqn{S_{m} = \frac{\sum_{f = 1}^N a_{mf}}{N}}
#'
#'Since this is analyzing the soundscape saturaion of a single file, no normality tests will be done.
#'
#'@references Burivalova, Z., Towsey, M., Boucher, T., Truskinger, A., Apelis, C., Roe, P., & Game, E. T. (2018). Using soundscapes to detect variable degrees of human influence on tropical forests in Papua New Guinea. Conservation Biology, 32(1), 205-215. https://doi.org/10.1111/cobi.12968
#'
#' @examples
#'
#' oldpar <- par(no.readonly = TRUE)
#'
#' ### Generating an artificial audio for the example
#' ## For this example we'll generate a sweep in a noisy soundscape
#' library(tuneR)
#'
#' # Define parameters for the artificial audio
#' samprate <- 12050
#' dur <- 59
#' n <- samprate * dur
#'
#' # White noise
#' set.seed(413)
#' noise <- rnorm(n)
#'
#' # Linear fade-out envelope
#' fade <- seq(1, 0, length.out = n)
#'
#' # Apply fade
#' signal <- noise * fade
#'
#' # Create Wave object
#' wave <- Wave(
#'   left = signal,
#'   samp.rate = samprate,
#'   bit = 16
#' )
#'
#' # Running singleSat() on the artificial audio
#' sat <- singleSat(wave, timeBin = 10)
#'
#' # Now we can plot the results
#' # In the left we have a periodogram and in the right saturaion values
#' # along one minute
#' par(mfrow = c(1,2))
#' image(periodogram(wave, width = 8192, normalize = FALSE), xlab = "Time (s)",
#' ylab = "Frequency (hz)", axes = FALSE)
#' axis(1, labels = seq(0,60, 10), at = seq(0,7e5,length.out = 7))
#' axis(2)
#' plot(sat, xlab = "Time (s)", ylab = "Soundscape Saturation (%)",
#' type = "b", pch = 16, axes = FALSE)
#' axis(1, labels = paste0(c("0-10","10-20","20-30","30-40","40-50","50-59"),
#' "s"), at = 1:6)
#' axis(2)
#'
#' par(oldpar)
#'
#' \donttest{
#' # Getting audiofile from the online Zenodo library
#' dir <- paste(tempdir(), "forExample", sep = "/")
#' dir.create(dir)
#' rec <- paste0("GAL24576_20250401_", sprintf("%06d", 0),".wav")
#' recDir <- paste(dir,rec , sep = "/")
#' url <- paste0("https://zenodo.org/records/17575795/files/", rec, "?download=1")
#'
#' # Downloading the file, might take some time denpending on your internet
#' download.file(url, destfile = recDir, mode = "wb")
#'
#' # Now we calculate soundscape saturation for both sides of the recording
#' sat <- singleSat(recDir)
#'
#' # Printing the results
#' print(sat)
#'
#' barplot(sat, col = c("darkgreen", "red"),
#'        names.arg = c("Left", "Right"), ylab = "Soundscape Saturation (%)")
#'
#' unlink(dir, recursive = TRUE)
#' }
singleSat <- function(soundfile,
                      channel = "stereo",
                      timeBin = 60,
                      dbThreshold = -90,
                      targetSampRate = NULL,
                      wl = 512,
                      window = signal::hamming(wl),
                      overlap = ceiling(length(window) / 2),
                      histbreaks = "FD",
                      powthr = 10,
                      bgnthr = 0.8,
                      beta = TRUE) {
  if (!(channel %in% c("left", "right", "stereo", "mono")))
    stop("channel must be 'stereo', 'mono', 'left', or 'right'")

  if (!is.numeric(timeBin))
    stop("timeBin must be numeric")

  if (!is.numeric(dbThreshold))
    stop("timeBin must be numeric")

  if (!is.null(targetSampRate)) {
    if (!is.numeric(targetSampRate))
      stop("targetSampRate must be either NULL or a numeric value")
  }

  halfWl <- round(wl / 2)

  BGNPOW <- bgNoise(
    soundfile,
    timeBin = timeBin,
    targetSampRate = targetSampRate,
    window = window,
    overlap = overlap,
    channel = channel,
    dbThreshold = dbThreshold,
    wl = wl,
    histbreaks = histbreaks
  )

  nBins <- length(BGNPOW$timeBins)

  if (BGNPOW$channel == "stereo") {
    BGN <- cbind(BGNPOW$values$left$BGN, BGNPOW$values$right$BGN)
    names <- paste0(rep(c("left", "right"), each = nBins), seq(nBins))
  } else {
    BGN <- BGNPOW$values[[BGNPOW$channel]]$BGN
    names <- paste0(rep(BGNPOW$channel, nBins), seq(nBins))
  }

  if (BGNPOW$channel == "stereo") {
    POW <- cbind(BGNPOW$values$left$POW, BGNPOW$values$right$POW)
  } else {
    POW <- BGNPOW$values[[BGNPOW$channel]]$POW
  }

  if (beta) {
    BGNQ <- quantile(unlist(BGN), bgnthr)

    singSat <- colMeans(BGN > BGNQ | POW > powthr)

  } else {
    singSat <- sapply(1:ncol(BGN), function(t) {
      sum(BGN[, t] > quantile(BGN[, t], bgnthr) |
            POW[, t] > powthr) / halfWl

    })

  }

  names(singSat) <- names

  return(singSat)

}
