#' @title Acoustic Activity Matrix
#'
#' @description Calculate the Acoustic Activity Matrix using the methodology proposed in Burivalova 2018
#'
#' @param soundfile tuneR Wave object or path to a valid audio
#' @param channel channel where the saturation values will be extract from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`.
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Defaults to `-90`, as set by Towsey 2017
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `signal::hanning(wl)` to use hanning instead. Defaults to `signal::hammning(wl)`
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges"`, `"scott"` or any numeric value (foe example = `100`). Defaults to `"FD"`
#' @param powthr single numeric value to calculate the activity matrix for soundscape power (in dB). Detauls to `10`
#' @param bgnthr single numeric value to calculate the activity matrix for background noise (in %). Detauls to `0.8`
#' @param beta how BGN thresholds are calculated. If `TRUE`, BGN thresholds are calculated using all recordings combined. If FALSE, BGN thresholds are calculated separately for each recording. Defaults to `TRUE`
#'
#' @returns This function returns a 0 and 1 matrix containing the activity for all time bins of the inputted file. The matrix's number of rows will equal to half the set window lenght (`wl`) and number of columns will equal the number of bins. Cells with the value of 1 represent the acoustically active frequency of a bin.
#'
#' @details To calculate the activity matrix, we use the methodology proposed by Burivalova 2018. We begin by applying the following formula to each time bin of the recording:
#'
#'\deqn{a_{mf} = 1\  if (BGN_{mf} > \theta_{1})\  or\  (POW_{mf} > \theta_{2});\  otherwise,\  a_{mf} = 0,}
#'
#'Where \eqn{\theta_{1}} equals the threshold of BGN values and \eqn{\theta_{2}} equals the threshold of dB values. We set 1 to active and 0 to inactive frequency windows.
#'
#'@references Burivalova, Z., Towsey, M., Boucher, T., Truskinger, A., Apelis, C., Roe, P., & Game, E. T. (2018). Using soundscapes to detect variable degrees of human influence on tropical forests in Papua New Guinea. Conservation Biology, 32(1), 205-215. https://doi.org/10.1111/cobi.12968
#'
#'@export
#'@importFrom methods is
#'@importFrom methods slot
#'@importFrom stats IQR
#'@importFrom stats quantile
#'@importFrom stats setNames
#'@importFrom stats shapiro.test
#'@importFrom nortest ad.test
#'
#' @examples
#' if (require("ggplot2")) {
#' ### Generating an artificial audio for the example
#' ## For this example we'll generate a sweep in a noisy soundscape
#' library(tuneR)
#' library(ggplot2)
#'
#' # Define parameters for the artificial audio
#' samprate <- 12050
#' dur <- 60
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
#' time <- 10
#' sat <- activity(wave, timeBin = time)
#'
#' # Now we can plot the results
#' satDim <- dim(sat)
#' numericTime <- seq(0, dur, by = time)
#' labels <- paste0(numericTime[-length(numericTime)], "-", numericTime[-1], "s")
#'
#' satDF <- data.frame(BIN = rep(paste0("BIN", seq(satDim[2])), each = satDim[1]),
#'                     WIN = rep(seq(satDim[1]), satDim[2]),
#'                     ACT = factor(c(sat), levels = c(0,1)))
#'
#' ggplot(satDF, aes(x = BIN, y = WIN, fill = ACT)) +
#'   geom_tile() +
#'   theme_bw() +
#'   scale_fill_manual(values = c("white", "black")) +
#'   scale_y_continuous(expand = c(0,0)) +
#'   scale_x_discrete(expand = c(0,0), labels = labels) +
#'   labs(x = "Time Bin", y = "Spectral Window") +
#'   guides(fill = guide_legend(title = "Activity"))
#'
#' }
activity <- function(soundfile,
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

    singSat <- BGN > BGNQ | POW > powthr

  } else {
    singSat <- sapply(1:ncol(BGN), function(t) {
      BGN[, t] > quantile(BGN[, t], bgnthr) | POW[, t] > powthr

    })

  }

  return(singSat * 1)

}
