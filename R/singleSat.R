#' @title Single Soundscape Saturation Index
#'
#' @param soundfile tuneR Wave object or path to a valid audio
#' @param channel channel where the background noise values will be extract from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`.
#' @param timeBin size (in seconds) of the time bin. Defaults to `60`.
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Defaults to `-90`, as set by Towsey 2017.
#' @param targetSampRate sample rate of the audios. Defaults to `NULL` to not change the sample rate. This argument is only used to down sample the audio.
#' @param wl window length of the spectrogram. Defaults to `512`.
#' @param window window used to smooth the spectrogram. Defaults to `signal::hammning(wl)`. Switch to `signal::hanning(wl)` if to use hanning instead.
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges`", `"scott"` and `100`. Defaults to `"FD"`.
#' <br>Can also be set to any number to limit or increase the amount of breaks.
#' @param powthr a single value to evaluate the activity matrix for Soundscape Power (in %dB). Defaults to `10`.
#' @param bgnthr a single value to evaluate the activity matrix for Background Noise (in %). Defaults to `0.8`
#'
#' @export
#' @returns A data frame containing the saturation values for all time bins of the inputed file
#' @details  Soundscape Saturation (`SAT`) is a measure of the proportion of frequency bins that are acoustically active in a determined window of time. It was developed by Burivalova et al. 2017 as an index to test the acoustic niche hypothesis.
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
#'@references Burivalova, Z., Towsey, M., Boucher, T., Truskinger, A., Apelis, C., Roe, P., & Game, E. T. (2017). Using soundscapes to detect variable degrees of human influence on tropical forests in Papua New Guinea. Conservation Biology, 32(1), 205-215. https://doi.org/10.1111/cobi.12968
#'
#' @examples
#'
#'### Generating an artificial audio for the example
#' ## For this example we'll generate a sweep in a noisy soundscape
#' library(tuneR)
#'
#' # Define the audio sample rate, duration and number of samples
#' samprate <- 12050
#' dur <- 59
#'
#' # Create a time vector
#' t <- seq(0, dur, by = 1/samprate)
#'
#' # It starts at the frequency 50hz and goes all the way up to 4000hz
#' # The sweep is exponnential
#' freqT <- 50 * (4000/50)^(t/dur)
#'
#' # Generate the signal
#' signal1 <- sin(10 * pi * cumsum(freqT) / samprate)
#' # We create an envelope to give the sweep a fade away
#' envelope <- exp(-4 * t / dur)
#' # Generating low noise for the audio
#' set.seed(413)
#' noise <- rnorm(length(t), sd = 0.3)
#' # Adding everything together in our signal
#' signal <- signal1 * envelope + noise
#'
#' # Normalize to 16-bit WAV range to create our Wave object
#' signalNorm <- signal / max(abs(signal))
#' wave_obj <- Wave(left = signalNorm, samp.rate = samprate, bit = 16)
#'
#' # Now we calculate soundscape saturation for our audio
#' # Here we are using timeBin = 10 so we get Soundscape Saturation values
#' # every 10 seconds on the audio
#' SAT <- singleSat(wave_obj, timeBin = 10)
#'
#' # Now we can plot the results
#' # In the left we have a periodogram and in the right saturaion values
#' # along one minute
#' par(mfrow = c(1,2))
#' image(periodogram(wave_obj, width = 64), xlab = "Time (s)",
#' ylab = "Frequency (hz)", log = "y", axes = FALSE)
#' axis(1, labels = seq(0,60, 10), at = seq(0,7e5,length.out = 7))
#' axis(2)
#' plot(SAT$left, xlab = "Time (s)", ylab = "Soundscape Saturation (%)",
#' type = "b", pch = 16, axes = FALSE)
#' axis(1, labels = paste0(c("0-10","10-20","20-30","30-40","40-50","50-59"),
#' "s"), at = 1:6)
#' axis(2)
#'
#'\donttest{
#'
#'oldpar <- par(no.readonly = TRUE)
#'
#'# Getting audiofile from the online Zenodo library
#' dir <- tempdir()
#' rec <- paste0("GAL24576_20250401_", sprintf("%06d", 0),".wav")
#' recDir <- paste(dir,rec , sep = "/")
#' url <- paste0("https://zenodo.org/records/17575795/files/", rec, "?download=1")
#'
#' # Downloading the file, might take some time denpending on your internet
#' download.file(url, destfile = recDir, mode = "wb")
#'
#' # Now we calculate soundscape saturation for both sides of the recording
#' sat <- singleSat(recDir, wl = 256)
#'
#' # Printing the results
#' print(sat)
#'
#' barplot(c(sat$left, sat$right), col = c("darkgreen", "red"),
#'        names.arg = c("Left", "Right"), ylab = "Soundscape Saturation (%)")
#'
#' unlink(recDir)
#' par(oldpar)
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
                      bgnthr = 0.8) {
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

  if (all(c("left", "right") %in% names(BGNPOW))) {
    BGNsaturation <- sapply(c("left", "right"), function(side) {
      list(apply(BGNPOW[[side]]$BGN, 2, function(BGN) {
        Q <- quantile(BGN, bgnthr)
        BGN > Q
      }))
    })

    POWsaturation <- sapply(c("left", "right"), function(side) {
      list(apply(BGNPOW[[side]]$POW, 2, function(POW) {
        POW > powthr
      }))
    })

    singSat <- as.data.frame(do.call(cbind, sapply(c("left", "right"), function(side) {
      list(sapply(1:length(BGNPOW$timeBins), function(i) {
        sum(BGNsaturation[[side]][, i] |
              POWsaturation[[side]][, i]) / halfWl
      }))
    })))

  } else if ("mono" %in% names(BGNPOW)) {
    BGNsaturation <- list(mono = apply(BGNPOW$mono$BGN, 2, function(BGN) {
      Q <- quantile(BGN, bgnthr)
      BGN > Q
    }))


    POWsaturation <- list(mono = apply(BGNPOW$mono$POW, 2, function(POW) {
      POW > powthr
    }))

    singSat <- data.frame("mono" = do.call(cbind, list(sapply(1:length(BGNPOW$timeBins), function(i) {
      sum(BGNsaturation$mono[, i] |
            POWsaturation$mono[, i]) / halfWl
    }))))

  } else {
    realChannel <- c("left", "right")[c("left", "right") %in% names(BGNPOW)]

    BGNsaturation <- list(apply(BGNPOW[[realChannel]]$BGN, 2, function(BGN) {
      Q <- quantile(BGN, bgnthr)
      BGN > Q
    })) |>
      setNames(realChannel)

    POWsaturation <- list(apply(BGNPOW[[realChannel]]$POW, 2, function(POW) {
      POW > powthr
    })) |>
      setNames(realChannel)

    singSat <- data.frame(do.call(cbind, list(sapply(1:length(BGNPOW$timeBins), function(i) {
      sum(BGNsaturation[[realChannel]][, i] |
            POWsaturation[[realChannel]][, i]) / halfWl
    })))) |>
      setNames(realChannel)

  }

  singSat

}
