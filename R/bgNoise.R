#' @title Background Noise and Soundscape Power Index
#'
#' @description Calculate the Background Noise and Soundscape Power values of a single audio using the methodology proposed in Towsey 2017
#'
#' @param soundfile tuneR Wave object or path to a valid audio file
#' @param channel channel where the metric values will be extract from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Defaults to `-90`, as set by Towsey 2017
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `signal::hanning(wl)` to use hanning instead. Defaults to `signal::hammning(wl)`
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges"`, `"scott"` or any numeric value (foe example = `100`). Defaults to `"FD"`
#'
#' @returns This function returns a list containing five objects. The first object (values) contain the values of BGN and POW. The second object (timeBins) contains the durations of each time bin analysed. The third object (sampRate) contains the audio's sampling rate. The fourth and last object (channel) contains the channels used for the calculation of the metric.
#'
#' @details Background Noise (`BGN`) is an acoustic metric that measures the most common continuous baseline level of acoustic energy in a frequency window and in a time bin. It was developed by Towsey 2017 using the Lamel et al 1981 algorithm.
#' The metric is calculated by taking the modal value of intensity values in temporal bin c in frequency window f of a reconding:
#'
#'\deqn{BGN_{f} = mode(dB_{cf})}
#'
#'Soundscape Power represents a measure of signal-to-noise ratio. It measures the relation of BGN to the loudest intensities in temporal bin c in frequency window f:
#'
#'\deqn{POW_{f} = max(dB_{cf}) - BGN_{cf}}
#'
#' This mean we'll have a value of BGN and POW to each frequency window of a recording.
#'
#' @references
#' Towsey, M. W. (2017). The calculation of acoustic indices derived from long-duration recordings of the natural environment. In eprints.qut.edu.au. https://eprints.qut.edu.au/110634/
#' <br>Lamel, L., Rabiner, L., Rosenberg, A., & Wilpon, J. (1981). An improved endpoint detector for isolated word recognition. IEEE Transactions on Acoustics, Speech, and Signal Processing, 29(4), 777-785 https://doi.org/10.1109/TASSP.1981.1163642
#'
#'@export
#'@importFrom signal specgram
#'@importFrom tuneR readWave
#'@importFrom tuneR readMP3
#'@importFrom tuneR downsample
#'
#' @examples
#' ### For our main example we'll create an artificial audio with
#' ### white noise to test its Background Noise
#' # We'll use the package tuneR
#' library(tuneR)
#'
#' # Define the audio sample rate, duration and number of samples
#' sampRate <- 12050
#' dur <- 59
#' sampN <- sampRate * dur
#'
#' # Then we Ggenerate the white noise for our audio and apply FFT
#' set.seed(413)
#' ruido <- rnorm(sampN)
#' spec <- fft(ruido)
#'
#' # Now we create a random spectral envelope for the audio and apply the spectral envelope
#' nPoints <- 50
#' env <- runif(nPoints)
#' env <- approx(env, n=nPoints)$y
#' specMod <- spec * env
#'
#' # Now we invert the FFT back to into a time waveform and normalize and convert to Wave
#' out <- Re(fft(specMod, inverse=TRUE)) / sampN
#' wave <- Wave(left = out, samp.rate = sampRate, bit = 16)
#' wave <- normalize(wave, unit = "16")
#'
#' # Here's our artificial audio
#' wave
#'
#' # Running the bgNoise function with all the default arguments
#' bgn <- bgNoise(wave)
#'
#' # Print the results
#' head(bgn$values$mono$BGN)
#' head(bgn$values$mono$POW)
#'
#' # Plotting background noise and soundscape profile for the first minute of the recording
#' par(mfrow = c(2,2))
#' plot(x = bgn$values$mono$BGN$BGN1, y = seq(1,bgn$sampRate, length.out = 256),
#'      xlab = "Background Noise (dB)", ylab = "Frequency (hz)",
#'      main = "BGN by Frequency",
#'      type = "l")
#' plot(x = bgn$values$mono$POW$POW1, y = seq(1,bgn$sampRate, length.out = 256),
#'      xlab = "Soundscape Power (dB)", ylab = "Frequency (hz)",
#'      main = "POW by Frequency",
#'      type = "l")
#' plot(bgn$values$mono$BGN$BGN1~bgn$values$mono$POW$POW1, pch = 16,
#'      xlab = "Soundscape Power (dB)", ylab = "Background Noise (dB)",
#'      main = "BGN~POW")
#' hist(bgn$values$mono$BGN$BGN1, main = "Histogram of BGN distribution",
#'       xlab = "Background Noise (BGN)")
#'
#'\donttest{
#' oldpar <- par(no.readonly = TRUE)
#' ### This is a secondary example using audio from a real soundscape
#' ### These audios are originated from the Escutadô Project
#' # Getting audiofile from the online Zenodo library
#' dir <- paste(tempdir(), "forExample", sep = "/")
#' dir.create(dir)
#' rec <- paste0("GAL24576_20250401_", sprintf("%06d", 0), ".wav")
#' recDir <- paste(dir, rec , sep = "/")
#' url <- paste0("https://zenodo.org/records/17575795/files/",
#'               rec,
#'               "?download=1")
#'
#' # Downloading the file, might take some time denpending on your internet
#' download.file(url, destfile = recDir, mode = "wb")
#'
#' # Running the bgNoise function with all the default arguments
#' bgn <- bgNoise(recDir)
#'
#' # Print the results
#' head(bgn$values$left$BGN)
#' head(bgn$values$left$POW)
#'
#' # Plotting background noise and soundscape profile for the first minute of the recording
#' par(mfrow = c(2, 2))
#' plot(
#'   x = bgn$values$left$BGN$BGN1,
#'   y = seq(1, bgn$sampRate, length.out = 256),
#'   xlab = "Background Noise (dB)",
#'   ylab = "Frequency (hz)",
#'   main = "BGN by Frequency",
#'   type = "l"
#' )
#' plot(
#'   x = bgn$values$left$POW$POW1,
#'   y = seq(1, bgn$sampRate, length.out = 256),
#'   xlab = "Soundscape Power (dB)",
#'   ylab = "Frequency (hz)",
#'   main = "POW by Frequency",
#'   type = "l"
#' )
#' plot(
#'   bgn$values$left$BGN$BGN1 ~ bgn$values$left$POW$POW1,
#'   pch = 16,
#'   xlab = "Soundscape Power (dB)",
#'   ylab = "Background Noise (dB)",
#'   main = "BGN~POW in left ear"
#' )
#' plot(
#'   bgn$values$right$BGN$BGN1 ~ bgn$values$right$POW$POW1,
#'   pch = 16,
#'   xlab = "Soundscape Power (dB)",
#'   ylab = "Background Noise (dB)",
#'   main = "BGN~POW in right ear"
#' )
#'
#' unlink(dir, recursive = TRUE)
#' par(oldpar)
#'}
bgNoise <- function(soundfile,
                    channel = "stereo",
                    timeBin = 60,
                    dbThreshold = -90,
                    targetSampRate = NULL,
                    wl = 512,
                    window = signal::hamming(wl),
                    overlap = ceiling(length(window) / 2),
                    histbreaks = "FD") {
  if (!channel %in% c("left", "right", "stereo", "mono")) {
    stop("Please provide a valid channel: 'left', 'right', 'stereo', or 'mono'.")
  }

  if(!(histbreaks %in% c("FD", "Sturges", "scott"))) {
    if(!is.numeric(histbreaks)) {
      stop("histbreaks must be 'FD', 'Sturges', 'scott' or a numeric value")
    }
  }

  audio <- if (is.character(soundfile)) {
    fileExt <- tolower(tools::file_ext(soundfile))
    if (fileExt %in% c("mp3", "wav")) {
      if (fileExt == "mp3") {
        tuneR::readMP3(soundfile)
      } else {
        tuneR::readWave(soundfile)
      }
    } else {
      stop("The audio file must be in MP3 or WAV format.")
    }
  } else {
    soundfile
  }

  if(channel == "mono" && audio@stereo) {
    audio <- tuneR::mono(audio, which = "both")
  }

  if (channel == "stereo" && !audio@stereo) {
    message("Audio is not stereo, defaulting to left channel.")
    channel <- "mono"
  }

  if (!is.null(targetSampRate)) {
    audio <- tuneR::downsample(audio, targetSampRate)
  }

  return(
    processChannel(
      audio,
      channel = channel,
      timeBin = timeBin,
      wl = wl,
      overlap = overlap,
      dbThreshold = dbThreshold,
      window = window,
      histbreaks = histbreaks
    )
  )

}
