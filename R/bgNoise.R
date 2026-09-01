#' @title Background Noise and Soundscape Power Index
#'
#' @description Calculate the Background Noise and Soundscape Power values of a single audio using the methodology proposed in Towsey 2017
#'
#' @param soundfile wav package numeric matrix, tuneR package Wave object or path to a `.wav` file
#' @param channel channel where the metric values will be extracted from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Set to `NULL` to leave db values unrestricted Defaults to `-90`, as set by Towsey 2017
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `signal::hanning(wl)` to use hanning instead. Defaults to `signal::hamming(wl)`
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges`", `"scott"` and `100`. Defaults to `"FD"`.
#' <br>Can also be set to any numerical value to limit or increase the amount of breaks.
#' @param DCfix if the DC offset should be removed before the metrics are calculated. Defaults to `TRUE`
#'
#' @returns This function returns a [noise.matrix-class] object
#'
#' @details Background Noise (`BGN`) is an acoustic metric that estimates the dominant baseline level of acoustic energy within a frequency window and time bin. It was described by Towsey (2017) based on the approach of Lamel et al. (1981).
#'
#' For each frequency window \eqn{f} and time bin \eqn{c}, `BGN` is defined as the modal value of the intensity distribution (in dB), representing the most frequently occurring sound level:
#'
#' \deqn{BGN_f = \mathrm{mode}(dB_{c,f})}
#'
#' This value approximates the continuous background component of the soundscape, filtering out transient acoustic events such as bird calls or other short-duration signals.
#'
#' Soundscape Power (`POW`) quantifies the contrast between this baseline level and the strongest acoustic events within the same frequency window and time bin. It is defined as:
#'
#' \deqn{POW_f = \max(dB_{c,f}) - BGN_f}
#'
#' where \eqn{\max(dB_{c,f})} is the maximum intensity observed. `POW` can be interpreted as a proxy for signal-to-noise ratio, with higher values indicating stronger or more prominent acoustic events relative to the background level.
#'
#' @seealso [ACIspec()] to calculate the Acoustic Complexity Index and [ENTspec()] to calculate Spectral Entropy from a single audio file. Also, check [activity()] and [singleSat()], which use this same Background Noise and Soundscape Power calculation to determine acoustic activity and saturation.
#'
#' @references
#' Towsey, M. W. (2017). The calculation of acoustic indices derived from long-duration recordings of the natural environment. In eprints.qut.edu.au. https://eprints.qut.edu.au/110634/
#' <br>Lamel, L., Rabiner, L., Rosenberg, A., & Wilpon, J. (1981). An improved endpoint detector for isolated word recognition. IEEE Transactions on Acoustics, Speech, and Signal Processing, 29(4), 777-785 https://doi.org/10.1109/TASSP.1981.1163642
#'
#'@export
#'@importFrom signal specgram
#'@importFrom tuneR readWave
#'@importFrom tuneR downsample
#'@importFrom wav read_wav
#'@importFrom grDevices nclass.FD
#'@importFrom grDevices nclass.Sturges
#'@importFrom grDevices nclass.scott
#'
#' @examples
#' ### For our main example we'll create an artificial audio with
#' ### white noise to test its Background Noise
#' # We'll use the package tuneR
#' library(tuneR)
#'
#' # Define the audio sample rate, duration and number of samples
#' samprate = 12050
#' dur = 60
#' n = samprate * dur
#'
#' # Then we generate white noise
#' set.seed(413)
#' noise = rnorm(n)
#'
#' # Linear fade-out envelope
#' fade = seq(1, 0, length.out = n)
#'
#' # Apply fade
#' signal = noise * fade
#'
#' wave = Wave(left = signal, right = signal,
#'             samp.rate = samprate,
#'             bit = 16)
#'
#' # Heres our artificial audio
#'
#' wave
#'
#' # Running the bgNoise function with all the default arguments
#' bgn = bgNoise(wave)
#'
#' # See the results
#' bgn
#'
#' # Plot background noise and soundscape power
#' plot(bgn)
#'
#'\donttest{
#' ### This is a secondary example using audio from a real soundscape
#' ### These audios are originated from the Escutadô Project, a project
#' ### that records the soundscapes of the brazilian semiarid
#' # Getting audiofile from the online Zenodo library
#' dir = paste(tempdir(), "forExample", sep = "/")
#' dir.create(dir)
#' rec = paste0("GAL24576_20250401_", sprintf("%06d", 0), ".wav")
#' recDir = paste(dir, rec , sep = "/")
#' url = paste0("https://zenodo.org/records/17575795/files/",
#'               rec,
#'               "?download=1")
#'
#' # Downloading the file, might take some time denpending on your internet
#' download.file(url, destfile = recDir, mode = "wb")
#'
#' # Running the bgNoise function with all the default arguments
#' bgn = bgNoise(recDir)
#'
#' # Here's the result
#' bgn
#'
#' # Plot background noise and soundscape power values
#' plot(bgn)
#'
#' # Plot the two indices against each other
#' plot(bgn@values$left$BGN$BGN1, bgn@values$left$POW$POW1,
#'      xlab = "BGN (dB)", ylab = "POW (dB)", pch = 16)
#'
#' # Now lets test and plot their correlation
#' BGNPOWlm = lm(bgn@values$left$BGN$BGN1~bgn@values$left$POW$POW1)
#' summary(BGNPOWlm)
#' abline(lm(bgn@values$left$BGN$BGN1~bgn@values$left$POW$POW1), col = "red")
#'}
bgNoise = function(soundfile,
                    channel = "stereo",
                    timeBin = 60,
                    dbThreshold = -90,
                    targetSampRate = NULL,
                    wl = 512,
                    window = signal::hamming(wl),
                    overlap = ceiling(length(window) / 2),
                    histbreaks = "FD",
                    DCfix = TRUE) {
  argHandler(
    FUN = "bgNoise",
    soundfile = soundfile,
    channel = channel,
    timeBin = timeBin,
    dbThreshold= dbThreshold,
    targetSampRate = targetSampRate,
    wl = wl,
    window = window,
    overlap = overlap,
    histbreaks = histbreaks,
    DCfix = DCfix
  )

  audio = typeof(soundfile)

  if (audio == "character") {
    if (tolower(tools::file_ext(soundfile)) == "wav") {
      soundfile = wav::read_wav(soundfile)
    } else {
      stop("The audio file must be in the WAV format.")
    }
  } else if (audio == "S4") {
    tempSamp = soundfile@samp.rate
    if (soundfile@stereo) {
      soundfile = matrix(c(soundfile@left, soundfile@right),
                         nrow = 2,
                         byrow = TRUE)
    } else {
      soundfile = matrix(soundfile@left, nrow = 1, byrow = TRUE)
    }
    attr(soundfile, "sample_rate") = tempSamp
  }

  savedAttr = attributes(soundfile)

  if (channel == "mono" && savedAttr$dim[1] > 1) {
    sampRate = attr(soundfile, "sample_rate")
    soundfile = matrix((soundfile[1, ] + soundfile[2, ]) / 2, nrow = 1)
    attr(soundfile, "sample_rate") = sampRate
    savedAttr$dim = dim(soundfile)
  }
  if (channel == "stereo" && nrow(soundfile) == 1) {
    message("Audio is not stereo, defaulting to left channel.")
    channel = "mono"
  }
  if (!is.null(targetSampRate)) {
    audioLen = length(soundfile[1, ])
    keepIdx  = round(seq(1, audioLen, length = targetSampRate * audioLen / savedAttr$sample_rate))
    soundfile = soundfile[1:savedAttr$dim[1], keepIdx, drop = FALSE]
    attr(soundfile, "sample_rate") = targetSampRate
  }

  BGNexp = processChannel.BGN(
    soundfile,
    samp.rate = attr(soundfile, "sample_rate"),
    channel = channel,
    timeBin = timeBin,
    wl = wl,
    overlap = overlap,
    dbThreshold = dbThreshold,
    window = window,
    histbreaks = histbreaks,
    DCfix = DCfix,
    noiseOBJ = new("noise.matrix")
  )

  if (BGNexp@channel == "stereo") {
    BGNexp@wl = nrow(BGNexp@values$left$BGN)

  } else {
    BGNexp@wl = nrow(BGNexp@values[[channel]]$BGN)

  }

  return(BGNexp)

}
