#' @title Spectral Temporal Entropy Index
#'
#' @description Calculate the Temporal Entropy values of a single audio using the methodology proposed in Towsey, et al. 2017
#'
#' @param soundfile wav package numeric matrix, tuneR package Wave object or path to a `.wav` file
#' @param channel channel where the metric values will be extracted from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `signal::hanning(wl)` to use hanning instead. Defaults to `signal::hamming(wl)`
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#'
#' @returns This function returns a [noise.matrix-class] object.
#'
#' @details The Temporal Entropy (`ENT`) quantifies how concentrated or dispersed acoustic energy is over time within each frequency bin. Unlike indices that track frame-to-frame change (e.g. `ACI`), `ENT` treats the distribution of energy across an entire time bin as a probability mass function and measures its Shannon entropy, capturing whether energy is spread evenly through time (high entropy, low concentration) or concentrated into brief pulses (low entropy, high concentration).
#'
#' In `Ruido`, `ENT` is computed independently within each time bin, using every time step of the spectrogram in that bin.
#'
#' For a given frequency bin \eqn{f}, the squared amplitude values \eqn{I_t} across all time steps \eqn{t} within the time bin are normalized to unit area, producing a probability mass function:
#'
#' \deqn{pmf_{f,t} = \frac{I_t^2}{\sum_{t = 1}^{N} I_t^2}}
#'
#' where \eqn{N} is the number of time steps in the time bin. The Shannon entropy of this distribution is then calculated as:
#'
#' \deqn{H[f] = \frac{-\sum_{t = 1}^{N} pmf_{f,t} \times \log_2(pmf_{f,t})}{\log_2 N}}
#'
#' To express the result as an intuitive measure of energy concentration rather than dispersion, `ENT` is calculated as the complement of \eqn{H}:
#'
#' \deqn{ENT[f] = 1 - H[f]}
#'
#' The result is a frequency-resolved representation of `ENT` for each time bin, rather than a single scalar value for the entire recording. Values close to `1` indicate energy concentrated in few time steps (e.g. transient calls or pulses), while values close to `0` indicate energy spread evenly across the time bin (e.g. steady background noise).
#'
#' @seealso [ACIspec()] to calculate the Acoustic Complexity Index and [bgNoise()] to calculate Background Noise and Soundscape Power.
#'
#' @references
#' Towsey, M. W. (2017). The calculation of acoustic indices derived from long-duration recordings of the natural environment. In eprints.qut.edu.au. https://eprints.qut.edu.au/110634/
#'
#'@export
#'@importFrom signal specgram
#'@importFrom tuneR readWave
#'@importFrom tuneR downsample
#'@importFrom wav read_wav
#'
#' @examples
#' \donttest{
#' ### This is an example using audio from a real soundscape
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
#' # Running the ENTspec function with all the default arguments
#' ent = ENTspec(recDir)
#'
#' # Here's the result
#' ent
#'
#' # Plot ENT values
#' plot(ent)
#'}
#'
ENTspec = function(soundfile,
                   channel = "stereo",
                   timeBin = 60,
                   targetSampRate = NULL,
                   wl = 512,
                   window = signal::hamming(wl),
                   overlap = ceiling(length(window) / 2)) {
  argHandler(
    FUN = "ENTspec",
    soundfile = soundfile,
    channel = channel,
    timeBin = timeBin,
    targetSampRate = targetSampRate,
    wl = wl,
    window = window,
    overlap = overlap
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

  ENTexp = processChannel.ENT(
    soundfile,
    samp.rate = attr(soundfile, "sample_rate"),
    channel = channel,
    timeBin = timeBin,
    wl = wl,
    overlap = overlap,
    window = window,
    noiseOBJ = new("noise.matrix")
  )

  if (ENTexp@channel == "stereo") {
    ENTexp@wl = nrow(ENTexp@values$left$ENT)

  } else {
    ENTexp@wl = nrow(ENTexp@values[[channel]]$ENT)

  }

  return(ENTexp)

}
