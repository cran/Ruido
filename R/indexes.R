#' @title Spectral Acoustic Complexity Index
#'
#' @description Calculate the Acoustic Complexity Index values of a single audio using the methodology proposed in Pieretti, et al. 2011
#'
#' @param soundfile wav package numeric matrix, tuneR package Wave object or path to a `.wav` file
#' @param channel channel where the metric values will be extracted from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param j size (in seconds) of the cluster interval. Set to `NULL` to use the entire bin as a single cluster. Defaults to `5`
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `hanning(wl)` to use hanning instead. Defaults to `hamming(wl)`
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#'
#' @returns This function returns a [noise.matrix-class] object.
#'
#' @details The Acoustic Complexity Index (`ACI`) quantifies the average proportional change in spectral amplitude between adjacent time steps  across frequency bins. Because biological sounds, particularly bird vocalizations, often exhibit rapid and irregular amplitude fluctuations over time, `ACI` captures this temporal variability as a proxy for acoustic complexity.
#'
#' In `Ruido`, `ACI` is computed independently within each time bin. Within a time bin, the signal is further subdivided into smaller temporal segments, here referred to as cluster intervals \eqn{j}.
#'
#' For a given frequency bin \eqn{f_l}, acoustic intensity values \eqn{I_k} are evaluated across consecutive time steps \eqn{k} within each cluster interval \eqn{j}. The absolute differences between adjacent time steps are calculated as:
#'
#' \deqn{d_k = |I_k - I_{k+1}|}
#'
#' These differences are summed within each cluster interval:
#'
#' \deqn{D_j = \sum_{k = 1}^{N} d_k}
#'
#' where \eqn{N} is the number of time steps \eqn{\Delta t_k} in interval \eqn{j}. The `ACI` for each cluster interval is then:
#'
#' \deqn{ACI_j = \frac{D_j}{\sum_{k = 1}^{N} I_k}}
#'
#' where \eqn{\sum_{k = 1}^{N} I_k} is the total acoustic intensity within the same interval.
#'
#' For each frequency bin \eqn{f_l}, `ACI` values are summed across all cluster intervals within the time bin:
#'
#' \deqn{ACI_{f_l} = \sum_{j = 1}^{m} ACI_j}
#'
#' where \eqn{m} is the number of cluster intervals in the time bin.
#'
#' The result is a frequency-resolved representation of `ACI` for each time bin, rather than a single scalar value for the entire recording.
#'
#' In the original formulation (Pieretti et al., 2011), `ACI` is further summed across all frequency bins:
#'
#' \deqn{ACI_{tot} = \sum_{l = 1}^{q} ACI_{f_l}}
#'
#' where \eqn{q} is the total number of frequency bins. This final aggregation step is not performed in this package.
#'
#' @seealso [ENTspec()] to calculate Spectral Entropy and [bgNoise()] to calculate Background Noise and Soundscape Power.
#'
#' @references
#' Pieretti, N., Farina, A., & Morri, D. (2011). A new methodology to infer the singing activity of an avian community: The Acoustic Complexity Index (ACI). Ecological Indicators, 11(3), 868–873. https://doi.org/10.1016/j.ecolind.2010.11.005
#'
#'@export
#'@importFrom tuneR readWave
#'@importFrom tuneR downsample
#'@importFrom wav read_wav
#'
#' @examples
#' \donttest{
#' ### This is an secondary example using audio from a real soundscape
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
#' # Running the ACIspec function with all the default arguments
#' aci = ACIspec(recDir)
#'
#' # Here's the result
#' aci
#'
#' # Plot ACI values
#' plot(aci)
#'}
#'
ACIspec = function(soundfile,
                   channel = "stereo",
                   timeBin = 60,
                   j = 5,
                   targetSampRate = NULL,
                   wl = 512,
                   window = hamming(wl),
                   overlap = ceiling(length(window) / 2)) {
  argHandler(
    FUN = "ACIspec",
    soundfile = soundfile,
    channel = channel,
    timeBin = timeBin,
    j = j,
    targetSampRate = targetSampRate,
    wl = wl,
    window = window,
    overlap = overlap
  )

  if (!is.null(timeBin)) {
    j = if (is.null(j))
      timeBin
    else
      min(j, timeBin)
  }

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

  ACIexp = processChannel.ACI(
    channelData = soundfile,
    samp.rate = attr(soundfile, "sample_rate"),
    channel = channel,
    timeBin = timeBin,
    j = j,
    wl = wl,
    overlap = overlap,
    window = window,
    noiseOBJ = new("noise.matrix")
  )

  if (ACIexp@channel == "stereo") {
    ACIexp@wl = nrow(ACIexp@values$left$ACI)

  } else {
    ACIexp@wl = nrow(ACIexp@values[[channel]]$ACI)

  }

  return(ACIexp)

}

#' @title Acoustic Activity Matrix
#'
#' @description Calculate the Acoustic Activity Matrix using the methodology proposed in Burivalova 2018
#'
#' @param soundfile wav package numeric matrix, tuneR package Wave object or path to a `.wav` file
#' @param channel channel where the saturation values will be extracted from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`.
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Set to `NULL` to leave db values unrestricted Defaults to `-90`, as set by Towsey 2017
#' @param dbRef reference value used to compute decibel values from the spectrogram. Defaults to `10`
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `hanning(wl)` to use hanning instead. Defaults to `hamming(wl)`
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges"`, `"scott"` and `100`. Defaults to `"FD"`.
#' <br>Can also be set to any numerical value to limit or increase the amount of breaks.
#' @param DCfix if the DC offset should be removed before the metrics are calculated. Defaults to `TRUE`
#' @param powthr single numeric value to calculate the activity matrix for soundscape power (in dB). Defaults to `10`
#' @param bgnthr single numeric value to calculate the activity matrix for background noise (in %). Defaults to `0.8`
#' @param beta how BGN thresholds are calculated. If `TRUE`, BGN thresholds are calculated using all recordings combined. If FALSE, BGN thresholds are calculated separately for each recording. Defaults to `TRUE`
#'
#' @returns This function returns a 0 and 1 matrix containing the activity for all time bins of the inputted file. The matrix's number of rows will equal to half the set window length (`wl`) and number of columns will equal the number of bins. Cells with the value of 1 represent the acoustically active frequency of a bin.
#'
#' @details To calculate the activity matrix, we use the methodology proposed by Burivalova 2018. We begin by applying the following formula to each time bin of the recording:
#'
#' \deqn{a_{m,f} = \begin{cases} 1, & \text{if } BGN_{m,f} > \theta_1 \ \text{ or } POW_{m,f} > \theta_2 \\ 0, & \text{otherwise} \end{cases}}
#'
#' where \eqn{\theta} is a user-defined threshold applied uniformly to both `BGN` and `POW`. We set 1 to active and 0 to inactive frequency windows.
#'
#' @seealso [multActivity()] to run this over multiple audio files and [singleSat()] to get the full saturation values this activity matrix is derived from.
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
#' library(ggplot2)
#' # We are going to load a sample noise.matrix object to demonstrate the basic usage of singleSat()
#' # To understand about the origin of this noise.matrix, check: ?sampleBGN
#' data("sampleBGN")
#'
#' # View the sample noise.matrix object
#' sampleBGN
#'
#' # Run the function
#' act = activity(sampleBGN)
#'
#' # Now we can plot the results for the left channel
#' actLeft = act@values$left$ACT
#' actDim = dim(actLeft)
#' numericTime = seq(0, sum(sampleBGN@timeBins), by = sampleBGN@timeBins[1])
#' labels = paste0(numericTime[-length(numericTime)], "-", numericTime[-1], "s")
#'
#' actDF = data.frame(BIN = rep(paste0("BIN", seq(actDim[2])), each = actDim[1]),
#'                     WIN = rep(seq(actDim[1]), actDim[2]),
#'                     ACT = factor(unlist(actLeft), levels = c(0, 1)))
#'
#' ggplot(actDF, aes(x = BIN, y = WIN, fill = ACT)) +
#'   geom_tile() +
#'   theme_bw() +
#'   scale_fill_manual(values = c("white", "black")) +
#'   scale_y_continuous(expand = c(0,0)) +
#'   scale_x_discrete(expand = c(0,0), labels = labels) +
#'   labs(x = "Time Bin", y = "Spectral Window") +
#'   guides(fill = guide_legend(title = "Activity"))
#'
#' }
activity = function(soundfile,
                    channel = "stereo",
                    timeBin = 60,
                    dbThreshold = -90,
                    dbRef = 10,
                    targetSampRate = NULL,
                    wl = 512,
                    window = hamming(wl),
                    overlap = ceiling(length(window) / 2),
                    histbreaks = "FD",
                    DCfix = TRUE,
                    powthr = 10,
                    bgnthr = 0.8,
                    beta = TRUE) {
  argHandler(
    FUN = "activity",
    soundfile = soundfile,
    channel = channel,
    timeBin = timeBin,
    dbThreshold = dbThreshold,
    dbRef = dbRef,
    targetSampRate = targetSampRate,
    wl = wl,
    window = window,
    overlap = overlap,
    histbreaks = histbreaks,
    DCfix = DCfix,
    powthr = powthr,
    bgnthr = bgnthr,
    beta = beta
  )

  BGNPOW = if (is(soundfile, "noise.matrix")) {
    soundfile
  } else {
    bgNoise.(
      soundfile,
      timeBin = timeBin,
      targetSampRate = targetSampRate,
      window = window,
      overlap = overlap,
      channel = channel,
      dbThreshold = dbThreshold,
      dbRef = dbRef,
      wl = wl,
      histbreaks = histbreaks,
      DCfix = DCfix,
      noiseOBJ = new("noise.matrix")
    )
  }

  nBins = length(BGNPOW@timeBins)

  if (beta) {
    if (BGNPOW@channel == "stereo") {
      BGN = cbind(BGNPOW@values$left$BGN, BGNPOW@values$right$BGN)
    } else {
      BGN = BGNPOW@values[[BGNPOW@channel]]$BGN
    }
  }

  BGNPOW@values = lapply(BGNPOW@values, function(ch) {
    chBGN = as.matrix(ch$BGN)
    chPOW = as.matrix(ch$POW)

    if (beta) {
      BGNQ = quantile(unlist(BGN), bgnthr)
      result = (chBGN > BGNQ | chPOW > powthr) * 1
    } else {
      result = (
        ch$BGN > matrix(
          matrixStats::colQuantiles(chBGN, probs = bgnthr),
          nrow = nrow(chBGN),
          ncol = ncol(chBGN),
          byrow = TRUE
        ) | chPOW > powthr
      ) * 1
    }
    colnames(result) = paste0("ACT", 1:nBins)
    list("ACT" = as.data.frame(result))
  })

  BGNPOW@index = "ACT"
  if (BGNPOW@channel == "stereo") {
    BGNPOW@wl = nrow(BGNPOW@values$left$ACT)
  } else {
    BGNPOW@wl = nrow(BGNPOW@values[[BGNPOW@channel]]$ACT)
  }

  return(BGNPOW)
}

#' @title Summarized Background Noise and Soundscape Power Index
#'
#' @description Calculate the summarized Background Noise and Soundscape Power values of a single audio using the methodology proposed in Towsey 2017
#'
#' @param soundfile wav package numeric matrix, tuneR package Wave object or path to a `.wav` file
#' @param channel channel where the metric values will be extracted from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Set to `NULL` to leave db values unrestricted Defaults to `-90`, as set by Towsey 2017
#' @param dbRef reference value used to compute decibel values from the spectrogram. Defaults to `10`
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl length of each waveform frame in samples. The waveform is divided into non-overlapping frames of this length before the envelope is calculated. Defaults to `512`
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges"`, `"scott"` and `100`. Defaults to `"FD"`.
#' <br>Can also be set to any numerical value to limit or increase the amount of breaks.
#' @param DCfix if the DC offset should be removed before the metrics are calculated. Defaults to `TRUE`
#'
#' @returns A list containing the `BGN` and `POW` values calculated for each time bin and selected channel.
#'
#' @details Background Noise (`BGN`) is an acoustic metric that estimates the dominant background sound level within a time bin from the waveform amplitude envelope. Following the approach described by Towsey (2017), each time bin is divided into non-overlapping waveform frames, and the maximum absolute amplitude of each frame is used to construct the waveform envelope.
#'
#' For each waveform frame \eqn{i}, the envelope amplitude \eqn{A_i} is converted to decibels as:
#'
#' \deqn{dB_i = dbRef . \log_{10}(A_i)}
#'
#' where \eqn{A_i} is the maximum absolute amplitude within the frame. Incomplete frames at the end of a time bin are discarded.
#'
#' The resulting dB values are grouped into histogram bins, and `BGN` is estimated from the modal histogram bin, representing the most frequently occurring sound level within the time bin:
#'
#' \deqn{BGN = \mathrm{mode}(dB_i)}
#'
#' `BGN` therefore represents an estimate of the dominant background level of the recording segment, rather than the background level of a particular frequency band.
#'
#' Soundscape Power (`POW`) quantifies the difference between the maximum observed sound level and the estimated background level within the same time bin. It is defined as:
#'
#' \deqn{POW = \max(dB_i) - BGN}
#'
#' Higher `POW` values indicate a greater difference between the strongest acoustic event and the estimated background level.
#'
#' @seealso [bgNoise()] to calculate Spectral Background Noise and Soundscape Power.
#'
#' @references
#' Towsey, M. W. (2017). The calculation of acoustic indices derived from long-duration recordings of the natural environment. In eprints.qut.edu.au. https://eprints.qut.edu.au/110634/
#' <br>Lamel, L., Rabiner, L., Rosenberg, A., & Wilpon, J. (1981). An improved endpoint detector for isolated word recognition. \emph{IEEE Transactions on Acoustics, Speech, and Signal Processing}, 29(4), 777-785 https://doi.org/10.1109/TASSP.1981.1163642
#'
#'@export
#'@importFrom grDevices nclass.FD
#'@importFrom grDevices nclass.Sturges
#'@importFrom grDevices nclass.scott
#'@importFrom matrixStats colMaxs
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
#' BGN = bgn(wave)
#'
#' # See the results
#' BGN
bgn = function(soundfile,
               channel = "stereo",
               timeBin = 60,
               dbThreshold = -90,
               dbRef = 10,
               targetSampRate = NULL,
               wl = 512,
               histbreaks = "FD",
               DCfix = TRUE) {
  argHandler(
    FUN = "bgn",
    soundfile = soundfile,
    channel = channel,
    timeBin = timeBin,
    wl = wl,
    dbThreshold = dbThreshold,
    dbRef = dbRef,
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
      soundfile = rbind(soundfile@left, soundfile@right)
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

  allSamples = if (is.null(timeBin)) {
    data.frame(b = 1, e = length(soundfile[1, ]))
  } else {
    getSampleBins(length(soundfile[1, ]),
                  attr(soundfile, "sample_rate"),
                  timeBin)
  }

  channelData = switch(
    channel,
    "stereo" = list(left  = soundfile[1, ], right = soundfile[2, ]),
    "mono"   = list(mono  = soundfile[1, ]),
    "left"   = list(left  = soundfile[1, ]),
    "right"  = list(right = soundfile[2, ])
  )

  hBreak = hBreaks(histbreaks)

  lapply(channelData, function(x) {
    if (DCfix) {
      x = x - mean(x)
    }
    x = abs(x)
    apply(allSamples, 1, function(y) {
      samples = x[y[1]:y[2]]
      mat = matrix(samples[seq_len(floor(length(samples) / wl) * wl)], nrow = wl)
      db = 10 * log10(matrixStats::colMaxs(mat))

      if (!is.null(dbThreshold)) {
        db[db < dbThreshold] = dbThreshold
      }

      dbMin = min(db)
      dbMax = max(db)
      numBins = hBreak(db)
      breaks = seq(dbMin, dbMax, length.out = numBins + 1)
      modalBin = which.max(tabulate(findInterval(x = db, vec = breaks)))
      modalIntensity = dbMin + modalBin * (breaks[2] - breaks[1])
      c(BGN = modalIntensity, POW = dbMax - modalIntensity)
    })

  })

}

#' @title Spectral Background Noise and Soundscape Power Index
#'
#' @description Calculate the spectral Background Noise and Soundscape Power values of a single audio using the methodology proposed in Towsey 2017
#'
#' @param soundfile wav package numeric matrix, tuneR package Wave object or path to a `.wav` file
#' @param channel channel where the metric values will be extracted from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Set to `NULL` to leave db values unrestricted Defaults to `-90`, as set by Towsey 2017
#' @param dbRef reference value used to convert spectrogram amplitude values to decibels. Defaults to `10`
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `hanning(wl)` to use hanning instead. Defaults to `hamming(wl)`
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges"`, `"scott"` and `100`. Defaults to `"FD"`.
#' <br>Can also be set to any numerical value to limit or increase the amount of breaks.
#' @param DCfix if the DC offset should be removed before the metrics are calculated. Defaults to `TRUE`
#'
#' @returns This function returns a [noise.matrix-class] object
#'
#' @details Background Noise (`BGN`) is an acoustic metric that estimates the dominant baseline level of acoustic energy within a frequency window and time bin. It was described by Towsey (2017) based on the approach of Lamel et al. (1981).
#'
#' Spectrogram amplitude values are converted to decibels according to:
#'
#' \deqn{dB_i = dbRef . \log_{10}(A_i)}
#'
#' where \eqn{A_i} is the amplitude of the \eqn{i}-th spectrogram value and \eqn{dbRef} is the coefficient used for the conversion. By default, \eqn{dbRef = 10}.
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
#' @seealso [bgn()] to calculate the summarized Background Noise and Soundscape Power. [ACIspec()] to calculate the Acoustic Complexity Index and [ENTspec()] to calculate Spectral Entropy from a single audio file. Also, check [activity()] and [singleSat()], which use this same Background Noise and Soundscape Power calculation to determine acoustic activity and saturation.
#'
#' @references
#' Towsey, M. W. (2017). The calculation of acoustic indices derived from long-duration recordings of the natural environment. In eprints.qut.edu.au. https://eprints.qut.edu.au/110634/
#' <br>Lamel, L., Rabiner, L., Rosenberg, A., & Wilpon, J. (1981). An improved endpoint detector for isolated word recognition. \emph{IEEE Transactions on Acoustics, Speech, and Signal Processing}, 29(4), 777-785 https://doi.org/10.1109/TASSP.1981.1163642
#'
#'@export
#'@importFrom stats mvfft
#'@useDynLib Ruido, .registration = TRUE
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
#' BGN = bgNoise(wave)
#'
#' # See the results
#' BGN
#'
#' # Plot background noise and soundscape power
#' plot(BGN)
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
#' BGN = bgNoise(recDir)
#'
#' # Here's the result
#' BGN
#'
#' # Plot background noise and soundscape power values
#' plot(BGN)
#'
#' # Plot the two indices against each other
#' plot(BGN@values$left$BGN$BGN1, BGN@values$left$POW$POW1,
#'      xlab = "BGN (dB)", ylab = "POW (dB)", pch = 16)
#'
#' # Now lets test and plot their correlation
#' BGNPOWlm = lm(BGN@values$left$BGN$BGN1~BGN@values$left$POW$POW1)
#' summary(BGNPOWlm)
#' abline(lm(BGN@values$left$BGN$BGN1~BGN@values$left$POW$POW1), col = "red")
#'}
bgNoise = function(soundfile,
                   channel = "stereo",
                   timeBin = 60,
                   dbThreshold = -90,
                   dbRef = 10,
                   targetSampRate = NULL,
                   wl = 512,
                   window = hamming(wl),
                   overlap = ceiling(length(window) / 2),
                   histbreaks = "FD",
                   DCfix = TRUE) {
  argHandler(
    FUN = "bgNoise",
    soundfile = soundfile,
    channel = channel,
    timeBin = timeBin,
    dbThreshold = dbThreshold,
    dbRef = dbRef,
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
      soundfile = rbind(soundfile@left, soundfile@right)
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
    channelData = soundfile,
    samp.rate = attr(soundfile, "sample_rate"),
    channel = channel,
    timeBin = timeBin,
    wl = wl,
    overlap = overlap,
    dbThreshold = dbThreshold,
    dbRef = dbRef,
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

#' @title Spectral Temporal Entropy Index
#'
#' @description Calculate the Temporal Entropy values of a single audio using the methodology proposed in Towsey, et al. 2017
#'
#' @param soundfile wav package numeric matrix, tuneR package Wave object or path to a `.wav` file
#' @param channel channel where the metric values will be extracted from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `hanning(wl)` to use hanning instead. Defaults to `hamming(wl)`
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
                   window = hamming(wl),
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
    channelData = soundfile,
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

#' @title Spectral Acoustic Event Index
#'
#' @description Calculate the number of acoustic events per second in each frequency bin of a single audio using the methodology described by Towsey (2017).
#'
#' @param soundfile wav package numeric matrix, tuneR package Wave object or path to a `.wav` file
#' @param channel channel where the metric values will be extracted from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Set to `NULL` to leave db values unrestricted Defaults to `-90`, as set by Towsey 2017
#' @param dbRef reference value used to compute decibel values from the spectrogram. Defaults to `10`
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `hanning(wl)` to use hanning instead. Defaults to `hamming(wl)`
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges"`, `"scott"` and `100`. Defaults to `"FD"`.
#' <br>Can also be set to any numerical value to limit or increase the amount of breaks.
#' @param DCfix if the DC offset should be removed before the metrics are calculated. Defaults to `TRUE`
#'
#' @returns This function returns a [noise.matrix-class] object containing the number of acoustic events per second in each frequency bin.
#'
#' @details
#' Spectral Acoustic Events (EVNsp) quantifies the number of acoustic events per second within each frequency bin. An acoustic event is defined as an upward crossing of a 3-dB threshold above the estimated background noise level.
#'
#' The spectrogram is first converted to decibels and the spectral Background Noise (BGN) is calculated independently for each frequency bin and time bin. The resulting background level is then subtracted from the corresponding frequency bin of the spectrogram to produce a noise-reduced spectrogram:
#'
#' \deqn{D_{f,t} = dB_{f,t} - BGN_f}
#'
#' where \eqn{D_{f,t}} is the noise-reduced intensity for frequency bin \eqn{f} at spectral frame \eqn{t}, \eqn{dB_{f,t}} is the corresponding spectrogram intensity, and \eqn{BGN_f} is the estimated background level for that frequency bin.
#'
#' An acoustic event is counted whenever the noise-reduced intensity crosses the 3-dB threshold from below to at or above the threshold:
#'
#' \deqn{EVN_{f} = \sum_t I(D_{f,t-1} < 3 ;\mathrm{dB} ;\mathrm{and}; D_{f,t} \geq 3 ;\mathrm{dB})}
#'
#' where \eqn{I} is an indicator function. Thus, consecutive spectral frames remaining above the threshold are considered part of the same event and are not counted repeatedly. Values are divided by the bin size in seconds.
#'
#' Event counts are calculated independently for each frequency bin. The resulting values represent the number of detected acoustic events per time bin and are expressed as events per second.
#'
#' @seealso [bgNoise()] to calculate spectral Background Noise and Soundscape Power, [bgn()] to calculate summarized Background Noise and Soundscape Power, and [ACIspec()] to calculate the Acoustic Complexity Index from a single audio file.
#'
#' @references
#' Towsey, M. W. (2017). The calculation of acoustic indices derived from long-duration recordings of the natural environment. In eprints.qut.edu.au. https://eprints.qut.edu.au/110634/
#' <br>Lamel, L., Rabiner, L., Rosenberg, A., & Wilpon, J. (1981). An improved endpoint detector for isolated word recognition. \emph{IEEE Transactions on Acoustics, Speech, and Signal Processing}, 29(4), 777-785 https://doi.org/10.1109/TASSP.1981.1163642
#'
#' @export
#'
#' @examples
#' \donttest{
#' ### These audios are originated from the Escutadô Project, a project
#' ### that records the soundscapes of the brazilian semiarid
#' # Getting audiofile from the online Zenodo library
#' dir = paste(tempdir(), "forExample", sep = "/")
#' dir.create(dir)
#' rec = paste0("GAL24576_20250401_", sprintf("%06d", 0), ".wav")
#' recDir = paste(dir, rec , sep = "/")
#' url = paste0("https://zenodo.org/records/17575795/files/",
#'              rec,
#'              "?download=1")
#'
#' # Downloading the file, might take some time denpending on your internet
#' download.file(url, destfile = recDir, mode = "wb")
#'
#' # Running the EVNspec function with all the default arguments
#' evn = EVNspec(recDir)
#'
#' # Here's the result
#' evn
#'
#' # Plot EVN values
#' plot(evn)
#' }
EVNspec = function(soundfile,
                   channel = "stereo",
                   timeBin = 60,
                   dbThreshold = -90,
                   dbRef = 10,
                   targetSampRate = NULL,
                   wl = 512,
                   window = hamming(wl),
                   overlap = ceiling(length(window) / 2),
                   histbreaks = "FD",
                   DCfix = TRUE) {
  argHandler(
    FUN = "EVNspec",
    soundfile = soundfile,
    channel = channel,
    timeBin = timeBin,
    dbThreshold = dbThreshold,
    dbRef = dbRef,
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
      soundfile = rbind(soundfile@left, soundfile@right)
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

  EVNexp = processChannel.EVN(
    channelData = soundfile,
    samp.rate = attr(soundfile, "sample_rate"),
    channel = channel,
    timeBin = timeBin,
    wl = wl,
    overlap = overlap,
    dbThreshold = dbThreshold,
    dbRef = dbRef,
    window = window,
    histbreaks = histbreaks,
    DCfix = DCfix,
    noiseOBJ = new("noise.matrix")
  )

  if (EVNexp@channel == "stereo") {
    EVNexp@wl = nrow(EVNexp@values$left$EVN)
  } else {
    EVNexp@wl = nrow(EVNexp@values[[channel]]$EVN)
  }

  return(EVNexp)

}

#' @title Single Soundscape Saturation Index
#'
#' @param soundfile wav package numeric matrix, tuneR package Wave object, Ruido noise.matrix object or path to a `.wav` file
#' @param channel channel where the background noise values will be extracted from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`.
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Set to `NULL` to leave db values unrestricted Defaults to `-90`, as set by Towsey 2017
#' @param dbRef reference value used to compute decibel values from the spectrogram.
#' @param targetSampRate sample rate of the audios. Defaults to `NULL` to not change the sample rate. This argument is only used to down sample the audio.
#' @param wl window length of the spectrogram. Defaults to `512`.
#' @param window window used to smooth the spectrogram. Defaults to `hamming(wl)`. Switch to `hanning(wl)` if to use hanning instead.
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges"`, `"scott"` and `100`. Defaults to `"FD"`.
#' <br>Can also be set to any numerical value to limit or increase the amount of breaks.
#' @param DCfix if the DC offset should be removed before the metrics are calculated. Defaults to `TRUE`
#' @param powthr a single value to evaluate the activity matrix for Soundscape Power (in %dB). Defaults to `10`.
#' @param bgnthr a single value to evaluate the activity matrix for Background Noise (in %). Defaults to `0.8`
#' @param beta how BGN thresholds are calculated. If TRUE, BGN thresholds are computed using all recordings combined.
#'
#' @export
#' @returns A list containing the saturation values for all time bins of the inputted file
#' @details Soundscape Saturation (`SAT`) quantifies the proportion of frequency bins that are acoustically active within a given time bin. It wasproposed by Burivalova et al. (2018) as a metric to evaluate the acoustic niche hypothesis.
#'
#' For each time bin \eqn{m}, an activity matrix \eqn{a_{m,f}} is first constructed across frequency bins \eqn{f}. A frequency bin is considered active if either its background level (`BGN`) or its soundscape power (`POW`) exceeds a defined threshold:
#'
#' \deqn{a_{m,f} = \begin{cases} 1, & \text{if } BGN_{m,f} > \theta_1 \ \text{ or } POW_{m,f} > \theta_2 \\ 0, & \text{otherwise} \end{cases}}
#'
#' where \eqn{\theta} is a user-defined threshold applied uniformly to both `BGN` and `POW`.
#'
#' Soundscape saturation for time bin \eqn{m} is then calculated as the proportion of active frequency bins:
#'
#' \deqn{S_m = \frac{\sum_{f = 1}^{N} a_{m,f}}{N}}
#'
#' where \eqn{N} is the total number of frequency bins. Higher values of `SAT` indicate a greater fraction of the frequency spectrum being occupied by acoustic activity.
#'
#' @seealso [soundSat()] and [soundMat()] to work with multiple audio files and [activity()] to get only the activity matrix.
#'
#'@references Burivalova, Z., Towsey, M., Boucher, T., Truskinger, A., Apelis, C., Roe, P., & Game, E. T. (2018). Using soundscapes to detect variable degrees of human influence on tropical forests in Papua New Guinea. Conservation Biology, 32(1), 205-215. https://doi.org/10.1111/cobi.12968
#'
#' @examples
#' # First example: Using a Ruido noise.matrix object
#' # We are going to load a sample noise.matrix object to demonstrate the basic usage of singleSat()
#' # To understand about the origin of this noise.matrix, check: ?sampleBGN
#' data("sampleBGN")
#'
#' # View the sample noise.matrix object
#' sampleBGN
#'
#' # Run the function
#' SAT = singleSat(sampleBGN)
#'
#' # View the results
#' SAT
#'
#' # Now lets plot our results to see the dynamics of soundscape saturation by minute
#' maxV = max(unlist(SAT))
#' minV = min(unlist(SAT))
#'
#' plot(x = c(1, 3), y = c(minV, maxV), type = "n",
#'      xlab = "Minute", ylab = "Soundscape Saturation (%)", xaxt = "n")
#' lines(x = 1:3, SAT$left, col = "#1ECBE1", type = "b", pch = 16)
#' axis(1, at = 1:3)
#' lines(x = 1:3, SAT$right, col = "#E1341E", type = "b", pch = 16)
#' legend("topright", legend = c("Left", "Right"), col = c("#1ECBE1", "#E1341E"), lty = 1, pch = 16)
#'
#' # Second example: Using a tuneR Wave-class object
#' # Lets produce an artificial audio with the tuneR package to demonstrate that
#' # the function can also read Wave-class objects (This is the same object used in
#' # the example of bgNoise!)
#'
#' library(tuneR)
#'
#' oldpar = par(no.readonly = TRUE)
#'
#' # Define parameters for the artificial audio
#' samprate = 12050
#' dur = 59
#' n = samprate * dur
#'
#' # White noise
#' set.seed(413)
#' noise = rnorm(n)
#'
#' # Linear fade-out envelope
#' fade = seq(1, 0, length.out = n)
#'
#' # Apply fade
#' signal = noise * fade
#'
#' # Create Wave object
#' wave = Wave(
#'   left = signal,
#'   samp.rate = samprate,
#'   bit = 16
#' )
#'
#' # Running singleSat() on the artificial audio
#' sat = singleSat(wave, timeBin = 10)
#'
#' # Now we can plot the results
#' # In the left we have a periodogram and in the right saturation values
#' # along one minute
#' par(mfrow = c(1,2))
#' image(periodogram(wave, width = 8192, normalize = FALSE), xlab = "Time (s)",
#' ylab = "Frequency (hz)", axes = FALSE)
#' axis(1, labels = seq(0,60, 10), at = seq(0,7e5,length.out = 7))
#' axis(2)
#' plot(sat$mono, xlab = "Time (s)", ylab = "Soundscape Saturation (%)",
#' type = "b", pch = 16, axes = FALSE)
#' axis(1, labels = paste0(c("0-10","10-20","20-30","30-40","40-50","50-59"),
#' "s"), at = 1:6)
#' axis(2)
#'
#' par(oldpar)
#'
#' # Third example: Calculating activity beforehand
#' ## You can use a noise.matrix with activity values to calculate saturation.
#'
#' data("sampleBGN")
#' sampleBGN
#' singleSat(activity(sampleBGN))
#'
#' \donttest{
#' # Fourth example: Reading a file directly
#' # Lets begin by loading an audio from the online Zenodo library and
#' # read it directly with the function
#' # Getting audiofile from the online Zenodo library
#' dir = paste(tempdir(), "forExample", sep = "/")
#' dir.create(dir)
#' rec = paste0("GAL24576_20250401_", sprintf("%06d", 0),".wav")
#' recDir = paste(dir,rec , sep = "/")
#' url = paste0("https://zenodo.org/records/17575795/files/", rec, "?download=1")
#'
#' # Downloading the file, might take some time denpending on your internet
#' download.file(url, destfile = recDir, mode = "wb")
#'
#' # Now we calculate soundscape saturation for both sides of the recording
#' sat = singleSat(recDir)
#'
#' # Printing the results
#' print(sat)
#'
#' barplot(unlist(sat), col = c("darkgreen", "red"),
#'        names.arg = c("Left", "Right"), ylab = "Soundscape Saturation (%)")
#'
#' unlink(dir, recursive = TRUE)
#' }
singleSat = function(soundfile,
                     channel = "stereo",
                     timeBin = 60,
                     dbThreshold = -90,
                     dbRef = 10,
                     targetSampRate = NULL,
                     wl = 512,
                     window = hamming(wl),
                     overlap = ceiling(length(window) / 2),
                     histbreaks = "FD",
                     DCfix = TRUE,
                     powthr = 10,
                     bgnthr = 0.8,
                     beta = TRUE) {
  argHandler(
    FUN = "singleSat",
    soundfile = soundfile,
    channel = channel,
    timeBin = timeBin,
    dbThreshold = dbThreshold,
    dbRef = dbRef,
    targetSampRate = targetSampRate,
    wl = wl,
    window = window,
    overlap = overlap,
    histbreaks = histbreaks,
    DCfix = DCfix,
    powthr = powthr,
    bgnthr = bgnthr,
    beta = beta
  )

  halfWl = as.integer(wl / 2)

  BGNPOW = if (is(soundfile, "noise.matrix")) {
    soundfile
  } else {
    bgNoise.(
      soundfile,
      timeBin = timeBin,
      targetSampRate = targetSampRate,
      window = window,
      overlap = overlap,
      channel = channel,
      dbThreshold = dbThreshold,
      dbRef = dbRef,
      wl = wl,
      histbreaks = histbreaks,
      DCfix = DCfix
    )
  }

  nBins = length(BGNPOW@timeBins)

  if (length(BGNPOW@index) == 1) {
    if (BGNPOW@channel == "stereo") {
      names = paste0(rep(c("left", "right"), each = nBins), seq(nBins))
      singSat = c(matrixStats::colMeans2(as.matrix(BGNPOW@values$left$ACT)),
                  matrixStats::colMeans2(as.matrix(BGNPOW@values$right$ACT)))
    } else {
      names = paste0(rep(BGNPOW@channel, nBins), seq(nBins))
      singSat = matrixStats::colMeans2(as.matrix(BGNPOW@values[[BGNPOW@channel]]$ACT))
    }

  } else {
    if (BGNPOW@channel == "stereo") {
      BGN = as.matrix(cbind(BGNPOW@values$left$BGN, BGNPOW@values$right$BGN))
      POW = as.matrix(cbind(BGNPOW@values$left$POW, BGNPOW@values$right$POW))
      names = paste0(rep(c("left", "right"), each = nBins), seq(nBins))
    } else {
      BGN = as.matrix(BGNPOW@values[[BGNPOW@channel]]$BGN)
      POW = as.matrix(BGNPOW@values[[BGNPOW@channel]]$POW)
      names = paste0(rep(BGNPOW@channel, nBins), seq(nBins))
    }

    if (beta) {
      BGNQ = quantile(unlist(BGN), bgnthr)

      singSat = matrixStats::colMeans2(BGN > BGNQ | POW > powthr)

    } else {
      singSat = matrixStats::colSums2((
        BGN > matrix(
          matrixStats::colQuantiles(BGN, probs = bgnthr),
          nrow = nrow(BGN),
          ncol = ncol(BGN),
          byrow = TRUE
        )
      ) | POW > powthr) / halfWl

    }
  }

  names(singSat) = names

  if (BGNPOW@channel == "stereo") {
    return(list(left  = singSat[seq(nBins)], right = singSat[seq(nBins + 1, nBins * 2)]))
  } else {
    return(setNames(list(singSat), BGNPOW@channel))
  }

}
