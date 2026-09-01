# processChannel.BGN ----------------------------------------------------------
## This version of processChannel handles the calculation of the BGN and POW index.
## This is called by bgNoise and bgNoise.
## This should not be seen or used by the user, as it is only intended for internal use by other functions.

processChannel.BGN = function(channelData,
                               samp.rate,
                               channel,
                               timeBin,
                               wl,
                               overlap,
                               dbThreshold,
                               window,
                               histbreaks,
                               DCfix,
                               noiseOBJ) {
  allSamples = if (is.null(timeBin)) {
    data.frame(b = 1, e = length(channelData[1, ]))
  } else {
    getSampleBins(length(channelData[1, ]), samp.rate, timeBin)
  }

  frameBin = nrow(allSamples)

  channelData = switch(
    channel,
    "stereo" = list(left  = channelData[1, ], right = channelData[2, ]),
    "mono"   = list(mono  = channelData[1, ]),
    "left"   = list(left  = channelData[1, ]),
    "right"  = list(right = channelData[2, ])
  )

  noiseOBJ@values = lapply(channelData, function(x) {
    if (DCfix) {
      x = x - mean(x)
    }

    tempHolder = apply(allSamples, 1, function(y) {
      list(signal::specgram(
        x = x[y[1]:y[2]],
        n = wl,
        Fs = samp.rate,
        window = window,
        overlap = overlap
      )$S)
    })

    BGNPOWdf = data.frame(do.call(cbind, lapply(lapply(tempHolder, function(singleBin) {
      spectS = abs(singleBin[[1]])

      spectS = 10 * log10(spectS / max(spectS))

      if (!is.null(dbThreshold)) {
        spectS[spectS < dbThreshold] = dbThreshold
      }

      apply(spectS, 1, function(z) {
        dbMax = max(z)
        dbMin = min(z)

        num_bins = if (is.numeric(histbreaks)) {
          histbreaks
        } else {
          switch(
            histbreaks,
            "FD"      = nclass.FD(z),
            "Sturges" = nclass.Sturges(z),
            "scott"   = nclass.scott(z)
          )
        }
        breaks   = seq(dbMin, dbMax, length.out = num_bins + 1)
        modalBin = which.max(tabulate(findInterval(x = z, vec = breaks)))
        modal_intensity = dbMin + modalBin * (breaks[2] - breaks[1])

        c(BGN = modal_intensity, POW = dbMax - modal_intensity)
      })

    }), function(x)
      data.frame(t(x)))))

    colnames(BGNPOWdf) = paste0(rep(c("BGN", "POW"), frameBin), rep(1:frameBin, each = 2))

    BGN = data.frame(BGNPOWdf[, grepl("BGN", colnames(BGNPOWdf)), drop = FALSE])
    POW = data.frame(BGNPOWdf[, grepl("POW", colnames(BGNPOWdf)), drop = FALSE])

    return(list(BGN = BGN, POW = POW))

  })

  noiseOBJ@index = c("BGN", "POW")
  noiseOBJ@timeBins = setNames(round((allSamples$e - allSamples$b) / samp.rate), paste0("BIN", seq(frameBin)))
  noiseOBJ@sampRate = samp.rate
  noiseOBJ@channel = channel

  return(noiseOBJ)

}

# processChannel.ACI ----------------------------------------------------------
## This version of processChannel handles the calculation of the ACI index.
## This is called by ACIspec()
## This should not be seen or used by the user, as it is only intended for internal use by other functions.

processChannel.ACI = function(channelData,
                               samp.rate,
                               channel,
                               timeBin,
                               j,
                               wl,
                               overlap,
                               window,
                               noiseOBJ) {
  allSamples = if (is.null(timeBin)) {
    data.frame(b = 1, e = length(channelData[1, ]))
  } else {
    getSampleBins(length(channelData[1, ]), samp.rate, timeBin)
  }

  frameBin = nrow(allSamples)

  channelData = switch(
    channel,
    "stereo" = list(left  = channelData[1, ], right = channelData[2, ]),
    "mono"   = list(mono  = channelData[1, ]),
    "left"   = list(left  = channelData[1, ]),
    "right"  = list(right = channelData[2, ])
  )

  noiseOBJ@values = lapply(channelData, function(x) {
    tempHolder = apply(allSamples, 1, function(y) {
      list(signal::specgram(
        x = x[y[1]:y[2]],
        n = wl,
        Fs = samp.rate,
        window = window,
        overlap = overlap
      )$S)
    })

    ACIdf = data.frame(do.call(cbind, lapply(tempHolder, function(singleBin) {
      spectS = abs(singleBin[[1]])

      specDim = dim(spectS)
      duration = length(spectS) / samp.rate

      jump = as.integer(j / (duration / specDim[2]))
      forJ = ceiling(specDim[2] / jump)

      ACIvals = vector("numeric", forJ)
      ACIvect = vector("numeric", specDim[1])

      for (s in 1:specDim[1]) {
        for (t in 1:forJ) {
          timeMin = (t - 1) * jump + 1
          timeMax = min(t * jump, specDim[2])

          D = sum(abs(diff(spectS[s, timeMin:(timeMax)])))

          ACIvals[t] = D / sum(spectS[s, timeMin:timeMax])
        }

        ACIvect[s] = sum(ACIvals)

      }

      return(ACIvect)

    })))

    colnames(ACIdf) = paste0("ACI", rep(1:frameBin))
    return(list(ACI = ACIdf))

  })

  noiseOBJ@index = "ACI"
  noiseOBJ@timeBins = setNames(round((allSamples$e - allSamples$b) / samp.rate), paste0("BIN", seq(frameBin)))
  noiseOBJ@sampRate = samp.rate
  noiseOBJ@channel = channel

  return(noiseOBJ)

}


# processChannel.ENT ------------------------------------------------------
## This version of processChannel handles the calculation of the ENT index.
## This is called by ENTspec()
## This should not be seen or used by the user, as it is only intended for internal use by other functions.

processChannel.ENT = function(channelData,
                               samp.rate,
                               channel,
                               timeBin,
                               wl,
                               overlap,
                               window,
                               noiseOBJ) {
  allSamples = if (is.null(timeBin)) {
    data.frame(b = 1, e = length(channelData[1, ]))
  } else {
    getSampleBins(length(channelData[1, ]), samp.rate, timeBin)
  }

  frameBin = nrow(allSamples)

  channelData = switch(
    channel,
    "stereo" = list(left  = channelData[1, ], right = channelData[2, ]),
    "mono"   = list(mono  = channelData[1, ]),
    "left"   = list(left  = channelData[1, ]),
    "right"  = list(right = channelData[2, ])
  )

  noiseOBJ@values = lapply(channelData, function(x) {
    tempHolder = apply(allSamples, 1, function(y) {
      list(signal::specgram(
        x = x[y[1]:y[2]],
        n = wl,
        Fs = samp.rate,
        window = window,
        overlap = overlap
      )$S)
    })

    ENTdf = data.frame(do.call(cbind, lapply(tempHolder, function(singleBin) {
      spectS = abs(singleBin[[1]])
      amp2    = spectS ** 2

      rowTot = rowSums(amp2)
      rowTot[rowTot == 0] = 1 ## Safeguard!
      pmf = amp2 / rowTot

      term = ifelse(pmf > 0, pmf * log2(pmf), 0)
      H = -rowSums(term) / log2(dim(amp2)[2])

      return(1 - H)
    })))

    colnames(ENTdf) = paste0("ENT", rep(1:frameBin))
    return(list(ENT = ENTdf))

  })

  noiseOBJ@index = "ENT"
  noiseOBJ@timeBins = setNames(round((allSamples$e - allSamples$b) / samp.rate), paste0("BIN", seq(frameBin)))
  noiseOBJ@sampRate = samp.rate
  noiseOBJ@channel = channel

  return(noiseOBJ)

}

# bgnoise alts ------------------------------------------------------------
## bgNoise
## This is a hidden version of the bgNoise function.
## In this version, the check for the arguments is skipped and is handled
## by a higher function.

bgNoise. = function(soundfile,
                     channel = "stereo",
                     timeBin = 60,
                     dbThreshold = -90,
                     targetSampRate = NULL,
                     wl = 512,
                     window = signal::hamming(wl),
                     overlap = ceiling(length(window) / 2),
                     histbreaks = "FD",
                     DCfix = TRUE) {

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
    noiseOBJ = new("noise.matrix.internal")
  )

}

## bgNoise
## This is another hidden version of the bgNoise function.
## This version is used in functions that take a folder as input.
## It skips checks and reads the audio files directly.

bgNoise.. = function(soundfile,
                      channel = "stereo",
                      timeBin = 60,
                      dbThreshold = -90,
                      targetSampRate = NULL,
                      wl = 512,
                      window = signal::hamming(wl),
                      overlap = ceiling(length(window) / 2),
                      histbreaks = "FD",
                      DCfix = TRUE) {

  if (tolower(tools::file_ext(soundfile)) == "wav") {
    soundfile = wav::read_wav(soundfile)
  } else {
    stop("The audio file must be in the WAV format.")
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
    noiseOBJ = new("noise.matrix.internal")
  )

}

# argHandler --------------------------------------------------------------
## Function to check if the inputted arguments are supported.

argHandler = function(FUN, ...) {
  args = list(...)

  #### soundpath ----
  if ("soundpath" %in% names(args)) {
    if (!all(file.exists(args$soundpath)))
      stop("all provided soundpaths must be valid")
  }

  #### channel ----
  if (length(args$channel) != 1 ||
      !(args$channel %in% c("left", "right", "stereo", "mono")))
    stop(
      paste0(
        'channel = ',
        capture.output(dput(args$channel)),
        '\nchannel must be set to either "stereo", "mono", "left", or "right"'
      )
    )

  #### timeBin ----
  if (!is.null(args$timeBin) &&
      (!is.numeric(args$timeBin) ||
       length(args$timeBin) != 1 || args$timeBin < 0)) {
    stop(
      paste0(
        'timeBin = ',
        capture.output(dput(args$timeBin)),
        '\ntimeBin must be NULL or a single non-negative number'
      ),
      call. = FALSE
    )
  }

  #### j ----
  if ("j" %in% names(args)) {
    if (!is.null(args$j) &&
        (!is.numeric(args$j) ||
         length(args$j) != 1 || args$j < 0)) {
      stop(
        paste0(
          'j = ',
          capture.output(dput(args$j)),
          '\nj must be NULL or a single non-negative number'
        ),
        call. = FALSE
      )
    }
  }

  #### dbThreshold ----
  if ("dbThreshold" %in% names(args)) {
    if (!is.null(args$dbThreshold) &&
        (
          !is.numeric(args$dbThreshold) ||
          length(args$dbThreshold) != 1 || args$dbThreshold >= 0
        )) {
      stop(
        paste0(
          'dbThreshold = ',
          capture.output(dput(args$dbThreshold)),
          '\ndbThreshold must be a single negative number'
        ),
        call. = FALSE
      )
    }
  }

  #### targetSampRate ----
  if (!is.null(args$targetSampRate) &&
      (
        !is.numeric(args$targetSampRate) ||
        length(args$targetSampRate) != 1 ||
        args$targetSampRate < 0
      )) {
    stop(
      paste0(
        'targetSampRate = ',
        capture.output(dput(args$targetSampRate)),
        '\ntargetSampRate must be NULL or a single non-negative number'
      ),
      call. = FALSE
    )
  }

  #### wl ----
  if (!is.numeric(args$wl) || length(args$wl) != 1 || args$wl < 0)
    stop(
      paste0(
        'wl = ',
        capture.output(dput(args$wl)),
        '\nwl must be a single non-negative number'
      ),
      call. = FALSE
    )

  #### window ----
  if (!is.numeric(args$window) ||
      length(args$window) != args$wl) {
    stop(
      paste0(
        "On window = ... \nPlease set window to signal::hamming(wl) or signal::hanning(wl)"
      ),
      call. = FALSE
    )
  }

  #### overlap ----
  if (!is.numeric(args$overlap) ||
      length(args$overlap) != 1 || args$overlap < 0)
    stop(
      paste0(
        'overlap = ',
        capture.output(dput(args$overlap)),
        '\noverlap must be a single non-negative number'
      ),
      call. = FALSE
    )

  #### histbreaks ----
  if ("histbreaks" %in% names(args)) {
    if (length(args$histbreaks) != 1 ||
        !is.numeric(args$histbreaks) &&
        !(args$histbreaks %in% c("FD", "Sturges", "scott")))
      stop(
        paste0(
          'histbreaks = ',
          capture.output(dput(args$histbreaks)),
          "\nhistbreaks must be 'FD', 'Sturges', 'scott' or a single non-negative number"
        ),
        call. = FALSE
      )
  }

  #### DCfix ----
  if ("DCfix" %in% names(args)) {
    if (length(args$DCfix) != 1 || !is.logical(args$DCfix) || is.na(args$DCfix)) {
      stop(paste0(
        'DCfix = ',
        capture.output(dput(args$DCfix)),
        "\nDCfix must be either TRUE or FALSE"
      ),
      call. = FALSE)
    }
  }

  #### powthr ----
  if ("powthr" %in% names(args)) {
    if (FUN %in% c("soundSat", "soundMat")) {
      if (!is.numeric(args$powthr))
        stop(paste0(
          'powthr = ',
          capture.output(dput(args$powthr)),
          "\npowthr must be numeric"
        ),
        call. = FALSE)
      if (length(args$powthr) != 3)
        stop(
          paste0(
            'powthr = ',
            capture.output(dput(args$powthr)),
            "\nFor ",
            FUN,
            "() powthr must have length 3"
          ),
          call. = FALSE
        )
      if (!all(args$powthr > 0))
        stop(
          paste0(
            'powthr = ',
            capture.output(dput(args$powthr)),
            "\nAll powthr values must be positive"
          ),
          call. = FALSE
        )
      if (args$powthr[1] >= args$powthr[2])
        stop(
          paste0(
            'On powthr = ',
            capture.output(dput(args$powthr)),
            '\nThe first value of powthr must be lower than the second\nTry: c(',
            paste(args$powthr[2], args$powthr[1], args$powthr[3], sep = ", "),
            ")"
          ),
          call. = FALSE
        )
    } else {
      if (!is.numeric(args$powthr))
        stop(paste0(
          'powthr = ',
          capture.output(dput(args$powthr)),
          "\npowthr must be numeric"
        ),
        call. = FALSE)
      if (length(args$powthr) != 1)
        stop(
          paste0(
            'powthr = ',
            capture.output(dput(args$powthr)),
            "\nFor ",
            FUN,
            "() powthr must a single value"
          ),
          call. = FALSE
        )
      if (!all(args$powthr > 0))
        stop(paste0(
          'powthr = ',
          capture.output(dput(args$powthr)),
          "\npowthr must be positive"
        ),
        call. = FALSE)
    }
  }

  #### bgnthr ----
  if ("bgnthr" %in% names(args)) {
    if (FUN %in% c("soundSat", "soundMat")) {
      if (!is.numeric(args$bgnthr))
        stop(paste0(
          'bgnthr = ',
          capture.output(dput(args$bgnthr)),
          "\nbgnthr must be numeric"
        ),
        call. = FALSE)
      if (length(args$bgnthr) != 3)
        stop(
          paste0(
            'bgnthr = ',
            capture.output(dput(args$bgnthr)),
            "\nFor ",
            FUN,
            "() bgnthr must have length 3"
          ),
          call. = FALSE
        )
      if (!all(args$bgnthr > 0))
        stop(
          paste0(
            'bgnthr = ',
            capture.output(dput(args$bgnthr)),
            "\nAll bgnthr values must be positive"
          ),
          call. = FALSE
        )
      if (args$bgnthr[1] >= args$bgnthr[2])
        stop(
          paste0(
            'On bgnthr = ',
            capture.output(dput(args$bgnthr)),
            '\nThe first value of bgnthr must be lower than the second\nTry: c(',
            paste(args$bgnthr[2], args$bgnthr[1], args$bgnthr[3], sep = ", "),
            ")"
          ),
          call. = FALSE
        )
    } else {
      if (!is.numeric(args$bgnthr))
        stop(paste0(
          'bgnthr = ',
          capture.output(dput(args$bgnthr)),
          "\nbgnthr must be numeric"
        ),
        call. = FALSE)
      if (length(args$bgnthr) != 1)
        stop(
          paste0(
            'bgnthr = ',
            capture.output(dput(args$bgnthr)),
            "\nFor ",
            FUN,
            "() bgnthr must have a single value"
          ),
          call. = FALSE
        )
      if (!all(args$bgnthr > 0))
        stop(paste0(
          'bgnthr = ',
          capture.output(dput(args$bgnthr)),
          "\nbgnthr must be positive"
        ),
        call. = FALSE)
    }
  }

  #### beta ----
  if ("beta" %in% names(args)) {
    if (!is.logical(args$beta) || length(args$beta) != 1) {
      stop(paste0(
        'beta = ',
        capture.output(dput(args$beta)),
        "\nbeta must be either TRUE or FALSE"
      ),
      call. = FALSE)
    }

  }

  #### backup ----
  if ("backup" %in% names(args)) {
    if (!is.null(args$backup) && !dir.exists(args$backup))
      stop(
        paste0(
          'backup = ',
          capture.output(dput(args$backup)),
          "\nPlease provide a valid directory for backup."
        ),
        call. = FALSE
      )
  }

  invisible(NULL)

}

# normHandler -------------------------------------------------------------
## This function deals with the normality tests

normHandler = function(normality) {
  if (normality == "shapiro.test") {
    answernorm = readline(
      "
      If you are working with a large dataset, then shapiro.test will most likely result in an error.
      Do you wish to use Anderson-Darling test instead? (Y/N).
      "
    )

    if (answernorm == "Y") {
      normality = "ad.test"
    } else if (answernorm == "N") {
      message("Using shapiro.test to test normality.")
    } else {
      stop("Please answer with Y or N next time.", call. = FALSE)
    }

  } else if (normality == "ks.test") {
    answernorm = readline(
      "ks.test is not supported since many combinations may have identifical values.
      Type N to ignore this warning.
      However, we recommend choosing one of these tests:
      a ad.test
      b cvm.test
      c lillie.test
      d pearson.test
      e sf.test
      (Type the letter to choose)
      "
    )

    normality = switch(
      answernorm,
      "a" = "ad.test",
      "b" = "cvm.test",
      "c" = "lillie.test",
      "d" = "pearson.test",
      "e" = "sf.test",
      "N" = "ks.test",
      "STOP"
    )

    if (normality == "STOP") {
      stop("Please pick a letter next time.", call. = FALSE)
    }

  }

  return(normality)

}

# getSampleBins -----------------------------------------------------------
## This is a helper function to dynamically transform seconds to samples

getSampleBins = function(samples, samp.rate, binSize) {
  b = seq(1, samples, by = samp.rate * binSize)
  e = pmin(b + samp.rate * binSize - 1, samples)

  keepThese = ((samp.rate * binSize) * 0.1) < e - b ## This is so we can keep only bins that are at least 10% the size of the audio's samp.rate
  data.frame(b, e)[keepThese, ]

}

#' Plot noise.matrix objects
#'
#' @param x an `noise.matrix` object
#' @param channel channel or channels to be ploted. By default, this set to `x@channel`, but can be changed to `left` or `right` if `x@channel` = `stereo`
#' @param bin temporal bin to be plotted. Defaults to `1`
#' @param index a character vector of length 1 or 2 with indexes to be plotted. Available indices are: `c("BGN", "POW")`, `"BGN"`, `"POW"` and `"ACI"`. Defaults to the indexes listed in `x@index`
#' @param nbreaks amount of breaks of the y axis. Defaults to `5`
#' @param yunit frequency unit to be used in plot. Available units are: `"hz"` and `"khz"`. Defaults to `"hz"`
#' @param main title for the plot. Set two strings if you are plotting and stereo noise.matrix. If set to `NULL`, default title will be `left/right/mono channel`
#' @param xlab label for the x-axis. Changes depending on `x@index`
#' @param ylab label for the y-axis. Defaults to `"Frequency"`
#' @param col plotting color for de indices. Defaults to `c("blue","red")`
#' @param type desired plot type. For details see [base::plot]
#' @param draw0 if a stripped line should be drawn at 0. Defaults to `TRUE`
#' @param box if a box should be drawn around the plot. Defaults to `TRUE`
#' @param axes if axes should be drawn. Defaults to `TRUE`
#' @param annotate if bin information should be added to the plot. Defaults to `TRUE`
#' @param ... further [graphical parameters] passed down to plot
#'
#' @details This is a method to quickly plot the results of [bgNoise]. This calls the helper function `plotBGN`, which is not meant to be used or seen by the user.
#'
#' @importFrom grDevices xy.coords
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics plot.new
#' @importFrom graphics plot.window
#' @importFrom graphics plot.xy
#' @importFrom methods new
#' @importFrom utils head
#'
plotNOISE = function(x,
                      channel,
                      bin,
                      index,
                      nbreaks,
                      yunit,
                      main,
                      xlab,
                      ylab,
                      col,
                      type,
                      draw0,
                      box,
                      axes,
                      annotate,
                      ...) {
  channels = if (channel == "stereo") {
    c("left", "right")
  } else {
    channel
  }

  sampRate = x@sampRate
  sampDivide = sampRate / nbreaks
  sampleStep = seq(1, sampRate, length.out = x@wl)
  roundAt = floor(sampDivide / 1000) * 1000

  if (yunit == "khz") {
    sampRate = sampRate / 1000
    sampDivide = sampDivide / 1000
    sampleStep = sampleStep / 1000
    roundAt = roundAt / 1000
  }

  geometry = if (channel == "stereo") {
    c(1, 2)

  } else {
    c(1, 1)

  }

  opar = par(no.readonly = TRUE)
  on.exit(par(opar))

  par(mfrow = geometry)

  chLoop = 0

  for (ch in channels) {
    chLoop = chLoop + 1

    values = sapply(index, function(ind) {
      x@values[[ch]][[ind]][, bin]

    })

    minV = min(values)
    maxV = max(values)

    plot.new()
    plot.window(c(minV, maxV), c(0, sampRate), yaxs = "i", ...)

    if (box)
      box()

    if (axes) {
      axis(side = 1, ...)
      axis(side = 2,
           at = seq(0 , sampRate, by = roundAt),
           ...)
    }

    if (draw0) {
      if (minV < 0) {
        abline(v = 0, lty = 2)

      }
    }

    adj = if (length(index) > 1) {
      c(((minV / 2) - minV) / (maxV - minV), ((maxV / 2) - minV) / (maxV - minV))
    } else {
      0.5
    }

    for (ind in 1:ncol(values)) {
      plot.xy(xy.coords(x = values[, ind], y = sampleStep),
              type = type,
              col = col[ind],
              ...)

      mtext(
        index[ind],
        side = 3,
        adj = adj[ind],
        line = 0.5,
        col = col[ind]
      )

    }

    title = if (length(main) == 1 && main == "channel") {
      paste(ch, "channel")
    } else {
      main[chLoop]
    }

    titleY = paste(ylab, ifelse(yunit == "khz", "(kHz)", "(Hz)"))

    title(main = title, ylab = titleY, ...)
    mtext(xlab, side = 1, line = 2, ...)

    if (annotate) {
      mtext(
        paste0(
          "Bin: ",
          bin,
          " | Bin Duration: ",
          x@timeBins[bin],
          "s | Channel: ",
          ch,
          " | Sampling Rate: ",
          sampRate,
          ifelse(yunit == "khz", "kHz", "Hz"),
          " | Window Length: ",
          x@wl
        ),
        side = 1,
        line = 3,
        ...
      )

    }

  }

}
