.importAudio = function(soundfile, channel, targetSampRate) {

  audio = typeof(soundfile)

  if (audio == "character") {
    path = soundfile
    if (tolower(tools::file_ext(soundfile)) == "wav") {
      soundfile = wav::read_wav(soundfile)
    } else {
      stop("The audio file must be in the WAV format.")
    }
  } else if (audio == "S4") {
    path = "soundfile"
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

  attr(soundfile, "channel") = channel
  attr(soundfile, "path") = path

  return(soundfile)

}

### Called by: ACTspec and SAT
.BGNspec = function(soundfile,
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
                    noiseOBJ = new("noise.matrix")) {

  channelData = .importAudio(soundfile = soundfile,
                             channel = channel,
                             targetSampRate = targetSampRate)

  BGNexp = processChannel.BGN(
    channelData = channelData,
    samp.rate = attr(channelData, "sample_rate"),
    channel = attr(channelData, "channel"),
    timeBin = timeBin,
    wl = wl,
    overlap = overlap,
    dbThreshold = dbThreshold,
    dbRef = dbRef,
    window = window,
    histbreaks = histbreaks,
    DCfix = DCfix,
    noiseOBJ = noiseOBJ
  )

  BGNexp@path = attr(channelData, "path")

  return(BGNexp)

}

### Called by: soundSAT, soundMAT, soundACT and satBackup
..BGNspec = function(soundfile,
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

  if (tolower(tools::file_ext(soundfile)) == "wav") {
    channelData = wav::read_wav(soundfile)
  } else {
    stop("The audio file must be in the WAV format.")
  }

  savedAttr = attributes(channelData)

  if (channel == "mono" && savedAttr$dim[1] > 1) {
    sampRate = attr(channelData, "sample_rate")
    channelData = matrix((channelData[1, ] + channelData[2, ]) / 2, nrow = 1)
    attr(channelData, "sample_rate") = sampRate
    savedAttr$dim = dim(channelData)
  }
  if (channel == "stereo" && nrow(channelData) == 1) {
    message("Audio is not stereo, defaulting to left channel.")
    channel = "mono"
  }
  if (!is.null(targetSampRate)) {
    audioLen = length(channelData[1, ])
    keepIdx  = round(seq(1, audioLen, length = targetSampRate * audioLen / savedAttr$sample_rate))
    channelData = channelData[1:savedAttr$dim[1], keepIdx, drop = FALSE]
    attr(channelData, "sample_rate") = targetSampRate
  }

  BGNexp = processChannel.BGN(
    channelData,
    samp.rate = attr(channelData, "sample_rate"),
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

  BGNexp@wl = if (BGNexp@channel == "stereo") {
    nrow(BGNexp@values$left$BGN)
  } else {
    nrow(BGNexp@values[[channel]]$BGN)
  }


  BGNexp@path = soundfile

  return(BGNexp)

}

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

  #### dbRef ----
  if ("dbRef" %in% names(args)) {
    if (!is.numeric(args$dbRef) ||
        length(args$dbRef) != 1 || args$dbRef < 0)
      stop(
        paste0(
          'dbRef = ',
          capture.output(dput(args$dbRef)),
          '\ndbRef must be a single non-negative number'
        ),
        call. = FALSE
      )
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
    stop(paste0(
      'wl = ',
      capture.output(dput(args$wl)),
      '\nwl must be a single non-negative number'
    ),
    call. = FALSE)

  #### window ----
  if ("window" %in% names(args)) {
    if (!is.numeric(args$window) ||
        length(args$window) != args$wl) {
      stop(paste0(
        "On window = ... \nPlease set window to hamming(wl) or hanning(wl)"
      ),
      call. = FALSE)
    }
  }

  #### overlap ----
  if ("overlap" %in% names(args)) {
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
  }

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
    if (length(args$DCfix) != 1 ||
        !is.logical(args$DCfix) || is.na(args$DCfix)) {
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
    if (FUN %in% c("soundSAT", "soundMAT")) {
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
    if (FUN %in% c("soundSAT", "soundMAT")) {
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

  #### saveNoise ----
  if ("saveNoise" %in% names(args)) {
    if (!is.logical(args$saveNoise) || length(args$saveNoise) != 1) {
      stop(paste0(
        'saveNoise = ',
        capture.output(dput(args$saveNoise)),
        "\nsaveNoise must be either TRUE or FALSE"
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

normHandler = function(normality, SIZE) {
  minN = c(
    shapiro.test = 3,
    sf.test = 5,
    lillie.test = 5,
    ad.test = 8,
    cvm.test = 8,
    pearson.test = 2
  )
  maxN = c(
    shapiro.test = 5000,
    sf.test = 5000,
    lillie.test = Inf,
    ad.test = Inf,
    cvm.test = Inf,
    pearson.test = Inf
  )

  default = if (SIZE > 5000)
    "ad.test"
  else
    "shapiro.test"

  if (normality == "auto")
    return(default)

  if (isTRUE(SIZE >= minN[normality] &&
             SIZE <= maxN[normality]))
    return(normality)

  warning("'",
          normality,
          "' doesn't fit n = ",
          SIZE,
          ". Using '",
          default,
          "' instead.")
  default
}

.spect = function(x, n, window, overlap) {
  winSize = length(window)

  if (length(x) > winSize) {
    offset = seq.int(1, length(x) - winSize, by = (winSize - overlap))
  } else {
    offset = 1
  }

  S = .Call(
    "spectMat",
    as.numeric(x),
    as.integer(n),
    as.numeric(window),
    as.integer(offset),
    PACKAGE = "Ruido"
  )

  S = mvfft(S)

  keepThese = if (n %% 2 == 1) {
    (n + 1) / 2
  } else {
    n / 2
  }

  S[1:keepThese, , drop = FALSE]
}

hBreaks.MAT = function(histbreaks) {
  if (is.numeric(histbreaks)) {
    function(z) {
      rowsNumeric(z, n = histbreaks)
    }
  } else {
    switch(
      histbreaks,
      "FD"      = function(z) {
        rowsFD(z, digits = 5)
      }
      ,
      "Sturges" = function(z) {
        rowsSturges(z)
      }
      ,
      "scott"   = function(z)
      {
        rowsscott(z)
      }
    )

  }
}

hBreaks = function(histbreaks) {
  if (is.numeric(histbreaks)) {
    function(x) {
      histbreaks + 1
    }
  } else {
    switch(
      histbreaks,
      "FD"      = function(z) {
        nclass.FD(z) + 1
      }
      ,
      "Sturges" = function(z) {
        nclass.Sturges(z) + 1
      }
      ,
      "scott"   = function(z)
      {
        nclass.scott(z) + 1
      }
    )

  }
}

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
#' @details This is a method to quickly plot the results of [BGNspec]. This calls the helper function `plotBGN`, which is not meant to be used or seen by the user.
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

rowsFD = function(z, digits = 5) {
  h = 2 * matrixStats::rowIQRs(z. <- signif(z, digits = digits))

  if (any(h == 0)) {
    who = which(h == 0)

    al = 1 / 4
    alMin = 1 / 512

    while (length(who) > 0 &&
           (al = al / 2) >= alMin) {
      q = matrixStats::rowQuantiles(z.[who, , drop = FALSE], probs = c(al, 1 - al)
                                    , drop = FALSE)
      h[who] = (q[, 2] - q[, 1]) / (1 - 2 * al)
      who = who[h[who] == 0]
    }
  }

  if (any(h == 0)) {
    who = which(h == 0)
    h[who] = 3.5 * matrixStats::rowSds(z.[who, , drop = FALSE])
  }

  out = ceiling(matrixStats::rowDiffs(matrixStats::rowRanges(z)) / h  * ncol(z)^(1 / 3)) + 1
  out[h < 1] = 1

  return(c(out))
}

rowsSturges = function(z) {
  rep(ceiling(log(ncol(z), base = 2) + 1), nrow(z))
}

rowsscott = function(z) {
  h = 3.5 * sqrt(matrixStats::rowVars(z)) * ncol(z)^(-1 / 3)
  positive = which(h > 0)
  negative = which(h < 0)
  h[negative] = 1
  h[positive] = pmax(1, ceiling(matrixStats::rowDiffs(matrixStats::rowRanges(z)) / h)) + 1

  return(h)
}

rowsNumeric = function(z, n) {
  rep(n + 1, nrow(z))
}
