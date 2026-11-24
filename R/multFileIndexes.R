#' @title Soundscape Saturation Index
#'
#' @description Calculate Soundscape Saturation for a combination of recordings using the methodology proposed in Burivalova 2018.
#'
#' @param soundpath single or multiple directories to your `.wav` audio files
#' @param channel channel where the saturation values will be extracted from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`.
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Set to `NULL` to leave db values unrestricted Defaults to `-90`, as set by Towsey 2017
#' @param dbRef reference value used to compute decibel values from the spectrogram.
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `hanning(wl)` to use hanning instead. Defaults to `hamming(wl)`
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges"`, `"scott"` and `100`. Defaults to `"FD"`.
#' <br>Can also be set to any numerical value to limit or increase the amount of breaks.
#' @param DCfix if the DC offset should be removed before the metrics are calculated. Defaults to `TRUE`
#' @param powthr numeric vector of length three containing the the range of thresholds used to evaluate the Soundscape Power of the  Activity Matrix (in dB). The values correspond to the minimum threshold, maximum threshold and step size respectively.
#' <br> Defaults to `c(5, 20, 1)`, which evaluates thresholds from 5 dB to 20 dB in increments of 1 dB
#' @param bgnthr numeric vector of length three containing the the range of thresholds used to evaluate the Background Noise of the  Activity Matrix (in %). The values correspond to the minimum threshold, maximum threshold and step size respectively.
#' <br> Defaults to `c(0.5, 0.9, 0.05)`, which evaluates thresholds from 50% to 90% in increments of 5%
#' @param normality character string with the normality test used to determine which threshold combination has the most normal distribution of values. One of `"auto"`, `"shapiro.test"`, `"sf.test"`, `"lillie.test"`, `"ad.test"`, `"cvm.test"` or `"pearson.test"`. Defaults to `"auto"`.
#' <br>`"ks.test"` is not available. Check details for more information.
#' @param beta how BGN thresholds are calculated. If `TRUE`, BGN thresholds are calculated using all recordings combined. If FALSE, BGN thresholds are calculated separately for each recording. Defaults to `TRUE`
#' @param saveNoise whether the noise.matrix objects generated internally should also be included in the output. Defaults to `TRUE`
#' @param backup path to save the backup. Defaults to `NULL`
#'
#' @returns A list containing five objects. The first and second objects (powthresh and bgnthresh) are the threshold values that yielded the most normal distribution of saturation values using the normality test set by the user. The third (normality) contains the statistics values of the normality test that yielded the most normal distribution. The fourth object (values) contains a data.frame with the values of saturation for each bin of each recording and the size of the bin in seconds. The fifth contains a data.frame with errors that occurred with specific files during the function.
#'
#' @details  Soundscape Saturation (`SAT`) quantifies the proportion of frequency bins that are acoustically active within a given time bin. It wasproposed by Burivalova et al. (2018) as a metric to evaluate the acoustic niche hypothesis.
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
#' After computing \eqn{S_m}, the function evaluates all tested threshold values and selects the one that yields the most normally distributed set of saturation values across time bins. Normality is assessed to identify a threshold that best stabilizes the distribution of `SAT`.
#'
#' If `backup` is set to a valid directory, a file named `"SATBACKUP.rds"` is automatically saved after every batch of five processed files. This file stores the current processing state and allows interrupted runs (e.g., due to manual termination, session crashes, or system shutdowns) to be resumed using [satBackup()].
#'
#' To resume processing, pass the saved file (e.g., `"path/SATBACKUP.rds"`) to [satBackup()]. Once a backup has been created, all original arguments and file paths must remain unchanged, unless they are explicitly modified within the saved `.rds` object.
#'
#' For the `normality` argument, each test only works within a range of sample sizes (3-5000 for `"shapiro.test`", 5-5000 for `"sf.test"`, 5 or more for `"lillie.test"`, 8 or more for `"ad.test"` and `"cvm.test"`, 2 or more for `"pearson.test"`). If the chosen test does not fit the number of values, a warning is issued and `"shapiro.test"` is used (or `"ad.test"` if there are more than 5000 values). With `"auto"`, `"shapiro.test"` is used for up to 5000 values and `"ad.test"` for more.
#'
#' @seealso [satBackup()] if you are working with bigger datasets and may need to generate backups.
#'
#'@references Burivalova, Z., Towsey, M., Boucher, T., Truskinger, A., Apelis, C., Roe, P., & Game, E. T. (2018). Using soundscapes to detect variable degrees of human influence on tropical forests in Papua New Guinea. Conservation Biology, 32(1), 205-215. https://doi.org/10.1111/cobi.12968
#'
#'@export
#'@importFrom methods is
#'@importFrom stats quantile
#'@importFrom stats setNames
#'@importFrom stats shapiro.test
#'@importFrom nortest ad.test
#'@importFrom nortest sf.test
#'@importFrom nortest lillie.test
#'@importFrom nortest cvm.test
#'@importFrom nortest pearson.test
#'@importFrom utils capture.output
#'
#' @examples
#' \donttest{
#' ### Downloading audiofiles from public Zenodo library
#' dir = paste(tempdir(), "forExample", sep = "/")
#' dir.create(dir)
#' recName = paste0("GAL24576_20250401_", sprintf("%06d", seq(0, 200000, by = 50000)),".wav")
#' recDir = paste(dir, recName, sep = "/")
#'
#' for(rec in recDir) {
#'  print(rec)
#'  url = paste0("https://zenodo.org/records/17575795/files/", basename(rec), "?download=1")
#'  download.file(url, destfile = rec, mode = "wb")
#' }
#'
#' ### Running the function
#' sat = soundSAT(dir)
#'
#' ### Preparing the plot
#' timeSplit = strsplit(sat$values$AUDIO, "_")
#' sides = sat$values$CHANNEL
#' date = sapply(timeSplit, function(x)
#'   x[2])
#' time = sapply(timeSplit, function(x)
#'   substr(x[3],1,6))
#' datePos = paste(substr(date,1,4), substr(date,5,6), substr(date,7,8), sep = "-")
#' timePos = paste(substr(time,1,2), substr(time,3,4), substr(time,5,6), sep = ":")
#' dateTime = as.POSIXct(paste(datePos, timePos), format = "%Y-%m-%d %H:%M:%OS")
#' leftEar = data.frame(SAT = sat$values$SAT[sides == "left"], HOUR = dateTime[sides == "left"])
#' rightEar = data.frame(SAT = sat$values$SAT[sides == "right"], HOUR = dateTime[sides == "right"])
#'
#' ### Plotting results
#'
#' plot(SAT~HOUR, data = leftEar, ylim = c(range(sat$values$SAT)),
#' col = "darkgreen", pch = 16,
#'      ylab = "Soundscape Saturation (%)", xlab = "Time of Day", type = "b")
#' points(SAT~HOUR, data = rightEar, ylim = c(range(sat$values$SAT)),
#' col = "red", pch = 16, type = "b")
#' legend("bottomright", legend = c("Left Ear", "Right Ear"),
#'        col = c("darkgreen", "red"), lty = 1)
#'
#' unlink(dir, recursive = TRUE)
#' }
soundSAT = function(soundpath,
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
                    powthr = c(5, 20, 1),
                    bgnthr = c(0.5, 0.9, 0.05),
                    normality = "auto",
                    beta = TRUE,
                    saveNoise = FALSE,
                    backup = NULL) {
  argHandler(
    FUN = "soundSAT",
    soundpath = soundpath,
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
    normality = normality,
    beta = beta,
    saveNoise = saveNoise,
    backup = backup
  )

  powthreshold = seq(powthr[1], powthr[2], powthr[3])
  names(powthreshold) = powthreshold
  bgnthreshold = seq(bgnthr[1], bgnthr[2], bgnthr[3])

  soundpath = normalizePath(soundpath, winslash = "/", mustWork = TRUE)
  soundfiles = list.files(soundpath, full.names = TRUE, recursive = TRUE)
  soundfiles = soundfiles[tolower(tools::file_ext(soundfiles)) == "wav"]

  thresholdCombinations = setNames(expand.grid(powthreshold, bgnthreshold),
                                   c("powthreshold", "bgnthreshold"))

  combinations = paste(thresholdCombinations[, 1], thresholdCombinations[, 2], sep = "/")

  message(
    paste(
      "Calculating saturation values for",
      length(soundfiles),
      "recordings using",
      length(combinations),
      "threshold combinations"
    )
  )

  halfWl = as.integer(wl / 2)

  if (!is.null(backup)) {
    relativeFiles = substring(soundfiles, nchar(ifelse(
      endsWith(soundpath, "/"),
      soundpath,
      paste0(soundpath, "/")
    )) + 1)
    folders = unique(dirname(relativeFiles))
    folders = folders[folders != "."]

    backupDir = file.path(normalizePath(backup, winslash = "/", mustWork = TRUE),
                          "SATbackup")
    dataDir = file.path(backupDir, "data")
    dir.create(dataDir, recursive = TRUE, showWarnings = FALSE)

    message(
      paste0(
        'Backup folder generated at ', backupDir,
        '\nCall satBackup("',backupDir, '") if this function is interrupted\n'
      )
    )

    for (fold in folders) {
      dir.create(paste0(dataDir, "/", fold),
                 recursive = TRUE,
                 showWarnings = FALSE)
    }

    ARGS = list(
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
      normality = normality,
      beta = beta,
      saveNoise = saveNoise,
      od = soundpath,
      backup = backupDir,
      files = soundfiles,
      relativeFiles = relativeFiles
    )
    saveRDS(ARGS, paste0(backupDir, "/args.rds"))
  }

  nFiles = length(soundfiles)

  SATdf = lapply(seq(nFiles), function(soundfile) {
    sPath = soundfiles[[soundfile]]
    result = tryCatch(
      ..BGNspec(
        sPath,
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
      ),
      error = function(e) {
        e$message = paste0(e$message, " (file: ", sPath, ")")
        e
      }
    )

    if (!is.null(backup)) {
      saveRDS(result, file = file.path(dataDir, paste0(relativeFiles[[soundfile]], ".rds")))
    }

    message(
      "\r(",
      basename(soundfiles[soundfile]),
      ") ",
      soundfile,
      " out of ",
      nFiles,
      " recordings concluded."
    )

    return(result)
  })

  whichError = sapply(SATdf, function(x) {
    is(x, "error")
  })

  ERRORS = SATdf[whichError]
  indexes = SATdf[!whichError]

  BGN = as.matrix(do.call(cbind, sapply(indexes, function(x) {
    if (x@channel == "stereo") {
      cbind(x@values$left$BGN, x@values$right$BGN)
    } else {
      x@values[[x@channel]]$BGN
    }
  })))

  POW = as.matrix(do.call(cbind, sapply(indexes, function(x) {
    if (x@channel == "stereo") {
      cbind(x@values$left$POW, x@values$right$POW)
    } else {
      x@values[[x@channel]]$POW
    }
  })))

  INFO = lapply(indexes, function(x) {
    nBins = length(x@timeBins)
    if (x@channel == "stereo") {
      list(
        rep(x@timeBins, 2),
        rep(x@sampRate, length(x@timeBins) * 2),
        rep(1:length(x@timeBins), 2),
        rep(c("left", "right"), each = nBins)
      )
    } else {
      list(
        x@timeBins,
        rep(x@sampRate, length(x@timeBins)),
        1:length(x@timeBins),
        rep(x@channel, nBins)
      )
    }
  })

  paths = unlist(sapply(indexes, function(x) {
    if (x@channel == "stereo") {
      rep(x@path, length(x@timeBins) * 2)
    } else {
      rep(x@path, length(x@timeBins))
    }
  }))

  SATinfo = data.frame(
    PATH = dirname(paths),
    AUDIO = basename(paths),
    CHANNEL = c(unlist(sapply(INFO, function(x) {
      x[[4]]
    }))),
    DURATION = c(unlist(sapply(INFO, function(x) {
      x[[1]]
    }))),
    BIN = c(unlist(sapply(INFO, function(x) {
      x[[3]]
    }))),
    SAMPRATE = c(unlist(sapply(INFO, function(x) {
      x[[2]]
    })))
  )

  if (beta) {
    BGNQ = quantile(unlist(BGN), probs = seq(bgnthr[1], bgnthr[2], bgnthr[3])) |>
      setNames(bgnthreshold)

    ACTdf = mapply(
      function(bgnthresh, powthresh) {
        list(BGN > BGNQ[as.character(bgnthresh)] |
               POW > powthresh)

      },
      thresholdCombinations$bgnthreshold,
      thresholdCombinations$powthreshold
    )

  } else {
    ACTdf = mapply(
      function(bgnthresh, powthresh) {
        list((
          BGN > matrix(
            matrixStats::colQuantiles(BGN, probs = bgnthresh),
            nrow = nrow(BGN),
            ncol = ncol(BGN),
            byrow = TRUE
          )
        ) | POW > powthresh)
      },
      thresholdCombinations$bgnthreshold,
      thresholdCombinations$powthreshold
    )

  }

  SATmat = lapply(ACTdf, function(ACT) {
    matrixStats::colSums2(ACT) / halfWl
  })

  export = list(values = data.frame(), errors = list())

  if (!is.null(normality)) {
    normality = normHandler(normality, SIZE = length(SATmat[[1]]))

    normal = lapply(SATmat, function(Q) {
      if (length(unique(Q)) != 1) {
        do.call(normality, list(Q))$statistic
      } else {
        NA
      }
    })


    if (normality %in% c("sf.test", "shapiro.test")) {
      normalID = which.max(normal)
      thresholds = unlist(strsplit(combinations[normalID], split = "/"))
    } else {
      normalID = which.min(normal)
      thresholds = unlist(strsplit(combinations[normalID], split = "/"))
    }

    normOUT = normal[normalID]

    normname = switch(
      normality,
      "shapiro.test" = "Shapiro-Wilk",
      "sf.test" = "Shapiro-Francia",
      "ad.test" = "Anderson-Darling",
      "cvm.test" = "Cram\u00e9r-von Mises",
      "lillie.test" = "Lilliefors",
      "pearson.test" = "Pearson chi-square"
    )
    normstat = switch(
      normality,
      "shapiro.test" = "W",
      "sf.test" = "W'",
      "ad.test" = "A",
      "cvm.test" = "W\u00b2",
      "lillie.test" = "D",
      "pearson.test" = "X\u00b2"
    )

    message(
      "\n           Soundscape Saturation Results\n\n",
      "POW Threshold = ",
      as.numeric(thresholds[1]),
      " dB        ",
      "BGN Threshold = ",
      as.numeric(thresholds[2]) * 100,
      "%\n",
      normname,
      " Test Statistic (",
      normstat ,
      ") = ",
      unlist(normOUT),
      "\n "
    )

    SATinfo$SAT = SATmat[[normalID]]

    export[["normality"]] = list()
    export[["powthresh"]] = as.numeric(thresholds[1])
    export[["bgnthresh"]] = as.numeric(thresholds[2]) * 100
    export[["normality"]]["statistic"] = normOUT
    export[["normality"]]["test"] = normality

  } else {
    SATmat = do.call(cbind, SATmat)
    colnames(SATmat) = combinations
    SATinfo = cbind(SATinfo, SATmat)

    export[["powthresh"]] = powthreshold
    export[["bgnthresh"]] = bgnthreshold * 100
  }

  export[["values"]] = SATinfo
  export[["errors"]] = ERRORS

  if (saveNoise) {
    names(SATdf) = soundfiles
    export[["noise.matrix"]] = SATdf
    names(ACTdf) = combinations
    export[["act"]] = ACTdf
  }

  return(export)

}

#' @title Backup for soundSAT
#'
#' @param backup path to the `.rds` file create by the backup of soundSAT
#'
#' @description
#' This function offers a way to continue an unfinished process of the [soundSAT()] function through a backup file.
#' Arguments can't be inputted nor changed since the function will automatically load them from the `.rds` file. However you may manually change them by editing the file (not recommended).
#'
#' @returns
#' This functions returns the same output of [soundSAT()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # It's impossible to demonstrate this function's intended use due to it's nature
#' # However, here is how this function is used:
#' ## This example will load an entire day of audios to your computer, so beware.
#'
#' ### Downloading audiofiles from public Zenodo library
#' dir = paste(tempdir(), "forExample", sep = "/")
#' dir.create(dir)
#' recName = paste0("GAL24576_20250401_", sprintf("%06d", seq(0, 230000, by = 10000)),".wav")
#' recDir = paste(dir, recName, sep = "/")
#'
#' for(rec in recName) {
#'   print(rec)
#'   url = paste0("https://zenodo.org/records/17575795/files/", rec, "?download=1")
#'   download.file(url, destfile = paste(dir, rec, sep = "/"), mode = "wb")
#' }
#'
#' sat = soundSAT(dir, backup = dir)
#'
#' # Now pretend the process was interrupted (manually/your R crashed/your computer turned off)
#' # We get the backup file
#'
#' list.files(dir)
#' backupDir = paste(dir, "SATBACKUP.rds", sep = "/")
#'
#' # To recall the backup you simply:
#'
#' satB = satBackup(backupDir)
#'
#' head(satB$values)
#'
#' unlink(dir, recursive = TRUE)
#' }
satBackup = function(backup) {
  backupDir = normalizePath(backup, winslash = "/", mustWork = TRUE)
  dataDir = file.path(backupDir, "data")

  ogARGS = readRDS(file.path(backupDir, "args.rds"))

  powthreshold = seq(ogARGS$powthr[1], ogARGS$powthr[2], ogARGS$powthr[3])
  names(powthreshold) = powthreshold
  bgnthreshold = seq(ogARGS$bgnthr[1], ogARGS$bgnthr[2], ogARGS$bgnthr[3])

  soundpath = ogARGS$od
  soundfiles = ogARGS$files

  Nsoundfiles = list.files(soundpath, full.names = TRUE, recursive = TRUE)
  Nsoundfiles = soundfiles[tolower(tools::file_ext(soundfiles)) == "wav"]

  soundfiles = if (all(soundfiles == Nsoundfiles)) {
    soundfiles
  } else {
    Nsoundfiles
  }

  Dsoundfiles = list.files(file.path(ogARGS$backup, "data"),
                           pattern = "rds",
                           recursive = TRUE)
  nFilesLeft = which(!(ogARGS$relativeFiles %in% tools::file_path_sans_ext(Dsoundfiles)))
  soundfiles = soundfiles[!(ogARGS$relativeFiles %in% tools::file_path_sans_ext(Dsoundfiles))]

  rm(Nsoundfiles)

  thresholdCombinations = setNames(expand.grid(powthreshold, bgnthreshold),
                                   c("powthreshold", "bgnthreshold"))

  combinations = paste(thresholdCombinations[, 1], thresholdCombinations[, 2], sep = "/")

  message(
    paste(
      "Calculating saturation values for",
      length(soundfiles),
      "recordings using",
      length(combinations),
      "threshold combinations"
    )
  )

  halfWl = as.integer(ogARGS$wl / 2)

  nFiles = length(soundfiles)

  SATdf = lapply(seq(nFiles), function(soundfile) {
    sPath = soundfiles[[soundfile]]
    result = tryCatch(
      ..BGNspec(
        sPath,
        timeBin = ogARGS$timeBin,
        targetSampRate = ogARGS$targetSampRate,
        window = ogARGS$window,
        overlap = ogARGS$overlap,
        channel = ogARGS$channel,
        dbThreshold = ogARGS$dbThreshold,
        dbRef = ogARGS$dbRef,
        wl = ogARGS$wl,
        histbreaks = ogARGS$histbreaks,
        DCfix = ogARGS$DCfix
      ),
      error = function(e) {
        e$message = paste0(e$message, " (file: ", sPath, ")")
        e
      }
    )

    saveRDS(result, file = file.path(dataDir, paste0(ogARGS$relativeFiles[[nFilesLeft[soundfile]]], ".rds")))

    message(
      "\r(",
      basename(soundfiles[soundfile]),
      ") ",
      soundfile,
      " out of ",
      nFiles,
      " recordings concluded."
    )

    return(result)
  })

  SATdf = c(lapply(file.path(dataDir, Dsoundfiles), function(soundfile) {
    readRDS(soundfile)
  }), SATdf)

  whichError = sapply(SATdf, function(x) {
    is(x, "error")
  })

  ERRORS = SATdf[whichError]
  indexes = SATdf[!whichError]

  BGN = as.matrix(do.call(cbind, sapply(indexes, function(x) {
    if (x@channel == "stereo") {
      cbind(x@values$left$BGN, x@values$right$BGN)
    } else {
      x@values[[x@channel]]$BGN
    }
  })))

  POW = as.matrix(do.call(cbind, sapply(indexes, function(x) {
    if (x@channel == "stereo") {
      cbind(x@values$left$POW, x@values$right$POW)
    } else {
      x@values[[x@channel]]$POW
    }
  })))

  INFO = lapply(indexes, function(x) {
    nBins = length(x@timeBins)
    if (x@channel == "stereo") {
      list(
        rep(x@timeBins, 2),
        rep(x@sampRate, length(x@timeBins) * 2),
        rep(1:length(x@timeBins), 2),
        rep(c("left", "right"), each = nBins)
      )
    } else {
      list(
        x@timeBins,
        rep(x@sampRate, length(x@timeBins)),
        1:length(x@timeBins),
        rep(x@channel, nBins)
      )
    }
  })

  paths = unlist(sapply(indexes, function(x) {
    if (x@channel == "stereo") {
      rep(x@path, length(x@timeBins) * 2)
    } else {
      rep(x@path, length(x@timeBins))
    }
  }))

  SATinfo = data.frame(
    PATH = dirname(paths),
    AUDIO = basename(paths),
    CHANNEL = c(unlist(sapply(INFO, function(x) {
      x[[4]]
    }))),
    DURATION = c(unlist(sapply(INFO, function(x) {
      x[[1]]
    }))),
    BIN = c(unlist(sapply(INFO, function(x) {
      x[[3]]
    }))),
    SAMPRATE = c(unlist(sapply(INFO, function(x) {
      x[[2]]
    })))
  )

  if (ogARGS$beta) {
    BGNQ = quantile(unlist(BGN), probs = seq(ogARGS$bgnthr[1], ogARGS$bgnthr[2], ogARGS$bgnthr[3])) |>
      setNames(bgnthreshold)

    ACTdf = mapply(
      function(bgnthresh, powthresh) {
        list(BGN > BGNQ[as.character(bgnthresh)] |
               POW > powthresh)

      },
      thresholdCombinations$bgnthreshold,
      thresholdCombinations$powthreshold
    )

  } else {
    ACTdf = mapply(
      function(bgnthresh, powthresh) {
        list((
          BGN > matrix(
            matrixStats::colQuantiles(BGN, probs = bgnthresh),
            nrow = nrow(BGN),
            ncol = ncol(BGN),
            byrow = TRUE
          )
        ) | POW > powthresh)
      },
      thresholdCombinations$bgnthreshold,
      thresholdCombinations$powthreshold
    )

  }

  SATmat = lapply(ACTdf, function(ACT) {
    matrixStats::colSums2(ACT) / halfWl
  })

  export = list(values = data.frame(), errors = list())

  if (!is.null(ogARGS$normality)) {
    normality = normHandler(ogARGS$normality, SIZE = length(SATmat[[1]]))

    normal = lapply(SATmat, function(Q) {
      if (length(unique(Q)) != 1) {
        do.call(normality, list(Q))$statistic
      } else {
        NA
      }
    })


    if (normality %in% c("sf.test", "shapiro.test")) {
      normalID = which.max(normal)
      thresholds = unlist(strsplit(combinations[normalID], split = "/"))
    } else {
      normalID = which.min(normal)
      thresholds = unlist(strsplit(combinations[normalID], split = "/"))
    }

    normOUT = normal[normalID]

    normname = switch(
      normality,
      "shapiro.test" = "Shapiro-Wilk",
      "sf.test" = "Shapiro-Francia",
      "ad.test" = "Anderson-Darling",
      "cvm.test" = "Cram\u00e9r-von Mises",
      "lillie.test" = "Lilliefors",
      "pearson.test" = "Pearson chi-square"
    )
    normstat = switch(
      normality,
      "shapiro.test" = "W",
      "sf.test" = "W'",
      "ad.test" = "A",
      "cvm.test" = "W\u00b2",
      "lillie.test" = "D",
      "pearson.test" = "X\u00b2"
    )

    message(
      "\n           Soundscape Saturation Results\n\n",
      "POW Threshold = ",
      as.numeric(thresholds[1]),
      " dB        ",
      "BGN Threshold = ",
      as.numeric(thresholds[2]) * 100,
      "%\n",
      normname,
      " Test Statistic (",
      normstat ,
      ") = ",
      unlist(normOUT),
      "\n "
    )

    SATinfo$SAT = SATmat[[normalID]]

    export[["normality"]] = list()
    export[["powthresh"]] = as.numeric(thresholds[1])
    export[["bgnthresh"]] = as.numeric(thresholds[2]) * 100
    export[["normality"]]["statistic"] = normOUT
    export[["normality"]]["test"] = normality

  } else {
    SATmat = do.call(cbind, SATmat)
    colnames(SATmat) = combinations
    SATinfo = cbind(SATinfo, SATmat)

    export[["powthresh"]] = powthreshold
    export[["bgnthresh"]] = bgnthreshold * 100
  }

  export[["values"]] = SATinfo
  export[["errors"]] = ERRORS

  if (ogARGS$saveNoise) {
    names(SATdf) = soundfiles
    export[["noise.matrix"]] = SATdf
    names(ACTdf) = combinations
    export[["act"]] = ACTdf
  }

  return(export)

}
