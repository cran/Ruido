#' @title Soundscape Saturation Index
#'
#' @description Calculate Soundscape Saturation for a combination of recordings using Burivalova 2018 methodology
#'
#' @param soundpath single directory or multiple directory to audio files. The directory must lead to a single folder or a combination of folders.
#' @param channel channel where the saturation values will be extract from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`.
#' @param timeBin size (in seconds) of the time bin. Defaults to `60`.
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Defaults to `-90`, as set by Towsey.
#' @param targetSampRate sample rate of the audios. Defaults to `NULL` to not change the sample rate. This argument is only used to down sample the audio.
#' @param wl window length of the spectrogram. Defaults to `512`.
#' @param window window used to smooth the spectrogram. Defaults to `signal::hammning(wl)`. Switch to `signal::hanning(wl)` if to use hanning instead.
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges"`, `"scott"` and `100`. Defaults to `"FD"`.
#' <br>Can also be set to any number to limit or increase the amount of breaks.
#' @param powthr vector of values to evaluate the activity matrix for Soundscape Power (in dB). The first value corresponds to the lowest dB value and the second corresponds to the highest, the third value is the step.
#' <br>Defaults to `c(5, 20, 1)`, which means the first thresholds starts at 5dB and jumps a whole number till 20dB.
#' @param bgnthr vector of values to evaluate the activity matrix for Background Noise (in %). The first value corresponds to the lowest quantile value and the second corresponds to the highest, the third value is the step.
#' <br>Defaults to `c(0.5, 0.9, 0.05)`, which means the first thresholds starts at 50% and jumps 5% till 90%.
#' @param normality normality test to determine which threshold combination has the most normal distribution of values. We recommend to pick any test from the `nortest` package. Input the test as a character. Defaults to `"ad.test"`.
#' <br>`"ks.test"` is not available. `"shapiro.test"` can be used, however we recommend using only when analyzing very few recordings.
#' @param backup path to backup the saturation values in case the computer is turned off during processing or in case you cannot be sure the computer will be on for the entire process. Defaults to `NULL`.
#' <br>The backup will be updated every 5 recordings processed.
#'
#' @returns A list containing five objects. The first and second objects (powthresh and bgnthresh) are the threshold values that yielded the most normal distribution of saturation values using the normality test set by the user. The third (normality) contains the statitics values of the normality test that yielded the most normal distribution. The fourth object (values) contains a data.frame with the the values of saturation for each bin of each recording and the size of the bin in seconds. The fifth contains a data.frame with errors that occurred with specific files during the function.
#'
#' @details Soundscape Saturation (`SAT`) is a measure of the proportion of frequency bins that are acoustically active in a determined window of time. It was developed by Burivalova et al. 2017 as an index to test the acoustic niche hypothesis.
#' To calculate this function, first we need to generate an activity matrix for each time bin of your recording with the following formula:
#'
#'\deqn{a_{mf} = 1\  if (BGN_{mf} > \theta_{1})\  or\  (POW_{mf} > \theta_{2});\  otherwise,\  a_{mf} = 0,}
#'
#'Where \eqn{\theta_{1}} is the threshold of BGN values and \eqn{\theta_{2}} is a threshold of dB values.
#'Since we define a interval for both the threshold, this means that an activity matrix will be generated for each bin of each recording.
#'For each combination of threshold a SAT measure will be taken with the following formula:
#'
#'\deqn{S_{m} = \frac{\sum_{f = 1}^N a_{mf}}{N}}
#'
#'After these equations are done, we check every threshold combination for normality and pick the combination that yields the most normal distribution of saturation values.
#'
#'If you set a path for the "path" argument, a single rds file will be written in your path. This objects can be loaded again through the "satBack" function to continue the calculation of saturation in case the process is suddenly stopped.
#'
#'@references Burivalova, Z., Towsey, M., Boucher, T., Truskinger, A., Apelis, C., Roe, P., & Game, E. T. (2017). Using soundscapes to detect variable degrees of human influence on tropical forests in Papua New Guinea. Conservation Biology, 32(1), 205-215. https://doi.org/10.1111/cobi.12968
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
#' \donttest{
#' ### Downloading audiofiles from public Zenodo library
#' dir <- tempdir()
#' recName <- paste0("GAL24576_20250401_", sprintf("%06d", seq(0, 200000, by = 50000)),".wav")
#' recDir <- paste(dir, recName, sep = "/")
#'
#' for(rec in recName) {
#'  print(rec)
#'  url <- paste0("https://zenodo.org/records/17575795/files/", rec, "?download=1")
#'  download.file(url, destfile = paste(dir, rec, sep = "/"), mode = "wb")
#' }
#'
#' ### Running the function
#' sat <- soundSat(dir, wl = 256)
#'
#' ### Preparing the plot
#' timeSplit <- strsplit(sat$values$AUDIO, "_")
#' sides <- ifelse(grepl("left", sat$values$BIN), "left", "right")
#' date <- sapply(timeSplit, function(x)
#'   x[2])
#' time <- sapply(timeSplit, function(x)
#'   substr(x[3],1,6))
#' datePos <- paste(substr(date,1,4), substr(date,5,6), substr(date,7,8), sep = "-")
#' timePos <- paste(substr(time,1,2), substr(time,3,4), substr(time,5,6), sep = ":")
#' dateTime <- as.POSIXct(paste(datePos, timePos), format = "%Y-%m-%d %H:%M:%OS")
#' leftEar <- data.frame(SAT = sat$values$SAT[sides == "left"], HOUR = dateTime[sides == "left"])
#' rightEar <- data.frame(SAT = sat$values$SAT[sides == "right"], HOUR = dateTime[sides == "right"])
#'
#' ### Plotting results
#'
#' plot(SAT~HOUR, data = leftEar, ylim = c(range(sat$values$SAT)),
#' col = "darkgreen", pch = 16,
#'      ylab = "Soundscape Saturation (%)", xlab = "Time of Day", type = "b")
#' points(SAT~HOUR, data = rightEar, ylim = c(range(sat$values$SAT)),
#' col = "red", pch = 16, type = "b")
#' legend("topright", legend = c("Left Ear", "Right Ear"),
#'        col = c("darkgreen", "red"), lty = 1)
#'
#' unlink(recDir)
#' }
soundSat <- function(soundpath,
                     channel = "stereo",
                     timeBin = 60,
                     dbThreshold = -90,
                     targetSampRate = NULL,
                     wl = 512,
                     window = signal::hamming(wl),
                     overlap = ceiling(length(window) / 2),
                     histbreaks = "FD",
                     powthr = c(5, 20, 1),
                     bgnthr = c(0.5, 0.9, 0.05),
                     normality = "ad.test",
                     backup = NULL) {
  if (all(!dir.exists(soundpath)))
    stop("all provided soundpaths must be valid.")

  if (!is.null(backup) && !dir.exists(backup))
    stop("please provide a valid folder for backup.")

  soundfiles <- list.files(soundpath, full.names = TRUE, recursive = TRUE)
  soundfiles <- soundfiles[tolower(tools::file_ext(soundfiles)) %in% c("mp3", "wav")]

  if (length(soundfiles) < 3)
    stop("please provide at least 3 recordings!")

  if (normality == "shapiro.test") {
    answernorm <- readline(
      "
      If you are working with a large dataset, then shapiro.test will most likely result in an error.
      Do you wish to use Anderson-Darling test instead? (Y/N).
      "
    )

    if (answernorm == "Y") {
      normality <- "ad.test"
    } else if (answernorm == "N") {
      message("Using shapiro.test to test normality.")
    } else {
      stop("Please answer with Y or N next time.")
    }

  } else if (normality == "ks.test") {
    answernorm <- readline(
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

    normality <- switch(
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
      stop("Please pick a letter next time.")
    }

  }

  powthreshold <- seq(powthr[1], powthr[2], powthr[3])
  names(powthreshold) <- powthreshold
  bgnthreshold <- seq(bgnthr[1], bgnthr[2], bgnthr[3])

  thresholdCombinations <- setNames(expand.grid(powthreshold, bgnthreshold),
                                    c("powthreshold", "bgnthreshold"))

  combinations <- paste(thresholdCombinations[, 1], thresholdCombinations[, 2], sep = "/")

  message(
    paste(
      "Calculating saturation values for",
      length(soundfiles),
      "recordings using",
      length(combinations),
      "threshold combinations"
    )
  )

  halfWl <- round(wl / 2)

  SATdf <- list()

  if (!is.null(backup)) {
    SATdf[["ogARGS"]] <- list(
      channel = channel,
      timeBin = timeBin,
      dbThreshold = dbThreshold,
      targetSampRate = targetSampRate,
      wl = wl,
      window = window,
      overlap = overlap,
      histbreaks = histbreaks,
      powthr = powthr,
      bgnthr = bgnthr,
      normality = normality,
      type = "soundSat"
    )
  }

  fiveSteps <- 1

  for (soundfile in soundfiles) {
    fiveSteps <- fiveSteps + 1

    BGNPOW <- tryCatch(
      bgNoise(
        soundfile,
        timeBin = timeBin,
        targetSampRate = targetSampRate,
        window = window,
        overlap = overlap,
        channel = channel,
        dbThreshold = dbThreshold,
        wl = wl,
        histbreaks = histbreaks
      ),
      error = function(e)
        e
    )

    SATdf[[soundfile]] <- if (is(BGNPOW, "error")) {
      warning("\n",
              basename(soundfile),
              "is not valid!\nError:",
              BGNPOW$message,
              "\n")

      BGNPOW

    } else {
      if (BGNPOW$channel == "stereo") {
        BGNQleft <- apply(BGNPOW$left$BGN, 2, function(n)
          setNames(
            quantile(n, probs = seq(bgnthr[1], bgnthr[2], bgnthr[3])),
            seq(bgnthr[1], bgnthr[2], bgnthr[3])
          ))

        BGNQright <- apply(BGNPOW$right$BGN, 2, function(n)
          setNames(
            quantile(n, probs = seq(bgnthr[1], bgnthr[2], bgnthr[3])),
            seq(bgnthr[1], bgnthr[2], bgnthr[3])
          ))

        BGNsaturation <- list(
          left = sapply(colnames(BGNQleft), function(Q) {
            list(sapply(BGNQleft[, Q], function(P)
              P < BGNPOW$left$BGN[, Q]))
          }),
          right = sapply(colnames(BGNQright), function(Q) {
            list(sapply(BGNQright[, Q], function(P)
              P < BGNPOW$right$BGN[, Q]))
          })
        )

        POWsaturation <- sapply(c("left", "right"), function(side) {
          list(sapply(colnames(BGNPOW[[side]][["POW"]]), function(Q) {
            list(sapply(powthreshold, function(P)
              P < BGNPOW[[side]][["POW"]][, Q]))
          }))
        })

        singsat <- do.call(rbind, sapply(c("left", "right"), function(side) {
          list(
            mapply(
              function(bgnthresh, powthresh) {
                sapply(1:length(BGNPOW$timeBins), function(i) {
                  sum(BGNsaturation[[side]][[paste0("BGN", i)]][, paste(bgnthresh)] |
                        POWsaturation[[side]][[paste0("POW", i)]][, paste(powthresh)]) / halfWl
                })
              },
              thresholdCombinations$bgnthreshold,
              thresholdCombinations$powthreshold
            )
          )

        }))

        binsUnique <- paste(rep(c('left', 'right'), each = nrow(singsat) /
                                  2), seq(nrow(singsat) / 2), sep = "_")

        DURATION <- rep(BGNPOW$timeBins, 2)
        SAMPRATE <- BGNPOW$sampRate

      } else if (BGNPOW$channel == "mono") {
        BGNQ <- apply(BGNPOW$mono$BGN, 2, function(n)
          setNames(
            quantile(n, probs = seq(bgnthr[1], bgnthr[2], bgnthr[3])),
            seq(bgnthr[1], bgnthr[2], bgnthr[3])
          ))

        BGNsaturation <- list(mono = sapply(colnames(BGNQ), function(Q) {
          list(sapply(BGNQ[, Q], function(P)
            P < BGNPOW$mono$BGN[, Q]))
        }))

        POWsaturation <- list(mono = sapply(colnames(BGNPOW$mono$POW), function(Q) {
          list(sapply(powthreshold, function(P)
            P < BGNPOW$mono$POW[, Q]))
        }))


        singsat <- mapply(
          function(bgnthresh, powthresh) {
            sapply(1:length(BGNPOW$timeBins), function(i) {
              sum(BGNsaturation$mono[[paste0("BGN", i)]][, paste(bgnthresh)] |
                    POWsaturation$mono[[paste0("POW", i)]][, paste(powthresh)]) / halfWl
            })
          },
          thresholdCombinations$bgnthreshold,
          thresholdCombinations$powthreshold
        )

        binsUnique <- paste("mono", seq(nrow(singsat)), sep = "_")

        DURATION <- BGNPOW$timeBins
        SAMPRATE <- BGNPOW$sampRate

      } else {
        realChannel <- c("left", "right")[c("left", "right") %in% names(BGNPOW)]

        BGNQ <- apply(BGNPOW[[realChannel]][["BGN"]], 2, function(n)
          setNames(
            quantile(n, probs = seq(bgnthr[1], bgnthr[2], bgnthr[3])),
            seq(bgnthr[1], bgnthr[2], bgnthr[3])
          ))

        BGNsaturation <- setNames(list(sapply(colnames(BGNQ), function(Q) {
          list(sapply(BGNQ[, Q], function(P)
            P < BGNPOW[[realChannel]][["BGN"]][, Q]))
        })), realChannel)

        POWsaturation <- setNames(list(sapply(colnames(BGNPOW[[realChannel]][['POW']]), function(Q) {
          list(sapply(powthreshold, function(P)
            P < BGNPOW[[realChannel]][['POW']][, Q]))
        })), realChannel)


        singsat <- mapply(
          function(bgnthresh, powthresh) {
            sapply(1:length(BGNPOW$timeBins), function(i) {
              sum(BGNsaturation[[realChannel]][[paste0("BGN", i)]][, paste(bgnthresh)] |
                    POWsaturation[[realChannel]][[paste0("POW", i)]][, paste(powthresh)]) / halfWl
            })
          },
          thresholdCombinations$bgnthreshold,
          thresholdCombinations$powthreshold
        )

        binsUnique <- paste(realChannel, seq(nrow(singsat)), sep = "_")

        DURATION <- BGNPOW$timeBins
        SAMPRATE <- BGNPOW$sampRate

      }

      message(
        "\r(",
        basename(soundfile),
        ") ",
        match(soundfile, soundfiles),
        " out of ",
        length(soundfiles),
        " recordinds concluded!",
        sep = ""
      )

      gc()

      list(
        SAT = singsat,
        DUR = DURATION,
        SMP = SAMPRATE,
        BIN = binsUnique,
        NAME = soundfile
      )

    }

    if (!is.null(backup) && fiveSteps %% 5 == 0) {
      saveRDS(SATdf, file = paste0(backup, "/SATBACKUP.RData"))
    }

  }

  if (!is.null(backup)) {
    SATdf["ogARGS"] <- NULL
    file.remove(paste0(backup, "/SATBACKUP.RData"))
  }

  which.error <- sapply(SATdf, function(x)
    is(x, "error"))
  ERRORS <- SATdf[which.error]
  DURATIONS <- c(unlist(sapply(SATdf[!which.error], function(x)
    x[["DUR"]])))
  SAMPRATES <- c(unlist(sapply(SATdf[!which.error], function(x)
    rep(x[["SMP"]], length(
      x[["BIN"]]
    )))))
  PATHS <- c(unlist(sapply(SATdf[!which.error], function(x)
    rep(x[["NAME"]], length(
      x[["BIN"]]
    )))))
  BINS <- c(unlist(sapply(SATdf[!which.error], function(x)
    x[["BIN"]])))
  SATdf <- do.call(rbind, lapply(SATdf[!which.error], function(x)
    x[["SAT"]]))

  colnames(SATdf) <- combinations

  normal <- apply(SATdf, 2, function(Q) {
    if (length(unique(Q)) != 1) {
      do.call(normality, list(Q))$statistic
    } else {
      NA
    }

  })

  if (normality %in% c("sf.test", "shapiro.test")) {
    thresholds <- unlist(strsplit(names(which.max(normal)), split = "/"))
    normOUT <- max(normal, na.rm = TRUE)
  } else {
    thresholds <- unlist(strsplit(names(which.min(normal)), split = "/"))
    normOUT <- min(normal, na.rm = TRUE)
  }

  normname <- switch(
    normality,
    "shapiro.test" = "Shapiro-Wilk",
    "sf.test" = "Shapiro-Francia",
    "ad.test" = "Anderson-Darling",
    "cvm.test" = "Cram\u00e9r-von Mises",
    "lillie.test" = "Lilliefors",
    "pearson.test" = "Pearson chi-square"
  )
  normstat <- switch(
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
    normOUT,
    "\n ",
    sep = ""
  )

  export <- list(
    powthresh = numeric(0),
    bgntresh = numeric(0),
    normality = numeric(0),
    values = data.frame(),
    errors = data.frame()
  )

  export["powthresh"] <- as.numeric(thresholds[1])
  export["bgntresh"] <- as.numeric(thresholds[2]) * 100
  export["normality"] <- normOUT
  export[["values"]] <- data.frame(
    PATH = dirname(PATHS),
    AUDIO = basename(PATHS),
    BIN = BINS,
    DURATION = DURATIONS,
    SAMPRATE = SAMPRATES,
    SAT = SATdf[, which(normOUT == normal)],
    row.names = NULL
  )
  export[["errors"]] <- data.frame(file = soundfiles[which.error], do.call(rbind, ERRORS))

  return(export)

}
