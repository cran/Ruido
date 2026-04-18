#' @title Soundscape Saturation Matrix
#'
#' @description Get the Soundscape Saturation matrix with all threshold combinations instead of the combination with the most normal distribution.
#'
#' @param soundpath single or multiple directories to your audio files
#' @param channel channel where the saturation values will be extract from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`.
##' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Set to `NULL` to leave db values unrestricted Defaults to `-90`, as set by Towsey 2017
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `signal::hanning(wl)` to use hanning instead. Defaults to `signal::hammning(wl)`
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges`", `"scott"` and `100`. Defaults to `"FD"`.
#' <br>Can also be set to any numerical value to limit or increase the amount of breaks.
#' @param DCfix if the DC offset should be removed before the metrics are calculated. Defaults to `TRUE`
#' @param powthr numeric vector of length three containing the the range of thresholds used to evaluate the Soundscape Power of the  Activity Matrix (in dB). The values correspond to the minimum threshold, maximum threshold and step size respectively.
#' <br> Defaults to `c(5, 20, 1)`, which evaluates thresholds from 5 dB to 20 dB in increments of 1 dB
#' @param bgnthr numeric vector of length three containing the the range of thresholds used to evaluate the Background Noise of the  Activity Matrix (in %). The values correspond to the minimum threshold, maximum threshold and step size respectively.
#' <br> Defaults to `c(0.5, 0.9, 0.05)`, which evaluates thresholds from 50% to 90% in increments of 5%
#' @param beta how BGN thresholds are calculated. If `TRUE`, BGN thresholds are calculated using all recordings combined. If FALSE, BGN thresholds are calculated separately for each recording. Defaults to `TRUE`
#' @param backup path to save the backup. Defaults to `NULL`
#'
#' @returns A list containing three objects. The first (info) contains the following variables from every audio file: PATH, AUDIO, CHANNEL, DURATION, BIN, SAMPRATE. The second (values) contains saturation values from all possible threshold combinations. The third (errors) contains the error messages and the paths to the files that returned an error during processing.
#'
#' @details Check [soundSat()] to see how the indices are calculated.
#'
#' If `backup` is set to a valid directory, a file named `"SATBACKUP.RData"` is saved after every batch of five processed files. If the function execution is interrupted (e.g., manual termination, an R session crash, or a system shutdown), this backup file can be passed to `satBackup()` (e.g., as `~path/SATBACKUP.RData`) to resume the original process. Once a backup is created, all arguments and file paths must remain unchanged, unless they are manually modified within the `.RData` object.
#'
#' @seealso [soundSat()] to get only the threshold with the most normal distribution and [multActivity()] to generate only activity matrices. Also, check [satBackup()] if you are working with larger datasets and want some safety.
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
#' \donttest{
#' oldpar <- par(no.readonly = TRUE)
#' ### Downloading audiofiles from public Zenodo library
#' dir <- paste(tempdir(), "forExample", sep = "/")
#' dir.create(dir)
#' recName <- paste0("GAL24576_20250401_", sprintf("%06d", seq(0, 200000, by = 50000)), ".wav")
#' recDir <- paste(dir, recName, sep = "/")
#'
#' for (rec in recName) {
#'   print(rec)
#'   url <- paste0("https://zenodo.org/records/17575795/files/",
#'                 rec,
#'                 "?download=1")
#'   download.file(url, destfile = paste(dir, rec, sep = "/"), mode = "wb")
#' }
#'
#' ### Running the function
#' sat <- soundMat(dir)
#'
#' ### Plotting results
#' sides <- sat$info$CHANNEL
#'
#' thresholds <- colnames(sat$values)
#' split <- strsplit(thresholds, "/")
#'
#' shapNorm <- apply(sat$values, 2, function(x)
#'
#'   if (var(x) == 0) {
#'     0
#'   } else {
#'     shapiro.test(x)$statistic
#'   })
#'
#' shapPos <- which.max(shapNorm)
#'
#' par(mfrow = c(3, 2))
#'
#' plot(
#'   sat$values[sides == "left", 1],
#'   main = paste0("POW = ", split[[1]][1], "dB | BGN = ", split[[1]][2], "%"),
#'   type = "b",
#'   ylim = c(0,1),
#'   xlab = "Time Index", ylab = "Soundsacpe Saturation (%)", col = "goldenrod"
#' )
#' points(sat$values[sides == "right", 1], col = "maroon", type = "b")
#'
#' hist(sat$values[,1], main = paste("Histogram of POW = ", split[[1]][1],
#' "dB | BGN = ", split[[1]][2], "%"), xlab = "Soundscape Saturation (%)")
#'
#' plot(
#' sat$values[sides == "left", 144],
#' main = paste0("POW = ", split[[144]][1], "dB | BGN = ", split[[144]][2], "%"),
#' type = "b",
#' ylim = c(0,1),
#' xlab = "Time Index", ylab = "Soundsacpe Saturation (%)", col = "goldenrod"
#' )
#' points(sat$values[sides == "right", 144], col = "maroon", type = "b")
#'
#' hist(sat$values[,144], main = paste("Histogram of POW = ", split[[144]][1],
#' "dB | BGN = ", split[[144]][2], "%"), xlab = "Soundscape Saturation (%)")
#'
#' plot(
#'   sat$values[sides == "left", shapPos],
#'   main = paste0(
#'     "POW = ",
#'     split[[shapPos]][1],
#'     "dB | BGN = ",
#'     split[[shapPos]][2],
#'     "%",
#'     "\nshapiro.test. statistic (W): ",
#'     which.max(shapNorm)
#'   ),
#'   type = "b",
#'   ylim = c(0,1),
#'   xlab = "Time Index", ylab = "Soundsacpe Saturation (%)", col = "goldenrod"
#' )
#' points(sat$values[sides == "right", shapPos], col = "maroon", type = "b")
#' hist(sat$values[,shapPos], main = paste("Histogram of POW = ",
#' split[[shapPos]][1], "dB | BGN = ", split[[shapPos]][2], "%"),
#' xlab = "Soundscape Saturation (%)")
#'
#' unlink(dir, recursive = TRUE)
#' par(oldpar)
#' }
soundMat <- function(soundpath,
                     channel = "stereo",
                     timeBin = 60,
                     dbThreshold = -90,
                     targetSampRate = NULL,
                     wl = 512,
                     window = signal::hamming(wl),
                     overlap = ceiling(length(window) / 2),
                     histbreaks = "FD",
                     DCfix = TRUE,
                     powthr = c(5, 20, 1),
                     bgnthr = c(0.5, 0.9, 0.05),
                     beta = TRUE,
                     backup = NULL) {

  argHandler(FUN = "soundMat", soundpath, channel, timeBin, dbThreshold, targetSampRate, wl,
             window, overlap, histbreaks, DCfix, powthr, bgnthr, beta, backup)

  powthreshold <- seq(powthr[1], powthr[2], powthr[3])
  names(powthreshold) <- powthreshold
  bgnthreshold <- seq(bgnthr[1], bgnthr[2], bgnthr[3])

  soundfiles <- list.files(soundpath, full.names = TRUE, recursive = TRUE)
  soundfiles <- soundfiles[tolower(tools::file_ext(soundfiles)) %in% c("mp3", "wav")]

  if (length(soundfiles) < 3)
    stop("please provide at least 3 recordings!")

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
      DCfix = DCfix,
      powthr = powthr,
      bgnthr = bgnthr,
      beta = beta,
      type = "soundMat",
      od = soundpath,
      nFiles = length(soundfiles),
      concluded = 0
    )
  }

  nFiles <- length(soundfiles)
  SATdf[["indexes"]] <- vector("list", nFiles)

  for (soundfile in 1:nFiles) {
    gc()

    sPath <- soundfiles[[soundfile]]

    SATdf[["indexes"]][[soundfile]] <- tryCatch(
      bgNoise.(
        sPath,
        timeBin = timeBin,
        targetSampRate = targetSampRate,
        window = window,
        overlap = overlap,
        channel = channel,
        dbThreshold = dbThreshold,
        wl = wl,
        histbreaks = histbreaks,
        DCfix = DCfix
      ),
      error = function(e)
        e
    )

    SATdf[["indexes"]][[soundfile]]@path <- sPath

    message(
      "\r(",
      basename(soundfiles[soundfile]),
      ") ",
      match(soundfiles[soundfile], soundfiles),
      " out of ",
      length(soundfiles),
      " recordings concluded!",
      sep = ""
    )

    if (!is.null(backup) && soundfile %% 5 == 1) {
      SATdf$ogARGS$concluded <- soundfile

      saveRDS(SATdf, file = paste0(backup, "/SATBACKUP.RData"))

    }

  }

  whichError <- sapply(SATdf[["indexes"]], function(x) {
    is(x, "error")
  })

  ERRORS <- SATdf$indexes[whichError]
  indexes <- SATdf$indexes[!whichError]

  BGN <- do.call(cbind, sapply(indexes, function(x) {
    if (x@channel == "stereo") {
      cbind(x@values$left$BGN, x@values$right$BGN)
    } else {
      x@values[[x@channel]]$BGN
    }
  }))

  POW <- do.call(cbind, sapply(indexes, function(x) {
    if (x@channel == "mono") {
      cbind(x@values$left$POW, x@values$right$POW)
    } else {
      x@values[[x@channel]]$POW
    }
  }))

  INFO <- lapply(indexes, function(x) {
    nBins <- length(x@timeBins)
    if (x@channel == "stereo") {
      list(
        rep(x@timeBins, each = 2),
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

  paths <- unlist(sapply(indexes, function(x) {
    if (x@channel == "stereo") {
      rep(x@path, length(x@timeBins) * 2)
    } else {
      rep(x@path, length(x@timeBins))
    }
  }))

  SATinfo <- data.frame(
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

  dimBGN <- dim(BGN)

  if (beta) {
    BGNQ <- quantile(unlist(BGN), probs = seq(bgnthr[1], bgnthr[2], bgnthr[3])) |>
      setNames(bgnthreshold)

    SATmat <- mapply(
      function(bgnthresh, powthresh) {
        sapply(1:ncol(BGN), function(t) {
          sum(BGN[, t] > BGNQ[names(BGNQ) == bgnthresh] |
                POW[, t] > powthresh) / halfWl

        })

      },
      thresholdCombinations$bgnthreshold,
      thresholdCombinations$powthreshold
    )

  } else {
    SATmat <- mapply(
      function(bgnthresh, powthresh) {
        sapply(1:ncol(BGN), function(t) {
          sum(BGN[, t] > quantile(BGN[, t], bgnthresh) |
                POW[, t] > powthresh) / halfWl

        })

      },
      thresholdCombinations$bgnthreshold,
      thresholdCombinations$powthreshold
    )

  }

  colnames(SATmat) <- combinations

  if (!is.null(backup)) {
    SATdf["ogARGS"] <- NULL
    file.remove(paste0(backup, "/SATBACKUP.RData"))
  }

  export <- list(info = data.frame(),
                 values = matrix(),
                 errors = list())

  export[["info"]] <- SATinfo
  export[["values"]] <- SATmat
  export[["errors"]] <- ERRORS

  return(export)

}
