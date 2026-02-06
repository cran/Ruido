#' @title Multiple Acoustic Activity Matrix
#'
#' @description Calculate the Acoustic Activity Matrix used in the the calculation of Soundscape Saturation using Burivalova 2018 methodology for a set of recordings
#'
#' @param soundpath single or multiple directories to your audio files
#' @param channel channel where the saturation values will be extract from. Available channels are: `"stereo"`, `"mono"`, `"left"` or `"right"`. Defaults to `"stereo"`.
#' @param timeBin size (in seconds) of the time bin. Set to `NULL` to use the entire audio as a single bin. Defaults to `60`
#' @param dbThreshold minimum allowed value of dB for the spectrograms. Defaults to `-90`, as set by Towsey 2017
#' @param targetSampRate desired sample rate of the audios.  This argument is only used to down sample the audio. If `NULL`, then audio's sample rate remains the same. Defaults to `NULL`
#' @param wl window length of the spectrogram. Defaults to `512`
#' @param window window used to smooth the spectrogram. Switch to `signal::hanning(wl)` to use hanning instead. Defaults to `signal::hammning(wl)`
#' @param overlap overlap between the spectrogram windows. Defaults to `wl/2` (half the window length)
#' @param histbreaks breaks used to calculate Background Noise. Available breaks are: `"FD"`, `"Sturges"`, `"scott"` or any numeric value (foe example = `100`). Defaults to `"FD"`
#' @param powthr single numeric value to calculate the activity matrix for soundscape power (in dB). Detauls to `10`
#' @param bgnthr single numeric value to calculate the activity matrix for background noise (in %). Detauls to `0.8`
#' @param beta how BGN thresholds are calculated. If `TRUE`, BGN thresholds are calculated using all recordings combined. If FALSE, BGN thresholds are calculated separately for each recording. Defaults to `TRUE`
#' @param backup directory to save the backup. Defaults to `NULL`
#'
#' @returns A list containing five objects. The first and second objects (powthresh and bgnthresh) are the threshold values inputted as arguments into the function. The third (info) contains the following variables from every audio file: PATH, AUDIO, CHANNEL, DURATION, BIN, SAMPRATE.. The fourth object (values) contains a matrix with the the values of activity for each bin of each recording and the size of the bin in seconds. The fifth contains a list with errors that occurred with specific files during the function.
#'
#' @details We generate an activity matrix using Burivalova 2018 methodology. For each time bin of the recording we apply the following formula:
#'
#'\deqn{a_{mf} = 1\  if (BGN_{mf} > \theta_{1})\  or\  (POW_{mf} > \theta_{2});\  otherwise,\  a_{mf} = 0,}
#'
#'Where \eqn{\theta_{1}} is the threshold of BGN values and \eqn{\theta_{2}} is a threshold of dB values. 1 = active and 0 = inactive.
#'
#' If `backup` is set to a valid directory, a file named `"SATBACKUP.RData"` is saved after every batch of five processed files. If the function execution is interrupted (e.g., manual termination, an R session crash, or a system shutdown), this backup file can be passed to `satBackup()` (e.g., `~path/SATBACKUP.RData`) to resume the original process. Once a backup is created, all arguments and file paths must remain unchanged, unless they are manually modified within the `.RData` object.
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
#'@examples
#' \donttest{
#' if (require("ggplot2") & require("patchwork")) {
#'   ### Generating an artificial audio for the example
#'   ## For this example we'll generate a sweep in a noisy soundscape
#'   library(ggplot2)
#'   library(patchwork)
#'
#'   ### Downloading audiofiles from public Zenodo library
#'   dir <- paste0(tempdir(), "/forExamples")
#'   dir.create(dir)
#'   recName <- paste0("GAL24576_20250401_", sprintf("%06d", seq(0, 200000, by = 50000)), ".wav")
#'   recDir <- paste(dir, recName, sep = "/")
#'
#'   for (rec in recName) {
#'     print(rec)
#'     url <- paste0("https://zenodo.org/records/17575795/files/",
#'                   rec,
#'                   "?download=1")
#'     download.file(url,
#'                   destfile = paste(dir, rec, sep = "/"),
#'                   mode = "wb")
#'   }
#'
#'   time <- sapply(strsplit(recName, "_"), function(x)
#'     paste(substr(x[3], 1, 2), substr(x[3], 3, 4), substr(x[3], 5, 6), sep = ":"))
#'   date <- sapply(strsplit(recName, "_"), function(x)
#'     paste(substr(x[2], 1, 4), substr(x[2], 5, 6), substr(x[2], 7, 8), sep = "-"))
#'
#'   dateTime <- as.POSIXct(paste(date, time))
#'
#'   timeLabels <- time[c(1, 7, 13, 19, 24)]
#'   timeBreaks <- as.character(dateTime[c(1, 7, 13, 19, 24)])
#'
#'   breaks <- round(c(1, cumsum(rep(256 / 6, 6))))
#'
#'   ### Running the function
#'   act <- multActivity(dir)
#'
#'   plotN <- 1
#'
#'   sDim <- dim(act$values)
#'
#'   sampRate <- act$info$SAMPRATE[1]
#'   kHz <- cumsum(c(0, rep(sampRate / 6, 6))) / 1000
#'
#'   plotList <- list()
#'
#'   for (cha in c("left", "right")) {
#'     actCurrent <- act$values[, act$info$CHANNEL == cha]
#'     actCurrentDF <- data.frame(
#'       TIME = as.character(rep(dateTime, each = sDim[1])),
#'       SPEC = rep(seq(sDim[1]), sDim[2]),
#'       VAL = factor(c(unlist(actCurrent)), levels = c(0, 1))
#'     )
#'
#'     plotList[[plotN]] <- ggplot(actCurrentDF, aes(x = TIME, y = SPEC, fill = VAL)) +
#'       geom_tile() +
#'       theme_classic() +
#'       scale_y_continuous(expand = c(NA, NA),
#'                          labels = kHz,
#'                          breaks = breaks) +
#'       scale_x_discrete(expand = c(0, 0),
#'                        labels = time) +
#'       scale_fill_manual(values = c("white", "black"),
#'                         labels = c("Inactive", "Active")) +
#'       guides(fill = guide_legend(title = "Acoustic Activity")) +
#'       labs(
#'         x = "Time of Day",
#'         y = "Frequency (kHz)",
#'         title = paste("Acoustic Activity in the", cha, "channel")
#'       )
#'
#'     plotN <- plotN + 1
#'
#'   }
#'
#'   plotList[[1]] + plotList[[2]] + plot_layout(guide = "collect")
#'
#'   unlink(recDir)
#'   unlink(dir)
#'
#'   }
#' }
multActivity <- function(soundpath,
                         channel = "stereo",
                         timeBin = 60,
                         dbThreshold = -90,
                         targetSampRate = NULL,
                         wl = 512,
                         window = signal::hamming(wl),
                         overlap = ceiling(length(window) / 2),
                         histbreaks = "FD",
                         powthr = 10,
                         bgnthr = 0.8,
                         beta = TRUE,
                         backup = NULL) {
  if (all(!dir.exists(soundpath)))
    stop("all provided soundpaths must be valid")

  if (!is.null(backup) && !dir.exists(backup))
    stop("please provide a valid folder for backup")

  if (!(channel %in% c("left", "right", "stereo", "mono")))
    stop("channel must be 'stereo', 'mono', 'left', or 'right'")

  if (!is.numeric(timeBin))
    stop("timeBin must be numeric")

  if (!is.numeric(dbThreshold))
    stop("timeBin must be numeric")

  if (!is.null(targetSampRate)) {
    if (!is.numeric(targetSampRate))
      stop("targetSampRate must be either NULL or a numeric value")
  }

  soundfiles <- list.files(soundpath, full.names = TRUE, recursive = TRUE)
  soundfiles <- soundfiles[tolower(tools::file_ext(soundfiles)) %in% c("mp3", "wav")]

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
      beta = beta,
      type = "multActivity",
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
      bgNoise(
        sPath,
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

    SATdf[["indexes"]][[soundfile]][["path"]] <- sPath

    message(
      "\r(",
      basename(soundfiles[soundfile]),
      ") ",
      match(soundfiles[soundfile], soundfiles),
      " out of ",
      length(soundfiles),
      " recordinds concluded!",
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
    if (x$channel == "stereo") {
      cbind(x$values$left$BGN, x$values$right$BGN)
    } else {
      x$values[[x$channel]]$BGN
    }
  }))

  POW <- do.call(cbind, sapply(indexes, function(x) {
    if (x$channel == "stereo") {
      cbind(x$values$left$POW, x$values$right$POW)
    } else {
      x$values[[x$channel]]$POW
    }
  }))

  INFO <- lapply(indexes, function(x) {
    nBins <- length(x$timeBins)
    if (x$channel == "stereo") {
      list(
        rep(x$timeBins, each = 2),
        rep(x$sampRate, length(x$timeBins) * 2),
        rep(1:length(x$timeBins), 2),
        rep(c("left", "right"), each = nBins)
      )
    } else {
      list(
        x$timeBins,
        rep(x$sampRate, length(x$timeBins)),
        1:length(x$timeBins),
        rep(x$channel, nBins)
      )
    }
  })

  paths <- unlist(sapply(indexes, function(x) {
    if (x$channel == "stereo") {
      rep(x$path, length(x$timeBins) * 2)
    } else {
      rep(x$path, length(x$timeBins))
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
    BGNQ <- quantile(unlist(BGN), probs = bgnthr) |>
      setNames(bgnthr)

    SATmat <- BGN > BGNQ |
      POW > powthr

  } else {
    SATmat <- sapply(1:ncol(BGN), function(t) {
      BGN[, t] > quantile(BGN[, t], bgnthr) |
        POW[, t] > powthr
    })
  }

  if (!is.null(backup)) {
    SATdf["ogARGS"] <- NULL
    file.remove(paste0(backup, "/SATBACKUP.RData"))
  }

  export <- list(
    powthresh = numeric(0),
    bgntresh = numeric(0),
    info = data.frame(),
    values = matrix(),
    errors = list()
  )

  export["powthresh"] <- powthr
  export["bgntresh"] <- bgnthr * 100
  export[["info"]] <- SATinfo
  export[["values"]] <- SATmat * 1
  export[["errors"]] <- ERRORS

  return(export)

}
