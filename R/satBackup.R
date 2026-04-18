#' @title Backup for Ruido's functions
#'
#' @param backup path to the `.RData` file create by the backup of soundSat, soundMat or multActivity
#'
#' @description
#' This function offers a way to continue an unfinished process of the [soundSat()], [soundMat()] or [multActivity()] functions through a backup file.
#' Arguments can't be inputted nor changed since the function will automatically load them from the `.RData` file. However you may manually change them by editing the file (not recommended).
#'
#' @returns
#' This functions returns the same output of [soundSat()], [soundMat()] or [multActivity()]
#'
#' @export
#' @importFrom stats window
#'
#' @examples
#' \dontrun{
#' # It's impossible to demonstrate this function's intended use due to it's nature
#' # However, here is how this function is used:
#' ## This example will load an entire day of audios to your computer, so beware.
#'
#' ### Downloading audiofiles from public Zenodo library
#' dir <- paste(tempdir(), "forExample", sep = "/")
#' dir.create(dir)
#' recName <- paste0("GAL24576_20250401_", sprintf("%06d", seq(0, 230000, by = 10000)),".wav")
#' recDir <- paste(dir, recName, sep = "/")
#'
#' for(rec in recName) {
#'   print(rec)
#'   url <- paste0("https://zenodo.org/records/17575795/files/", rec, "?download=1")
#'   download.file(url, destfile = paste(dir, rec, sep = "/"), mode = "wb")
#' }
#'
#' sat <- soundSat(dir, backup = dir)
#'
#' # Now pretend the process was interrupted (manually/your R crashed/your computer turned off)
#' # We get the backup file
#'
#' list.files(dir)
#' backupDir <- paste(dir, "SATBACKUP.RData", sep = "/")
#'
#' # To recall the backup you simply:
#'
#' satB <- satBackup(backupDir)
#'
#' head(satB$values)
#'
#' unlink(dir, recursive = TRUE)
#' }
satBackup <- function(backup) {
  SATdf <- readRDS(backup)

  list2env(SATdf$ogARGS, envir = environment())

  soundfiles <- list.files(od, full.names = TRUE, recursive = TRUE)
  soundfiles <- soundfiles[tolower(tools::file_ext(soundfiles)) %in% c("mp3", "wav")]

  if (type != "multActivity") {
    powthreshold <- seq(powthr[1], powthr[2], powthr[3])
    names(powthreshold) <- powthreshold
    bgnthreshold <- seq(bgnthr[1], bgnthr[2], bgnthr[3])

    thresholdCombinations <- setNames(expand.grid(powthreshold, bgnthreshold),
                                      c("powthreshold", "bgnthreshold"))

    combinations <- paste(thresholdCombinations[, 1], thresholdCombinations[, 2], sep = "/")
  }

  halfWl <- wl / 2

  if (concluded == nFiles) {
    message("All files have already been processed!")

  } else {
    for (soundfile in concluded:nFiles) {
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

        saveRDS(SATdf, file = backup)
      }

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
    if (x@channel == "stereo") {
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
    if (type != "multActivity") {
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
      BGNQ <- quantile(unlist(BGN), probs = bgnthr) |>
        setNames(bgnthr)

      SATmat <- BGN > BGNQ |
        POW > powthr

    }

  } else {
    if (type != "multActivity") {
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

    } else {
      SATmat <- sapply(1:ncol(BGN), function(t) {
        BGN[, t] > quantile(BGN[, t], bgnthr) |
          POW[, t] > powthr
      })

    }

  }

  if (type != "multActivity") {
    colnames(SATmat) <- combinations

  }

  if (type == "soundSat") {
    normal <- apply(SATmat, 2, function(Q) {
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

    SATinfo$SAT <- SATmat[, which(normal == normOUT)]

    export <- list(
      powthresh = numeric(0),
      bgnthresh = numeric(0),
      normality = list(),
      values = data.frame(),
      errors = list()
    )

    export["powthresh"] <- as.numeric(thresholds[1])
    export["bgnthresh"] <- as.numeric(thresholds[2]) * 100
    export[["normality"]]["test"] <- normality
    export[["normality"]]["statistic"] <- normOUT
    export[["values"]] <- SATinfo
    export[["errors"]] <- ERRORS

  } else if (type == "soundMat") {
    export <- list(info = data.frame(),
                   values = matrix(),
                   errors = list())

    export[["info"]] <- SATinfo
    export[["values"]] <- SATmat
    export[["errors"]] <- ERRORS

    return(export)

  } else if (type == "multActivity") {

    export <- list(
      powthresh = numeric(0),
      bgnthresh = numeric(0),
      info = data.frame(),
      values = matrix(),
      errors = list()
    )

    export["powthresh"] <- powthr
    export["bgnthresh"] <- bgnthr * 100
    export[["info"]] <- SATinfo
    export[["values"]] <- SATmat * 1
    export[["errors"]] <- ERRORS

  }

  if (!is.null(backup)) {
    file.remove(backup)
  }

  return(export)

}
