#' Backup for Soundscape Saturation Index
#'
#' @param backupPath path you set in your `backup` in the `soundSat()` function. Audiofiles already finished will be drawn from this path
#' @param od path or paths containing your original audiofiles
#'
#' @description
#' This function is a way to continue an unfinished process of the `soundSat()` or `soundMat()` functions through a backup file.
#' Arguments can't be inputted nor changed since the function will automatically load them from the original `soundSat()` run.
#'
#' @returns
#' A list containing five objects. The first and second objects (powthresh and bgnthresh) are the threshold values that yielded the most normal distribution of saturation values using the normality test set by the user. The third (normality) contains the statitics values of the normality test that yielded the most normal distribution. The fourth object (values) contains a data.frame with the the values of saturation for each bin of each recording and the size of the bin in seconds. The fifth contains a data.frame with errors that occurred with specific files during the function.
#'
#' @export
#' @importFrom stats window
#'
#' @examples
#' \dontrun{
#' # It's impossible to create a functioning example since you would have to manually stop the process
#' # However, here is how this function is used:
#' ## This example will load an entire day of audios to your computer, so beware.
#'
#' ### Downloading audiofiles from public Zenodo library
#' dir <- tempdir()
#' recName <- paste0("GAL24576_20250401_", sprintf("%06d", seq(0, 230000, by = 10000)),".wav")
#' recDir <- paste(dir, recName, sep = "/")
#'
#' for(rec in recName) {
#'   print(rec)
#'   url <- paste0("https://zenodo.org/records/17575795/files/", rec, "?download=1")
#'   download.file(url, destfile = paste(dir, rec, sep = "/"), mode = "wb")
#' }
#'
#' sat <- soundSat(dir, backup = dir, wl = 256)
#'
#' # Now pretend the process was interrupted (manually/your R crashed/your computer turned off)
#' # To recall the backup you simply:
#'
#' satB <- satBackup(dir, dir)
#'
#' unlink(recDir)
#' }
satBackup <- function(backupPath, od) {
  backfile <- paste0(backupPath, "/SATBACKUP.RData")
  SATdf <- readRDS(backfile)
  originalfiles <- list.files(od, full.names = TRUE, recursive = TRUE)
  originalfiles <- originalfiles[tools::file_ext(originalfiles) %in% c("mp3", "wav")]

  remainingfiles <- originalfiles[!(basename(originalfiles) %in% basename(names(SATdf)))]

  list2env(SATdf$ogARGS, envir = environment())

  powthreshold <- seq(powthr[1], powthr[2], powthr[3])
  names(powthreshold) <- powthreshold
  bgnthreshold <- seq(bgnthr[1], bgnthr[2], bgnthr[3])

  thresholdCombinations <- setNames(expand.grid(powthreshold, bgnthreshold),
                                    c("powthreshold", "bgnthreshold"))

  combinations <- paste(thresholdCombinations[, 1], thresholdCombinations[, 2], sep = "/")

  halfWl <- wl / 2

  if (length(remainingfiles) == 0) {
    message("All files have already been processed!")

    SATdf <- lapply(SATdf, readRDS)

  } else {
    fiveSteps <- 1

    for (soundfile in remainingfiles) {
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
        message("\n",
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
          match(soundfile, remainingfiles),
          " out of ",
          length(remainingfiles),
          " remaining recordinds concluded!",
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

      saveRDS(SATdf, file = paste0(backupPath, "/SATBACKUP.RData"))

    }
  }

  SATdf["ogARGS"] <- NULL
  file.remove(paste0(backupPath, "/SATBACKUP.RData"))

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

  if(type == "soundSat") {

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
      SAT = SATdf[, which(normOUT == normal)]
    )
    export[["errors"]] <- data.frame(file = remainingfiles[which.error], do.call(rbind, ERRORS))

  } else {

    export <- list(
      powthreshs = numeric(0),
      bgntreshs = numeric(0),
      info = data.frame(),
      matrix = data.frame(),
      errors = data.frame()
    )

    export[["powthreshs"]] <- powthreshold
    export[["bgntreshs"]] <- bgnthreshold
    export[["info"]] <- data.frame(NAME = basename(PATHS), BINS, SAMPRATES)
    export[["matrix"]] <- SATdf
    export[["errors"]] <- data.frame(file = remainingfiles[which.error], do.call(rbind, ERRORS))

  }

  return(export)

}
