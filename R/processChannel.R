processChannel <- function(channelData,
                            channel,
                            timeBin,
                            wl,
                            overlap,
                            dbThreshold,
                            window,
                            histbreaks) {

  samp.rate <- channelData@samp.rate

  allSamples <- if(is.null(timeBin)) {
    data.frame(b = 1, e = length(channelData))
  } else {
    getSampleBins(length(channelData), samp.rate, timeBin)
  }

  frameBin <- nrow(allSamples)

  channelData <- switch(
    channel,
    "stereo" = list("left" = channelData@left, "right" = channelData@right),
    "mono" = list(mono = channelData@left),
    setNames(list(slot(channelData, channel)), channel)
  )

  BGNexp <- list(values = lapply(channelData, function(x) {

    offset <- x - mean(x)

    tempHolder <- apply(allSamples, 1, function(y) {
      list(
        signal::specgram(
          x = offset[y[1]:y[2]],
          n = wl,
          Fs = samp.rate,
          window = window,
          overlap = overlap
        )$S
      )
    })

    BGNPOWdf <- data.frame(do.call(cbind, lapply(lapply(tempHolder, function(singleBin) {
      spectS <- abs(singleBin[[1]])

      spectS <- 10 * log10(spectS / max(spectS))

      spectS[spectS < dbThreshold] <- dbThreshold

      apply(spectS, 1, function(z) {
        dbMax <- max(z)
        dbMin <- min(z)

        num_bins <- ifelse(is.numeric(histbreaks), histbreaks, eval(parse(
          text = paste0("nclass.", histbreaks, "(z)")
        )))

        modal_intensity <- dbMin + ((which.max(tabulate(
          findInterval(
            x = z,
            vec = seq(dbMin, dbMax, length.out = num_bins)
          )
        ))) * 2 * IQR(z) / length(z)^(1 / 3))

        c(BGN = modal_intensity, POW = dbMax - modal_intensity)
      })

    }), function(x)
      data.frame(t(x)))))

    colnames(BGNPOWdf) <- paste0(rep(c("BGN", "POW"), frameBin), rep(1:frameBin, each = 2))

    BGN <- data.frame(BGNPOWdf[, grepl("BGN", colnames(BGNPOWdf)), drop = FALSE])
    POW <- data.frame(BGNPOWdf[, grepl("POW", colnames(BGNPOWdf)), drop = FALSE])

    return(list(BGN = BGN, POW = POW))

  }))

  BGNexp[["timeBins"]] <- setNames(round((allSamples$e - allSamples$b) / samp.rate),
                                     paste0("BIN", seq(frameBin)))
  BGNexp[["sampRate"]] <- samp.rate
  BGNexp[["channel"]] <- channel

  return(BGNexp)

}
