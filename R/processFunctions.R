processChannel.BGN = function(channelData,
                              samp.rate,
                              channel,
                              timeBin,
                              wl,
                              overlap,
                              dbThreshold,
                              dbRef,
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

  hBreak = hBreaks.MAT(histbreaks)

  if (DCfix) {
    channelData = channelData - matrixStats::rowMeans2(channelData)
  }

  bVec = allSamples$b
  eVec = allSamples$e

  noiseOBJ@values = lapply(1:nrow(channelData), function(x) {
    BGNPOWdf = lapply(seq_along(bVec), function(i) {
      tempHolder = abs(.spect(
        x = channelData[x,bVec[i]:eVec[i]],
        n = wl,
        window = window,
        overlap = overlap
      ))
      singleBin = dbRef * log10(tempHolder / max(tempHolder))
      if (!is.null(dbThreshold)) {
        singleBin[singleBin < dbThreshold] = dbThreshold
      }

      dbDif = matrixStats::rowRanges(singleBin)
      dbMax = dbDif[, 2]
      dbMin = dbDif[, 1]

      numBins = hBreak(singleBin)
      binWidth = (dbMax - dbMin) / numBins
      bins = floor((singleBin - dbMin) / binWidth) + 1L

      modalBin = vapply(seq_len(wl / 2), function(i) {
        which.max(tabulate(bins[i, ], nbins = numBins[i]))
      }, integer(1))
      modalIntensity = dbMin + (modalBin - 0.5) * binWidth
      rbind(BGN = modalIntensity, POW = dbMax - modalIntensity)
    })

    return(list(
      BGN = as.data.frame(do.call(
        cbind, lapply(BGNPOWdf, function(df)
          df[1, ])
      )) |>
        setNames(paste0("BGN", 1:frameBin)),
      POW = as.data.frame(do.call(
        cbind, lapply(BGNPOWdf, function(df)
          df[2, ])
      )) |>
        setNames(paste0("POW", 1:frameBin))
    ))

  })

  names(noiseOBJ@values) = switch(
      channel,
      "stereo" = c("left", "right"),
      "mono"   = "mono",
      "left"   = "left",
      "right"  = "right"
    )

  noiseOBJ@index = c("BGN", "POW")
  noiseOBJ@timeBins = setNames(round((allSamples$e - allSamples$b) / samp.rate), paste0("BIN", seq(frameBin)))
  noiseOBJ@sampRate = samp.rate
  noiseOBJ@channel = channel

  return(noiseOBJ)

}

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

  bVec = allSamples$b
  eVec = allSamples$e

  noiseOBJ@values = lapply(1:nrow(channelData), function(x) {
    ACIdf = lapply(seq_along(bVec), function(i) {
      tempHolder = abs(.spect(
        x = channelData[x,bVec[i]:eVec[i]],
        n = wl,
        window = window,
        overlap = overlap
      ))

      specDim = dim(tempHolder)

      jump = if (is.null(j)) {
        ncol(tempHolder)
      } else {
        max(1L, as.integer(j * samp.rate / (wl - overlap)))
      }

      forJ = ceiling(specDim[2] / jump)

      timeMin = (1:forJ - 1) * jump + 1
      timeMax = pmin((1:forJ * jump), specDim[2])

      vapply(seq(specDim[1]), function(s) {
        matrixStats::sum2(vapply(seq(timeMin), function(t) {
          matrixStats::sum2(abs(matrixStats::diff2(tempHolder[s, timeMin[t]:timeMax[t]]))) /
            matrixStats::sum2(tempHolder[s, timeMin[t]:timeMax[t]])
        }, numeric(1)))
      }, numeric(1))
    })

    return(list(ACI = data.frame(ACIdf) |>
                  setNames(paste0(
                    "ACI", 1:frameBin
                  ))))

  })

  names(noiseOBJ@values) = switch(
    channel,
    "stereo" = c("left", "right"),
    "mono"   = "mono",
    "left"   = "left",
    "right"  = "right"
  )

  noiseOBJ@index = "ACI"
  noiseOBJ@timeBins = setNames(round((allSamples$e - allSamples$b) / samp.rate), paste0("BIN", seq(frameBin)))
  noiseOBJ@sampRate = samp.rate
  noiseOBJ@channel = channel

  return(noiseOBJ)

}


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

  bVec = allSamples$b
  eVec = allSamples$e

  noiseOBJ@values = lapply(1:nrow(channelData), function(x) {
    ENTdf = lapply(seq_along(bVec), function(i) {
      tempHolder = .spect(
        x = channelData[x,bVec[i]:eVec[i]],
        n = wl,
        window = window,
        overlap = overlap
      )

      spectS = abs(tempHolder)
      amp2 = spectS ** 2

      rowTot = matrixStats::rowSums2(amp2)
      rowTot[rowTot == 0] = 1 ## Safeguard!
      pmf = amp2 / rowTot

      term = ifelse(pmf > 0, pmf * log2(pmf), 0)
      H = -matrixStats::rowSums2(term) / log2(dim(amp2)[2])

      return(1 - H)
    })

    return(list(ENT = data.frame(ENTdf) |>
                  setNames(paste0(
                    "ENT", 1:frameBin
                  ))))

  })

  names(noiseOBJ@values) = switch(
    channel,
    "stereo" = c("left", "right"),
    "mono"   = "mono",
    "left"   = "left",
    "right"  = "right"
  )

  noiseOBJ@index = "ENT"
  noiseOBJ@timeBins = setNames(round((allSamples$e - allSamples$b) / samp.rate), paste0("BIN", seq(frameBin)))
  noiseOBJ@sampRate = samp.rate
  noiseOBJ@channel = channel

  return(noiseOBJ)

}

processChannel.EVN = function(channelData,
                              samp.rate,
                              channel,
                              timeBin,
                              wl,
                              overlap,
                              dbThreshold,
                              dbRef,
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

  hBreak = hBreaks.MAT(histbreaks)

  if (DCfix) {
    channelData = channelData - matrixStats::rowMeans2(channelData)
  }

  bVec = allSamples$b
  eVec = allSamples$e

  noiseOBJ@values = lapply(1:nrow(channelData), function(x) {
    EVNdf = lapply(seq_along(bVec), function(i) {
      tempHolder = abs(.spect(
        x = channelData[x,bVec[i]:eVec[i]],
        n = wl,
        window = window,
        overlap = overlap
      ))
      singleBin = dbRef * log10(tempHolder / max(tempHolder))
      if (!is.null(dbThreshold)) {
        singleBin[singleBin < dbThreshold] = dbThreshold
      }

      dbDif = matrixStats::rowRanges(singleBin)
      dbMax = dbDif[, 2]
      dbMin = dbDif[, 1]

      numBins = hBreak(singleBin)
      binWidth = (dbMax - dbMin) / numBins
      bins = floor((singleBin - dbMin) / binWidth) + 1L

      modalBin = vapply(seq_len(wl / 2), function(i) {
        which.max(tabulate(bins[i, ], nbins = numBins[i]))
      }, integer(1))
      modalIntensity = dbMin + (modalBin - 0.5) * binWidth
      noiseReduced = (singleBin - modalIntensity) >= 3

      matrixStats::rowSums2(noiseReduced[, -1, drop = FALSE] &
                              !noiseReduced[, -ncol(noiseReduced), drop = FALSE]) / ((ncol(singleBin) - 1) * (wl - overlap) / samp.rate)
    })

    return(list(EVN = data.frame(EVNdf) |>
                  setNames(paste0(
                    "EVN", 1:frameBin
                  ))))

  })

  names(noiseOBJ@values) = switch(
    channel,
    "stereo" = c("left", "right"),
    "mono"   = "mono",
    "left"   = "left",
    "right"  = "right"
  )

  noiseOBJ@index = c("EVN")
  noiseOBJ@timeBins = setNames(round((allSamples$e - allSamples$b) / samp.rate), paste0("BIN", seq(frameBin)))

  noiseOBJ@sampRate = samp.rate
  noiseOBJ@channel = channel

  return(noiseOBJ)

}
