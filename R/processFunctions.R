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

  channelData = switch(
    channel,
    "stereo" = list(left  = channelData[1, ], right = channelData[2, ]),
    "mono"   = list(mono  = channelData[1, ]),
    "left"   = list(left  = channelData[1, ]),
    "right"  = list(right = channelData[2, ])
  )

  if (DCfix) {
    channelData = lapply(channelData, function(x)
      x - mean(x))
  }

  noiseOBJ@values = lapply(channelData, function(x) {
    bVec = allSamples$b
    eVec = allSamples$e

    tempHolder = lapply(seq_along(bVec), function(i) {
      abs(.spect(
        x = x[bVec[i]:eVec[i]],
        n = wl,
        window = window,
        overlap = overlap
      ))
    })

    globalMax = do.call(max, tempHolder)

    BGNPOWdf = lapply(tempHolder, function(singleBin) {
      singleBin = dbRef * log10(singleBin / globalMax)

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

  channelData = switch(
    channel,
    "stereo" = list(left  = channelData[1, ], right = channelData[2, ]),
    "mono"   = list(mono  = channelData[1, ]),
    "left"   = list(left  = channelData[1, ]),
    "right"  = list(right = channelData[2, ])
  )

  noiseOBJ@values = lapply(channelData, function(x) {
    tempHolder = apply(allSamples, 1, function(y) {
      list(.spect(
        x = x[y[1]:y[2]],
        n = wl,
        window = window,
        overlap = overlap
      ))
    })

    ACIdf = lapply(tempHolder, function(singleBin) {
      spectS = abs(singleBin[[1]])

      specDim = dim(spectS)
      duration = length(spectS) / samp.rate

      jump = if (is.null(j))
        ncol(spectS)
      else
        max(1L, as.integer(j * samp.rate / (wl - overlap)))
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

      ACIvect

    })

    return(list(ACI = data.frame(ACIdf) |>
                  setNames(paste0(
                    "ACI", 1:frameBin
                  ))))

  })

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

  channelData = switch(
    channel,
    "stereo" = list(left  = channelData[1, ], right = channelData[2, ]),
    "mono"   = list(mono  = channelData[1, ]),
    "left"   = list(left  = channelData[1, ]),
    "right"  = list(right = channelData[2, ])
  )

  noiseOBJ@values = lapply(channelData, function(x) {
    tempHolder = apply(allSamples, 1, function(y) {
      list(.spect(
        x = x[y[1]:y[2]],
        n = wl,
        window = window,
        overlap = overlap
      ))
    })

    ENTdf = lapply(tempHolder, function(singleBin) {
      spectS = abs(singleBin[[1]])
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

  channelData = switch(
    channel,
    "stereo" = list(left  = channelData[1, ], right = channelData[2, ]),
    "mono"   = list(mono  = channelData[1, ]),
    "left"   = list(left  = channelData[1, ]),
    "right"  = list(right = channelData[2, ])
  )

  if (DCfix) {
    channelData = lapply(channelData, function(x)
      x - mean(x))
  }

  noiseOBJ@values = lapply(channelData, function(x) {
    bVec = allSamples$b
    eVec = allSamples$e

    tempHolder = lapply(seq_along(bVec), function(i) {
      abs(.spect(
        x = x[bVec[i]:eVec[i]],
        n = wl,
        window = window,
        overlap = overlap
      ))
    })

    globalMax = do.call(max, tempHolder)

    EVNdf = lapply(tempHolder, function(singleBin) {
      singleBin = dbRef * log10(singleBin / globalMax)

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

  noiseOBJ@index = c("EVN")
  noiseOBJ@timeBins = setNames(round((allSamples$e - allSamples$b) / samp.rate), paste0("BIN", seq(frameBin)))

  noiseOBJ@sampRate = samp.rate
  noiseOBJ@channel = channel

  return(noiseOBJ)

}
