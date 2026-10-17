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

  noiseOBJ@values = lapply(channelData, function(x) {
    if (DCfix) {
      x = x - mean(x)
    }

    bVec = allSamples$b
    eVec = allSamples$e

    tempHolder = lapply(seq_along(bVec), function(i) {
      .spect(
        x = x[bVec[i]:eVec[i]],
        n = wl,
        window = window,
        overlap = overlap
      )
    })

    BGNPOWdf = lapply(tempHolder, function(singleBin) {
      spectS = abs(singleBin)

      spectS = dbRef * log10(spectS / max(spectS))

      if (!is.null(dbThreshold)) {
        spectS[spectS < dbThreshold] = dbThreshold
      }

      dbDif = matrixStats::rowRanges(spectS)
      dbMax = dbDif[, 2]
      dbMin = dbDif[, 1]

      numBins = hBreak(spectS)
      binWidth = (dbMax - dbMin) / numBins

      modalBin = vapply(seq_len(wl / 2), function(i) {
        bins = floor((spectS[i, ] - dbMin[i]) / binWidth[i]) + 1L
        which.max(tabulate(bins, nbins = numBins[i]))
      }, integer(1))
      modalIntensity = dbMin + modalBin * binWidth
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

      jump = as.integer(j / (duration / specDim[2]))
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

  noiseOBJ@values = lapply(channelData, function(x) {
    if (DCfix) {
      x = x - mean(x)
    }

    tempHolder = lapply(seq_along(allSamples$b), function(i) {
      .spect(
        x = x[allSamples$b[i]:allSamples$e[i]],
        n = wl,
        window = window,
        overlap = overlap
      )
    })

    EVNdf = lapply(tempHolder, function(singleBin) {
      spectS = abs(singleBin)

      spectS = dbRef * log10(spectS / max(spectS))

      if (!is.null(dbThreshold)) {
        spectS[spectS < dbThreshold] = dbThreshold
      }

      dbDif = matrixStats::rowRanges(spectS)
      dbMax = dbDif[, 2]
      dbMin = dbDif[, 1]

      numBins = hBreak(spectS)
      binWidth = (dbMax - dbMin) / numBins

      modalBin = vapply(seq_len(wl / 2), function(i) {
        bins = floor((spectS[i, ] - dbMin[i]) / binWidth[i]) + 1
        which.max(tabulate(bins, nbins = numBins[i]))
      }, integer(1))

      modalIntensity = dbMin + modalBin * binWidth
      noiseReduced = spectS - modalIntensity
      matrixStats::rowSums2(diff(noiseReduced >= 3) == 1) / ((ncol(spectS) - 1) * (wl - overlap) / samp.rate)

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
