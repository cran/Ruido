# Tell R CMD check that these arguments are imported externally in satBackup
utils::globalVariables(c(
  'powthr',
  'bgnthr',
  'wl',
  'timeBin',
  'targetSampRate',
  'overlap',
  'channel',
  'dbThreshold',
  'histbreaks',
  'normality',
  'type'
))
