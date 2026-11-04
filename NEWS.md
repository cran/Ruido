# Ruido 1.2.2

- Standardized function names across the package. Previous function names have been removed and replaced with the following:
  - activity -> ACTspec()
  - bgNoise() -> BGNspec()
  - bgn -> BGN()
  - singleSat() -> SAT()
  - soundSat() -> soundSAT()
  - soundMat() -> soundMAT()
  - multActivity() -> soundACT()
- Standardized the `noise.matrix` S4 class to always include `@path` and `@wl` slots
- Reverted the global maximum normalization introduced in version 1.2.1. Spectrogram frames are once again normalized independently within each temporal bin, consistent with the original methodology
- Bug fixes:
  - Fixed histogram bin calculation in `BGN()` to prevent `BGN` values from being overestimated by half a bin
  - Fixed an issue in `BGN()` that caused the lowest histogram bin to be omitted during `BGN` value tabulation
  - Fixed an issue in `BGN()` where `dbRef` was hardcoded to `-10`, preventing the user-specified value from being properly applied
- Updated documentation to use standardized index names (e.g., BGNsp, POWsp, etc.)
- Expanded function examples with additional audio signals generated using the tuneR package

# Ruido 1.2.1

- Added the `saveNoise` argument to batch processing functions (`soundSat()`, `soundMat()` and `multActivity()`), allowing internally generated `noise.matrix` objects to be included in the output. Defaults to `FALSE` to preserve the output structure of previous versions. This argument is also inherited by `satBackup()`
- Updated batch processing functions to no longer automatically delete backup files after processing, leaving their management to the user
- Improved the processing speed of `activity()` and `singleSat()` through the use of `matrixStats` functions and vectorized operations
- Moved DC offset correction outside the time-bin processing loop, improving both performance and implementation of the algorithm
- Bug fixes:
  - Fixed `EVNspec()` returning one fewer row than expected in `noise.matrix@values`
  - Fixed histogram bin calculation to prevent BGN values from being overestimated by half a bin
  - Fixed decibel conversion so that spectrogram frames are no longer normalized independently by their individual maximum amplitudes. A single global maximum is now used for the entire spectrogram
- Documentation updates:
  - Corrected the height description in `sampleBGN()` from 4'11" to 1.5 meters
  - Fixed documentation typos
  - Fixed function masking in `bgNoise()` examples
  - Added the missing `@slot` documentation for `noise.matrix`

# Ruido 1.2.0

- Added `matrixStats` dependency for efficient matrix operations
- Removed `signal` dependency; spectrogram and window generation are now implemented natively
- Added two new acoustic index functions:
  - `bgn()` for Summarized Background Noise and Soundscape Power
  - `EVNspec()` for the Spectral Acoustic Event Index
- Added the `dbRef` argument to functions that perform decibel conversion, allowing users to specify the decibel reference. The default is `10`, so existing results remain unchanged
- Reworked `normHandler()` to eliminate user prompts and automatically select an appropriate `normality` test based on the number of files
- Improved `bgNoise()` processing performance by:
  - Reducing unnecessary data conversions
  - Using the `matrixStats` package to optimize matrix calculations
  - Making object formatting more efficient when the input is an S4 object
  - Implementing C code for more efficient spectrogram matrix calculation
  - Vectorizing histogram breaks calculations
- Updated and expanded function documentation, including:
  - Adding the decibel conversion formula to `bgNoise()`.
  - Specifying "Spectral" in the titles of relevant spectral acoustic indices
  - Documenting the normality tests accepted by the `normality` argument

# Ruido 1.1.1

- Updated `activity()` to output a `noise.matrix` object
- Updated `singleSat()` to accept `noise.matrix` objects produced by `activity()`
- Added examples and documentation clarifying the use of `activity()` together with `singleSat()`
- Fixed a bug in the `plot()` method for `noise.matrix` objects, where the `annotate` argument was being ignored

# Ruido 1.1.0

- Updated package sub-title and description to better reflect current and future development focus
- Added two new acoustic index functions:
  - `ACIspec()` for Spectral Acoustic Complexity Index
  - `ENTspec()` for Spectral Temporal Entropy Index
- Removed native `.mp3` file support. Users can convert `.mp3` to `.wav` or load as `tuneR` Wave objects beforehand
- Migrated primary WAV file reading from `tuneR::readWave()` to `wav::read_wav()` for improved performance
  - Full backward compatibility maintained: `tuneR` Wave objects continue to work as input
- Improved error handling across all functions for clearer user feedback
- Fixed logical errors in channel handling for stereo/mono audio processing
- Removed redundant code and unused variables
- Cleaned up batch processing functions (removed explicit `gc()` calls since R manages memory by itself)
- Enhanced internal helper functions to support new and future acoustic indices
- Updated documentation throughout
  
# Ruido 1.0.3

- Added `noise.matrix` object and associated methods to improve function output handling
- Added two internal helper functions for error handling and normality testing:
  - `argHandler()`
  - `normHandler()`
- Added the `DCfix` argument across all functions to allow users to optionally skip DC offset removal
- Added sample dataset `sampleBGN` for testing
- Added package `testhat` to package suggests
- Updated the `timeBin` argument to actually accept `NULL` (#10)
- Updated `singleSat()` output to be consistent when working with both stereo and mono files (#11)
- Moved all auxiliary functions to a single file

# Ruido 1.0.2

- Added functions:
 - `activity()`
 - `multActivity()`
- Updated the output of already existing functions to include channel, normality statistics and path to audio. New functions also follow the same output pattern
- Added extra examples to show more of the package's utilities. Also updated existing examples to accommodate the changes made to the outputs
- Greatly reduced the number of lines by improving how audios are processed and managed inside the functions
- Improved error handling for all functions
- Overall documentation improvements

# Ruido 1.0.1

- `bgnoise()`, `soundSat()` and `soundMat()` now have better mono audio files handling
- In `soundSat()` and `soundMat()` warnings messages will no longer be included in the output
- `satBackup()` can now continue an unfinished process of `soundMat()`
- Removed the `spelling` package from suggests

# Ruido 1.0.0

* This is the initial CRAN submission.

- Functions:
  - `soundSat()`
  - `singleSat()`
  - `bgnNoise()`
  - `satBackup()`
  - `soundMat()`
