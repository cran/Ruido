#' Calculate A Hann Window
#'
#' @description Generates a Hann window of a specified length.
#'
#' @param n single positive integer greater than 1 specifying the length of the window
#'
#' @returns A numeric vector containing the Hann window.
#'
#' @details The Hann window is a cosine window commonly used to reduce spectral leakage before Fourier transformation.
#'
#' The window is calculated following the formula described in Harris, 1978:
#' \deqn{
#' w_k = 0.5 - 0.5 \cos\left(\frac{2\pi k}{n - 1}\right)
#' }{
#' w[k] = 0.5 - 0.5 cos(2 pi k / (n - 1))
#' }
#' for \eqn{k = 0, \ldots, n - 1}.
#'
#' @export
#'
#' @references
#' Harris, F. J. (1978). On the use of windows for harmonic analysis with the discrete Fourier transform. \emph{Proceedings of the IEEE}, 66(1), 51–83. https://doi.org/10.1109/proc.1978.10837
#'
#' @examples
#' hannWindow = hanning(256)
#' plot(hannWindow, type = "l")
hanning = function(n) {
  if (length(n) != 1 || n <= 1 || n != as.integer(n)) {
    stop("'n' must be a single integer greater than 1")
  }

  n = n - 1
  0.5 - 0.5 * cos(2 * pi * 0:n / n)

}

#' Calculate A Hamming Window
#'
#' @description Generates a Hamming window of a specified length.
#'
#' @param n single positive integer greater than 1 specifying the length of the window
#'
#' @returns A numeric vector containing the Hamming window.
#'
#' @details The Hamming window is a cosine window commonly used to reduce spectral leakage before Fourier transformation.
#'
#' The window is calculated following the formula described in Harris, 1978:
#' \deqn{
#' w_k = 0.54 - 0.46 \cos\left(\frac{2\pi k}{n - 1}\right)
#' }{
#' w[k] = 0.54 - 0.46 cos(2 pi k / (n - 1))
#' }
#' for \eqn{k = 0, \ldots, n - 1}.
#'
#' @export
#'
#' @references
#' Harris, F. J. (1978). On the use of windows for harmonic analysis with the discrete Fourier transform. \emph{Proceedings of the IEEE}, 66(1), 51–83. https://doi.org/10.1109/proc.1978.10837
#'
#' @examples
#' hammingWindow = hamming(256)
#' plot(hammingWindow, type = "l")
hamming = function(n) {
  if (length(n) != 1 || n <= 1 || n != as.integer(n)) {
    stop("'n' must be a single integer greater than 1")
  }

  n = n - 1
  0.54 - 0.46 * cos(2 * pi * 0:n / n)

}
