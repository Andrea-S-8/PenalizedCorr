######
#' @title bandtap
#'
#' @description Calculated the corrected (partial) autocorrelation function
#' @param acf incoming acf value
#' @param n incoming length of x from acf funtion
#' @param ... additional arguments passed to plotting.
#' 
#' @references Gallagher, C., Killick, R., Tan, X. (2024+) Penalized M-estimation 
#' for Autocorrelation. \emph{Submitted.}
#' McMurray & Politis (2010) forecast package for CI
#####

bandtap <- function (acf, n, ...){
  
  n <- length(x)
  lh <- 2*sqrt(log10(n)/n)
  #compares the absolute value of each acf with the threshold.
  thresh = abs(acf) < lh
  rle = rle(thresh)
  #finds position in the run-length encoding where a value repeats 5 times
  rle5 = which(rle==5)
  r <- which(rle[rle5] == TRUE)
  w <- c(rep(1,r), (2-((r+1):2*r)/r), rep(0, n-2*r))
  return (acf*w)
}
