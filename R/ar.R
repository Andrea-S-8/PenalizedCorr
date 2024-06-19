######
#' @title ar
#'
#' @description Fit an autoregressive time series model to the data using penalized correlation estimation, by default selecting the order by AIC.
#'series = deparse1(substitute(x)), 
#' @param x a univariate numeric time series.
#' @param aic 'logical', If 'TRUE' then the Akaike Information Criterion is used to choose the order of the autoregressive model. If FALSE, the model of order order.max is fitted.
#' @param order.max maximum lag at which to calculate the penalized/sample PACF. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.
#' @param method character string specifying the method to fit the regular autoregressive time series model. Default to "yule-walker".
#' @param na.action function to be called to handle missing values. 'na.pass' can be used.
#' @param series names for the series. Defaults to \code{deparse1(substitute(x))}.
#' @param ... additional arguments for \code{ar} function.
#'
#' @return An object of penalized ar model fit with the following elements:
#' \describe{
#' \item{\code{order}}{The order of the fitted model. This is chosen by minimizing the AIC if 'AIC=TRUE', otherwise it is 'order.max'.}
#' \item{\code{ar}}{Estimated penalized autorregression coefficients for the fitted model.}
#' \item{\code{aic}}{The vector of AIC values.}
#' \item{\code{n.used}}{The number of observations in the time series.}
#' \item{\code{order.max}}{The value of the 'order.max' argument.}
#' \item{\code{call}}{The matched call.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' # Example for penalized ar model fit
#' ar(data)
#' }
#' @export
#' @importFrom stats na.fail
#' @importFrom stats ar.yw
#' @importFrom stats ar.burg
#' @importFrom stats ar.ols
#' @importFrom stats ar.mle
#####

ar <-
function (x, aic = TRUE, order.max = NULL, method = c("penyw","yule-walker", 
    "burg", "ols", "mle", "yw"), na.action = na.fail, series = deparse1(substitute(x)), 
    ...) 
{
    res <- switch(match.arg(method), penyw = ar.penyw(x,aic=aic,
        order.max=order.max,na.action=na.action,series=series,,...), yw=,
        `yule-walker` = ar.yw(x, 
        aic = aic, order.max = order.max, na.action = na.action, 
        series = series, ...), burg = ar.burg(x, aic = aic, order.max = order.max, 
        na.action = na.action, series = series, ...), ols = ar.ols(x, 
        aic = aic, order.max = order.max, na.action = na.action, 
        series = series, ...), mle = ar.mle(x, aic = aic, order.max = order.max, 
        na.action = na.action, series = series, ...))
    res$call <- match.call()
    res
}
