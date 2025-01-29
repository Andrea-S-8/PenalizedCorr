######
#' @title Fit Autoregressive Models to Time Series
#'
#' @description Fit an autoregressive time series model to the data, by default selecting the complexity by AIC.
#' Updates the \code{stats::ar} function to include penalized acf/pacf estimation in the Yule-Walker method as the default.
#'
#' @param x a univariate or multivariate numeric time series.
#' @param aic 'logical', If \code{TRUE} then the Akaike Information Criterion is used to choose the order of the autoregressive model. If \code{FALSE}, the model of order \code{order.max} is fitted.
#' @param order.max maximum order (or order) of model to fit.  Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)}
#'  where \eqn{N} is the number of non-missing observations except for \code{method="mle"} where it is the minimum of this quantity and 12.  \code{nser} is the number of series.
#' @param method character string specifying the method to fit the model. Must be one of the strings in the default argument (the first few characters are sufficient).  Default to "penyw".
#' @param na.action function to be called to handle missing values. Currently, via \code{na.action=na.pass}, only (penalized) Yule-Walker methods can handle missing values which must be consistent within a time point; either all variables must be missing or none.
#' @param series names for the series. Defaults to \code{deparse1(substitute(x))}.
#' @param lh vector of length 1 (value used for all lags), or length lag.max. Default uses formula in the description.
#' @param ... additional arguments for specific methods.
#'
#' @details
#' For definiteness, note that the AR coefficients have the sign in
#' \eqn{x_t - \mu = a_1(x_{t-1}-\mu)+\ldots+a_p(x_{t-p}-\mu)+e_t}
#' 
#' \code{ar} is just a wrapper for the functions \code{ar.penyw}, \code{ar.yw}, \code{ar.burg}, \code{ar.ols} and \code{ar.mle}.
#' 
#' Order selection is done by AIC if \code{aic} is true. This is problematic, as of the methods here only \code{ar.mle} performs true maximum likelihood estimation. The AIC is computed as if the variance estimate were the MLE, omitting the determinant term from the likelihood. Note that this is not the same as the Gaussian likelihood evaluated at the estimated parameter values. In \code{ar.yw} and \code{ar.penyw} the variance matrix of the innovations is computed from the fitted coefficients and the autocovariance of \eqn{x}.
#' 
#' \code{ar.burg} allows two methods to estimate the innovations variance and hence AIC. Method 1 is to use the update given by the Levinson-Durbin recursion (Brockwell and Davis, 1991, (8.2.6) on page 242), and follows S-PLUS. Method 2 is the mean of the sum of squares of the forward and backward prediction errors (as in Brockwell and Davis, 1996, page 145). Percival and Walden (1998) discuss both. In the multivariate case the estimated coefficients will depend (slightly) on the variance estimation method.
#' 
#' Remember that \code{ar} includes, by default, a constant in the model, by removing the overall mean of \eqn{x} before fitting the AR model, or (\code{ar.mle}) estimating a constant to subtract.

#' 
#' @return An object of class \code{ar} with the following elements:
#' \describe{
#' \item{\code{order}}{The order of the fitted model. This is chosen by minimizing the AIC if \code{aic=TRUE}, otherwise it is \code{order.max}.}
#' \item{\code{ar}}{Estimated autorregression coefficients for the fitted model.}
#' \item{\code{var.pred}}{The prediction variance: an estimate of the portion of the variance of the time series that is not explained by the autoregressive model.}
#' \item{\code{x.mean}}{The estimated mean of the series used in fitting and for use in prediction.}
#' \item{\code{aic}}{The differences in AIC between each model and the best-fitting model.  Note that the latter can have an AIC of \code{-Inf}.}
#' \item{\code{n.used}}{The number of observations in the time series, including missing.}
#' \item{\code{n.obs}}{The number of non-missing observationsin the time series.}
#' \item{\code{order.max}}{The value of the \code{order.max} argument.}
#' \item{\code{partialacf}}{The estimate of the partial autocorrelation function up to lag \code{order.max}.}
#' \item{\code{resid}}{Residuals from the fitted model, conditioning on the first \code{order} observations.  The first \code{order} residuals are set to \code{NA}.  If \eqn{x} is a time series, so is \code{resid}.}
#' \item{\code{method}}{Value of the \code{method} argument.}
#' \item{\code{series}}{The name(s) of the time series.}
#' \item{\code{frequency}}{The frequency of the time series.}
#' \item{\code{call}}{The matched call.}
#' }
#'
#' @references Gallagher, C., Killick, R., Tan, X. (2024+) Penalized M-estimation 
#' for Autocorrelation. \emph{Submitted.}
#' 
#' @examples
#' \dontrun{
#' set.seed(1234)
#' data <- arima.sim(n=100, model=list(ar=0.5))
#'
#' ar(data) # penalized ar model fit
#' ar(data,method="yw") # default stats::ar() estimate
#' 
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
    lh=NULL,...) 
{
    res <- switch(match.arg(method), penyw = ar.penyw(x,aic=aic,
        order.max=order.max,na.action=na.action,series=series,lh=lh,...), yw=,
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
