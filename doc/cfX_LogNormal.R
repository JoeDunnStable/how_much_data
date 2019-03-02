#' @title Characteristic function of Lognormal distribution
#'
#' @description
#' cfX_LogNormal(t,mu,sigma) Computes the characteristic function cf(t) of
#' the Lognormal distribution with parameters mu (real) and sigma > 0,
#' computed for real (vector) argument t, i.e.
#' cf(t) = cfX_LogNormal(t,mu,sigma);
#'
#' @importFrom stats integrate
#'
#' @details
#' In probability theory, a log-normal (or lognormal) distribution is a
#' continuous probability distribution of a random variable whose logarithm
#' is normally distributed. The lognormal distribution is defined for x in
#' (0,+inf) by its PDF/CDF/CF, as follows
#'  \eqn{pdf(x) = 1/(x*sigma*sqrt(2*pi))*exp(-(ln(x)-mu)^2/(2*sigma^2))
#'  cdf(x) = 1/2+1/2*erf((ln(x)-mu)/(sqrt(2)*sigma))
#'  cf(t)  = sum_0^infinity{(it)^n/n!*exp(n*mu + (n*sigma)^2/2)}.}
#' As noted, this representation is asymptotically divergent but sufficient
#' for numerical purposes.
#'
#' cfX_LogNormal is based on the standard integral representation of the
#' characteristic function of the lognormal distribution, i.e.
#'  cf(t) = Integral_0^inf exp(i*t*x) * PDF(x) dx.
#' By using the half-space Fourier integral transformation we get
#'  cf(t) = Integral_0^inf (i/t) * exp(-x) * PDF(i*x/t) dx.
#' If we define the integrand as funCF(t,x) = (i/t) * exp(-x) * PDF(i*x/t),
#' then by using a stabilizing transformation from [0,inf] to [0,1], we can
#' evaluate the CF by the following (well behaved) integral:
#'  cf(t) = Integral_0^1 2x/(1-x)^3 * funCF(t,(x/(1-x))^2) dx.
#'
#' cfX_LogNormal evaluates this integral by using the R built in
#' function 'integrate', with precission specified by tolerance tol (default
#' value is tol = 1e-6).
#'
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Log-normal_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param mu real, default value mu = 0
#' @param sigma > 0, default value sigma = 1
#' @param tol tolerance, default value tol = 1e-6
#'
#' @return characteristic function cf(t) of the Lognormal distribution computed for real (vector)
#'
#' @example R/Examples/example_cfX_LogNormal.R
#'
#' @export

cfX_LogNormal <- function(t, mu = 0, sigma = 1, tol = 1e-6) {
  reltol <- tol
  
  szt <- dim(t)
  t <- c(t)
  
  cf <- rep(1+0i,length(t))
  cfIm <- Im(cf)
  cfRe <- Re(cf)
  id <- t != 0
  integrand_Re <- function(x) Re(funCF(mu, sigma, t, (x/(1-x))^2) * (2*x/(1-x)^3))
  integrand_Im <- function(x) Im(funCF(mu, sigma, t, (x/(1-x))^2) * (2*x/(1-x)^3))
  cfRe[id] <- unlist(lapply(t[id], FUN = function(t) integrate(integrand_Re, 0, 1, rel.tol = reltol)$value))
  cfIm[id] <- unlist(lapply(t[id], FUN = function(t) integrate(integrand_Im, 0, 1, rel.tol = reltol)$value))
  cf <- array(cfRe + 1i*cfIm)
  
  return(cf)
}

funCF <- function(mu, sigma, t, x) {
  t  <- exp(mu)*t
  szx <- length(x)
  x  <- c(x)
  
  funPDF <- function(x, s) exp(-0.5*(log(x)/s)^2) / (sqrt(2*pi)*s*x)
  
  if (sigma >= 1/3) {
    t <- 1/t
    f <- (1i*t) * exp(-x) * funPDF(1i*x*t, sigma)
  } else {
    # Set optimum limits for small and large abs(t)
    small <- 7*sqrt(1/sigma)
    large <- 25*sqrt(1/sigma)

    if (abs(t) <= small) {
      f <- exp(1i*t*x) * funPDF(x,sigma)
    }
    
    if (abs(t) > small){
      f = exp(1i*x*sign(t))*funPDF(x/abs(t),sigma)/abs(t)
    }

  }
  return(f)
}