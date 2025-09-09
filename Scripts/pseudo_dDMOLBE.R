dDMOLBE <- function(x, mu=1, sigma=1, log=FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(x < 0))       stop(paste("x must be >=0", "\n", ""))
  k1 <- (1+x/mu) * exp(-x/mu)
  k2 <- (1+(x+1)/mu) * exp(-(x+1)/mu)
  p1 <- log(sigma) + log(k1 - k2)
  p2 <- -log(1-(1-sigma)*k1) - log(1-(1-sigma)*k2)
  res <- p1 + p2
  if(log)
    return(res)
  else
    return(exp(res))
}
#' @export
#' @rdname dDMOLBE
pDMOLBE <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(q < 0))       stop(paste("q must be >=0", "\n", ""))
  ly <- max(length(q), length(mu), length(sigma))
  q <- rep(q, length = ly)
  mu <- rep(mu, length = ly)
  sigma <- rep(sigma, length = ly)
  num <- 1 - (1+(q+1)/mu) * exp(-(q+1)/mu)
  den <- 1 - (1-sigma) * (1 + (q+1)/mu) * exp(-(q+1)/mu)
  cdf <- num / den
  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @importFrom stats runif
#' @export
#' @rdname dDMOLBE
rDMOLBE <- function(n, mu=1, sigma=1) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(n <= 0))      stop(paste("n must be a positive integer", "\n", ""))
  # Begin auxiliar function
  one_random_DMOLBE <- function(u, mu, sigma) {
    p <- dDMOLBE(x=0, mu=mu, sigma=sigma, log=FALSE)
    F <- p
    i <- 0
    while (u >= F) {
      i <- i + 1
      p <- dDMOLBE(x=i, mu=mu, sigma=sigma, log=FALSE)
      F <- F + p
    }
    return(i)
  }
  one_random_DMOLBE <- Vectorize(one_random_DMOLBE)
  # End auxiliar function
  one_random_DMOLBE(u=runif(n), mu, sigma)
}
#' @export
#' @rdname dDMOLBE
qDMOLBE <- function(p, mu = 1, sigma = 1, lower.tail = TRUE,
                    log.p = FALSE) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  # Begin auxiliar function
  one_quantile_DMOLBE <- function(p, mu, sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      prob <- dDMOLBE(x=0, mu=mu, sigma=sigma, log=FALSE)
      F <- prob
      i <- 0
      while (p >= F) {
        i <- i + 1
        prob <- dDMOLBE(x=i, mu=mu, sigma, log=FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_DMOLBE <- Vectorize(one_quantile_DMOLBE)
  # End auxiliar function
  one_quantile_DMOLBE(p=p, mu=mu, sigma=sigma)
}