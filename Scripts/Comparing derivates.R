
# Derivada con respecto a mu

dldm_manual = function(y, mu, sigma) {
  t1 <- (1 + y/mu) * exp(-y/mu)
  t2 <- (1 + (y+1)/mu) * exp(-(y+1)/mu)
  A  <- t1 - t2
  B1 <- 1 - (1 - sigma) * t1
  B2 <- 1 - (1 - sigma) * t2
  dA  <- (y^2 / mu^3) * exp(-y/mu) - ((y+1)^2 / mu^3) * exp(-(y+1)/mu)
  dB1 <- -(1 - sigma) * (y^2 / mu^3) * exp(-y/mu)
  dB2 <- -(1 - sigma) * ((y+1)^2 / mu^3) * exp(-(y+1)/mu)
  dldm <- (dA / A) - (dB1 / B1) - (dB2 / B2)

  dldm
}

dldm_compu = function(y, mu, sigma) {
  dm   <- gamlss::numeric.deriv(dDMOLBE(y, mu, sigma, log=TRUE),
                                theta="mu",
                                delta=0.00001)
  dldm <- as.vector(attr(dm, "gradient"))
  dldm
}

y <- 2:6
mu <- 2:6
sigma <- 2:6

dldm_manual(y, mu, sigma)
dldm_compu(y, mu, sigma)

library(microbenchmark)

res <- microbenchmark(dldm_manual(y, mu, sigma),
                      dldm_compu(y, mu, sigma),
                      times=100)
plot(res)

# Derivada con respecto a sigma

dldd_manual = function(y, mu, sigma) {
  t1 <- (1 + y/mu) * exp(-y/mu)
  t2 <- (1 + (y+1)/mu) * exp(-(y+1)/mu)
  B1 <- 1 - (1 - sigma) * t1
  B2 <- 1 - (1 - sigma) * t2
  dldd <- (1 / sigma) - (t1 / B1) - (t2 / B2)
  dldd
}

dldd_compu = function(y, mu, sigma) {
  dd   <- gamlss::numeric.deriv(dDMOLBE(y, mu, sigma, log=TRUE),
                                theta="sigma",
                                delta=0.00000001)
  dldd <- as.vector(attr(dd, "gradient"))
  dldd
}

y <- 2:6
mu <- 2:6
sigma <- 2:6

dldd_manual(y, mu, sigma)
dldd_compu(y, mu, sigma)

library(microbenchmark)

res <- microbenchmark(dldd_manual(y, mu, sigma),
                      dldd_compu(y, mu, sigma),
                      times=100)
plot(res)

