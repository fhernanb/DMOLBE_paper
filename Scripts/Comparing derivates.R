
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
                      dldm_compu(y, mu, sigma))
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

y <- 2:15
mu <- 2:15
sigma <- 2:15

dldd_manual(y, mu, sigma)
dldd_compu(y, mu, sigma)

library(microbenchmark)

res <- microbenchmark(dldd_manual(y, mu, sigma),
                      dldd_compu(y, mu, sigma))
plot(res)

# Derivada con respecto a mu dos veces

d2ldm2_manual = function(y, mu, sigma){
  ifelse(y == 0,
         -4*(sigma + 1)*mu/((mu - sigma - 1)^2*(mu + sigma + 1)^2),
         (y + 1)/(mu + sigma + 1)^2 + (1 - y)/(mu + sigma - 1)^2 - 1/mu^2
  )
}

d2ldm2_manual2 = function(y, mu, sigma){
  ifelse(y == 0,
         1/(mu+sigma+1)^2-1/(-mu+sigma+1)^2,
         (y+1)/(mu+sigma+1)^2+(1-y)/(mu+sigma-1)-1/mu^2
  )
}

d2ldm2_compu = function(y, mu, sigma) {
  dm   <- gamlss::numeric.deriv(dDMOLBE(y, mu, sigma, log=TRUE),
                                theta="mu",
                                delta=0.00001)
  dldm <- as.vector(attr(dm, "gradient"))

  d2ldm2 <- - dldm * dldm
  d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
  d2ldm2
}

d2ldm2_manual(y, mu, sigma)
d2ldm2_manual2(y, mu, sigma)
d2ldm2_compu(y, mu, sigma)    # Alert: d2ldm2 seems to be incorrect

# Derivada con respecto a mu y luego a sigma

d2ldmdd_manual = function(y, mu, sigma){
  ifelse(y == 0,
         2*((sigma + 1)^2 + mu^2)/((sigma - mu + 1)^2*(mu + sigma + 1)^2),
         2*((mu + sigma)*(mu + sigma - 2*y) + 1)/
           ((mu + sigma + 1)^2*(mu + sigma - 1)^2)
  )
}

d2ldmdd_compu = function(y, mu, sigma) {
  dm   <- gamlss::numeric.deriv(dDMOLBE(y, mu, sigma, log=TRUE),
                                theta="mu",
                                delta=0.00001)
  dldm <- as.vector(attr(dm, "gradient"))

  dd   <- gamlss::numeric.deriv(dDMOLBE(y, mu, sigma, log=TRUE),
                                theta="sigma",
                                delta=0.00001)
  dldd <- as.vector(attr(dd, "gradient"))

  d2ldmdd <- - dldm * dldd
  d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
  d2ldmdd
}

d2ldmdd_manual(y, mu, sigma)
d2ldmdd_compu(y, mu, sigma)   # Alert: d2ldmdd seems to be incorrect

# Derivada con respecto a sigma dos veces

d2ldd2_manual = function(y, mu, sigma){
  ifelse(y == 0,
         -4*mu*(sigma + 1)/((sigma - mu + 1)^2 * (sigma + mu + 1)^2),
         (y + 1)/(mu + sigma + 1)^2 + (1 - y)/(mu + sigma - 1)^2
  )
}

d2ldd2_compu  = function(y, mu, sigma) {
  dd   <- gamlss::numeric.deriv(dDMOLBE(y, mu, sigma, log=TRUE),
                                theta="sigma",
                                delta=0.00001)
  dldd <- as.vector(attr(dd, "gradient"))

  d2ldd2 <- - dldd * dldd
  d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
  d2ldd2
}

d2ldd2_manual(y, mu, sigma)
d2ldd2_compu(y, mu, sigma) # Alert: d2ldd2 seems to be incorrect


