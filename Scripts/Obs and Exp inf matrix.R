# Function to obtain the observed information matrix J

obs_inf_DMOLBE <- function(y, mu, sigma) {
  # Total sample size
  n <- length(y)
  
  # --- 1. Structural Components ---
  A <- (1 + y / mu) * exp(-y / mu)
  B <- (1 + (y + 1) / mu) * exp(-(y + 1) / mu)
  
  N <- A - B
  D1 <- 1 - (1 - sigma) * A
  D2 <- 1 - (1 - sigma) * B
  
  # --- 2. First Derivatives wrt mu ---
  dA <- (y^2 / mu^3) * exp(-y / mu)
  dB <- ((y + 1)^2 / mu^3) * exp(-(y + 1) / mu)
  
  dN  <- dA - dB
  dD1 <- -(1 - sigma) * dA
  dD2 <- -(1 - sigma) * dB
  
  # --- 3. Second Derivatives wrt mu ---
  d2A <- (y^2 / mu^4) * (y / mu - 3) * exp(-y / mu)
  d2B <- ((y + 1)^2 / mu^4) * ((y + 1) / mu - 3) * exp(-(y + 1) / mu)
  
  d2N  <- d2A - d2B
  d2D1 <- -(1 - sigma) * d2A
  d2D2 <- -(1 - sigma) * d2B
  
  # --- 4. Matrix Entries Evaluation (Element-wise over vector y) ---
  # J_mumu component elements
  term_N  <- (d2N * N - dN^2) / N^2
  term_D1 <- (d2D1 * D1 - dD1^2) / D1^2
  term_D2 <- (d2D2 * D2 - dD2^2) / D2^2
  J_mumu_elements <- -term_N + term_D1 + term_D2
  
  # J_sigmasigma component elements
  J_sig_elements <- (1 / sigma^2) - (A / D1)^2 - (B / D2)^2
  
  # J_musigma cross-component elements
  J_musig_elements <- (dA / D1^2) + (dB / D2^2)
  
  # --- 5. Total Observed Information Matrix (Summing over the sample) ---
  J_mumu   <- sum(J_mumu_elements)
  J_sigsig <- sum(J_sig_elements)
  J_musig  <- sum(J_musig_elements)
  
  # Constructing the 2x2 matrix
  info_matrix <- matrix(c(J_mumu,  J_musig, 
                          J_musig, J_sigsig), 
                        nrow = 2, ncol = 2)
  
  colnames(info_matrix) <- rownames(info_matrix) <- c("mu", "sigma")
  
  return(info_matrix)
}


exp_inf_DMOLBE <- function(n, mu, sigma, tol = 1e-10, max_iter = 10000) {
  
  # Initialize sums for the elements of the single-observation expected info
  E_mumu   <- 0
  E_sigsig <- 0
  E_musig  <- 0
  
  y <- 0
  converged <- FALSE
  
  while (!converged && y <= max_iter) {
    # 1. Compute the probability weight using your custom PMF function
    pmf_y <- dDMOLBE(y, mu = mu, sigma = sigma, log = FALSE)
    
    # If the probability is zero or negligible, skip to avoid NaN/Inf divisions
    if (pmf_y < 1e-15) {
      y <- y + 1
      next
    }
    
    # 2. Structural Components for the current 'y'
    A <- (1 + y / mu) * exp(-y / mu)
    B <- (1 + (y + 1) / mu) * exp(-(y + 1) / mu)
    
    # Check if components are structural zero to avoid division by zero
    N <- A - B
    if (N <= 0) N <- 1e-15 
    
    D1 <- 1 - (1 - sigma) * A
    D2 <- 1 - (1 - sigma) * B
    
    # 3. First Derivatives wrt mu
    dA <- (y^2 / mu^3) * exp(-y / mu)
    dB <- ((y + 1)^2 / mu^3) * exp(-(y + 1) / mu)
    
    dN  <- dA - dB
    dD1 <- -(1 - sigma) * dA
    dD2 <- -(1 - sigma) * dB
    
    # 4. Second Derivatives wrt mu
    d2A <- (y^2 / mu^4) * (y / mu - 3) * exp(-y / mu)
    d2B <- ((y + 1)^2 / mu^4) * ((y + 1) / mu - 3) * exp(-(y + 1) / mu)
    
    d2N  <- d2A - d2B
    d2D1 <- -(1 - sigma) * d2A
    d2D2 <- -(1 - sigma) * d2B
    
    # 5. Observed Information Matrix entries evaluated at the current 'y'
    term_N  <- (d2N * N - dN^2) / N^2
    term_D1 <- (d2D1 * D1 - dD1^2) / D1^2
    term_D2 <- (d2D2 * D2 - dD2^2) / D2^2
    
    J_mumu   <- -term_N + term_D1 + term_D2
    J_sigsig <- (1 / sigma^2) - (A / D1)^2 - (B / D2)^2
    J_musig  <- (dA / D1^2) + (dB / D2^2)
    
    # 6. Weight by PMF and add to the running expectation sum
    inc_mumu   <- J_mumu * pmf_y
    inc_sigsig <- J_sigsig * pmf_y
    inc_musig  <- J_musig * pmf_y
    
    E_mumu   <- E_mumu + inc_mumu
    E_sigsig <- E_sigsig + inc_sigsig
    E_musig  <- E_musig + inc_musig
    
    # 7. Check Convergence Criteria 
    # If the largest increment is below the tolerance threshold, stop iterating
    if (max(abs(c(inc_mumu, inc_sigsig, inc_musig))) < tol && y > 10) {
      converged <- TRUE
    }
    
    y <- y + 1
  }
  
  # Scale by sample size n to get the total Expected Information Matrix
  exp_matrix <- matrix(c(E_mumu,  E_musig, 
                         E_musig, E_sigsig), 
                       nrow = 2, ncol = 2) * n
  
  colnames(exp_matrix) <- rownames(exp_matrix) <- c("mu", "sigma")
  
  return(exp_matrix)
}

get_parameter_se <- function(model) {
  # Ensure the object is a gamlss model
  if (!inherits(model, "gamlss")) {
    stop("The input model must be a gamlss object.")
  }
  
  # 1. Extract estimated intercept coefficients on the link scale
  beta0  <- coef(model, what="mu")[1]
  gamma0 <- coef(model, what="sigma")[1]
  
  # 2. Back-transform to original scale point estimates
  mu_est    <- exp(beta0)
  sigma_est <- exp(gamma0)
  
  # 3. Extract the standard errors on the link scale from the summary/vcov
  # vcov(model) extracts the full asymptotic variance-covariance matrix
  # of the predictors on the link scale
  var_matrix <- vcov(model)
  
  # Extract specific link-scale standard errors via the square root of diagonal variances
  # Note: names match the internal gamlss structure for intercepts
  se_beta0  <- sqrt(var_matrix[1, 1])
  se_gamma0 <- sqrt(var_matrix[2, 2])
  
  # 4. Apply the Delta Method to find original scale standard errors
  se_mu    <- mu_est * se_beta0
  se_sigma <- sigma_est * se_gamma0
  
  # 5. Build a neat summary matrix
  result <- matrix(
    c(mu_est, se_mu, sigma_est, se_sigma),
    nrow=2, byrow=TRUE,
    dimnames=list(c("mu", "sigma"), c("Estimate", "Std. Error"))
  )
  
  return(result)
}

# Example data sample
library(DiscreteDists)

n <- 1500
mu <- 0.5
sigma <- 0.5

y <- rDMOLBE(n=n, mu=mu, sigma=sigma) # Hypothetical discrete data

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, sigma.fo=~1, family=DMOLBE)

get_parameter_se(mod1)

# Using fitdistrplus
library(fitdistrplus)
fit <- fitdist(data = y, distr = "DMOLBE",
               start = list(mu=1, sigma=1))
fit

llsurface(data=y, distr="DMOLBE", plot.arg=c("mu", "sigma"), 
          min.arg=c(0.2, 0.2), max.arg=c(0.8, 0.8), back.col = FALSE,
          main="log-likelihood for DMOLBE distribution")
points(x=mu, y=sigma, col="blue", pch=20)

# Extracting the fitted values for mu and sigma
# using the inverse link function
mu_est <- exp(coef(mod1, what="mu"))
sigma_est <- exp(coef(mod1, what="sigma"))

# Get the observed information matrix
J <- obs_inf_DMOLBE(y=y, mu=mu_est, sigma=sigma_est)
sqrt(diag(solve(J)))

I <- exp_inf_DMOLBE(mu=mu, sigma=sigma, n=n)
sqrt(diag(solve(I)))

# Simulation study --------------------------------------------------------

nrep <- 1500
results <- matrix(NA, ncol=2, nrow=nrep)

for (i in 1:nrep) {
  y <- rDMOLBE(n=n, mu=mu, sigma=sigma)
  mod_aux <- gamlss(y~1, sigma.fo=~1, family=DMOLBE)
  mu_est <- exp(coef(mod_aux, what="mu"))
  sigma_est <- exp(coef(mod_aux, what="sigma"))
  results[i, ] <- c(mu_est, sigma_est)
}

head(results)

apply(results, MARGIN=2, FUN=mean, na.rm=TRUE)
apply(results, MARGIN=2, FUN=sd, na.rm=TRUE)
sqrt(diag(solve(J)))
sqrt(diag(solve(I)))

