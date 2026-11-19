# Example 1 ---------------------------------------------------------------

# Data Set I (death due to coronavirus in China). The first data set is the number
# of deaths due to coronavirus in China from 23 January to 28 March.
# The data sets used in the paper was collected from 2020 year. The data set
# is reported in https://www.worldometers.info/coronavirus/country/china/.
# The data are:

y <- c(8, 16, 15, 24, 26, 26, 38, 43, 46, 45, 57, 64, 65, 73, 
       73, 86, 89, 97, 108, 97, 146, 121, 143, 142, 105, 98, 
       136, 114, 118, 109, 97, 150, 71, 52, 29, 44, 47, 35, 
       42, 31, 38, 31, 30, 28, 27, 22, 17, 22, 11, 7, 13, 10, 
       14, 13, 11, 8, 3, 7, 6, 9, 7, 4, 6, 5, 3, 5)

# Fitting the model
library(DiscreteDists)
library(gamlss)

mod1 <- gamlss(y~1, sigma.fo=~1, family=DMOLBE,
               control=gamlss.control(n.cyc=500))

summary(mod1)

# Extracting the fitted values for mu and sigma
# using the inverse link function
mu_hat <- exp(coef(mod1, what="mu"))
mu_hat
sigma_hat <- exp(coef(mod1, what="sigma"))
sigma_hat

# Create a grid of values for mu and sigma around the true values
mu_vals     <- seq(30, 50, length.out=30)
sigma_vals  <- seq(0.1, 1.5, length.out=30)

# Define the log-likelihood function (using sigma, not sigma^2)
log_lik_DMOLBE <- function(mu, sigma, data) {
  ll <- sum(dDMOLBE(x=data, mu=mu, sigma=sigma, log=TRUE))
  return(ll)
}

# Evaluate the log-likelihood on the grid
# We use Vectorize to easily fill the matrix over the grid elements
log_lik_vectorized <- function(m, s) log_lik_DMOLBE(m, s, y)
z_matrix <- outer(mu_vals, sigma_vals, Vectorize(log_lik_vectorized))

# Plot the 3D surface using persp
# Adjust theta and phi to rotate the 3D angle view
persp(
  x=mu_vals, 
  y=sigma_vals, 
  z=z_matrix,
  theta=35,          # Azimuthal direction (rotation around z-axis)
  phi=30,            # Colatitude direction (tilt angle)
  expand=0.6,        # Z-axis expansion factor
  col="lightblue",   # Color of the surface facets
  border="black",    # Color of the line grid
  ltheta=120,        # Light source direction
  shade=0.4,         # Shadow intensity for a 3D effect
  ticktype = "detailed", 
  xlab="mu",
  ylab="sigma",
  zlab="Log-Likelihood",
  main="Log-Likelihood Surface for DMOLBE distribution (mu, sigma)"
)

fit <- fitdist(data=y, distr="DMOLBE", 
               start = list(mu=2, sigma=3))
print(fit)

# Example 2 ---------------------------------------------------------------

# Data Set II (daily death due to coronavirus in Pakistan). The second data
# set is the daily deaths due to coronavirus in Pakistan from 18 March
# to 30 June. The data sets used in the paper was collected from 2020 year.
# The data is reported in
# https://www.worldometers.info/coronavirus/country/Pakistan.
# The data are:

y <- c(1, 6, 6, 4, 4, 4, 1, 20, 5, 2, 3, 15, 17, 7, 8, 25, 
       8, 25, 11, 25, 16, 16, 12, 11, 20, 31, 42, 32, 23, 
       17, 19, 38, 50, 21, 14, 37, 23, 47, 31, 24, 9, 64, 
       39, 30, 36, 46, 32, 50, 34, 32, 34, 30, 28, 35, 57, 
       78, 88, 60, 78, 67, 82, 68, 97, 67, 65, 105, 83, 
       101, 107, 88, 178, 110, 136, 118, 136, 153, 119, 
       89, 105, 60, 148, 59, 73, 83, 49, 137, 91)

# Fitting the model
library(DiscreteDists)
library(gamlss)

mod2 <- gamlss(y~1, sigma.fo=~1, family=DMOLBE,
               control=gamlss.control(n.cyc=500))

summary(mod2)

# Extracting the fitted values for mu and sigma
# using the inverse link function
mu_hat <- exp(coef(mod2, what="mu"))
mu_hat
sigma_hat <- exp(coef(mod2, what="sigma"))
sigma_hat

# Create a grid of values for mu and sigma around the true values
mu_vals     <- seq(30, 40, length.out=30)
sigma_vals  <- seq(0.1, 1.5, length.out=30)

# Define the log-likelihood function (using sigma, not sigma^2)
log_lik_DMOLBE <- function(mu, sigma, data) {
  ll <- sum(dDMOLBE(x=data, mu=mu, sigma=sigma, log=TRUE))
  return(ll)
}

# Evaluate the log-likelihood on the grid
# We use Vectorize to easily fill the matrix over the grid elements
log_lik_vectorized <- function(m, s) log_lik_DMOLBE(m, s, y)
z_matrix <- outer(mu_vals, sigma_vals, Vectorize(log_lik_vectorized))

# Plot the 3D surface using persp
# Adjust theta and phi to rotate the 3D angle view
persp(
  x=mu_vals, 
  y=sigma_vals, 
  z=z_matrix,
  theta=35,          # Azimuthal direction (rotation around z-axis)
  phi=30,            # Colatitude direction (tilt angle)
  expand=0.6,        # Z-axis expansion factor
  col="lightblue",   # Color of the surface facets
  border="black",    # Color of the line grid
  ltheta=120,        # Light source direction
  shade=0.4,         # Shadow intensity for a 3D effect
  ticktype = "detailed", 
  xlab="mu",
  ylab="sigma",
  zlab="Log-Likelihood",
  main="Log-Likelihood Surface for DMOLBE distribution (mu, sigma)"
)

fit <- fitdist(data=y, distr="DMOLBE", 
               start = list(mu=2, sigma=3))
print(fit)


