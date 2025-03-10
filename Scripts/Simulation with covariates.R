require(DiscreteDists)
require(gamlss)

# To perform the simulation -----------------------------------------------

library("parSim")

gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(2.21 + 2.91 * x1) # 39 approximately
  sigma <- exp(1.26 - 5.14 * x2) # 0.27 approximately
  y <- rDMOLBE(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

estim_mu_sigma_DMOLBE <- function(y) {
  mod <- nlminb(start=c(0, 0),
               objective=logLik_DMOLBE,
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}

estim_mu_sigma_DMOLBE(dat$y)

dat <- gendat(n=350)

mod <- gamlss(y~x1, sigma.fo=~x2, family=DMOLBE, data=dat,
              control=gamlss.control(n.cyc=500, trace=FALSE))

coef(mod, what="mu")
coef(mod, what="sigma")


parSim(
  ### SIMULATION CONDITIONS
  n = c(50, 100, 150, 200, 300),
  
  reps = 10,                      # repetitions
  write = TRUE,                    # Writing to a file
  name = "res_with_covariates_01",  # Name of the file
  nCores = 1,                      # Number of cores to use
  
  expression = {
    # True parameter values
    dat <- gendat(n=n)
    
    library(gamlss)
    mod <- gamlss(y~x1, sigma.fo=~x2, family=DMOLBE, data=dat,
                  control=gamlss.control(n.cyc=500, trace=FALSE))
    
    beta_0_hat  <- coef(mod, what="mu")[1]
    beta_1_hat  <- coef(mod, what="mu")[2]
    gamma_0_hat <- coef(mod, what="sigma")[1]
    gamma_1_hat <- coef(mod, what="sigma")[2]
    
    # Results list:
    Results <- list(
      beta_0_hat = beta_0_hat,
      beta_1_hat = beta_1_hat,
      gamma_0_hat = gamma_0_hat,
      gamma_1_hat = gamma_1_hat
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------

datos1 <- read.table("Simuls/res_with_covariates_01.txt", header=TRUE)

prop.table(table(datos1$error == ))

View(datos1)


datos <- rbind(datos1, datos2, datos3)



# To analize the results --------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

trim <- 0.2 # percentage of values to be trimmed

dat <- datos %>% group_by(n, mu, sigma, case) %>% 
  summarise(nobs = n(),
            mean_mu = mean(mu_hat, trim=trim, na.rm=TRUE),
            mean_si = mean(sigma_hat, trim=trim, na.rm=TRUE),
            mse_mu = mean((mu_hat - mu)^2, trim=trim, na.rm=TRUE),
            mse_si = mean((sigma_hat - sigma)^2, trim=trim, na.rm=TRUE),
            bias_mu = mean(mu_hat-mu, trim=trim, na.rm=TRUE),
            bias_si = mean(sigma_hat-sigma, trim=trim, na.rm=TRUE),
  )

dat

# Plots


