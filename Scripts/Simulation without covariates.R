require(DiscreteDists)
require(gamlss)

# To perform the simulation -----------------------------------------------

library("parSim")

parSim(
  ### SIMULATION CONDITIONS
  n = c(50, 100, 150, 200, 300),
  mu = c(0.5, 1, 1.5),
  sigma = c(0.5, 1.5),
  
  reps = 100,           # repetitions
  write = TRUE,          # Writing to a file
  name = "res02",  # Name of the file
  nCores = 1,            # Number of cores to use
  
  expression = {
    # True parameter values
    y <- rDMOLBE(n=n, mu, sigma)
    
    library(gamlss)
    mod <- gamlss(y~1, family=DMOLBE,
                  control=gamlss.control(n.cyc=1000, trace=TRUE))
    
    mu_hat    <- exp(coef(mod, what="mu"))
    sigma_hat <- exp(coef(mod, what="sigma"))
    
    # Results list:
    Results <- list(
      mu_hat = mu_hat,
      sigma_hat = sigma_hat
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------

datos1 <- read.table("Results_simuls/res01.txt", header=TRUE)
datos2 <- read.table("Results_simuls/res02.txt", header=TRUE)

datos <- rbind(datos1, datos2)

datos$case <- with(datos, ifelse(mu==0.5 & sigma==0.5, 1, 
                   ifelse(mu==0.5 & sigma==1.5, 2,
                          ifelse(mu==1.0 & sigma==0.5, 3,
                                 ifelse(mu==1.0 & sigma==1.5, 4, 
                                        ifelse(mu==1.5 & sigma==0.5, 5, 6))))))

datos$case <- as.factor(datos$case)

# To analize the results --------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

trim <- 0.2

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


p1 <- ggplot(dat, aes(x=n, y=bias_mu, colour=case)) +
  geom_line() + 
  ylab(expression(paste("Bias for ", mu)))

p2 <- ggplot(dat, aes(x=n, y=bias_si, colour=case)) +
  geom_line() + 
  ylab(expression(paste("Bias for ", sigma)))

p1 + p2

ggsave(filename="Figs/bias_simul1.pdf", width=12, height=6,
       plot=p1+p2)


p3 <- ggplot(dat, aes(x=n, y=mse_mu, colour=case)) +
  geom_line() + 
  ylab(expression(paste("MSE for ", mu)))

p4 <- ggplot(dat, aes(x=n, y=mse_si, colour=case)) +
  geom_line() + 
  ylab(expression(paste("MSE for ", sigma)))

p3 + p4

ggsave(filename="Figs/mse_simul1.pdf", width=12, height=6,
       plot=p3+p4)


