require(DiscreteDists)
require(gamlss)


# To perform the simulation -----------------------------------------------

library("parSim")

gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(-1.3 + 3.4 * x1) # 1.5 approximately
  sigma <- exp( 2.1 - 3.4 * x2) # 1.5 approximately
  y <- rDMOLBE(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

parSim(
  ### SIMULATION CONDITIONS
  n = seq(from=400, to=1000, by=100),
  
  reps = 800,                      # repetitions
  write = TRUE,                    # Writing to a file
  name = "Simuls/sim2_07",  # Name of the file
  nCores = 1,                      # Number of cores to use
  
  expression = {
    # True parameter values
    dat <- gendat(n=n)
    
    mod <- gamlss(y~x1, sigma.fo=~x2, family=DMOLBE, data=dat,
                  control=gamlss.control(n.cyc=1000, trace=FALSE))
    
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

archivos <- list.files(pattern = "^sim2.*\\.txt$", 
                       path="Simuls",
                       full.names = TRUE)
archivos

lista_datos <- lapply(archivos, read.table, header = TRUE, 
                      sep = "", stringsAsFactors = FALSE)
datos <- do.call(rbind, lista_datos)


prop.table(table(datos$error == TRUE))


# To analize the results --------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

trim <- 0.10

dat <- datos %>% group_by(n) %>% 
  summarise(nobs = n(),
            
            bias_b0 = mean(beta_0_hat - (-1.3), trim=trim, na.rm=TRUE),
            bias_b1 = mean(beta_1_hat - (3.4), trim=trim, na.rm=TRUE),
            bias_g0 = mean(gamma_0_hat - (2.1), trim=trim, na.rm=TRUE),
            bias_g1 = mean(gamma_1_hat - (-3.4), trim=trim, na.rm=TRUE),
            
            mse_b0 = mean((beta_0_hat - (-1.3))^2, trim=trim, na.rm=TRUE),
            mse_b1 = mean((beta_1_hat - (3.4))^2, trim=trim, na.rm=TRUE),
            mse_g0 = mean((gamma_0_hat - (2.1))^2, trim=trim, na.rm=TRUE),
            mse_g1 = mean((gamma_1_hat - (-3.4))^2, trim=trim, na.rm=TRUE)

  )

dat

# Legend and colores
leyenda <- c(expression(hat(beta)[0]), 
             expression(hat(beta)[1]), 
             expression(hat(gamma)[0]), 
             expression(hat(gamma)[1]))

colores <- c("#F8766D", "#00BA38", "#619CFF", "#F564E3")


d <- pivot_longer(data=dat, 
                  cols=c("bias_b0", "bias_b1", 
                         "bias_g0", "bias_g1"),
                  names_to="Estimator",
                  values_to="value")

# Plots
p1 <- ggplot(d, aes(x=n, y=value, colour=Estimator)) +
  geom_line() + 
  ylab("Bias") + 
  scale_color_manual(labels=leyenda,
                     values=colores)

p1

d <- pivot_longer(data=dat, 
                  cols=c("mse_b0", "mse_b1", 
                         "mse_g0", "mse_g1"),
                  names_to="Estimator",
                  values_to="value")

p2 <- ggplot(d, aes(x=n, y=value, colour=Estimator)) +
  geom_line() + 
  ylab("MSE") + 
  scale_color_manual(labels=leyenda,
                     values=colores)

p2

ggsave(filename="Figs/bias_mse_simul2.pdf", width=12, height=6,
       plot=p1+p2)



