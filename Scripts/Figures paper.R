library(DiscreteDists)

# Figure pmf --------------------------------------------------------------

pdf("Figs/pmf.pdf", height=12, width=12)

par(mfrow=c(3, 3))

# Part A
mu_A <- 0.5
mu <- mu_A
sigma <- 0.5
x <- 0:15
plot(x=x, y=dDMOLBE(x, mu, sigma), ylim=c(0, 0.8),
     xlab="y", ylab="P(Y=y)", type="h", lwd=2,
     main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))
  
mu <- mu_A
sigma <- 1
x <- 0:15
plot(x=x, y=dDMOLBE(x, mu, sigma), ylim=c(0, 0.8),
     xlab="y", ylab="P(Y=y)", type="h", lwd=2,
     main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))

mu <- mu_A
sigma <- 3
x <- 0:15
plot(x=x, y=dDMOLBE(x, mu, sigma), ylim=c(0, 0.8),
     xlab="y", ylab="P(Y=y)", type="h", lwd=2,
     main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))

# Part B
mu_B <- 1
mu <- mu_B
sigma <- 0.5
x <- 0:15
plot(x=x, y=dDMOLBE(x, mu, sigma), ylim=c(0, 0.8),
     xlab="y", ylab="P(Y=y)", type="h", lwd=2,
     main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))

mu <- mu_B
sigma <- 1
x <- 0:15
plot(x=x, y=dDMOLBE(x, mu, sigma), ylim=c(0, 0.8),
     xlab="y", ylab="P(Y=y)", type="h", lwd=2,
     main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))

mu <- mu_B
sigma <- 3
x <- 0:15
plot(x=x, y=dDMOLBE(x, mu, sigma), ylim=c(0, 0.8),
     xlab="y", ylab="P(Y=y)", type="h", lwd=2,
     main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))

# Part C
mu_C <- 3
mu <- mu_C
sigma <- 0.5
x <- 0:15
plot(x=x, y=dDMOLBE(x, mu, sigma), ylim=c(0, 0.8),
     xlab="y", ylab="P(Y=y)", type="h", lwd=2,
     main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))

mu <- mu_C
sigma <- 1
x <- 0:15
plot(x=x, y=dDMOLBE(x, mu, sigma), ylim=c(0, 0.8),
     xlab="y", ylab="P(Y=y)", type="h", lwd=2,
     main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))

mu <- mu_C
sigma <- 3
x <- 0:15
plot(x=x, y=dDMOLBE(x, mu, sigma), ylim=c(0, 0.8),
     xlab="y", ylab="P(Y=y)", type="h", lwd=2,
     main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))


dev.off()

