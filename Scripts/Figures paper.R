library(DiscreteDists)

# Figure pmf --------------------------------------------------------------

pdf("Figs/pmf.pdf", height=4, width=12)

par(mfrow=c(1, 3))

mu <- 0.5
sigma <- 0.5
x <- 0:10
plot(x=x, y=dDMOLBE(x, mu, sigma), 
       xlab="y", ylab="P(Y=y)", type="h", lwd=2,
       main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))
  
mu <- 5
sigma <- 1
x <- 0:25
plot(x=x, y=dDMOLBE(x, mu, sigma), 
     xlab="y", ylab="P(Y=y)", type="h", lwd=2,
     main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))

mu <- 5
sigma <- 15
x <- 0:50
plot(x=x, y=dDMOLBE(x, mu, sigma), 
     xlab="y", ylab="P(Y=y)", type="h", lwd=2,
     main=bquote(paste(mu, "=", .(mu), ", ", sigma, "=", .(sigma))))


dev.off()

