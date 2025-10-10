# This script contains the analysis of E(Y) and Var(Y)


# Main functions ----------------------------------------------------------

mean_DMOLBE <- function(mu, sigma, tol = 1e-10, max_iter = 1e6) {
  y <- 1
  mean_sum <- 0
  
  repeat {
    common <- (1 + y/mu) * exp(-y/mu)
    denom  <- 1 - (1 - sigma) * common
    term   <- sigma * common / denom
    
    mean_sum <- mean_sum + term
    
    if (abs(term) < tol) break
    
    y <- y + 1
    if (y > max_iter) {
      warning("Series did not converge in mean_DMOLBE.")
      break
    }
  }
  
  return(mean_sum)
}

mean_DMOLBE <- Vectorize(mean_DMOLBE)


var_DMOLBE <- function(mu, sigma, tol = 1e-10, max_iter = 1e6) {
  y <- 1
  mean_sum <- 0
  var_sum  <- 0
  
  repeat {
    common <- (1 + y/mu) * exp(-y/mu)
    denom  <- 1 - (1 - sigma) * common
    
    term_mean <- sigma * common / denom
    term_var  <- sigma * (2*y - 1) * common / denom
    
    mean_sum <- mean_sum + term_mean
    var_sum  <- var_sum + term_var
    
    if (abs(term_mean) < tol && abs(term_var) < tol) break
    
    y <- y + 1
    if (y > max_iter) {
      warning("Series did not converge in var_DMOLBE.")
      break
    }
  }
  
  variance <- var_sum - mean_sum^2
  return(variance)
}

var_DMOLBE <- Vectorize(var_DMOLBE)


# Comparing the theoretical values with simulated sample ------------------

mu    <- 1.3
sigma <- 1.5

mean_DMOLBE(mu=mu, sigma=sigma)
 var_DMOLBE(mu=mu, sigma=sigma)

library(DiscreteDists)
x <- rDMOLBE(n=10000, mu=mu, sigma=sigma)
mean(x)
var(x)

# Exploring E(X) and Var(X) for combinations ------------------------------

# Ratio function Var(Y)/E(Y)
ratio <- function(mu, sigma) {
  ratio <- var_DMOLBE(mu=mu, sigma=sigma) / mean_DMOLBE(mu=mu, sigma=sigma)
  ratio
}

# Possible values for mu and sigma
mus    <- seq(from=0.1, to=5, length.out=100)
sigmas <- seq(from=0.1, to=5, length.out=100)

# Exploring the mean ------------------------------------------------------

res <- outer(X=mus, Y=sigmas, FUN=mean_DMOLBE)
colnames(res) <- mus
rownames(res) <- sigmas
res

# Plot
pdf("Figs/EY.pdf", height=4, width=6)

# Contour plot
filled.contour(x=mus, y=sigmas, z=res, nlevels=50,
               main="Contour plot for E(Y)",
               xlab=expression(mu), ylab=expression(sigma),
               color = terrain.colors,
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(list(x=mus, y=sigmas, z=res), add = TRUE)
               })

dev.off()

# Exploring the variance --------------------------------------------------

res <- outer(X=mus, Y=sigmas, FUN=var_DMOLBE)
colnames(res) <- mus
rownames(res) <- sigmas
res

# Plot
pdf("Figs/VarY.pdf", height=4, width=6)

# Contour plot
filled.contour(x=mus, y=sigmas, z=res, nlevels=50,
               main="Contour plot for Var(Y)",
               xlab=expression(mu), ylab=expression(sigma),
               color = terrain.colors,
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(list(x=mus, y=sigmas, z=res), add = TRUE)
               })

dev.off()


# Exploring the ratio -----------------------------------------------------

res <- outer(X=mus, Y=sigmas, FUN=ratio)
colnames(res) <- mus
rownames(res) <- sigmas
res

# Plot
pdf("Figs/ratio.pdf", height=4, width=6)

# Contour plot
filled.contour(x=mus, y=sigmas, z=res, nlevels=20,
               main="Contour plot for the ratio Var(Y)/E(Y)",
               xlab=expression(mu), ylab=expression(sigma),
               color = terrain.colors,
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(mus, sigmas, res, levels=1, add=TRUE, 
                         col="white", lwd=3)
               })

text(x=2, y=2, "Overdispersion \n zone", 
     col="black", cex=1.2)

text(x=-0.1, y=4, "Underdispersion \n zone", 
     col="black", cex=1.2, pos=4)

dev.off()


