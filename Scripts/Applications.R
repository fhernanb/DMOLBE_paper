
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

mod1 <- gamlss(y~1, sigma.fo=~1, family=DMOLBE)

summary(mod1)

# Extracting the fitted values for mu and sigma
# using the inverse link function
mu_hat <- exp(coef(mod1, what="mu"))
mu_hat
sigma_hat <- exp(coef(mod1, what="sigma"))
sigma_hat

# Some measures
logLik(mod1)
AIC(mod1)
AIC(mod1, k=log(length(y)))

# Competitive models
mod12 <- gamlss(y~1, sigma.fo=~1, family=DBH)
inv_logit <- function(x) exp(x) / (1+exp(x))
inv_logit(coef(mod12, what="mu"))
logLik(mod12)
AIC(mod12)

mod13 <- gamlss(y~1, sigma.fo=~1, family=PO)
exp(coef(mod13, what="mu"))
logLik(mod13)
AIC(mod13)

mod14 <- gamlss(y~1, sigma.fo=~1, family=NBI)
exp(coef(mod14, what="mu"))
exp(coef(mod14, what="sigma"))
logLik(mod14)
AIC(mod14)

# Comparisons

pdf("Figs/cdf_app1.pdf", height=6, width=6)

emp_cumulative <- ecdf(y)
plot(emp_cumulative, 
     xlab="y", ylab="CDF",
     main="CDF for application 1")

prob_hat <- pDMOLBE(q=sort(y), mu=mu_hat, sigma=sigma_hat)
points(x=sort(y), y=prob_hat, col="orange", pch=19)

prob_hat <- pDBH(q=sort(y), mu=1)
points(x=sort(y), y=prob_hat, col="tomato", pch=19)

prob_hat <- pPO(q=sort(y), mu=49.74242)
points(x=sort(y), y=prob_hat, col="blue3", pch=19)

legend("bottomright", bty="n", pch=19,
       legend=c("Empirical", "DMOLBE", "DHB", "POISSON"),
       lty=1, col=c("black", "orange", "tomato", "blue3"))

dev.off()




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
               control=gamlss.control(n.cyc=30, trace=TRUE))

summary(mod2)

# Extracting the fitted values for mu and sigma
# using the inverse link function
mu_hat <- exp(coef(mod2, what="mu"))
mu_hat
sigma_hat <- exp(coef(mod2, what="sigma"))
sigma_hat

# Some measures
logLik(mod2)
AIC(mod2)
AIC(mod2, k=log(length(y)))


# Competitive models
mod22 <- gamlss(y~1, sigma.fo=~1, family=DBH)
inv_logit <- function(x) exp(x) / (1+exp(x))
inv_logit(coef(mod22, what="mu"))
logLik(mod22)
AIC(mod22)

mod23 <- gamlss(y~1, sigma.fo=~1, family=PO)
exp(coef(mod23, what="mu"))
logLik(mod23)
AIC(mod23)

# Comparisons

pdf("Figs/cdf_app2.pdf", height=6, width=6)

emp_cumulative <- ecdf(y)
plot(emp_cumulative, 
     xlab="y", ylab="CDF",
     main="CDF for application 2")

prob_hat <- pDMOLBE(q=sort(y), mu=mu_hat, sigma=sigma_hat)
points(x=sort(y), y=prob_hat, col="orange", pch=19)

prob_hat <- pDBH(q=sort(y), mu=1)
points(x=sort(y), y=prob_hat, col="tomato", pch=19)

prob_hat <- pPO(q=sort(y), mu=50.057)
points(x=sort(y), y=prob_hat, col="blue3", pch=19)


legend("bottomright", bty="n", pch=19,
       legend=c("Empirical", "DMOLBE", "DHB", "POISSON"),
       lty=1, col=c("black", "orange", "tomato", "blue3"))

dev.off()




