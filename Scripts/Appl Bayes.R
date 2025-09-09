
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
library(bamlss)

mod1 <- bamlss(y~1, sigma.fo=~1, family=DMOLBE)

summary(mod1)

plot(mod1)

plot(mod1, which = "samples")

plot(mod1, which = c("hist-resid", "qq-resid"))

plot(b)
plot(b, which = 3:4)
plot(b, which = "samples")

