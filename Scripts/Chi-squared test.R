# Example 1 ---------------------------------------------------------------

# Data Set I (death due to coronavirus in China). The first data set is the number
# of deaths due to coronavirus in China from 23 January to 28 March.
# The data sets used in the paper was collected from 2020 year. The data set
# is reported in https://www.worldometers.info/coronavirus/country/china/.
# The data are:

y <- c(8, 16, 15, 24, 26, 26, 38, 43, 46, 45, 57, 64, 65, 73, 73, 86, 89, 97, 108, 97,
       146, 121, 143, 142, 105, 98, 136, 114, 118, 109, 97, 150, 71, 52, 29, 44, 47,
       35, 42, 31, 38, 31, 30, 28, 27, 22, 17, 22, 11, 7, 13, 10, 14, 13, 11, 8, 3, 7,
       6, 9, 7, 4, 6, 5, 3, 5)

library(fitdistrplus)

fit <- gof <- NULL
fit <- fitdist(y, "DMOLBE", start=list(mu=2, sigma=3))
fit
gof <- gofstat(fit, discrete=TRUE)
gof[c("chisq", "chisqpvalue")]

fit <- gof <- NULL
fit <- fitdist(y, "PO", start=list(mu=2))
fit
gof <- gofstat(fit, discrete=TRUE)
gof[c("chisq", "chisqpvalue")]

fit <- gof <- NULL
fit <- fitdist(y, "DBH", start=list(mu=0.5))
fit
gof <- gofstat(fit, discrete=TRUE)
gof[c("chisq", "chisqpvalue")]

fit <- gof <- NULL
fit <- fitdist(y, "COMPO2", start=list(mu=10, sigma=0))
fit
gof <- gofstat(fit, discrete=TRUE)
gof[c("chisq", "chisqpvalue")]


fit <- gof <- NULL
fit <- fitdist(y, "NBI", start=list(mu=0.5, sigma=0.5))
fit
gof <- gofstat(fit, discrete=TRUE)
gof[c("chisq", "chisqpvalue")]

library(DiscreteDists)
library(gamlss)

chisq_test(y, family=DMOLBE)
chisq_test(y, family=DBH)
chisq_test(y, family=PO)
chisq_test(y, family=NBI)

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

library(DiscreteDists)
library(gamlss)

chisq_test(y, family=DMOLBE)
chisq_test(y, family=DBH)
chisq_test(y, family=PO)
chisq_test(y, family=NBI)

