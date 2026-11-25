# Dataset

dat <- data.frame(
  attentions = c(
    2, 3, 9, 3, 22, 9, 3, 0, 2, 2, 2, 1, 9, 0, 2, 4, 2, 1, 4, 9, 1, 8, 4, 2, 3,
    3, 3, 1, 5, 4, 1, 4, 2, 6, 5, 2, 2, 1, 3, 2, 1, 13, 6, 13, 3, 3, 1, 4, 0, 1,
    4, 5, 7, 3, 3, 4, 1, 0, 0, 3, 5, 7, 2, 0, 2, 0, 11, 8, 10, 2, 6, 1, 1, 5, 3,
    11, 7, 3, 3, 3, 6, 18, 0, 1, 4, 9, 6, 2, 11, 1, 6, 9, 8, 2, 4, 4, 11, 0, 0, 1,
    7, 15, 2, 1, 2, 1, 3, 0, 9, 2, 2, 1, 3, 4, 2, 3, 5, 2, 5, 1, 0, 2, 3, 10, 11,
    5, 1, 1, 2, 5, 15, 2, 2, 0, 13, 0, 0, 4, 2
  ),
  age = c(
    65, 64, 63, 67, 80, 78, 35, 81, 37, 49, 24, 75, 56, 30, 48, 37, 24, 22, 77, 66,
    29, 60, 75, 22, 29, 76, 31, 48, 57, 71, 26, 59, 29, 71, 77, 32, 77, 49, 22, 51,
    70, 60, 76, 75, 22, 65, 53, 33, 44, 20, 67, 69, 78, 37, 45, 29, 27, 39, 40, 62,
    57, 85, 22, 22, 31, 37, 76, 80, 65, 69, 57, 54, 39, 73, 70, 46, 21, 67, 57, 53,
    53, 71, 40, 23, 60, 57, 68, 55, 80, 57, 46, 82, 47, 64, 34, 46, 54, 50, 44, 24,
    56, 71, 58, 51, 80, 31, 34, 34, 71, 39, 23, 39, 43, 50, 80, 41, 21, 42, 81, 37,
    45, 37, 62, 75, 48, 85, 42, 21, 34, 66, 68, 52, 73, 24, 83, 43, 75, 59, 41
  ),
  comorbidity = factor(
    c(
      1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1,
      1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0,
      0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1,
      1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0,
      1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1,
      0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1
    ),
    levels = c(0, 1)
  )
)

dat

# Exploratory analysis

pdf("Figs/app_hosp.pdf", height=5, width=11)
par(mfrow=c(1, 2))

with(dat, plot(
  x = age,
  y = attentions,
  pch = 20,
  xlab = "Age (years)",
  ylab = "Number of medical attentions per year",
  las = 1
))


with(dat, boxplot(
  attentions ~ comorbidity,
  pch = 20,
  names = c("Single condition (0)", "Multiple comorbidities (1)"),
  xlab = "Comorbidity status",
  ylab = "Number of medical attentions per year",
  col = "white",
  las = 1
))

dev.off()

# Fitting the models
library(gamlss)
library(DiscreteDists)

mod1 <- gamlss(attentions ~ age + comorbidity,
               sigma.formula = ~ comorbidity,
               family=DMOLBE, data=dat,
               control=gamlss.control(n.cyc=500, trace=FALSE))

mod1 <- gamlss(attentions ~ age + comorbidity,
               sigma.formula = ~ 1,
               family=DMOLBE, data=dat,
               control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod1)

mod2 <- gamlss(attentions ~ age + comorbidity,
               family=PO, data=dat)
mod3 <- gamlss(attentions ~ age + comorbidity,
               family=DBH, data=dat, 
               control=gamlss.control(n.cyc=700, trace=TRUE))
mod4 <- gamlss(attentions ~ age + comorbidity,
               sigma.formula = ~ comorbidity,
               family=COMPO2, data=dat,
               control=gamlss.control(n.cyc=700, trace=TRUE))
mod5 <- gamlss(attentions ~ age + comorbidity,
               sigma.formula = ~ comorbidity,
               family=NBI, data=dat)
mod6 <- gamlss(attentions ~ age + comorbidity,
               sigma.formula = ~ comorbidity,
               family=ZIP, data=dat)

# Metrics
AIC(mod1, mod2, mod3, mod4, mod5, mod6)

Rsq(mod1)
Rsq(mod2)
Rsq(mod3)
Rsq(mod4)
Rsq(mod5)
Rsq(mod6)

# Estimated mean for each model with simulation
library(model)
m1 <- est_param(mod1, m=10000, fun="mean")
m2 <- est_param(mod2, m=10000, fun="mean")
m3 <- est_param(mod3, m=10000, fun="mean")
m4 <- est_param(mod4, m=10000, fun="mean")
m5 <- est_param(mod5, m=10000, fun="mean")
m6 <- est_param(mod6, m=10000, fun="mean")
cor(dat$age, m1)
cor(dat$age, m2)
cor(dat$age, m3)
cor(dat$age, m4)
cor(dat$age, m5)
cor(dat$age, m6)

# worm plot
pdf("Figs/app_hosp_wp.pdf", height=8, width=13)
par(mfrow=c(2, 3))
wp(mod1)
title("DMOLBE")
wp(mod2)
title("PO")
wp(mod3)
title("DBH")
wp(mod4)
title("COMPO2")
wp(mod5)
title("NBI")
wp(mod6)
title("ZIP")
dev.off()

plot(mod1)
plot(mod2)
plot(mod3)
plot(mod4)
plot(mod5)
plot(mod6)

