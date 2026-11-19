
chisq_test <- function(y, family, min_expected = 5) {
  
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("Package 'gamlss' is required.")
  }
  if (!"package:gamlss" %in% search()) {
    suppressMessages(library(gamlss))
  }
  
  FA <- gamlss.dist::as.gamlss.family(family)
  fam_name <- FA$family[1]
  
  fit <- gamlss::histDist(y = y, family = FA, plot = FALSE)
  
  if (isFALSE(fit$converged)) {
    warning("Model did not converge for family: ", fam_name)
  }
  
  par_names  <- fit$parameters
  par_values <- sapply(par_names, function(p) fitted(fit, what = p)[1])
  names(par_values) <- par_names
  m <- length(par_names)   # number of estimated parameters -> needed for correct df
  
  d_fun_name <- paste0("d", fam_name)
  p_fun_name <- paste0("p", fam_name)
  
  if (!exists(d_fun_name, mode = "function")) {
    stop("Function '", d_fun_name, "' not found.")
  }
  if (!exists(p_fun_name, mode = "function")) {
    stop("Function '", p_fun_name, "' not found.")
  }
  
  d_fun <- get(d_fun_name, mode = "function")
  p_fun <- get(p_fun_name, mode = "function")
  
  n <- length(y)
  vmin <- min(y); vmax <- max(y)
  vals <- vmin:vmax
  dens_vals <- do.call(d_fun, c(list(vals), as.list(par_values)))
  
  obs_class <- c(); exp_class <- c()
  obs_acc <- 0; exp_acc <- 0
  
  for (i in seq_along(vals)) {
    obs_acc <- obs_acc + sum(y == vals[i])
    exp_acc <- exp_acc + n * dens_vals[i]
    if (exp_acc >= min_expected) {
      obs_class <- c(obs_class, obs_acc)
      exp_class <- c(exp_class, exp_acc)
      obs_acc <- 0; exp_acc <- 0
    }
  }
  
  if (length(obs_class) == 0) {
    stop("No class reached the minimum expected frequency (", min_expected,
         "). Sample size may be too small for this test.")
  }
  
  # merge the trailing incomplete class into the last completed one
  if (exp_acc > 0) {
    obs_class[length(obs_class)] <- obs_class[length(obs_class)] + obs_acc
    exp_class[length(exp_class)] <- exp_class[length(exp_class)] + exp_acc
  }
  
  # tail mass outside [vmin, vmax] -- MUST use the CDF, not the density
  p_below <- do.call(p_fun, c(list(vmin - 1), as.list(par_values)))
  p_above <- 1 - do.call(p_fun, c(list(vmax), as.list(par_values)))
  exp_class[1] <- exp_class[1] + n * p_below
  exp_class[length(exp_class)] <- exp_class[length(exp_class)] + n * p_above
  
  k  <- length(obs_class)
  df <- k - 1 - m   # k classes - 1 (sum = n constraint) - m (estimated parameters)
  
  if (df <= 0) {
    warning("Non-positive degrees of freedom (k = ", k, ", m = ", m,
            "). Test result is not valid; consider a larger sample, ",
            "a lower min_expected, or fewer estimated parameters.")
  }
  
  chi_stat <- sum((obs_class - exp_class)^2 / exp_class)
  p_value  <- if (df > 0) 1 - pchisq(chi_stat, df) else NA_real_
  
  result <- list(
    family     = fam_name,
    parameters = par_values,
    table      = data.frame(class = 1:k, observed = obs_class,
                            expected = round(exp_class, 2)),
    statistic  = chi_stat,
    df         = df,
    p.value    = p_value,
    converged  = fit$converged,
    AIC        = AIC(fit)
  )
  
  class(result) <- "chisq_gof"
  return(result)
}

print.chisq_gof <- function(x, ...) {
  cat("Chi-squared goodness-of-fit test\n")
  cat("Family:", x$family, "\n")
  cat("Number of classes (k):", nrow(x$table), "\n")
  cat("Estimated parameters:",
      paste(names(x$parameters), round(x$parameters, 4),
            sep = " = ", collapse = ", "), "\n")
  cat("Chi-squared statistic:", round(x$statistic, 3), "\n")
  cat("df:", x$df, "\n")
  cat("p-value:", signif(x$p.value, 4), "\n")
  invisible(x)
}

