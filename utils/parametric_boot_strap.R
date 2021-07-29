
parametric_boot_strap <- function(data,formula,N){
  mod <- lm(data = data, formula)
  residuals <- mod$residuals
  fitted <- mod$fitted.values
  best_fit <- mod$coefficients
  dat_new <- data
  samples <- best_fit
  for(i in 1:N){
    dat_new$y <- fitted + sample(mod$residuals,replace = T)
    mod_new <- lm(data = dat_new, formula)
    samples <- rbind(samples,mod_new$coefficients)
  }
  return(list(model = mod, best_fit = best_fit, samples = samples))
}

