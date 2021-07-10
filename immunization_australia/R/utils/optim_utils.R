# Utility functions to optimize vaccine distribution
# to minimize: Reff, deaths, hospitalizations or
# doses.

optimize_herd_immunity = function(params,
                                  uptake = 0.8) {
  doses_required = function(x) {
    if (any(x < 0) || any(x > 1)) {
      return(Inf)
    }
    
    target_vacc_prop = c(0,
                         0,
                         rep(x[1], 1),
                         rep(x[2], 3),
                         rep(x[3], 6),
                         rep(x[4], 4))
    actual_vacc_prop = sapply(target_vacc_prop,
                              function(x)
                                min(x, uptake))
    Reff = calc_Reff(params,
                     actual_vacc_prop)
    print(paste(c("Target: ", target_vacc_prop)))
    print(paste(c("Actual: ", actual_vacc_prop)))
    print(paste(c("Reff: ", Reff)))
    if (Reff < 1) {
      return(sum(target_vacc_prop * params$population_size))
    } else {
      return(Inf)
    }
  }
  
  optimal_solution = optim(rep(1, 4),
                           doses_required,
                           method = "SANN")
  
  return(optimal_solution)
}


