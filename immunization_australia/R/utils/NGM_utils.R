

# Calculate Next-Generation-Matrix
calc.NGM = function(base_params, 
                    vacc_params = list(coverage = rep(0,16),
                                       Vacq = rep(0,16),
                                       Vsev = rep(0, 16),
                                       Vtrans = rep(0,16),
                                       Vmor = rep(0,16))) {
  with(as.list(c(base_params, vacc_params)), {
    pre_u = scale_factor * 
      susceptibility * 
      (1 - coverage) * (1 - seropositivity) * population_size
    pre_v = scale_factor * 
      susceptibility * 
      (1 - Vacq ) * 
      coverage * (1 - seropositivity) * population_size
    post_u = (clinical_fraction * (preclinical_period + clinical_period) + 
                (1 - clinical_fraction) * rel_infectious * asymptomatic_period) /
      (population_size)
    post_v = ((1 - Vsev) * clinical_fraction * (preclinical_period + clinical_period) + 
                (1 - (1 - Vsev) * clinical_fraction) * rel_infectious * asymptomatic_period) *
      (1 - Vtrans) / (population_size)
    
    NGM_uu = diag(pre_u) %*% contact_matrix %*% diag(post_u)
    NGM_vu = diag(pre_v) %*% contact_matrix %*% diag(post_u)
    NGM_uv = diag(pre_u) %*% contact_matrix %*% diag(post_v)
    NGM_vv = diag(pre_v) %*% contact_matrix %*% diag(post_v)
    
    NGM = cbind(rbind(NGM_uu, NGM_vu),rbind(NGM_uv, NGM_vv))
    
    return(NGM)
  })
}




# Calculate the effective reproduction number
calc.Reff = function(base_params,
                     vacc_params = list(coverage = rep(0,16),
                                        Vacq = rep(0,16),
                                        Vsev = rep(0, 16),
                                        Vtrans = rep(0,16),
                                        Vmor = rep(0,16))) {
  NGM = calc.NGM(base_params, vacc_params)
  
  return(max(abs(eigen(NGM, only.values=TRUE)$values)))
}



# Calculate the scale factor (i.e., the (psuedo-)probability of transmission given contact)
calc.scale.factor = function(base_params, R0) {
  base_params$seropositivity = rep(0,16)
  R0_unscaled = calc.Reff(base_params)
  
  return(R0 / R0_unscaled)
}