# Utility functions to calculate the next-generation matrix (NGM), R0, Reff and HIT
source("utils/dose_utils.R")

# Calculate the next-generation matrix
calc_NGM = function(params, vacc_prop = rep(0,16)) {
  with(params,{
    
    pre = scale_factor * epi$susceptibility * (1 - seropositivity) * population_size
    
    post_vac = (epi$clinical_fraction * (period$preclinical + period$clinical) + (1 - epi$clinical_fraction) * rel_infectious * period$asymptomatic) *
                  (1 - vacc_prop)
    post_unvac = ( (1 - efficacy$Vsev) * epi$clinical_fraction * (period$preclinical + period$clinical) + 
                     (1 - (1 - efficacy$Vsev) * epi$clinical_fraction) * rel_infectious * period$asymptomatic ) *
                    (1 - efficacy$Vtrans) * (1 - efficacy$Vinf) * vacc_prop
    
    NGM = diag(pre) %*% contact_matrix %*% (diag(post_vac) + diag(post_unvac)) %*% diag(1 / population_size)
    return(NGM)
  })
}


# Calculate the effective reproduction number
calc_Reff = function(params, vacc_prop = rep(0, 16)) {
  NGM = calc_NGM(params, vacc_prop)
  
  return(max(abs(eigen(NGM, only.values=TRUE)$values)))
}

calc_Reff_demerged_R = function(params, vacc_prop = rep(0, 16)) {
  NGM = calc_NGM_demergered_R(params, vacc_prop)
  
  return(max(abs(eigen(NGM, only.values=TRUE)$values)))
}


# Calculate the scale factor (i.e., the (psuedo-)probability of transmission given contact)
calc_scale_factor = function(params, R0) {
  params$seropositivity = rep(0,16)
  R0_unscaled = calc_Reff(params)
  
  return(R0 / R0_unscaled)
}


# Calculate the herd immunity threshold
calc_HIT = function(params, R0 = calc_Reff(params)) {
  Reff = calc_Reff(params, rep(1, each = 16))
  
  return((R0 - 1) / (R0 - Reff))
}

# Calculate herd immunity threshold for different vaccination strategies
calc_strategy_HIT = function(params, 
                             strategy = "Untargeted",
                             uptake = 0.95,
                             eligibility_cutoff = 5) {
  
  if (eligibility_cutoff == 5) {
    max_coverage = calc_age_specific_coverage_over_5(target_coverage = 1, 
                                                     params = params, 
                                                     strategy = strategy, 
                                                     uptake = uptake)
    Rmin = calc_Reff(params, 
                     vacc_prop = max_coverage)
    print(Rmin)
    if (Rmin > 1){
      return(Inf)
    } else{
      f = function(x) {
        vacc_prop = calc_age_specific_coverage_over_5(target_coverage = x,
                                                      params = params,
                                                      strategy = strategy,
                                                      uptake = uptake)
        Reff = calc_Reff(params,
                         vacc_prop)
        return(Reff - 1)
      }
    }
  } else if (eligibility_cutoff == 10) {
    max_coverage = calc_age_specific_coverage_over_10(target_coverage = 1, 
                                                     params = params, 
                                                     strategy = strategy, 
                                                     uptake = uptake)
    Rmin = calc_Reff(params, 
                     vacc_prop = max_coverage)
    print(Rmin)
    if (Rmin > 1){
      return(Inf)
    } else{
      f = function(x) {
        vacc_prop = calc_age_specific_coverage_over_10(target_coverage = x,
                                                      params = params,
                                                      strategy = strategy,
                                                      uptake = uptake)
        Reff = calc_Reff(params,
                         vacc_prop)
        return(Reff - 1)
      }
    }
  } else if (eligibility_cutoff == 15) {
    max_coverage = calc_age_specific_coverage_over_15(target_coverage = 1, 
                                                      params = params, 
                                                      strategy = strategy, 
                                                      uptake = uptake)
    Rmin = calc_Reff(params, 
                     vacc_prop = max_coverage)
    print(Rmin)
    if (Rmin > 1){
      return(Inf)
    } else{
      f = function(x) {
        vacc_prop = calc_age_specific_coverage_over_15(target_coverage = x,
                                                       params = params,
                                                       strategy = strategy,
                                                       uptake = uptake)
        Reff = calc_Reff(params,
                         vacc_prop)
        return(Reff - 1)
      }
    }
  }
  
  
  herd_immunity_threshold = uniroot(f, c(0,1))$root
  return(herd_immunity_threshold)
}


calc_NGM_demergered_R = function(params, vacc_prop = rep(0, 16)) {
  with(params, {
    pre_u = scale_factor * 
      epi$susceptibility * 
      (1 - vacc_prop) * (1 - seropositivity) * population_size
    pre_v = scale_factor * 
      epi$susceptibility * 
      (1 - efficacy$Vinf ) * 
      vacc_prop * (1 - seropositivity) * population_size
    post_u = (epi$clinical_fraction * (period$preclinical + period$clinical) + 
              (1 - epi$clinical_fraction) * rel_infectious * period$asymptomatic) /
      (population_size)
    post_v = ((1 - efficacy$Vsev) * epi$clinical_fraction * (period$preclinical + period$clinical) + 
                (1 - (1 - efficacy$Vsev) * epi$clinical_fraction) * rel_infectious * period$asymptomatic) *
      (1 - efficacy$Vtrans) / (population_size)
    
    NGM_uu = diag(pre_u) %*% contact_matrix %*% diag(post_u)
    NGM_vu = diag(pre_v) %*% contact_matrix %*% diag(post_u)
    NGM_uv = diag(pre_u) %*% contact_matrix %*% diag(post_v)
    NGM_vv = diag(pre_v) %*% contact_matrix %*% diag(post_v)
    
    NGM = cbind(rbind(NGM_uu, NGM_vu),rbind(NGM_uv, NGM_vv))
    
    return(NGM)
  })
}
