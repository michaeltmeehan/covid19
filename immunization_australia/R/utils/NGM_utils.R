# Utility functions to calculate the next-generation matrix (NGM), R0, Reff and HIT


# Calculate the next-generation matrix
calc_NGM = function(params, vacc_prop = rep(0, 16)) {
  with(params, {
    pre = scale_factor * epi$susceptibility * (1 - efficacy$Vinf * vacc_prop) * (1 - seropositivity) * population_size
    post = ((1 - efficacy$Vsev * vacc_prop) * epi$clinical_fraction * (period$preclinical + period$clinical) + 
              (1 - (1 - efficacy$Vsev * vacc_prop) * epi$clinical_fraction) * rel_infectious * period$asymptomatic) *
          (1 - efficacy$Vtrans * vacc_prop) / population_size
    NGM = diag(pre) %*% contact_matrix %*% diag(post)
    
    return(NGM)
  })
}


# Calculate the final-size next-generation matrix
# calc_NGM_FS = function(params, vacc_prop = rep(0, 16)) {
#   with(params, {
#     NGM = calc_NGM(params, vacc_prop)
#     
#     return()
#   })
#   
# }


# Calculate the effective reproduction number
calc_Reff = function(params, vacc_prop = rep(0, 16)) {
  NGM = calc_NGM(params, vacc_prop)
  
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
    #check with Michael -- just dividing by pop size makes sense
    #rather than vacc_prop * pop_size etc?????
    NGM_uu = diag(pre_u) %*% contact_matrix %*% diag(post_u)
    NGM_vu = diag(pre_v) %*% contact_matrix %*% diag(post_u)
    NGM_uv = diag(pre_u) %*% contact_matrix %*% diag(post_v)
    NGM_vv = diag(pre_v) %*% contact_matrix %*% diag(post_v)
    
    NGM = cbind(rbind(NGM_uu, NGM_vu),rbind(NGM_uv, NGM_vv))
    
    return(NGM)
  })
}
