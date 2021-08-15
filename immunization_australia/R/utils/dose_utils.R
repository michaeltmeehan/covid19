


age_groups = c(
  "0-4",
  "5-9",
  "10-14",
  "15-19",
  "20-24",
  "25-29",
  "30-34",
  "35-39",
  "40-44",
  "45-49",
  "50-54",
  "55-59",
  "60-64",
  "65-69",
  "70-74",
  "75+"
)


calc.vacc.params = function(doses,
                            base_params) {
  
  # Extract list of vaccines included in schedule
  vaccines = colnames(doses)
  
  with(base_params, {
    # Collect the relevant vaccine efficacies
    efficacies = do.call(rbind, lapply(vaccines, function(x) {
      efficacy[efficacy$vaccine == x,]
    }))
    
    # Convert age-specific doses to coverage
    coverage = rowSums(doses) / population_size
    
    # Convert coverage to mean efficacies
    mean_efficacies = as.matrix(doses) %*% as.matrix(efficacies[, -1]) / rowSums(doses)
    
    # Clean up NaNs
    mean_efficacies[which(is.nan(mean_efficacies[, "Vacq"])),] = 0
    
    return(
      list(
        coverage = coverage,
        Vacq = mean_efficacies[, "Vacq"],
        Vsev = mean_efficacies[, "Vsev"],
        Vtrans = mean_efficacies[, "Vtrans"],
        Vmor = mean_efficacies[, "Vmor"]
      )
    )
  })
}




calc.age.specific.doses = function(
  base_params,
  coverage,
  program="Mix",
  strategy="Vulnerable",
  age_break=55,
  age_cutoff=15,
  uptake=0.9
) {
  
  idx_break = age_break / 5 + 1
  idx_cutoff = age_cutoff / 5 + 1
  
  doses = rep(0,16)
  
  with(base_params, {
    
    ## Calculate age-specific coverages
    
    # Determine size of eligible age groups
    unvaccinated = round(uptake * population_size)
    
    doses_remaining = coverage * sum(population_size)
    
    # print(doses_remaining)
    
    if (strategy == "Untargeted") {
      if (doses_remaining < sum(unvaccinated[idx_cutoff:16])) {
        doses[idx_cutoff:16] = floor(doses_remaining * population_size[idx_cutoff:16] / sum(population_size[idx_cutoff:16]))
        
      } else{
        doses[idx_cutoff:16] = unvaccinated[idx_cutoff:16]
        
      }
     # print(doses) 
    } else{
      if (strategy == "Transmitters") {
        priority_idx = idx_cutoff:(idx_break - 1)
        excess_idx = idx_break:16
      } else if (strategy == "Vulnerable") {
        priority_idx = idx_break:16
        excess_idx = idx_cutoff:(idx_break - 1)
      }
      
      if (doses_remaining < sum(unvaccinated[priority_idx])) {
        doses[priority_idx] = floor(doses_remaining * population_size[priority_idx] / sum(population_size[priority_idx]))
        
      } else {
        doses[priority_idx] = unvaccinated[priority_idx]
        
        doses_remaining = doses_remaining - sum(unvaccinated[priority_idx])
        
        if (doses_remaining < sum(unvaccinated[excess_idx])) {
          doses[excess_idx] = floor(doses_remaining * population_size[excess_idx] / sum(population_size[excess_idx]))
        } else {
          doses[excess_idx] = unvaccinated[excess_idx]
        }
      }
    }
      
    ## Assign vaccines to different age groups
    if (program == "Mix") {
      # doses = data.frame(Pfizer = c(doses[1:(idx_break-1)], rep(0, 16 - idx_break + 1)),
      #                    AstraZeneca = c(rep(0, idx_break-1), doses[idx_break:16]))
      doses = data.frame(Pfizer = c(doses[1:(idx_break)], rep(0, 16 - idx_break)),
                         AstraZeneca = c(rep(0, idx_break), doses[(idx_break+1):16]))
    } else if (program == "Pfizer") {
      doses = data.frame(Pfizer = doses)
    } else if (program == "AstraZeneca") {
      doses = data.frame(AstraZeneca = doses)
    }
    return(doses)
  })
}
