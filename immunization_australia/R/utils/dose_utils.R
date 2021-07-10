# Utility functions to calculate the age-specific
# vaccination coverage (accounting for strategy
# and uptake)



vaccination_strategies = c("Uniform",
                           "Vulnerable",
                           "Transmitters")


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


calc_age_specific_coverage = function(target_coverage,
                                      params,
                                      strategy = "Uniform",
                                      uptake = 1.0) {
  if (is.numeric(strategy)) {
    age_group_priorities = strategy
    names(age_group_priorities) = age_groups
  }
  else if (is.character(strategy) || is.factor(strategy)) {
    if (strategy == "Uniform") {
      # Redistribute doses across only >10 year olds
      effective_coverage = target_coverage * sum(params$population_size) / sum(params$population_size[3:16])
      age_specific_coverage = c(0, 0, rep(min(c(
        uptake, effective_coverage, 1
      )), 14))
      names(age_specific_coverage) = age_groups
      return(age_specific_coverage)
    }
    else if (strategy == "Vulnerable") {
      age_group_priorities = c(
        #"0-4" = Inf,
        #"5-9" = Inf,
        "10-14" = 3,
        "15-19" = 3,
        "20-24" = 3,
        "25-29" = 3,
        "30-34" = 2,
        "35-39" = 2,
        "40-44" = 2,
        "45-49" = 2,
        "50-54" = 2,
        "55-59" = 1,
        "60-64" = 1,
        "65-69" = 1,
        "70-74" = 1,
        "75+" = 1
      )
    }
    else if (strategy == "Transmitters") {
      age_group_priorities = c(
        #"0-4" = NA,
        #"5-9" = NA,
        "10-14" = 1,
        "15-19" = 1,
        "20-24" = 1,
        "25-29" = 1,
        "30-34" = 1,
        "35-39" = 1,
        "40-44" = 1,
        "45-49" = 1,
        "50-54" = 1,
        "55-59" = 2,
        "60-64" = 2,
        "65-69" = 2,
        "70-74" = 2,
        "75+" = 2
      )
    }
  }
  # Initialize age-specific coverage vector
  age_specific_coverage = rep(0, 16)
  names(age_specific_coverage) = age_groups
  
  # Initialize priority group
  priority = 0
  
  with(params,
       {
         # Step 1. Convert target population-level coverage to number of available doses
         doses_available = round(target_coverage * sum(population_size))
         
         # Whilst available doses remain work through remaining priority groups
         while (doses_available > 0 &&
                priority < max(age_group_priorities, na.rm = TRUE)) {
           # Set current priority group
           priority = priority + 1
           priority_age_groups = names(age_group_priorities[age_group_priorities == priority])
           priority_group_size = sum(population_size[priority_age_groups])
           
           # Calculate number of doses required to achieve full uptake in current priority group
           doses_needed = uptake * priority_group_size
           
           # If enough doses are available, set coverage to full uptake
           if (doses_available >= doses_needed) {
             # print("Enough doses available for current group")
             age_specific_coverage[priority_age_groups] = uptake
             
           }
           # If not, maximize coverage evenly across currently prioritized age groups using available doses
           else{
             # print("Not enough doses available for current group")
             age_specific_coverage[priority_age_groups] = doses_available / priority_group_size
           }
           
           # Calculate the number of doses remaining after latest distribution
           doses_available = doses_available - doses_needed
         }
         return(age_specific_coverage)
         
       })
}




calc_age_specific_coverage_children_included = function(target_coverage,
                                                        params,
                                                        strategy = "Uniform",
                                                        uptake = 1.0) {
  if (is.numeric(strategy)) {
    age_group_priorities = strategy
    names(age_group_priorities) = age_groups
  }
  else if (is.character(strategy) || is.factor(strategy)) {
    if (strategy == "Uniform") {
      # Redistribute doses across only >10 year olds
      age_specific_coverage = c(rep(min(c(
        uptake, target_coverage, 1
      )), 16))
      names(age_specific_coverage) = age_groups
      return(age_specific_coverage)
    }
    else if (strategy == "Vulnerable") {
      age_group_priorities = c(
        "0-4" = 3,
        "5-9" = 3,
        "10-14" = 3,
        "15-19" = 3,
        "20-24" = 3,
        "25-29" = 3,
        "30-34" = 2,
        "35-39" = 2,
        "40-44" = 2,
        "45-49" = 2,
        "50-54" = 2,
        "55-59" = 1,
        "60-64" = 1,
        "65-69" = 1,
        "70-74" = 1,
        "75+" = 1
      )
    }
    else if (strategy == "Transmitters") {
      age_group_priorities = c(
        "0-4" = 2,
        "5-9" = 1,
        "10-14" = 1,
        "15-19" = 1,
        "20-24" = 1,
        "25-29" = 1,
        "30-34" = 1,
        "35-39" = 1,
        "40-44" = 1,
        "45-49" = 1,
        "50-54" = 1,
        "55-59" = 2,
        "60-64" = 2,
        "65-69" = 2,
        "70-74" = 2,
        "75+" = 2
      )
    }
  }
  # Initialize age-specific coverage vector
  age_specific_coverage = rep(0, 16)
  names(age_specific_coverage) = age_groups
  
  # Initialize priority group
  priority = 0
  
  with(params,
       {
         # Step 1. Convert target population-level coverage to number of available doses
         doses_available = round(target_coverage * sum(population_size))
         
         # Whilst available doses remain work through remaining priority groups
         while (doses_available > 0 &&
                priority < max(age_group_priorities, na.rm = TRUE)) {
           # Set current priority group
           priority = priority + 1
           priority_age_groups = names(age_group_priorities[age_group_priorities == priority])
           priority_group_size = sum(population_size[priority_age_groups])
           
           # Calculate number of doses required to achieve full uptake in current priority group
           doses_needed = uptake * priority_group_size
           
           # If enough doses are available, set coverage to full uptake
           if (doses_available >= doses_needed) {
             # print("Enough doses available for current group")
             age_specific_coverage[priority_age_groups] = uptake
             
           }
           # If not, maximize coverage evenly across currently prioritized age groups using available doses
           else{
             # print("Not enough doses available for current group")
             age_specific_coverage[priority_age_groups] = doses_available / priority_group_size
           }
           
           # Calculate the number of doses remaining after latest distribution
           doses_available = doses_available - doses_needed
         }
         return(age_specific_coverage)
         
       })
}
