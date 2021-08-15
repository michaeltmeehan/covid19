# Utility functions to calculate the final size of an epidemic, the total number of hospitalizations and deaths

library(nleqslv)

source("./utils/NGM_utils.R")
source("./utils/dose_utils.R")


# Calculate the final size using the final size equation
calc.final.size =  function(base_params, vacc_params) {
  NGM = calc.NGM(base_params, vacc_params)
  with(c(base_params, vacc_params), {
    coverage <- coverage + rep(0.000001, 16)
    pop_size = unlist(c((1 - coverage) * population_size, coverage * population_size))
    #print(pop_size)
    f = function(x) {
      x - exp(- (diag(1/pop_size) %*% NGM %*% diag(pop_size)) %*% (1 - x))
    }
    J = function(x) {
      diag(32) - diag(as.vector(exp(-NGM %*% (1 - x)))) %*% NGM
    }
    
    final_size = 1 - nleqslv(rep(0.05, 32), 
                             f, 
                             jac = J, 
                             global="cline",
                             control = list(ftol=1e-15,
                                            xtol=1e-15))$x
    #print(sum(final_size*pop_size/(sum(pop_size))))
    return(sapply(final_size, function(x){max(x, 0)}))
  }
  )
}



# Calculate the number of infections, hospitalizations and deaths
calc.final.burden = function(base_params,
                             vacc_params) {
  final_size = calc.final.size(base_params,
                               vacc_params)
  
  with(c(base_params, vacc_params), {
    pop_size = unlist(c((1 - coverage) * population_size,
                        coverage * population_size
    ))
    
    Infections = final_size * pop_size
    Hospitalizations = hospitalization_rate * (1 - seropositivity) *
      c(final_size[1:16] * pop_size[1:16],
        (1 - Vsev) * (1 - Vmor) * final_size[17:32] * pop_size[17:32])
    Deaths = infection_fatality_rate * (1 - seropositivity) *
      c(final_size[1:16] * pop_size[1:16],
        (1 - Vsev) * (1 - Vmor) * final_size[17:32] * pop_size[17:32])
    YLL = Deaths * rep(life_expectancy - seq(2.5, 77.5, 5), times=2)
    
    return(as.data.frame(rbind(Infections,
                               Hospitalizations,
                               Deaths,
                               YLL)))
  })
  
}



calc.segregated.totals = function(burden) {
  
  Outcome = rownames(burden)
  rownames(burden) = NULL
  
  Unvaccinated = rowSums(burden[,1:16])
  Vaccinated = rowSums(burden[,17:32])
  Total = Unvaccinated + Vaccinated
  
  
  return(data.frame(Outcome = Outcome,
                    Unvaccinated = Unvaccinated,
                    Vaccinated = Vaccinated,
                    Total = Total))
}

