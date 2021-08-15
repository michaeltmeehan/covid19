# Utility functions for generating sampled parameter values

library(lhs)
library(EnvStats)

# Parameters varied through sensitivity analysis

## Parameter                /  Dimension  /  Distribution

#  Relative infectiousness  /       1     /   T(0.15, 0.25, 1.0)
#  Pfizer VE (overall)      /       1     /   T(0.782, 0.879, 0.932)
#  AstraZeneca VE (overall) /       1     /   T(0.289, 0.598, 0.773)
#  odds ratio hospital.     /       1     /   T(1.8, 2.08, 2.38)
#  odds ratio death         /       1     /   T(1.47, 2.32, 3.3)


# Read in baseline estimates of parameters
strain = "Delta"

epi = import.epi.data(strain)
efficacy = import.efficacy.data(strain)

baseline_params = get.params()

generate.samples = function(n_samples = 10) {
  samples = lhs::randomLHS(n = n_samples, k = 5)
  colnames(samples) = c(
    "rel_infectious",
    "Pfizer_VE",
    "AstraZeneca_VE",
    "or_hospitalization",
    "or_death"
  )
  
  # Transform sampled parameters using assumed disttribution
  samples[, "rel_infectious"] = qtri(samples[, "rel_infectious"],
                                     min = 0.15,
                                     max = 1.0,
                                     mode = 0.25)
  samples[, "Pfizer_VE"] = qtri(samples[, "Pfizer_VE"],
                                min = 0.782,
                                max = 0.932,
                                mode = 0.879)
  samples[, "AstraZeneca_VE"] = qtri(samples[, "AstraZeneca_VE"],
                                     min = 0.289,
                                     max = 0.773,
                                     mode = 0.598)
  samples[, "or_hospitalization"] = qtri(samples[, "or_hospitalization"],
                                         min = 1.8,
                                         max = 2.38,
                                         mode = 2.08)
  samples[, "or_death"] = qtri(samples[, "or_death"],
                               min = 1.47,
                               max = 3.3,
                               mode = 2.32)
  
  # Back-calculate Vacq and Vsev from VE
  Pfizer_Vacq = runif(n_samples, 
                      min = 0.25 * samples[, "Pfizer_VE"], 
                      max = samples[, "Pfizer_VE"])
  Pfizer_Vsev = 1 - (1 - samples[, "Pfizer_VE"]) / (1 - Pfizer_Vacq)
  
  # VE = 1 - (1 - Vacq) * (1 - Vsev)
  # 1 - VE = (1 - Vacq) * (1 - Vsev)
  # 1 - Vsev = (1 - VE) / (1 - Vacq)
  # Vsev = 1 - (1 - VE) / (1 - Vacq)
  
  AstraZeneca_Vacq = runif(n_samples, 
                           min = 0.25 * samples[, "AstraZeneca_VE"],
                           max = samples[, "AstraZeneca_VE"])
  AstraZeneca_Vsev = 1 - (1 - samples[, "AstraZeneca_VE"]) / (1 - AstraZeneca_Vacq)
  
  samples = cbind(samples,
                  Pfizer_Vacq,
                  Pfizer_Vsev,
                  AstraZeneca_Vacq,
                  AstraZeneca_Vsev)
  
  return(data.frame(samples))
}


merge.base.params = function(params) {
  # print(params)
  
  base_params = baseline_params
  
  base_params$rel_infectious = params["rel_infectious"]
  
  base_params$efficacy[1, c(2, 3)] =
    c(params["Pfizer_Vacq"],
      params["Pfizer_Vsev"])
  
  base_params$efficacy[2, c(2, 3)] =
    c(params["AstraZeneca_Vacq"],
      params["AstraZeneca_Vsev"])
  
  base_params$hospitalization_rate = transform.odds(epi$hospitalization_rate, params["or_hospitalization"])
  base_params$infection_fatality_rate = transform.odds(epi$infection_fatality_rate, params["or_death"])
  
  # base_params$hospitalization_rate = base_params$hospitalization_rate / base_params$rr_hospitalization * params["rr_hospitalization"]
  # base_params$infection_fatality_rate = base_params$infection_fatality_rate / base_params$rr_death * params["rr_death"]
  
  base_params$or_hospitalization = unname(params["or_hospitalization"])
  base_params$or_death = unname(params["or_death"])
  
  # print(base_params)
  
  # With parameters now modified, we need to recalibrate the scale factor
  # base_params$scale_factor = base_params$scale_factor / baseline_params$R0 * base_params$R0
  
  base_params$scale_factor = 1
  base_params$scale_factor = calc.scale.factor(base_params, params["R0"])
  
  return(base_params)
  
}



calc.prcc = function(out,
                     outcome = "Deaths",
                     program = "Mix",
                     strategy = "Vulnerable") {
  prcc = out %>%
    dplyr::filter(Outcome == outcome, Program == program, Strategy == Strategy) %>%
    dplyr::select(rel_infectious,
                  Pfizer_VE,
                  AstraZeneca_VE,
                  or_hospitalization,
                  or_death,
                  Total) %>%
    epiR::epi.prcc()
  
  return(cbind(
    Outcome = rep(outcome, nrow(prcc)),
    Program = rep(program, nrow(prcc)),
    Strategy = rep(strategy, nrow(prcc)),
    Parameter = c(
      "rel_infectious",
      "Pfizer_VE",
      "AstraZeneca_VE",
      "or_hospitalization",
      "or_death"
    ),
    prcc
  ))
  
}
