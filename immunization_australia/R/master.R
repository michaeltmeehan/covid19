
# Source all utility files
sapply(sapply(list.files("utils"), function(x)
  paste0("utils/", x)), source)


master.base = function(coverage,
                       R0 = 5,
                       program = "Mix",
                       strategy = "Untargeted",
                       age_cutoff = 15,
                       uptake = 0.9) {
  # Import parameters
  base_params = get.params(R0 = R0)
  
  # Generate dosing schedule
  doses = calc.age.specific.doses(base_params,
                                  coverage,
                                  program,
                                  strategy,
                                  age_break = 55,
                                  age_cutoff,
                                  uptake)
  
  # Calculate corresponding vaccination parameters
  vacc_params = calc.vacc.params(doses, base_params)
  
  # Calculate infections, hospitalizations, deaths and YLL
  burden = calc.final.burden(base_params, vacc_params)
  
  return(list(burden = burden,
              vacc_params = vacc_params))
  
}


master.outcomes = function(Coverage = seq(0, 1, 0.5),
                           R0 = c(3, 5, 7),
                           Program = c("Pfizer", "AstraZeneca", "Mix"),
                           Strategy = c("Untargeted", "Vulnerable", "Transmitters"),
                           Cutoff = c(15, 5),
                           Uptake = 0.9) {
  out = expand.grid(
    Coverage = Coverage,
    R0 = R0,
    Program = Program,
    Strategy = Strategy,
    Cutoff = Cutoff,
    Uptake = Uptake
  )
  
  fcn = function(coverage,
                 R0 = 5,
                 program = "Mix",
                 strategy = "Untargeted",
                 age_cutoff = 15,
                 uptake = 0.9) {
    burden = master.base(coverage,
                         R0,
                         program,
                         strategy,
                         age_cutoff,
                         uptake)$burden
    
    # Calculated total burden split by vaccination status
    segregated_totals = calc.segregated.totals(burden)
    
    # Combine calculated total with input parameters
    params = data.frame(
      Coverage = rep(coverage, nrow(segregated_totals)),
      R0 = rep(R0, nrow(segregated_totals)),
      Program = rep(program, nrow(segregated_totals)),
      Strategy = rep(strategy, nrow(segregated_totals)),
      Cutoff = rep(age_cutoff, nrow(segregated_totals)),
      Uptake = rep(uptake, nrow(segregated_totals))
    )
    
    return(cbind(params, segregated_totals))
  }
  
  out = bind_rows(with(out, {
    mapply(fcn,
           Coverage,
           R0,
           Program,
           Strategy,
           Cutoff,
           Uptake,
           SIMPLIFY = FALSE)
  }))
  
}



master.breakdown = function(Coverage = seq(0, 1, 0.2),
                            R0 = 5,
                            Program = "Mix",
                            Strategy = c("Untargeted", "Vulnerable", "Transmitters"),
                            Cutoff = 15,
                            Uptake = 0.9) {
  out = expand.grid(
    Coverage = Coverage,
    R0 = R0,
    Program = Program,
    Strategy = Strategy,
    Cutoff = Cutoff,
    Uptake = Uptake
  )
  
  
  fcn = function(coverage,
                 R0 = 5,
                 program = "Mix",
                 strategy = "Untargeted",
                 age_cutoff = 15,
                 uptake = 0.9) {
    res = master.base(coverage,
                      R0,
                      program,
                      strategy,
                      age_cutoff,
                      uptake)
    
    burden = res$burden[c("Hospitalizations", "Deaths"),]
    
    vacc_params = res$vacc_params
    
    # Calculated total burden split by vaccination status
    segregated_totals = calc.segregated.totals(burden)
    segregated_totals$Total = NULL
    
    # Add in events averted
    segregated_totals$Severity_averted = with(vacc_params, {
      rowSums(burden[, 17:32] * (1 / ((1 - Vsev) * (1 - Vmor)) - 1))
    })
    
    # Combine calculated total with input parameters
    params = data.frame(
      Coverage = rep(coverage, nrow(segregated_totals)),
      R0 = rep(R0, nrow(segregated_totals)),
      Program = rep(program, nrow(segregated_totals)),
      Strategy = rep(strategy, nrow(segregated_totals)),
      Cutoff = rep(age_cutoff, nrow(segregated_totals)),
      Uptake = rep(uptake, nrow(segregated_totals))
    )
    
    return(cbind(params, segregated_totals))
  }
  
  
  out = bind_rows(with(out, {
    mapply(fcn,
           Coverage,
           R0,
           Program,
           Strategy,
           Cutoff,
           Uptake,
           SIMPLIFY = FALSE)
  }))
  
  
  unmitigated = out %>% group_by(Strategy, Outcome) %>% summarize(Unmitigated = max(Unvaccinated))
  
  out = merge(x = out, y = unmitigated, all = TRUE) %>%
    dplyr::mutate(Infection_averted = Unmitigated - Unvaccinated - Vaccinated - Severity_averted)
  
  return(out)
}




master.sensitivity = function(n_samples = 1e1,
                              Coverage = seq(0, 1, 0.2),
                              R0 = c(3, 5, 7),
                              Program = c("Pfizer", "AstraZeneca", "Mix"),
                              Strategy = c("Untargeted", "Vulnerable", "Transmitters"),
                              Cutoff = c(15),
                              Uptake = 0.9) {
  # Generate set of sensitivity parameters
  print("Generating parameters...")
  
  samples = generate.samples(n_samples)
  
  print("Parameters generated.")
  
  # Create grid of desired outputs
  out = merge(
    x = samples,
    y = expand.grid(
      Coverage = Coverage,
      R0 = R0,
      Program = Program,
      Strategy = Strategy,
      Cutoff = Cutoff,
      Uptake = Uptake
    ),
    all = TRUE
  )
  
  print("Importing parameters...")
  
  base_params = apply(out[, c(
    "rel_infectious",
    "Pfizer_Vacq",
    "Pfizer_Vsev",
    "AstraZeneca_Vacq",
    "AstraZeneca_Vsev",
    "or_hospitalization",
    "or_death",
    "R0"
  )], 1, merge.base.params)
  
  print("Parameters imported.")
  
  print("Calculating doses...")
  doses = with(
    out,
    mapply(
      calc.age.specific.doses,
      base_params,
      Coverage,
      Program,
      Strategy,
      age_cutoff = Cutoff,
      uptake = Uptake,
      SIMPLIFY = FALSE
    )
  )
  print("Doses calculated.")
  
  vacc_params = mapply(calc.vacc.params, doses, base_params, SIMPLIFY =
                         FALSE)
  
  print("Calculating burden...")
  burden = mapply(calc.final.burden, base_params, vacc_params, SIMPLIFY = FALSE)
  
  print("Burden calculated.")
  
  segregated_totals = bind_rows(lapply(burden, calc.segregated.totals))
  
  out = cbind(out[rep(seq_len(nrow(out)), each = 4),],
              segregated_totals)
  
  return(out)
  
}





#### Outcome analysis ####
outcomes = c("Infections", "Hospitalizations", "Deaths", "YLL")

out_panel = master.outcomes(Coverage = seq(0, 1, 0.04))
plts_panel = sapply(outcomes, simplify =
                      FALSE, function(x)
                        plot.outcomes.panel(out_panel, outcome = x))
names(plts_panel) = outcomes
sapply(outcomes, function(x) ggsave(filename = paste0("../figures/outcomes/",
                                                      x,
                                                      ".uptake.9.png"),
                                    plot = plts_panel[[x]],
                                    dpi=600,
                                    height=6,
                                    width=8))





#### Uptake analysis ####
out_uptake = master.outcomes(Coverage = 1,
                             Uptake = seq(0,1,0.04),
                             Strategy = c("Untargeted"))

plts_uptake = sapply(outcomes, simplify =
                       FALSE, function(x)
                         plot.uptake.panel(out_uptake, outcome = x))
names(plts_uptake) = outcomes
sapply(outcomes, function(x) ggsave(filename = paste0("../figures/uptake/",
                                                      x,
                                                      ".uptake.9.png"),
                                    plot = plts_uptake[[x]],
                                    dpi=600,
                                    height=6,
                                    width=8))




#### Breakdown analysis ####
breakdown_grid = expand.grid(
  R0 = c(3, 5, 7),
  Cutoff = c(5, 15),
  Uptake = c(0.85, 0.9, 0.95)
)

out_breakdown = with(breakdown_grid,
{
  mapply(function(x, y, z)
  {
    master.breakdown(R0 = x,
                     Cutoff = y,
                     Uptake = z)
  }, R0, Cutoff, Uptake, SIMPLIFY = FALSE)
})
plts_breakdown = lapply(out_breakdown, plot.segragated.burden)
lbs = with(breakdown_grid, 
           mapply(function(x,y,z) paste0("R0.",x,".Cutoff.",y,".Uptake.",z), 
                  R0, Cutoff, Uptake))
names(plts_breakdown) = lbs

sapply(lbs, function(x) ggsave(filename = paste0("../figures/breakdown/",
                                                      x,
                                                      ".png"),
                                    plot = plts_breakdown[[x]],
                                    dpi=600,
                                    height=3,
                                    width=9))




#### Sensitivity analysis ####
out_sensitivity = master.sensitivity(n_samples = 500,
                                     Coverage = seq(0, 1, 0.04))

plts_sensitivity_cont = sapply(outcomes, simplify =
                            FALSE, function(x)
                              plot.sensitivity.panel.cont(out_sensitivity, outcome = x))
names(plts_sensitivity_cont) = outcomes
sapply(outcomes, function(x) ggsave(filename = paste0("../figures/sensitivity/continuous.",
                                                      x,
                                                      ".uptake.9.png"),
                                    plot = plts_sensitivity_cont[[x]],
                                    dpi=600,
                                    height=6,
                                    width=8))



plts_sensitivity_disc = sapply(outcomes, simplify =
                                 FALSE, function(x)
                                   plot.sensitivity.panel.disc(out_sensitivity, outcome = x))
names(plts_sensitivity_disc) = outcomes
sapply(outcomes, function(x) ggsave(filename = paste0("../figures/sensitivity/discrete.",
                                                      x,
                                                      ".uptake.9.png"),
                                    plot = plts_sensitivity_disc[[x]],
                                    dpi=600,
                                    height=6,
                                    width=8))

prccs = dplyr::bind_rows(with(
  expand.grid(
    Outcome = c("Infections",
                "Hospitalizations",
                "Deaths",
                "YLL"),
    Program = c("Pfizer",
                "AstraZeneca",
                "Mix")
  ),
  mapply(
    function(x, y, z)
      calc.prcc(out_sensitivity,
                outcome = x,
                program = y),
    Outcome,
    Program,
    SIMPLIFY = FALSE
  )
))


plts_prcc = plot.prccs(prccs)
ggsave("../figures/sensitivity/plts_prcc.png", plot=plts_prcc, dpi=600, width=5, height=7)



# Save workspace
save.image(paste0("../outputs/outputs_", format(Sys.time(), "%Y%m%d_%H%M%S_"), ".RData"))
