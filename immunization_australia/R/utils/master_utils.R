# Utility functions that act as wrappers around
# the utilities from other utils files
source("utils/dose_utils.R")
source("utils/final_size_utils.R")
source("utils/import_utils.R")
source("utils/NGM_utils.R")
source("utils/optim_utils.R")

library(scales)

calc_targets = function(params,
                        targets = c("Infections", "Hospitalizations", "Deaths"),
                        strategies = c("Untargeted", "Vulnerable", "Transmitters"),
                        target_coverages = seq(0, 1.0, 0.2),
                        vaccine_combinations = list(
                          "Pfizer" = c("BNT162b2", "BNT162b2"),
                          "AstraZeneca" = c("ChAdOx1", "ChAdOx1"),
                          "Mix" = c("BNT162b2", "ChAdOx1")
                        ),
                        eligibility_cutoff = 15,
                        uptake = 1.0) {
  out = expand.grid(
    Strategy = strategies,
    Vaccine = names(vaccine_combinations),
    Coverage = target_coverages
  )

  out$Infections = rep(NA, nrow(out))
  out$Hospitalizations = rep(NA, nrow(out))
  out$Deaths = rep(NA, nrow(out))

  out_lower = out
  out_upper = out

  params_central = params
  params_lower = params
  params_upper = params

  for (i in 1:nrow(out)) {
    # Calculate age-specific vaccination coverage
    if (eligibility_cutoff == 15) {
      vacc_prop = calc_age_specific_coverage_over_15(
        target_coverage = out$Coverage[i],
        params = params,
        strategy = out$Strategy[i],
        uptake = uptake
      )
    } else if (eligibility_cutoff == 10) {
      vacc_prop = calc_age_specific_coverage_over_10(
        target_coverage = out$Coverage[i],
        params = params,
        strategy = out$Strategy[i],
        uptake = uptake
      )
    } else if (eligibility_cutoff == 5) {
      vacc_prop = calc_age_specific_coverage_over_5(
        target_coverage = out$Coverage[i],
        params = params,
        strategy = out$Strategy[i],
        uptake = uptake
      )
    }
    # Determine current vaccine combination
    vacc_comb = vaccine_combinations[[out$Vaccine[i]]]

    # Update the vaccine efficacy given the current vaccine combination
    params_central$efficacy = import_efficacy_data(params$strain,
                                                   vaccine_1 = vacc_comb[1],
                                                   vaccine_2 = vacc_comb[2])

    # Incorporate uncertainty in vaccine efficacy
    # params_lower$efficacy = import_efficacy_data(
    #   params$strain,
    #   vaccine_1 = paste(vacc_comb[1], "(lower)", sep =
    #                       " "),
    #   vaccine_2 = paste(vacc_comb[2], "(lower)", sep =
    #                       " ")
    # )
    #
    # params_upper$efficacy = import_efficacy_data(
    #   params$strain,
    #   vaccine_1 = paste(vacc_comb[1], "(upper)", sep =
    #                       " "),
    #   vaccine_2 = paste(vacc_comb[2], "(upper)", sep =
    #                       " ")
    # )


    # Incorporate uncertainty in R0
    params_lower = get_params(R0 = 3,
                              vaccine_1 = vacc_comb[1],
                              vaccine_2 = vacc_comb[2])

    params_upper = get_params(R0 = 7,
                              vaccine_1 = vacc_comb[1],
                              vaccine_2 = vacc_comb[2])

    # Calculate the age-specific infections / hospitalizations / deaths
    final_burden = calc_final_burden(params_central, vacc_prop)

    out$Infections[i] = sum(final_burden$infections) * 1e5 / sum(params$population_size)
    out$Hospitalizations[i] = sum(final_burden$hospitalizations) * 1e5 / sum(params$population_size)
    out$Deaths[i] = sum(final_burden$deaths) * 1e5 / sum(params$population_size)


    final_burden_lower = calc_final_burden(params_lower, vacc_prop)

    out_lower$Infections[i] = sum(final_burden_lower$infections) * 1e5 / sum(params$population_size)
    out_lower$Hospitalizations[i] = sum(final_burden_lower$hospitalizations) * 1e5 / sum(params$population_size)
    out_lower$Deaths[i] = sum(final_burden_lower$deaths) * 1e5 / sum(params$population_size)


    final_burden_upper = calc_final_burden(params_upper, vacc_prop)

    out_upper$Infections[i] = sum(final_burden_upper$infections) * 1e5 / sum(params$population_size)
    out_upper$Hospitalizations[i] = sum(final_burden_upper$hospitalizations) * 1e5 / sum(params$population_size)
    out_upper$Deaths[i] = sum(final_burden_upper$deaths) * 1e5 / sum(params$population_size)

  }

  # Melt, and then stitch together lower, central and upper estimates
  out = out %>%
    dplyr::select(-Hospitalizations) %>%
    reshape2::melt(id.vars = c("Strategy", "Vaccine", "Coverage"),
                   variable.name = "Target")
  out_lower = out_lower %>%
    dplyr::select(-Hospitalizations) %>%
    reshape2::melt(id.vars = c("Strategy", "Vaccine", "Coverage"),
                   variable.name = "Target")
  out_upper = out_upper %>%
    dplyr::select(-Hospitalizations) %>%
    reshape2::melt(id.vars = c("Strategy", "Vaccine", "Coverage"),
                   variable.name = "Target")

  out$lower = out_lower$value
  out$upper = out_upper$value

  return(out)
}



calc_targets_vec = function(params,
                        R0s = c(3,5,7),
                        strategies = c("Untargeted", "Vulnerable", "Transmitters"),
                        target_coverages = seq(0, 1.0, 0.2),
                        vaccine_combinations = list(
                          "Pfizer" = c("BNT162b2", "BNT162b2"),
                          "AstraZeneca" = c("ChAdOx1", "ChAdOx1"),
                          "Mix" = c("BNT162b2", "ChAdOx1")
                        ),
                        eligibility_cutoffs = c(5,15),
                        uptake = 0.9) {
  out = expand.grid(
    R0 = R0s,
    Age = eligibility_cutoffs,
    Strategy = strategies,
    Vaccine = names(vaccine_combinations),
    Coverage = target_coverages
  )
  
  out$Infections = rep(NA, nrow(out))
  out$Hospitalizations = rep(NA, nrow(out))
  out$Deaths = rep(NA, nrow(out))
  
  params_temp = params
  
  for (i in 1:nrow(out)) {
    
    vacc_comb = vaccine_combinations[[out$Vaccine[i]]]
    
    params_temp = get_params(R0 = out$R0[i],
                             vaccine_1 = vacc_comb[1],
                             vaccine_2 = vacc_comb[2])
    
    lower_age_limit = out$Age[i]
    
    
    # Calculate age-specific vaccination coverage
    if (lower_age_limit == 15) {
      vacc_prop = calc_age_specific_coverage_over_15(
        target_coverage = out$Coverage[i],
        params = params_temp,
        strategy = out$Strategy[i],
        uptake = uptake
      )
    } else if (lower_age_limit == 10) {
      vacc_prop = calc_age_specific_coverage_over_10(
        target_coverage = out$Coverage[i],
        params = params_temp,
        strategy = out$Strategy[i],
        uptake = uptake
      )
    } else if (lower_age_limit == 5) {
      vacc_prop = calc_age_specific_coverage_over_5(
        target_coverage = out$Coverage[i],
        params = params_temp,
        strategy = out$Strategy[i],
        uptake = uptake
      )
    }
    # Determine current vaccine combination
    vacc_comb = vaccine_combinations[[out$Vaccine[i]]]
    
    # Calculate the age-specific infections / hospitalizations / deaths
    final_burden = calc_final_burden(params_temp, vacc_prop)
    
    out$Infections[i] = sum(final_burden$infections) * 1e5 / sum(params$population_size)
    out$Hospitalizations[i] = sum(final_burden$hospitalizations) * 1e5 / sum(params$population_size)
    out$Deaths[i] = sum(final_burden$deaths) * 1e5 / sum(params$population_size)
    
  }

  return(out)
}




calc_target_breakdown =  function(params,
                                  target = "Deaths",
                                  strategies = c("Untargeted", "Vulnerable", "Transmitters"),
                                  target_coverages = seq(0, 1.0, 0.2),
                                  vaccine_combination = c("BNT162b2", "ChAdOx1"),
                                  eligibility_cutoff = 15,
                                  uptake = 1.0) {
  out = expand.grid(Strategy = strategies,
                    Coverage = target_coverages)
  
  out$`Unvaccinated (realized)` = rep(NA, nrow(out))
  out$`Vaccinated (realized)` = rep(NA, nrow(out))
  out$`Severity (averted)` = rep(NA, nrow(out))
  out$`Infection (averted)` = rep(NA, nrow(out))
  
  unmitigated_burden = calc_final_burden(params, rep(0,16))
  unmitigated_deaths = sum(unmitigated_burden$deaths)
  
  for (i in 1:nrow(out)) {
    # Calculate age-specific vaccination coverage
    if (eligibility_cutoff == 15) {
      vacc_prop = calc_age_specific_coverage_over_15(
        target_coverage = out$Coverage[i],
        params = params,
        strategy = out$Strategy[i],
        uptake = uptake
      )
    } else if (eligibility_cutoff == 10) {
      vacc_prop = calc_age_specific_coverage_over_10(
        target_coverage = out$Coverage[i],
        params = params,
        strategy = out$Strategy[i],
        uptake = uptake
      )
    } else if (eligibility_cutoff == 5) {
      vacc_prop = calc_age_specific_coverage_over_5(
        target_coverage = out$Coverage[i],
        params = params,
        strategy = out$Strategy[i],
        uptake = uptake
      )
    }
    
    # Update the vaccine efficacy given the current vaccine combination
    params$efficacy = import_efficacy_data(params$strain,
                                           vaccine_1 = vaccine_combination[1],
                                           vaccine_2 = vaccine_combination[2])
    
    
    pop_size = unlist(c((1 - vacc_prop) * params$population_size,
                        vacc_prop * params$population_size))
    
    
    # Calculate the age-specific infections
    final_size = calc_final_size(params, vacc_prop)
    
    out$`Unvaccinated (realized)`[i] = sum(params$epi$infection_fatality_rate * (1 - params$seropositivity) *
      final_size[1:16] * pop_size[1:16])
    
    out$`Vaccinated (realized)`[i]  = sum(params$epi$infection_fatality_rate * (1 - params$seropositivity) * 
      (1 - params$efficacy$Vsev) * final_size[17:32] * pop_size[17:32])
    
    out$`Severity (averted)`[i] =  sum(params$epi$infection_fatality_rate * (1 - params$seropositivity) * 
                                      final_size[17:32] * pop_size[17:32]) - out$`Vaccinated (realized)`[i]
    
    out$`Infection (averted)`[i] = max(c(0,unmitigated_deaths - (out$`Unvaccinated (realized)`[i] + 
                                                       out$`Vaccinated (realized)`[i] +
                                                       out$`Severity (averted)`[i])))
    
  }
  
  return(out)
}



plot_targets = function(out) {
  zero_data = data.frame(Strategy = factor(rep("Untargeted", 6)),
                         Vaccine = factor(rep(c("Pfizer", "AstraZeneca","Mix"),each=2), levels=c("Pfizer", "AstraZeneca", "Mix")),
                         Coverage = rep(0,6),
                         Target = factor(rep(c("Infections", "Deaths"), 3), levels=c("Infections", "Deaths")),
                         value = rep(0, 6),
                         lower = rep(0,6),
                         upper = rep(0,6))
  plt = out %>%
    ggplot2::ggplot(aes(
      x = Coverage,
      y = value,
      col = Strategy,
      fill = Strategy
    )) +
    geom_line(lwd = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.15,
                colour = NA) +
    geom_point(data=zero_data, aes(x=Coverage, y=value), alpha=0) +
    facet_grid(Target ~ Vaccine, scales = "free") +
    ylab("Cases per 100,000 population") +
    xlab("Target Coverage") +
    scale_x_continuous(
      label = scales::percent,
      breaks = seq(0, 1, by = 0.2),
      limits = c(0, 1)
    ) +
    scale_y_continuous(label=scales::comma) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.major.x = element_blank(),
      axis.line = element_line(size = 0.4, colour = "grey10"),
      text = element_text(size = 12,  family = "serif"),
      legend.key = element_rect(fill = "white", colour = "white"),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 12, colour = 'black'),
      panel.spacing = unit(1, "lines")
    )
  
  print(plt)
}


plot_breakdown = function(out, target="Deaths") {
  out$Coverage = factor(paste(as.integer(out$Coverage * 100), "%", sep=""),
                           levels = c("0%", "20%", "40%", "60%", "80%", "100%"))
  out[,3:6] = out[,3:6] / rowSums(out[,3:6])
  out = reshape2::melt(out, id.vars = c("Strategy", "Coverage"),
                       variable.name = "Group")
  out$Group = factor(out$Group, levels = c("Infection (averted)",
                                           "Severity (averted)",
                                           "Unvaccinated (realized)",
                                           "Vaccinated (realized)"))
  plt = ggplot2::ggplot(out, aes(x=Coverage, y=value, fill=Group)) +
    geom_col() +
    ylab(paste("%", target, "(relative to unmitigated epidemic)")) +
    xlab("Target Coverage") +
    facet_grid(.~Strategy) +
    scale_fill_manual(values = c("deepskyblue4", "deepskyblue", "darkred", "red")) +
    scale_y_continuous(expand = expansion(mult = c(0.03, .0)),
                       label = scales::percent,
                       breaks = seq(0, 1, by = 0.2),
                       limits = c(0, 1.02)
                       ) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_line(colour = "grey90"),
      panel.grid.major.x = element_blank(),
      axis.line = element_line(size = 0.4, colour = "grey10"),
      text = element_text(size = 12,  family = "serif"),
      legend.key = element_rect(fill = "white", colour = "white"),
      legend.title = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 12, colour = 'black'),
      panel.spacing = unit(1, "lines")
    )
  print(plt)
}


generate_plots = function(cutoffs = c(5,10,15),
                          uptakes = c(0.9, 0.95, 1.0)) {
  # Loop over all input parameters
  for (age in cutoffs) {
    for (uptake in uptakes) {
      # Target plots
      params = get_params()
      
      out = calc_targets(
        params,
        target_coverages = seq(0, 1, 0.02),
        eligibility_cutoff = age,
        uptake = uptake
      )
      
      plt = plot_targets(out)
      
      ggsave(
        paste0(
          "../MJA_figures/targets",
          "_uptake.",
          uptake,
          "_cutoff.",
          age,
          ".png"
        ),
        plot = plt,
        dpi = 600,
        width = 12,
        height = 6
      )
      
      
      # Breakdown plots
      bout = calc_target_breakdown(params,
                                   eligibility_cutoff = age,
                                   uptake = uptake)
      
      bplt = plot_breakdown(bout)
      
      ggsave(
        paste0(
          "../MJA_figures/deaths_breakdown",
          "_uptake.",
          uptake,
          "_cutoff.",
          age,
          ".png"
        ),
        plot = bplt,
        dpi = 600,
        width = 12,
        height = 4
      )
      
    }
  }
}

plot_infections = function(out){

out$R0 = factor(out$R0)
# out$R0_lab = factor(sapply(levels(out$R0), function(x){paste(c("R[0] =", x))}))
out$`Eligibility cutoff` = factor(out$Age, levels=c("15", "5"))

R0.labs = sapply(levels(out$R0), function(x){paste("Reff =", x)})

ggplot2::ggplot(out, aes(x=Coverage, y=Infections, col=Strategy, lty=`Eligibility cutoff`)) +
  geom_line(lwd=0.5) +
  facet_grid(R0 ~ Vaccine,
             labeller=labeller(R0=R0.labs),
             # labeller=label_parsed
             ) +
  ylab("Infections per 100,000 population") +
  xlab("Target Coverage") +
  scale_x_continuous(
    label = scales::percent,
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  ) +
  scale_y_continuous(label=scales::comma) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.major.x = element_line(colour = "grey90"),
    axis.line = element_line(size = 0.4, colour = "grey10"),
    text = element_text(size = 12,  family = "serif"),
    legend.key = element_rect(fill = "white", colour = "white"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12, colour = 'black'),
    panel.spacing = unit(1, "lines")
  )
}


plot_deaths = function(out){
  
  out$R0 = factor(out$R0)
  # out$R0_lab = factor(sapply(levels(out$R0), function(x){paste(c("R[0] =", x))}))
  out$`Eligibility cutoff` = factor(out$Age, levels=c("15", "5"))
  
  R0.labs = sapply(levels(out$R0), function(x){paste("Reff =", x)})
  
  ggplot2::ggplot(out, aes(x=Coverage, y=Deaths, col=Strategy, lty=`Eligibility cutoff`)) +
    geom_line(lwd=0.5) +
    facet_grid(R0 ~ Vaccine,
               labeller=labeller(R0=R0.labs),
               # labeller=label_parsed
    ) +
    ylab("Deaths per 100,000 population") +
    xlab("Target Coverage") +
    scale_x_continuous(
      label = scales::percent,
      breaks = seq(0, 1, by = 0.2),
      limits = c(0, 1)
    ) +
    scale_y_continuous(label=scales::comma) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.major.x = element_line(colour = "grey90"),
      axis.line = element_line(size = 0.4, colour = "grey10"),
      text = element_text(size = 12,  family = "serif"),
      legend.key = element_rect(fill = "white", colour = "white"),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 12, colour = 'black'),
      panel.spacing = unit(1, "lines")
    )
}
