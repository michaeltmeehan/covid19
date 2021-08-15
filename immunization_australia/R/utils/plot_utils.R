

library(reshape2)
library(ggplot2)
library(scales)


# Plot outcome panels
plot.outcomes.panel = function(out,
                               outcome = "Deaths",
                               denominator = 1) {
  out$R0 = factor(out$R0)
  out$`Eligibility cutoff` = factor(out$Cutoff, levels = c("15", "5"))
  
  R0_labs = sapply(levels(out$R0), function(x) {
    paste("Reff =", x)
  })
  
  y_label = ifelse(outcome == "YLL", "Years of life lost", outcome)
  
  if (denominator > 1) {
    y_label = paste0(y_label, " per 100,000 population")
  }
  
  plt = out %>%
    dplyr::filter(Outcome == outcome) %>%
    ggplot2::ggplot(aes(
      x = Coverage,
      y = Total / denominator,
      col = Strategy,
      lty = `Eligibility cutoff`
    )) +
    geom_line(lwd = 0.5) +
    facet_grid(R0 ~ Program,
               labeller = labeller(R0 = R0_labs)) +
    ylab(y_label) +
    xlab("Target Coverage") +
    scale_x_continuous(
      label = scales::percent,
      breaks = seq(0, 1, by = 0.2),
      limits = c(0, 1)
    ) +
    scale_y_continuous(label = scales::comma) +
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
  
  return(plt)
}


# Plot outcome panels
plot.uptake.panel = function(out,
                             outcome = "Deaths",
                             denominator = 1
) {
  out$R0 = factor(out$R0)
  out$`Eligibility cutoff` = factor(out$Cutoff, levels = c("15", "5"))
  
  R0_labs = sapply(levels(out$R0), function(x) {
    paste("Reff =", x)
  })
  
  y_label = ifelse(outcome == "YLL", "Years of life lost", outcome)
  
  if (denominator > 1) {
    y_label = paste0(y_label, " per 100,000 population")
  }
  
  plt = out %>%
    dplyr::filter(Outcome == outcome) %>%
    ggplot2::ggplot(aes(
      x = Uptake,
      y = Total / denominator,
      lty = `Eligibility cutoff`
    )) +
    geom_line(lwd = 0.5) +
    facet_grid(R0 ~ Program,
               labeller = labeller(R0 = R0_labs)) +
    ylab(y_label) +
    xlab("Uptake") +
    scale_x_continuous(
      label = scales::percent,
      breaks = seq(0, 1, by = 0.2),
      limits = c(0, 1)
    ) +
    scale_y_continuous(label = scales::comma) +
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
  
  return(plt)
}




# Plot segragated burden
plot.segragated.burden = function(out) {
  plt = out %>%
    dplyr::select(-Unmitigated) %>%
    reshape2::melt(
      id.vars = c(
        "Strategy",
        "Outcome",
        "Coverage",
        "R0",
        "Program",
        "Cutoff",
        "Uptake"
      ),
      variable.name = "Group",
      value.name = "Cases"
    ) %>%
    dplyr::mutate(
      Group = dplyr::case_when(
        Group == "Unvaccinated" ~ "Unvaccinated (realized)",
        Group == "Vaccinated" ~ "Vaccinated (realized)",
        Group == "Severity_averted" ~ "Severity (averted)",
        Group == "Infection_averted" ~ "Infection (averted)"
      )
    ) %>%
    dplyr::filter(Outcome == "Deaths") %>%
    ggplot2::ggplot(aes(x = Coverage, y = Cases, fill = Group)) +
    geom_col() +
    facet_grid(Outcome ~ Strategy, scales = "free") +
    scale_x_continuous(label = scales::percent, breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(label = scales::comma) +
    scale_fill_manual(values = c("deepskyblue4", "deepskyblue", "darkred", "red")) +
    xlab("Target Coverage") +
    ylab("") +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.major.x = element_line(colour = "grey90"),
      axis.line = element_line(size = 0.4, colour = "grey10"),
      text = element_text(size = 12,  family = "serif"),
      legend.key = element_rect(fill = "white", colour = "white"),
      legend.title = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 12, colour = 'black'),
      panel.spacing = unit(1, "lines")
    )
  
  return(plt)
}


# Plot sensitivity outcome panels
plot.sensitivity.panel.cont = function(out,
                                  outcome = "Deaths",
                                  denominator = 1) {
  out$R0 = factor(out$R0)
  # out$`Eligibility cutoff` = factor(out$Cutoff, levels = c("15", "5"))
  
  R0_labs = sapply(levels(out$R0), function(x) {
    paste("Reff =", x)
  })
  
  y_label = ifelse(outcome == "YLL", "Years of life lost", outcome)
  
  plt = out %>%
    dplyr::filter(Outcome == outcome) %>%
    dplyr::select(-Unvaccinated, -Vaccinated) %>%
    dplyr::group_by(Coverage, R0, Program, Strategy, Cutoff, Uptake) %>%
    dplyr::summarize(lower = quantile(Total, prob=0.025),
                     median = quantile(Total, prob=0.5),
                     upper = quantile(Total, prob=0.975)) %>%
    ggplot2::ggplot(aes(
      x = Coverage,
      y = median / denominator,
      col = Strategy
    )) +
    geom_line(lwd = 0.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=Strategy), alpha=0.2) +
    facet_grid(R0 ~ Program,
               labeller = labeller(R0 = R0_labs)) +
    ylab(y_label) +
    xlab("Target Coverage") +
    scale_x_continuous(
      label = scales::percent,
      breaks = seq(0, 1, by = 0.2),
      limits = c(0, 1)
    ) +
    scale_y_continuous(label = scales::comma) +
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
  
  return(plt)
}



plot.sensitivity.panel.disc = function(out,
                                       outcome = "Deaths",
                                       denominator = 1) {
  out$R0 = factor(out$R0)
  # out$`Eligibility cutoff` = factor(out$Cutoff, levels = c("15", "5"))
  
  out$Coverage = factor(paste0(out$Coverage * 100, "%"),
                        levels = c("0%", "20%", "40%", "60%", "80%", "100%"))
  
  R0_labs = sapply(levels(out$R0), function(x) {
    paste("Reff =", x)
  })
  
  y_label = ifelse(outcome == "YLL", "Years of life lost", outcome)
  
  plt = out %>%
    dplyr::filter(Outcome == outcome) %>%
    ggplot2::ggplot(aes(
      x = Coverage,
      y = Total,
      fill = Strategy
    )) +
    geom_boxplot(lwd = 0.2, position = position_dodge(1)) +
    facet_grid(R0 ~ Program,
               labeller = labeller(R0 = R0_labs)) +
    ylab(y_label) +
    xlab("Target Coverage") +
    scale_x_discrete(expand = c(0.05, 0.05)) +
    scale_y_continuous(label = scales::comma) +
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
  
  return(plt)
}


plot.prccs = function(prccs) {
  
  prccs$Parameter = factor(prccs$Parameter,
                           levels = c(
                             "rel_infectious",
                             "Pfizer_VE",
                             "AstraZeneca_VE",
                             "or_hospitalization",
                             "or_death"
                           ))
  
  plt = prccs %>%
    ggplot2::ggplot(aes(x=Parameter, y=est, col=Parameter, shape=Parameter)) +
    geom_hline(yintercept = 0, lty=2, lwd=1) +
    # geom_errorbar(aes(ymin=lower, ymax=upper), lwd=1.5, width=0, position=position_dodge(0.7)) +
    geom_point(size=2, position=position_dodge(0.7)) +
    ylab("Partial Rank Correlation Coefficient") +
    facet_grid(Outcome ~ Program) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.major.x = element_line(colour = "grey90"),
      axis.line = element_line(size = 0.4, colour = "grey10"),
      text = element_text(size = 12,  family = "serif"),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1),
      legend.position = "none",
      legend.title = element_blank(),
      # legend.key = element_rect(fill = "white", colour = "white"),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 12, colour = 'black'),
      panel.spacing = unit(1, "lines")
    )
  
}
