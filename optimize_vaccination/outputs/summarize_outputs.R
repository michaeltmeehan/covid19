library(dplyr)
library(reshape2)

source("../R/utils/plot_utils.R")

##### Variable doses #####

summarize_var_results = function(res){
  colnames(res) = c("country", "R0", "HIT", "optimum", "0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+")
  
  optima = res %>%
    dplyr::group_by(country) %>%
    dplyr::filter(optimum == min(optimum)) %>%
    dplyr::ungroup()
  
  res_summary = res %>%
    reshape2::melt(id.vars=c("country","R0","HIT","optimum"),variable.name="age_group",value.name="coverage") %>%
    dplyr::group_by(country, age_group) %>%
    dplyr::summarise(q0 = min(coverage),
                     q025 = quantile(coverage, 0.025),
                     q25 = quantile(coverage, 0.25),
                     q50 = quantile(coverage,0.5),
                     q75 = quantile(coverage,0.75),
                     q975 = quantile(coverage,0.975),
                     q1 = max(coverage)) %>%
    dplyr::ungroup()
  
  return(list(optima=optima, summary=res_summary))
}

popdata = read.csv("../data/popdata/poptotal_summarized.csv")

res_examples = summarize_var_results(read.csv("../outputs/var_dose/example_countries_infection_constant_eff_res.csv"))

res_all = summarize_var_results(read.csv("../outputs/var_dose/all_countries_infection_constant_eff_res.csv"))

res_all$optima = left_join(res_all$optima, popdata[,c("country","median_age")])

ggplot(res_all$optima %>% filter(R0 < 5), aes(x=median_age, y=optimum/HIT, col=R0)) +
  geom_point() +
  xlab("Median age") +
  ylab("Optimized coverage reduction") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0.5, 1, by=0.1), limits=c(0.5,1)) +
  theme_bw() +
  theme(text = element_text(size=16,  family="serif"),
        panel.spacing = unit(2, "lines"),
        strip.background =element_blank(),
        strip.text = element_text(size=16))

plot_optimal_var_strategy_by_age("India",res_examples$optima, res_examples$summary)
plot_optimal_var_strategy_by_age("China",res_examples$optima, res_examples$summary, mark=TRUE)
plot_optimal_var_strategy_by_age("United Kingdom",res_examples$optima, res_examples$summary)



#### Fixed doses ####

summarize_fixed_results = function(res, trim=TRUE){
  colnames(res) = c("country", "R0", "dose", "Reff", "0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+")
  
  
  optima = res %>%
    dplyr::group_by(country, dose) %>%
    dplyr::filter(Reff == min(Reff)) %>%
    dplyr::ungroup()
  
  optima = optima[!duplicated(optima$Reff), ]
  
  if(trim){
    res = res %>%
      dplyr::group_by(country, dose) %>%
      dplyr::filter(Reff <= 1.01 * min(Reff)) %>%
      dplyr::ungroup()
  }
  
  summary = res %>%
    reshape2::melt(id.vars=c("country","R0","dose","Reff"),variable.name="age_group",value.name="coverage") %>%
    dplyr::group_by(country, age_group, dose) %>%
    #dplyr::group_by(country, age_group) %>%
    #dplyr::group_by(dose) %>%
    dplyr::summarise(q0 = min(coverage),
                     q025 = quantile(coverage, 0.025),
                     q25 = quantile(coverage, 0.25),
                     q50 = quantile(coverage,0.5),
                     q75 = quantile(coverage,0.75),
                     q975 = quantile(coverage,0.975),
                     q1 = max(coverage)) %>%
    dplyr::ungroup()
  
  return(list(optima=optima, summary=summary))
}


fixed_res_infection = summarize_fixed_results(read.csv("../outputs/fixed_dose/example_countries_infection_constant_eff_res.csv"))

fixed_res_disease = summarize_fixed_results(read.csv("../outputs/fixed_dose/example_countries_disease_constant_eff_res.csv"))

fixed_both = rbind(fixed_res_infection$optima, fixed_res_disease$optima)
fixed_both$vac_mode = relevel(factor(rep(c("Infection", "Disease"), each=nrow(fixed_res_infection$optima))), ref="Infection")

fixed_both = fixed_both[,c("country", "dose", "Reff", "vac_mode")]

uniform_vac = data.frame(country = rep(rep(c("India", "China", "United Kingdom"), each=6), times=2),
                         dose = rep(rep(seq(0,1,0.2),times=3),times=2),
                         Reff = c(2.200000000000005,
                                  1.8920000000000023,
                                  1.584000000000001,
                                  1.2760000000000031,
                                  0.9680000000000022,
                                  0.6600000000000004,
                                  2.5999999999999988,
                                  2.235999999999999,
                                  1.8720000000000003,
                                  1.5079999999999973,
                                  1.1439999999999997,
                                  0.7799999999999998,
                                  2.300000000000002,
                                  1.9780000000000006,
                                  1.6560000000000032,
                                  1.3340000000000019,
                                  1.0120000000000013,
                                  0.6900000000000013,
                                  2.200000000000005,
                                  2.1301629124567065,
                                  2.060427260452016,
                                  1.9908042613326287,
                                  1.921306793531393,
                                  1.8519497031827041,
                                  2.5999999999999988,
                                  2.478387609837972,
                                  2.357200280134312,
                                  2.2365094647713892,
                                  2.1164021292339537,
                                  1.9969846698224212,
                                  2.300000000000002,
                                  2.2140383340966734,
                                  2.1282220074275564,
                                  2.0425690075304748,
                                  1.9571003246909355,
                                  1.871840583016692),
                         vac_mode = rep(c("Infection","Disease"),each=18))

uniform_vac$Allocation = "Uniform"
fixed_both$Allocation = "Targeted"

fixed_both = rbind(fixed_both, uniform_vac)

fixed_both$country = factor(fixed_both$country, levels=c("India", "China", "United Kingdom"))

ggplot(fixed_both, aes(x=dose, y=Reff, lty = Allocation, col=vac_mode)) +
  geom_point(size=2) +
  geom_line(lwd=0.7) +
  geom_hline(yintercept = 1, col="black", lty=3) +
  xlab("Population-level coverage") +
  ylab(expression(R["eff"])) +
  labs(lty="Vaccine allocation", col="Vaccination mode") +
  scale_y_continuous(breaks=seq(0,3,by=0.5), limits=c(0,3)) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, by=0.2), limits=c(0,1)) +
  facet_grid(.~country) +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(size = 0.4, colour = "grey10"),
        text = element_text(size=12,  family="serif"),
        legend.key = element_rect(fill = "white", colour = "white"),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size = 12, colour = 'black'))




fixed_plot_age_dist = list()
fixed_plot_age_dist[[1]] = plot_optimal_fixed_strategy_by_age("India", fixed_res_infection$optima, fixed_res_infection$summary, error_bars=TRUE) + xlab("India") + ylab("")
fixed_plot_age_dist[[4]] = plot_optimal_fixed_strategy_by_age("India", fixed_res_disease$optima, fixed_res_disease$summary, error_bars=TRUE) + xlab("") + ylab("")
fixed_plot_age_dist[[2]] = plot_optimal_fixed_strategy_by_age("China", fixed_res_infection$optima, fixed_res_infection$summary, error_bars=TRUE) + xlab("China") + ylab("")
fixed_plot_age_dist[[5]] = plot_optimal_fixed_strategy_by_age("China", fixed_res_disease$optima, fixed_res_disease$summary, error_bars=TRUE) + xlab("") + ylab("")
fixed_plot_age_dist[[3]] = plot_optimal_fixed_strategy_by_age("United Kingdom", fixed_res_infection$optima, fixed_res_infection$summary, error_bars=TRUE) + xlab("United Kingdom") + ylab("")
fixed_plot_age_dist[[6]] = plot_optimal_fixed_strategy_by_age("United Kingdom", fixed_res_disease$optima, fixed_res_disease$summary, error_bars=TRUE) + xlab("") + ylab("")

fig = ggarrange(plotlist=fixed_plot_age_dist, nrow=2, ncol=3, common.legend = TRUE, labels = c("A","B","C","D","E","F"))
annotate_figure(fig,
                left = text_grob("Age-specific coverage", rot=90, family="serif", size=16),
                bottom = text_grob("Population-level coverage", family="serif", size=16))
