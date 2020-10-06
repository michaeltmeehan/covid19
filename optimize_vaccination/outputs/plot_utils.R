library(ggplot2)
library(gplots)
library(dplyr)
library(reshape2)


# plot_R0_calibration = function(df){
#   
#   plt = ggplot(df, aes(x=R0_cdc, y=R0_cal)) +
#     geom_point(alpha=0.7, col="royalblue2") +
#     geom_abline(slope=1, intercept=0, lty=2, col="royalblue2") +
#     xlab(expression(R["0,CDC"])) +
#     ylab(expression(R["0,Cal"])) +
#     ylim(c(0.5, 4.5)) +
#     xlim(c(0.5, 4.5)) +
#     theme_bw() +
#     theme(text = element_text(size=16,  family="serif"),
#           panel.spacing = unit(2, "lines"),
#           strip.background =element_blank(),
#           strip.text = element_text(size=16))
#   
#   return(plt)
#   
# }



# plot_Reff = function(df, countries, add_ref=TRUE){
#   
#   plt = df %>%
#     dplyr::filter(country %in% countries) %>%
#     ggplot(aes(x=dose, y=Reff, col=country)) +
#       geom_point(size=2) +
#       geom_line(lwd=0.5) +
#       geom_hline(yintercept=1, lty=2) +
#       xlab("Vaccination proportion") +
#       ylab(expression(R["eff"])) +
#       scale_y_continuous(breaks=seq(0,3,by=0.5), limits=c(0,3)) +
#       scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, by=0.2), limits=c(0,1)) +
#       theme(panel.background = element_rect(fill = "white", colour = "white"),
#             panel.grid.major = element_line(colour = "grey90"),
#             panel.grid.major.x = element_blank(),
#             axis.line = element_line(size = 0.4, colour = "grey10"),
#             text = element_text(size=12,  family="serif"),
#             legend.key = element_rect(fill = "white", colour = "white"),
#             strip.background =element_rect(fill="royalblue"),
#             strip.text = element_text(size = 10, colour = 'white'),
#             legend.title=element_blank())
#   
#   if (add_ref){
#     ref_data = data.frame(dose=seq(0,1,by=0.05))
#     ref_data$Reff = 2.43622 * (1 - 0.7 * ref_data$dose)
#     plt = plt + 
#       geom_point(data=ref_data, aes(x=dose, y=Reff), col="black", size=2) +
#       geom_line(data=ref_data, aes(x=dose, y=Reff), col="black", lty=2, lwd=0.5)
#   }
#   
#   return(plt)
# }
  
  

# plot_optimal_var_strategy_by_age = function(cname, optima, res_summary, error_bars = TRUE, mark=FALSE){
#   
#   optima_ = optima %>%
#     dplyr::filter(country==cname) %>%
#     reshape2::melt(id.vars=c("country","R0","HIT","optimum"),variable.name="age_group",value.name="coverage")
#   
#   plt = ggplot(optima_, aes(x=age_group, y=coverage)) +
#     geom_col(fill="royalblue2", alpha=0.9, col="royalblue2", lwd=1) +
#     geom_hline(yintercept=max(optima_$HIT), lty=2, lwd=1) +
#     geom_hline(yintercept=min(optima_$optimum), col="firebrick2", lty=2, lwd=1) +
#     ylab("Vaccination coverage") +
#     xlab("Age group (years)") +
#     scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, by=0.2), limits=c(0,1)) +
#     theme_bw() +
#     theme(text = element_text(size=16,  family="serif"),
#           panel.spacing = unit(2, "lines"),
#           strip.background =element_blank(),
#           strip.text = element_text(size=16),
#           legend.title = element_blank())
#   
#   if(mark){
#     plt = plt +
#       annotate(geom = "text", label="Uniform target", hjust=0.5, x=2, y = max(HIT) *1.05, family="serif") +
#       annotate(geom = "text", label="Optimized target", hjust=0.5, x=2, y = min(optimum) *1.08, col="firebrick2", family="serif")
#   }
#   
#   if(error_bars){
#     res_summary_ = res_summary %>%
#       dplyr::filter(country==cname)
#     plt = plt +
#       geom_errorbar(data=res_summary_, mapping=aes(x=age_group,ymin=q025, ymax=q975), inherit.aes = FALSE, width=0, col="dimgray") +
#       geom_point(data=res_summary_, aes(x=age_group,y=q50), col="dimgray")
#       
#   }
#   
#   return(plt)
#   
# }


# plot_optimal_var_results_all_countries = function(df){
#   
#   plt = ggplot(df, aes(x=age, y=optimum / HIT)) +
#     geom_point(alpha=0.7, col="royalblue2") +
#     ylab("Relative reduction") +
#     xlab("Median age (years)") +
#     scale_y_continuous(labels = scales::percent, breaks = seq(0.5, 1, by=0.1), limits=c(0.5,1)) +
#     theme_bw() +
#     theme(text = element_text(size=16,  family="serif"),
#           panel.spacing = unit(2, "lines"),
#           strip.background =element_blank(),
#           strip.text = element_text(size=16),
#           legend.title = element_blank())
#   
#   return(plt)
# }


plot_optimal_fixed_strategy_by_age = function(cname, optima, res_summary=NULL, error_bars = FALSE){
  
  optima_ = optima %>%
    dplyr::filter(country == cname) %>%
    dplyr::filter(dose %in% c(0.2, 0.4, 0.6, 0.8)) %>%
    reshape2::melt(id.vars=c("country","R0","dose", "Reff"), variable.name="age_group", value.name="coverage") %>%
    dplyr::mutate(dose = factor(dose, levels=c(0.2, 0.4, 0.6, 0.8), labels=c("20%", "40%", "60%", "80%")))

  if(!error_bars){
  plt = ggplot(optima_, aes(x=dose, y=coverage, fill=age_group, col=age_group)) +
    geom_col(position=position_dodge(0.7), width=0.7, alpha=0.9, lwd=0.5) +
    ylab("Age-specific coverage") +
    xlab("Population-level coverage") +
    labs(fill="Age group (years)", col="Age group (years)") +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, by=0.2), limits=c(0,1)) +
    theme_bw() +
    theme(text = element_text(size=12,  family="serif"),
          panel.spacing = unit(2, "lines"),
          strip.background =element_blank(),
          strip.text = element_text(size=12))
  }
  else{
    res_summary_ = res_summary %>%
      dplyr::filter(country == cname) %>%
      dplyr::filter(dose %in% c(0.2, 0.4, 0.6, 0.8)) %>%
      #reshape2::melt(id.vars=c("country", "age_group", "dose"), variable.name="quantile", value.name="coverage") %>%
      dplyr::mutate(dose = factor(dose, levels=c(0.2, 0.4, 0.6, 0.8), labels=c("20%", "40%", "60%", "80%")))

    res_summary_$optima = optima_$coverage
    
    plt = ggplot(res_summary_, aes(x=dose, fill=age_group, col=age_group)) +
      #geom_col(aes(y=optimum), position=position_dodge(0.7), width=0.7, alpha=0.9, lwd=1) +
      geom_col(data=optima_, aes(x=dose, y=coverage, fill=age_group, col=age_group), position=position_dodge(0.7), width=0.7, alpha=0.9, lwd=0.25) +
      geom_point(aes(y=q50), position=position_dodge(0.7), col="dimgray") +
      geom_errorbar(aes(ymin=q025, ymax=q975), position=position_dodge(0.7), width=0, col="dimgray") +
      ylab("Age-specific coverage") +
      xlab("Population-level coverage") +
      labs(fill="Age group (years)", col="Age group (years)") +
      scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, by=0.2), limits=c(0,1)) +
      theme_bw() +
      theme(text = element_text(size=14,  family="serif"),
            panel.spacing = unit(2, "lines"),
            strip.background =element_blank(),
            strip.text = element_text(size=14))
    
    
    
  }
  
  return(plt)
  
}

make_pop_pyramid = function(cname){
  
  popfemale = read.csv("popdata/popfemale.csv") %>% dplyr::filter(countryname == cname)
  popmale = read.csv("popdata/popmale.csv") %>% dplyr::filter(countryname == cname)
  
  popdata = rbind(popfemale, popmale)
  names(popdata) = c("is03c", "countryname", "year", "0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79","80-84","85-89","90-94","95-99","100+","total")
  popdata$Gender = as.factor(c("Female","Male"))
  
  popdata = popdata %>%
    dplyr::select(-c("is03c","countryname","year","total")) %>%
    reshape2::melt(id.vars=c("Gender"), variable.name="Age",value.name="population") %>%
    dplyr::mutate(Age = factor(Age, levels=Age, labels=Age),
                  Proportion = population / (popfemale$total + popmale$total),
                  Proportion = ifelse(Gender == "Male", -Proportion, Proportion))
  
  plt = ggplot(popdata, aes(x = Age, y = Proportion, fill = Gender)) + 
    geom_bar(subset = .(Gender == "Female"), stat = "identity") + 
    geom_bar(subset = .(Gender == "Male"), stat = "identity") + 
    scale_y_continuous(breaks = seq(-0.3, 0.3, 0.05),
                       labels = paste0(as.character(c(seq(30, 0, -5), seq(5, 30, 5))), "%"),
                       limits=c(-0.15, 0.15)) +
    coord_flip() + 
    scale_fill_brewer(palette = "Set1") + 
    theme_bw() +
    theme(text = element_text(size=12,  family="serif"))
  
  return(plt)
}

age_cuts = seq(0, 80, by=5) # Specify a vector of cut-off values for each age partition
# Stitch adjacent age cut-offs together to generate the interval labels
age_lbs = glue("{age_cuts[1:length(age_cuts) - 1]}-{age_cuts[2: length(age_cuts)]}")
age_lbs[length(age_lbs)] = "75+"

age_lbs = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+")

mat_plot = function(mat, x_lb, y_lb, rnames=age_lbs, cnames=age_lbs, val_lb = "Daily contacts", lim=NULL){
  
  data = expand.grid(x = rnames, y = cnames)
  data$vals = matrix(mat, ncol=1, nrow=nrow(mat) * ncol(mat))
  
  names(data) = c("x","y", "z")
  
  #print(data)
  
  plt = ggplot(data, aes(x=x, y=y, fill=z)) +
    geom_tile() +
    scale_fill_viridis_c(limits=lim) +
    xlab(x_lb) +
    ylab(y_lb) +
    labs(fill=val_lb) +
    theme(axis.text.x = element_text(angle = 75, hjust=1),
          text = element_text(size=12,  family="serif"),
          panel.spacing = unit(2, "lines"),
          strip.background =element_blank(),
          strip.text = element_text(size=12))
  
  return(plt)
  
  
}
