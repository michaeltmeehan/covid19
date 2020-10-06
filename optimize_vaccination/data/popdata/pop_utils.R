library(plyr)

read_country_data = function(cname){
  read.csv("C:/Users/jc213439/Dropbox/Optimize_vaccination/data/original/popdata.txt") %>%
    dplyr::filter(country == cname) %>% 
    dplyr::mutate(agegroup = cut(age0, 
                                 breaks = c(seq(-0.5, 74.5, by = 5), Inf),
                                 labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+")
    )
    ) %>%
    dplyr::select(agegroup, population) %>%
    dplyr::group_by(agegroup) %>%
    dplyr::summarize(population = sum(population)) %>%
    dplyr::ungroup() # %>%
    #dplyr::pull(population)
}


make_pop_pyramid = function(cname){
  cdata = read.csv("C:/Users/jc213439/Dropbox/Optimize_vaccination/data/original/popdata.txt") %>%
    dplyr::filter(country == cname) %>% 
    dplyr::mutate(agegroup = cut(age0, 
                                 breaks = c(seq(-0.5, 74.5, by = 5), Inf),
                                 labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+")
    )
    ) %>%
    dplyr::mutate(Population = ifelse(sex == "M", -population, population),
                  Gender = ifelse(sex=="M", "Male", "Female")) %>%
    dplyr::select(agegroup, Gender, Population)
  
  cdata$Age = factor(cdata$agegroup, levels=cdata$agegroup, labels=cdata$agegroup)
  cdata$Proportion = cdata$Population / sum(abs(cdata$Population))
  
  plt = ggplot(cdata, aes(x = Age, y = Proportion, fill = Gender)) + 
    geom_bar(subset = .(Gender == "Female"), stat = "identity") + 
    geom_bar(subset = .(Gender == "Male"), stat = "identity") + 
    scale_y_continuous(breaks = seq(-0.3, 0.3, 0.01),
                       labels = paste0(as.character(c(seq(30, 0, -1), seq(1, 30, 1))), "%"),
                       limits=c(-0.05, 0.05)) +
    coord_flip() + 
    scale_fill_brewer(palette = "Set1") + 
    theme_bw() +
    theme(text = element_text(size=16,  family="serif"))
  
  return(plt)
}


calc_median_age = function(cname){
  
  cdata = read.csv("C:/Users/jc213439/Dropbox/Optimize_vaccination/data/original/popdata.txt") %>%
    dplyr::filter(country == cname) %>% 
    dplyr::select(age0, population) %>%
    dplyr::group_by(age0) %>%
    dplyr::summarize(population = sum(population)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(age = age0 + 2.5) %>%
    dplyr::select(age, population)
  
  return(calc_median(cdata))
  
}

calc_median = function(freq_table){
  
  freq_table = freq_table %>%
                  dplyr::arrange(age) %>%
                  dplyr::mutate(cum_population = cumsum(population))
  
  tot_pop = sum(freq_table$population)
  
  n = ceiling(tot_pop / 2)
  
  for (idx in 1:nrow(freq_table)){
    median = freq_table$age[idx]
    if (freq_table$cum_population[idx] > n){
      break
    }
  }
  return(median)
}


calc_median_age = function(pop_vec){
  
  age_mid_points = seq(2.5, 102.5, 5)

  n = ceiling(pop_vec["total"] / 2)

  cum_pop = cumsum(as.numeric(pop_vec[1:(length(pop_vec) - 1)]))
  
  for (idx in 1:length(cum_pop)){
    median = age_mid_points[idx]
    if (cum_pop[idx] > n){
      break
    }
  } 
  return(median)
}
