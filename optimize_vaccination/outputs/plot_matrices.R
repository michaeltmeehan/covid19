
source("../utils/R/plot_utils.R")

countries = list("India", "China", "United Kingdom")
names(countries) = countries

import_contact_matrix = function(cname){
  
  contact_raw = read.csv(paste0("contactdata/survey/", cname, "/", cname, "_all.csv"))[,2:17]
  
  contact_raw = as.matrix(contact_raw)
  
  return(contact_raw)
  
}

contact_matrices = lapply(countries, import_contact_matrix)

contact_plots = lapply(contact_matrices, mat_plot, x_lb = "Age of contactor", y_lb = "Age of contactee")


import_NGM_matrix = function(cname){
  
  NGM_raw = read.csv(paste0("../outputs/country_data/NGM_", cname, ".csv"))
  
  NGM_raw = as.matrix(NGM_raw)
  
  return(NGM_raw)
  
}

NGMs = lapply(countries, import_NGM_matrix)

NGM_plots = lapply(NGMs, mat_plot, x_lb="Age of infector", y_lb="Age of infectee", val_lb="Infections", lim=c(0,1.5))
