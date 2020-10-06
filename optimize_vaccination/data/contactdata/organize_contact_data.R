
setwd("C:/Users/jc213439/Dropbox/Optimize_vaccination/data/")
source("contactdata/contact_utils.R")

all_matrices = readRDS("original/all_matrices.rds")
country_map = read.csv("country_map.csv")

for (idx in 1:nrow(country_map)){
  
  folder = paste0("contactdata/", country_map$pop_label[idx])
  
  # Make new directory for each country in contactdata folder
  dir.create(folder)
  
  # Extract contact matrices for current country
  contact_matrices = all_matrices[[country_map$contact_label[idx]]]
  
  contact_matrices[["all"]] = Reduce('+', contact_matrices)
  
  for(loc in names(contact_matrices)){
    write.csv(contact_matrices[[loc]], paste0(folder, "/", country_map$pop_label[idx], "_", loc,".csv"))
  }
  
}