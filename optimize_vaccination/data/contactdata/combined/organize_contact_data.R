library(openxlsx)

country_map = read.xlsx("country_map.xlsx")

complete_country_map = country_map[!is.na(country_map$contact_country), ]

load("original/prem_2020/contact_all.rdata")

contact_countries = names(contact_all)

for (cname in contact_countries){
  
  full_cname = complete_country_map$consensus_name[complete_country_map$contact_country == cname]
  
  folder = paste0("contactdata/prem_2020/", full_cname)
  
  # Make new directory for each country in contactdata folder
  dir.create(folder)
  
  write.csv(contact_all[[cname]], paste0(folder, "/", full_cname, "_all.csv"))
  
  
}