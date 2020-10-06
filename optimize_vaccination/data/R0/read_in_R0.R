

rt_1 = read.csv("original/lshtm_rt_estimates_1.txt")
rt_2 = read.csv("original/lshtm_rt_estimates_2.txt")

# Rename countries in rt_2 to match rt_1
rt_2$country[rt_2$country == "Bosnia & Herzegovina"] = "Bosnia and Herzegovina"
rt_2$country[rt_2$country == "United States" ] = "United States of America"
rt_2$country[rt_2$country == "CÃ´te dâ???TIvoire" ] = "Cote dIvoire"
rt_2$country[rt_2$country == "Guinea-Bissau" ] = "Guinea Bissau"
rt_2$country[rt_2$country == "Palestinian Territories" ] = "Palestine"

# Clean up other rt_2 names
rt_2$country[rt_2$country == "SÃ£o TomÃ© & PrÃ?ncipe" ] = "Sao Tome and Principe"
rt_2$country[rt_2$country == "Trinidad & Tobago" ] = "Trinidad and Tobago"
rt_2$country[rt_2$country == "Antigua & Barbuda" ] = "Antigua and Barbuda"
rt_2$country[rt_2$country == "Turks & Caicos Islands" ] = "Turks and Caicos Islands"

# Rename countries in rt_1 to match rt_2
rt_1$country[rt_1$country == "Congo"] = "Congo - Brazzaville"
rt_1$country[rt_1$country == "Democratic Republic of the Congo"] = "Congo - Kinshasa"


rt = rbind(rt_1, rt_2[,-3])



extract_R0 = function(cname, rt_table=rt, val = "upper_90"){
  
  if (cname %in% unique(rt_table$country)){
    return(rt_table %>% filter(country == cname) %>% select(all_of(val)) %>% max())
  }else{
    return(NA)
  }
}


R0 = sapply(unique(rt$country), extract_R0)

out = data.frame(country=names(R0), R0 = unname(R0))

write.csv(out, "R0/LSHTM_R0_estimates.csv", row.names=FALSE)
