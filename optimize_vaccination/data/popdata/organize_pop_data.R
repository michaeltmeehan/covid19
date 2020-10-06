library(ggplot2)
library(dplyr)



source("popdata/pop_utils.R")


poptotal = read.csv("popdata/poptotal.csv")

poptotal$median_age = apply(poptotal[,4:25], 1, calc_median_age)

poptotal = poptotal %>% 
  mutate("age75+" = age75 + age80 + age85 + age90 + age95 + age100) %>%
  select("iso3c", "country", "year", "median_age", "total",  "age0", "age5", "age10", "age15", "age20", "age25", "age30", "age35", "age40", "age45", "age50", "age55", "age60", "age65", "age70", "age75+")

write.csv(poptotal, "popdata/poptotal_summarized.csv", row.names=FALSE)

