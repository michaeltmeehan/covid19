# https://www.sharpsightlabs.com/blog/map-oil-production-country-r/

# load libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# load data
opti_df_constant <- read.csv("global_doses_infection_constant_rel_eff_50_rel_inf_filtered.csv")
opti_df_variable <- read.csv("global_doses_infection_variable_rel_eff_50_rel_inf_filtered.csv")

# format and combine constant and variable data
colnames(opti_df_constant)[4:5] <- paste0("constant_", colnames(opti_df_constant)[4:5])
colnames(opti_df_variable)[4:5] <- paste0("variable_", colnames(opti_df_variable)[4:5])

opti_df <- merge(opti_df_constant[,c("country", "constant_HIT", "constant_optimum")],
                 opti_df_variable[,c("country", "variable_HIT", "variable_optimum")],
                 by = "country")

# get world map data
map.world <- map_data('world')

# find mismatched country names
mismatches <- anti_join(opti_df, map.world, by = c('country' = 'region'))
unique(mismatches$country)

# list map.world country names
map.world %>%
  group_by(region) %>%
  summarise() %>%
  print(n = Inf)

# update mismatched country names in dataset
opti_updated <- opti_df %>% 
  mutate(country = recode(country
                          , 'Antigua and Barbuda' = "Antigua"
                          , 'Brunei Darussalam' = "Brunei"
                          , 'Cote dIvoire' = "Ivory Coast"
                          , 'Congo - Kinshasa' = "Democratic Republic of the Congo"
                          , 'Congo - Brazzaville' = "Republic of Congo"
                          , "Czechia" = "Czech Republic"
                          , 'United Kingdom' = "UK"
                          , 'Guinea Bissau' = "Guinea-Bissau"
                          , 'China, Hong Kong SAR' = "China"
                          , "Lao People's Democratic Republic" = "Laos"
                          , "North Macedonia" = "Macedonia"
                          , 'Eswatini' = "Swaziland"
                          , 'Trinidad and Tobago' = "Trinidad"
                          , "United Republic of Tanzania" = "Tanzania"
                          , "United States of America" = "USA"))

# fill Inf values with 1.2 (> max value)
opti_updated$constant_optimum[is.finite(opti_updated$constant_optimum) == FALSE] <- 1.2
opti_updated$variable_optimum[is.finite(opti_updated$variable_optimum) == FALSE] <- 1.2

# set all HIT > 1 to 1.2 (to bin all HITs above 100% for color breaks)
opti_updated$constant_HIT[opti_updated$constant_HIT >= 1] <- 1.2
opti_updated$variable_HIT[opti_updated$variable_HIT >= 1] <- 1.2

# create difference metrics and replace NAs with zero
opti_updated$constant_difference <- abs(opti_updated$constant_optimum - opti_updated$constant_HIT)
opti_updated$variable_difference <- abs(opti_updated$variable_optimum - opti_updated$variable_HIT)

# join data with world map data
coverage <- left_join(map.world, opti_updated, by = c('region' = 'country'))

# create mapping function
mapFun <- function(response, titleName, standardColors, standardLabels, standardBreaks, legend, legPos){
  map <- ggplot(coverage, aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill = response)) +
    scale_fill_gradientn(colours = standardColors
                         ,labels = standardLabels
                         ,breaks = standardBreaks
    ) +
    guides(fill = guide_legend(reverse = T)) +
    labs(fill = 'Population-level\ncoverage'
         ,x = NULL
         ,y = NULL
         # ,subtitle = titleName
         ) +
    theme(text = element_text(color = '#EEEEEE')
          ,plot.title = element_text(size = 28)
          ,plot.subtitle = element_text(size = 14)
          ,axis.ticks = element_blank()
          ,axis.text = element_blank()
          ,panel.grid = element_blank()
          ,panel.background = element_rect(fill = '#333333')
          ,plot.background = element_rect(fill = '#333333')
          ,legend.position = 'none'
    )
  if(legend == TRUE){
    map <- map + 
      theme(legend.position = legPos,
            ,legend.background = element_blank()
            ,legend.key = element_blank()
            )
  }
  map + 
}

standardColors <- c("#440154FF", "#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF", "#FDE725FF")
standardLabels <- c("<20 %", "20 - 39 %", "40 - 59 %", "60 - 79 %", "80 - 100 %", ">100 %\n(other strategies req'd)")
standardBreaks <- c(0.19999, 0.39999, 0.59999, 0.79999, 0.99999, 1.19999)
standardLegPosition <- c(.10,.37)

# create maps with uniform vaccination ---------------------------------------
CHITmap <- mapFun(coverage[, "constant_HIT"], titleName = "", standardColors, standardLabels, standardBreaks, legend = TRUE, legPos = standardLegPosition)
COmap = mapFun(response = coverage[,"constant_optimum"], titleName = "", standardColors, standardLabels, standardBreaks, legend = FALSE, legPos = standardLegPosition)
constantMapLong <- grid.arrange(CHITmap, COmap) #, CDiffMap
ggsave(file = 'constant_map_2long.pdf', constantMapLong)

# create maps with targeted vaccination ---------------------------------------
VHITmap <- mapFun(coverage[, "variable_HIT"], titleName = "", standardColors, standardLabels, standardBreaks, legend = TRUE, legPos = c(0.13, 0.40))
VOmap = mapFun(response = coverage[,"variable_optimum"], titleName = "", standardColors, standardLabels, standardBreaks, legend = FALSE)
varyingMapLong <- grid.arrange(VHITmap, VOmap) #, CDiffMap
ggsave(file = 'varying_map_2long.pdf', varyingMapLong)