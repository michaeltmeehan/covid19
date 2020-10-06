library(glue)

age_cuts = seq(0, 80, by=5) # Specify a vector of cut-off values for each age partition
# Stitch adjacent age cut-offs together to generate the interval labels
age_lbs = glue("{age_cuts[1:length(age_cuts) - 1]}-{age_cuts[2: length(age_cuts)]}")
age_lbs[length(age_lbs)] = "75+"

mat_plot = function(mat, x_lb, y_lb, rnames=age_lbs, cnames=age_lbs, val_lb = "Contacts / day"){
  
  data = expand.grid(x = rnames, y = cnames)
  data$vals = matrix(mat, ncol=1, nrow=length(mat))
  
  names(data) = c("x","y", "z")
  
  print(data)
  
  plt = ggplot(data, aes(x=x, y=y, fill=z)) +
    geom_tile() +
    scale_fill_viridis_c() +
    xlab(x_lb) +
    ylab(y_lb) +
    theme(axis.text.x = element_text(angle = 75, hjust=1),
          text = element_text(size=16,  family="serif"),
          panel.spacing = unit(2, "lines"),
          strip.background =element_blank(),
          strip.text = element_text(size=16))
  
  return(plt)

  
}
