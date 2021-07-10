# Utility functions for importing country-specific, strain-specific and vaccine-specific data

library(openxlsx)
library(tidyverse)

source("./utils/NGM_utils.R")


age_groups = c(
  "0-4",
  "5-9",
  "10-14",
  "15-19",
  "20-24",
  "25-29",
  "30-34",
  "35-39",
  "40-44",
  "45-49",
  "50-54",
  "55-59",
  "60-64",
  "65-69",
  "70-74",
  "75+"
)


#### Import country-specific data ####
# Import country population data
import_pop_data = function(country) {
  popdata = read.csv("../data/popdata/poptotal_summarized.csv")
  population_size = as.numeric(popdata[popdata$country == country, 6:21])
  names(population_size) = age_groups
  return(population_size)
}


# Import country contact data
import_contact_data = function(country) {
  contact_data = read.csv(paste0(
    "../data/contactdata/prem_2020/",
    country,
    "/",
    country,
    "_all.csv"
  ))
  contact_matrix = as.matrix(contact_data[,-1])
  rownames(contact_matrix) = age_groups
  colnames(contact_matrix) = age_groups
  return(contact_matrix)
}


# Import country seropositivity
import_seropositivity_data = function(country) {
  seropositivity_data = openxlsx::read.xlsx("../data/epidata/country_specific_seroprevalence.xlsx")
  seropositivity = as.numeric(seropositivity_data[seropositivity_data$country == country, 2:17])
  names(seropositivity) = age_groups
  return(seropositivity)
}




#### Import strain-specific data ####
import_epi_data = function(strain) {
  epi_data = openxlsx::read.xlsx("../data/epidata/strain_specific_epi_data.xlsx", sheet = strain)
  epi_data$age = age_groups
  return(epi_data)
}




#### Import strain- and vaccine-specific data ####
import_efficacy_data = function(strain, vaccine_1, vaccine_2 = NULL) {
  efficacy_data = openxlsx::read.xlsx("../data/epidata/strain_specific_vaccine_efficacy.xlsx",
                                      sheet = strain)
  
  if (is.null(vaccine_2)){
    all_ages_efficacy_data = efficacy_data[efficacy_data$vaccine == vaccine_1,]
    Vinf = rep(all_ages_efficacy_data$Vinf, 16)
    Vsev = rep(all_ages_efficacy_data$Vsev, 16)
    Vtrans = rep(all_ages_efficacy_data$Vtrans, 16)  
  } else{
    under_60_efficacy_data = efficacy_data[efficacy_data$vaccine == vaccine_1,]
    over_60_efficacy_data = efficacy_data[efficacy_data$vaccine == vaccine_2,]
    Vinf = c(rep(under_60_efficacy_data$Vinf, 12), rep(over_60_efficacy_data$Vinf, 4))
    Vsev = c(rep(under_60_efficacy_data$Vsev, 12), rep(over_60_efficacy_data$Vsev, 4))
    Vtrans = c(rep(under_60_efficacy_data$Vtrans, 12), rep(over_60_efficacy_data$Vtrans, 4))
  }
  
  return(list(
    "Vinf" = Vinf,
    "Vsev" = Vsev,
    "Vtrans" = Vtrans
  ))
  
}


# Import all parameters
get_params = function(country,
                      vaccine_1 = "BNT162b2",
                      vaccine_2 = NULL,
                      strain = "Delta",
                      R0 = 5,
                      rel_infectious = 0.25,
                      period = list("preclinical" = 2.1,
                                    "clinical" = 2.9,
                                    "asymptomatic" = 5.0),
                      seropositivity=NULL) {
  if (is.null(seropositivity)){
    seropositivity = import_seropositivity_data(country)
  }
  params = list(
    population_size = import_pop_data(country),
    contact_matrix = import_contact_data(country),
    epi = import_epi_data(strain),
    efficacy = import_efficacy_data(strain, vaccine_1, vaccine_2),
    rel_infectious = rel_infectious,
    R0 = R0,
    period = period,
    seropositivity = seropositivity,
    scale_factor = 1,
    strain = strain
  )
  params$scale_factor = calc_scale_factor(params, R0)
  return(params)
}
