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
import.pop.data = function(country) {
  popdata = read.csv("../data/popdata/poptotal_summarized.csv")
  population_size = as.numeric(popdata[popdata$country == country, 6:21])
  names(population_size) = age_groups
  return(population_size)
}


# Import life expectancy data
import.life.expectancy = function(country) {
  life_expectancy_data = read.csv("../data/popdata/life_expectancy.csv")
  life_expectancy = life_expectancy_data$life_expectancy
  names(life_expectancy) = age_groups
  return(life_expectancy)
}


# Import country contact data
import.contact.data = function(country) {
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
import.seropositivity.data = function(country) {
  seropositivity_data = openxlsx::read.xlsx("../data/epidata/country_specific_seroprevalence.xlsx")
  seropositivity = as.numeric(seropositivity_data[seropositivity_data$country == country, 2:17])
  names(seropositivity) = age_groups
  return(seropositivity)
}




#### Import strain-specific data ####

# Import susceptibility / severity data
import.epi.data = function(strain) {
  epi_data = openxlsx::read.xlsx("../data/epidata/strain_specific_epi_data.xlsx", sheet = strain)
  epi_data$age = age_groups
  return(epi_data)
}


# Transform risk of hospitalization and death for Delta
transform.odds = function(rate, odds_ratio) {
  
  odds = rate / (1 - rate)
  
  return(odds_ratio * odds / (1 + odds_ratio * odds))
  
}


# Import vaccine efficacy data
import.efficacy.data = function(strain) {
  
  efficacy_data = openxlsx::read.xlsx("../data/epidata/strain_specific_vaccine_efficacy.xlsx",
                                      sheet = strain)
  
  return(efficacy_data)
  
}




#### Master import ####

# Import all parameters
get.params = function(country="Australia",
                      strain = "Delta",
                      R0 = 5,
                      rel_infectious = 0.25,
                      preclinical_period = 2.1,
                      clinical_period = 2.9,
                      asymptomatic_period = 5.0,
                      or_hospitalization = 2.08, # Fisman et al. 95% CI: 1.8 - 2.38
                      or_death = 2.32, # Fisman et al. 95% CI: 1.47 - 3.30
                      seropositivity=NULL) {
  if (is.null(seropositivity)){
    seropositivity = import.seropositivity.data(country)
  }
  
  epi_data = import.epi.data(strain)
  
  params = list(
    population_size = import.pop.data(country),
    life_expectancy = import.life.expectancy(country),
    contact_matrix = import.contact.data(country),
    susceptibility = epi_data$susceptibility,
    clinical_fraction = epi_data$clinical_fraction,
    hospitalization_rate = transform.odds(epi_data$hospitalization_rate, or_hospitalization),
    infection_fatality_rate = transform.odds(epi_data$infection_fatality_rate, or_death),
    efficacy = import.efficacy.data(strain),
    rel_infectious = rel_infectious,
    R0 = R0,
    preclinical_period = preclinical_period,
    clinical_period = clinical_period,
    asymptomatic_period = asymptomatic_period,
    seropositivity = seropositivity,
    scale_factor = 1,
    strain = strain,
    or_hospitalization = or_hospitalization,
    or_death = or_death
  )
  params$scale_factor = calc.scale.factor(params, R0)
  return(params)
}
