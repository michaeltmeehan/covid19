# Optimizing age-specific COVID-19 vaccination

The contents of this repository contain all of the necessary input data and code required to reproduce the results of the research article: *Age-targeted dose allocation can halve COVID-19 vaccine requirements* available on [medRxiv](https://www.medrxiv.org/content/10.1101/2020.10.08.20208108v1).

## Input data

Country-specific demographic and contact data can be found in the data/popdata and data/contactdata folders respectively, alongside estimates of the basic reproduction number in each country (data/R0). The data folder also contains original data sources used to derive the input data, including the raw London School of Tropical Health Medicine time-series estimates of Reff, as well as the calibrated values of the age-dependent susceptibility and clinical fraction.

## Running simulations

The code used to calculate the optimal vaccine allocation policies in each country is written in the Julia programming language (v1.5.0). To run the code, open a Julia console and navigate to the `julia` directory. From here, execute the command: 

`include("run_all_countries.jl")`

This will run the optimization search routine for all 179 countries included in the analysis. Outputs are stored in a newly created directory `output_data` and each country's results can be found in the country-specific folder. If you only wish to run the analysis for a subset of countries, you can edit the `run_all_countries.jl` julia script at line 13.

The baseline parameters used in the optimization routine can be found in the `baseline_params.jl` script, and the necessary helper functions required to execute the optimization search (including the functions used to import the input data and calculate the next-generation matrix and reproduction numbers) are stored within the julia/utils folder.
