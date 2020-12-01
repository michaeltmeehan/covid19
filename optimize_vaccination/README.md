# Optimizing age-specific COVID-19 vaccination

The contents of this repository contain all of the necessary input data and code required to reproduce the results of the research article: *Age-targeted dose allocation can halve COVID-19 vaccine requirements* available on [medRxiv](https://www.medrxiv.org/content/10.1101/2020.10.08.20208108v1).

## Input data

Country-specific demographic and contact data can be found in the data/popdata and data/contactdata folders respectively, alongside estimates of the basic reproduction number in each country (data/R0). The data folder also contains original data sources used to derive the input data, including the raw London School of Tropical Health Medicine time-series estimates of Reff, as well as the calibrated values of the age-dependent susceptibility and clinical fraction.

## Running simulations

The code used to calculate the optimal vaccine allocation policies in each country is written in the Julia programming language (v1.5.0). To run the code, open a Julia console and navigate to the `julia` directory. From here, the master optimization function (`optima_master()`) can be imported by executing the command: 

`include("master.jl")`

Individual optimization searches for combinations of countries, optimization targets, parameters, etc. can be run by calling the `optim_master()` function and supplying several keyword arguments. For example, to optimize the hospitalizations and deaths in China, India and the United Kingdom, run the command

`optima_master(countries=["China", "India", "United Kingdom"], targets=["hospitalizations", "deaths"])

The defaults for each argument can be found in the master.jl script.

The results of the optimization search will be saved as .csv files in country-specific folders in the designated `output_dir`. To run more specialized searches for particular optimization targets you can use the optimization functions imported from the optim_utils.jl script.

The baseline parameters used in the optimization routine can be found in the `baseline_params.jl` script, and the necessary helper functions required to execute the optimization search (including the functions used to import the input data and calculate the next-generation matrix and reproduction numbers) are stored within the julia/utils folder.
