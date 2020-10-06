include("utils\\import_utils.jl")
include("utils\\NGM_utils.jl")
include("utils\\optim_utils.jl")
include("baseline_params.jl")

using DataFrames
using CSV

# Import list of countries
#countries = import_country_list().country
#countries = ["India", "China", "United Kingdom"]
#countries = ["China"]

#countries = ["New Zealand", "Nepal", "North Macedonia"]

countries = ["Canada"]

# Loop over relative efficacies by age

vac_mode = "disease"

rel_eff = 1.0

iterations = 10^5
nsol = 100


doses = 0.2:0.2:0.8

optima = optimize_fixed_doses(countries, re=rel_eff, doses, vac_mode=vac_mode, source="prem_2020", nsol=nsol, verbose=true)

res = format_fixed_optima(optima, countries, doses, re=rel_eff)

#CSV.write("..\\outputs\\fixed_dose\\example_countries_disease_constant_eff_res.csv", res)  

