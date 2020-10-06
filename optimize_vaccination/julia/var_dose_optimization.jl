include("utils\\import_utils.jl")
include("utils\\NGM_utils.jl")
include("utils\\optim_utils.jl")
include("baseline_params.jl")

using DataFrames
using CSV

# Import list of countries
countries = import_country_list()
#countries = ["Zimbabwe", "China", "United Kingdom"]
#countries = ["India", "China", "United Kingdom"]
#countries = ["India", "United Kingdom"]

#countries = ["New Zealand", "Nepal", "North Macedonia"]
#countries = ["Guinea Bissau"]

countries = ["Canada"]

# Loop over relative efficacies by age

vac_mode = "infection"

rel_eff = 1.0

iterations = 10^5
nsol = 100

optima = optimize_var_doses(countries, re=rel_eff, source="prem_2020", keep="all", vac_mode=vac_mode, scale=0.5, nsol=nsol, iterations=iterations, verbose=true)

res = format_var_optima(optima, countries)

#res_summary = summarize_var_optima(optima, countries)

#if re == 1
#	fdir = string("..\\outputs\\", mode, "\\constant_efficacy\\variable_dose\\")
#elseif re == 0.5
#	fdir = string("..\\outputs\\", mode, "\\variable_efficacy\\variable_dose\\")
#end

#CSV.write(string(fdir, "all_results.csv"), res)
#CSV.write(string(fdir, "sum_results.csv"), res_summary)
