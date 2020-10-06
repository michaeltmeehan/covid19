include("utils\\import_utils.jl")
include("utils\\NGM_utils.jl")
include("utils\\optim_utils.jl")
include("baseline_params.jl")

using DataFrames
using CSV

# Set contact data source
source = "combined"

# Import list of countries
countries = import_country_list()
#countries = ["India", "China", "United Kingdom"]


vac_modes = ["infection", "disease"]

rel_eff = [1.0, 0.5]

iterations = 10^4
nsol = 50

doses = 0.2:0.2:0.8

# Loop over all countries
for country in countries

	# Create new output directory for country
	out_dir = string("..\\outputs\\all_countries\\", country)
	#out_dir = string("..\\outputs\\sensitivity\\uniform_both\\all_countries\\", country)
	#out_dir = string("..\\outputs\\sensitivity\\uniform_susceptibility\\all_countries\\", country)
	#out_dir = string("..\\outputs\\sensitivity\\uniform_clinical_fraction\\all_countries\\", country)
	if ~isdir(out_dir)
		mkdir(out_dir)
	end
	
	# Calculate and save next-generation matrix
	NGM_country = calc_NGM(import_params(country, source=source))
	CSV.write(string(out_dir, "\\", country, "_NGM.csv"), convert(DataFrame, NGM_country))
	
	# Fixed disease optima
	for mode in vac_modes
		
		for re in rel_eff
		
			optima = optimize_fixed_doses([country], re=re, doses, vac_mode=mode, source=source, nsol=nsol, verbose=true)
			optima = format_fixed_optima(optima, [country], doses, vac_mode=mode, re=re)
			
			# Calculate Reff reduction under random vaccine allocation for comparison
			uniform = map(0:0.2:1.0) do d 
						 calc_Reff(fill(d, 8), merge(import_params(country, source=source), (;re=re)), vac_mode=mode) 
					  end 			
			uniform = convert(DataFrame, [0:0.2:1.0  uniform])
			names!(uniform, Symbol.(["dose", "Reff"]))
			
			out_file = string(out_dir, "\\", country, "_fixed_doses_", mode, "_")
			
			if re == 1.0
				CSV.write(string(out_file, "constant_eff_optima.csv"), optima)
				CSV.write(string(out_file, "constant_eff_uniform.csv"), uniform)
			elseif re == 0.5
				CSV.write(string(out_file, "variable_eff_optima.csv"), optima)
				CSV.write(string(out_file, "variable_eff_uniform.csv"), uniform)
			end

		end
		
	end
	
	
	# Variable disease optima
	for re in rel_eff

		optima = optimize_var_doses([country], re=re, source=source, keep="all", vac_mode="infection", scale=0.5, nsol=20, iterations=iterations, verbose=true)
		res = format_var_optima(optima, [country])
		
		out_file = string(out_dir, "\\", country, "_optimal_variable_doses_infection_")
			
		if re == 1.0
			CSV.write(string(out_file, "constant_eff.csv"), res)
		elseif re == 0.5
			CSV.write(string(out_file, "variable_eff.csv"), res)
		end
		
	end
		
end