include("utils\\import_utils.jl")
include("utils\\NGM_utils.jl")
include("utils\\final_size_utils.jl")
include("utils\\optim_utils.jl")
include("baseline_params.jl")

using DataFrames
using CSV



function optim_master(;countries = import_country_list(),
					   targets = ["Reff", "doses", "hospitalizations", "deaths"],
					   vac_modes = ["infection", "disease"],
					   rel_eff = [1.0, 0.5],
					   rel_inf = [0.25, 0.5, 0.75],
					   doses = 0.2:0.2:0.8,
					   output_dir = "..\\outputs\\all_countries\\",
					   iterations = 10^3,
					   nsol = 100,
					   source = "combined",
					   writeNGM = true,
					   include_baseline = true)


	for country in countries

		# Create new output directory for country
		out_dir = string(output_dir, country)
		if ~isdir(out_dir)
			mkdir(out_dir)
		end
	
		if writeNGM
			# Calculate and save next-generation matrix
			NGM_country = calc_NGM(import_params(country, source=source))
			CSV.write(string(out_dir, "\\", country, "_NGM.csv"), convert(DataFrame, NGM_country))
		end
		
		
		for mode in vac_modes
			
			if include_baseline
			
			for re in rel_eff
			
				if "Reff" in targets
					
					optima = optimize_Reff([country], re=re, f=0.5, doses, vac_mode=mode, source=source, nsol=nsol, verbose=true, iterations=iterations)
					optima = format_Reff_optima(optima, [country], doses, vac_mode=mode, re=re, f=0.5)
			
					out_file = string(out_dir, "\\", country, "_Reff_", mode, "_")
			
					if re == 1.0
						CSV.write(string(out_file, "constant_rel_eff_50_rel_inf.csv"), optima)
					elseif re == 0.5
						CSV.write(string(out_file, "variable_rel_eff_50_rel_inf.csv"), optima)
					end
					
				end
				
				
				if "hospitalizations" in targets
					
					optima = optimize_hospitalizations([country], re=re, f=0.5, doses, vac_mode=mode, source=source, nsol=nsol, verbose=true, iterations=iterations)
					optima = format_hospitalizations_optima(optima, [country], doses, vac_mode=mode, re=re, f=0.5)
			
					out_file = string(out_dir, "\\", country, "_hospitalizations_", mode, "_")
			
					if re == 1.0
						CSV.write(string(out_file, "constant_rel_eff_50_rel_inf.csv"), optima)
					elseif re == 0.5
						CSV.write(string(out_file, "variable_rel_eff_50_rel_inf.csv"), optima)
					end
					
				end
				
				if "deaths" in targets
					
					optima = optimize_deaths([country], re=re, f=0.5, doses, vac_mode=mode, source=source, nsol=nsol, verbose=true, iterations=iterations)
					optima = format_deaths_optima(optima, [country], doses, vac_mode=mode, re=re, f=0.5)
			
					out_file = string(out_dir, "\\", country, "_deaths_", mode, "_")
			
					if re == 1.0
						CSV.write(string(out_file, "constant_rel_eff_50_rel_inf.csv"), optima)
					elseif re == 0.5
						CSV.write(string(out_file, "variable_rel_eff_50_rel_inf.csv"), optima)
					end
					
				end
				
				if "doses" in targets
					
					optima = optimize_doses([country], re=re, vac_mode=mode, source=source, nsol=nsol, verbose=true, iterations=iterations)
					optima = format_doses_optima(optima, [country], re=re)
			
					out_file = string(out_dir, "\\", country, "_doses_", mode, "_")
			
					if re == 1.0
						CSV.write(string(out_file, "constant_rel_eff_50_rel_inf.csv"), optima)
					elseif re == 0.5
						CSV.write(string(out_file, "variable_rel_eff_50_rel_inf.csv"), optima)
					end
					
				end
				
			end
			
			end
			
			for ri in filter!(x->x≠0.5,rel_inf)
			
				if "Reff" in targets
					
					optima = optimize_Reff([country], re=re, f=ri, doses, vac_mode=mode, source=source, nsol=nsol, verbose=true, iterations=iterations)
					optima = format_Reff_optima(optima, [country], doses, vac_mode=mode, re=re, f=ri)
			
					out_file = string(out_dir, "\\", country, "_Reff_", mode, "_")
			
					if ri == 0.25
						CSV.write(string(out_file, "constant_rel_eff_", 25, "_rel_inf.csv"), optima)
					elseif ri == 0.75
						CSV.write(string(out_file, "constant_rel_eff_", 75, "_rel_inf.csv"), optima)
					end
					
				end
				
				if "hospitalizations" in targets
					
					optima = optimize_hospitalizations([country], re=re, f=ri, doses, vac_mode=mode, source=source, nsol=nsol, verbose=true, iterations=iterations)
					optima = format_hospitalizations_optima(optima, [country], doses, vac_mode=mode, re=re, f=ri)
			
					out_file = string(out_dir, "\\", country, "_hospitalizations_", mode, "_")
			
					if ri == 0.25
						CSV.write(string(out_file, "constant_rel_eff_", 25, "_rel_inf.csv"), optima)
					elseif ri == 0.75
						CSV.write(string(out_file, "constant_rel_eff_", 75, "_rel_inf.csv"), optima)
					end
					
				end
				
				if "deaths" in targets
					
					optima = optimize_deaths([country], re=re, f=ri, doses, vac_mode=mode, source=source, nsol=nsol, verbose=true, iterations=iterations)
					optima = format_deaths_optima(optima, [country], doses, vac_mode=mode, re=re, f=ri)
			
					out_file = string(out_dir, "\\", country, "_deaths_", mode, "_")
			
					if ri == 0.25
						CSV.write(string(out_file, "constant_rel_eff_", 25, "_rel_inf.csv"), optima)
					elseif ri == 0.75
						CSV.write(string(out_file, "constant_rel_eff_", 75, "_rel_inf.csv"), optima)
					end
					
				end
			
			end
			
		end
		
	end

end
