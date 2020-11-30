using CSV
using DataFrames

include("NGM_utils.jl")
include("..\\baseline_params.jl")

function import_country_list()

	countries = CSV.read("..\\data\\country_list.csv")
	
	return countries.consensus_name

end

function import_pop_data(country)

	popdata = CSV.read("..\\data\\popdata\\poptotal_summarized.csv")
	
	return popdata[popdata.country .== country, :]

end


#function import_median_age(;countries=nothing)

#	median_age_table = CSV.read("..\\data\\popdata\\median_age.csv")
	
#	if isnothing(countries)
#		countries = median_age_table.country
#	end
	
#	return median_age_table[in(countries).(median_age_table.country), :age ]
	
	#if all
	#	return median_age_table
	#else
	#	return median_age_table[median_age_table.country .== country, :age]
	#end

#end


function import_contact_data(country; loc="all", as_array=true, source="prem_2020")

	contact_matrix = CSV.read("..\\data\\contactdata\\" * source * "\\" * country * "\\" * country * "_" * loc * ".csv")
	
	if as_array
		contact_matrix = convert(Matrix, contact_matrix[:, 2:end])
	end
	
end


function import_R0_data(; all=false, country=nothing)

	#R0_table = CSV.read("..\\data\\R0\\cdc_R0_estimates.csv")
	#R0_table = CSV.read("..\\data\\R0\\LSHTM_R0_estimates.csv")
	R0_table = CSV.read("..\\data\\R0\\combined_R0_estimates.csv")
	
	if isnothing(country)
		all = true
	end
	
	if all
		return R0_table
	else
		return R0_table[R0_table.country .== country, :R0]
	end
end


function import_all_data(country; source="combined")

	# Import population data and extract values
	pop_data = import_pop_data(country)
	S = vec(convert(Array, pop_data[6:21]))
	N0 = S
	median_age = pop_data.median_age
	
	
	# Import contact data
	c = import_contact_data(country; source=source)
	
	# Import CDC R0 estimates
	R0 = import_R0_data(country=country)
	
	
	return (;S,N0,c,R0,median_age)

end


function import_params(country; source="combined", sf_mode="var", f=0.5)

	psub_country = import_all_data(country; source=source)
	
	p_country = merge(p_base, psub_country, (;f=f))
	
	if sf_mode == "fixed"
		sf = 0.06530640076170585 # Consensus fit
		#sf = 1.3605945118527203	# Parametric fit
	elseif sf_mode == "var"
		sf = calc_scale_factor(p_country, psub_country.R0)
	end
	
	p_country = merge(p_country, (;sf,S = (1 .- p_country.pre) .* p_country.S))

end