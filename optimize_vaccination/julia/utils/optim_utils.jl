using Optim
using DataFrames
using Statistics

include("NGM_utils.jl")
include("final_size_utils.jl")



function Reff_cost(x, p, doses; vac_mode="infection")
	
	#x = [calc_rem_doses(x, p.S, doses) x']
	x = [calc_rem_doses(x, p.N0, doses) x']

	if (any(x .< 0) || any(x .> 1))
		return Inf
	end

	return calc_Reff(x, p; vac_mode)

end



function doses_cost(x, p; vac_mode="infection")

	if (any(x .< 0) || any(x .> 1))
		return Inf
	end

	Reff = calc_Reff(x, p; vac_mode)
	
	if Reff < 1
		return sum(repeat(x, inner=2) .* p.S) / sum(p.S)
	else
		return Inf
	end
end



function hospitalizations_cost(x, p, doses; vac_mode="infection")

	#x = [calc_rem_doses(x, p.S, doses) x']
	x = [calc_rem_doses(x, p.N0, doses) x']
	
	if (any(x.<0) || any(x.>1))
		return Inf
	end
	
	# Check that Reff > 1
	#Reff = calc_Reff(x, p; vac_mode)
	#if Reff < 1
	#	return 0
	#else
		return sum(calc_final_hospitalizations(x, p; vac_mode))
	#end
end



function deaths_cost(x, p, doses; vac_mode="infection")

	#x = [calc_rem_doses(x, p.S, doses) x']
	x = [calc_rem_doses(x, p.N0, doses) x']
	
	if (any(x.<0) || any(x.>1))
		return Inf
	end
	
	# Check that Reff > 1
	#Reff = calc_Reff(x, p; vac_mode)
	#if Reff < 1
	#	return 0
	#else
		return sum(calc_final_mortality(x, p; vac_mode))
	#end
end



function calc_rem_doses(x, S, doses)

	rem_doses = doses * sum(S) - sum(repeat(x, inner=2) .* S[3:16])

	return rem_doses / sum(S[1:2])

end



function optimize_vaccination(cost_fcn, x_init, p; keep="all", scale=0.1, nsol=1, iterations=1000)

	objective_fcn(x) = cost_fcn(x, p)
	
	optima = fill(NaN, nsol, length(x_init) + 1)
	
	for i = 1:nsol
	
		optimum = optimize(objective_fcn, x_init, SimulatedAnnealing(neighbor=(x,y)->my_neighbor!(x,y,scale), keep_best=true), Optim.Options(iterations=iterations))
		#optimum = optimize(objective_fcn, zeros(8), ones(8), x_init, SAMIN(rt=0.3), Optim.Options(iterations=iterations))
		
		println("optimum: ", optimum)
		
		optima[i, :] = [optimum.minimum optimum.minimizer']
		
	end

	optima = optima[sortperm(optima[:, 1]), :]

	if keep == "best"
		return optima[1,:]
	elseif keep == "all"
		return optima
	end
end


function my_neighbor!(x::AbstractArray{T}, x_proposal::AbstractArray, scale) where T
    @assert size(x) == size(x_proposal)
    for i in 1:length(x)
        #@inbounds x_proposal[i] = x[i] + T(scale * randn()) # workaround because all types might not have randn
		@inbounds x_proposal[i] = minimum([1, maximum([0, x[i] + T(scale * randn())])])
    end
    return
end



function optimize_doses(countries; re=1.0, f=0.5, keep="all", vac_mode="infection", source="combined", sf_mode="var", scale=0.1, nsol=10, iterations = 10000, verbose=false)

	optima = fill(NaN, nsol, length(countries), 11)
	
	for i in 1:length(countries)
		
		p_country = import_params(countries[i], source=source, sf_mode=sf_mode, f=f)
		
		p_country = merge(p_country, (;re))
	
		HIT_country = calc_HIT(p_country)
		
		Reff_test = calc_Reff(fill(HIT_country,8), p_country)
		
		if Reff_test < 1
			x_init = fill(HIT_country, 8)
		else
			x_init = fill(0.999, 8)
		end
	
		optima_country = optimize_vaccination((x,p)->doses_cost(x,p;vac_mode), x_init, p_country; keep=keep, scale=scale, nsol=nsol, iterations=iterations)

		optima[:, i, :] = [fill(calc_R0(p_country), nsol) fill(HIT_country, nsol) optima_country]
	
		if verbose
			println("Country: ", countries[i])
			println("Optimum: ", optima[1, i, :])
		end
	
	end
	
	if keep == "best"
		return optima[1, :, :]
	elseif keep == "all"
		return optima
	end
end



function format_doses_optima(optima, countries; re=1.0)

	res = convert(DataFrame, reshape(optima, (size(optima)[1] * size(optima)[2], size(optima)[3])))

	insert!(res, 1, repeat(countries, inner=size(optima)[1]), :country)

	colnames = ["country","R0","HIT","optimum","0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+"]

	names!(res, Symbol.(colnames))

	return res

end



function optimize_Reff(countries, doses; re=1.0, f=0.5, source="combined", keep="all", scale=0.1, vac_mode="infection", sf_mode="var", nsol=10, iterations=10000, verbose=false)

	optima = fill(NaN, nsol, length(countries), length(doses), 10)
	
	for i in 1:length(countries)
	
		p_country = import_params(countries[i], source=source, sf_mode=sf_mode, f=f)
		
		p_country = merge(p_country, (;re=re))	
				
		for j in 1:length(doses)
		
			Reff_country = calc_Reff(fill(doses[j],8), p_country, vac_mode=vac_mode)
		
			optima_country = optimize_vaccination((x,p)->Reff_cost(x,p,doses[j];vac_mode), fill(doses[j], 7), p_country; keep=keep, scale=scale, nsol=nsol, iterations=iterations)
					
			rem_doses = map(1:nsol) do row
							#calc_rem_doses(optima_country[row, 2:end], p_country.S, doses[j])
							calc_rem_doses(optima_country[row, 2:end], p_country.N0, doses[j])
						end
		
			optima[:, i, j, :] = [fill(Reff_country, nsol) optima_country[:,1] rem_doses optima_country[:,2:end] ]	
		
			if verbose
				println("Country: ", countries[i], "; dose: ", doses[j])
				println("Optimum: ", optima[:, i, j, :])
			end
			
		end
		
	end

	return optima
	
end



function format_Reff_optima(optima, countries, doses; vac_mode="infection", re=1.0, f=0.5)

	res = convert(DataFrame, reshape(optima, (size(optima)[1] * size(optima)[2] * size(optima)[3], size(optima)[4])))

	insert!(res, 2, repeat(doses, inner=size(optima)[1] * size(optima)[2]), :dose)

	insert!(res, 1, repeat(countries, inner=size(optima)[1], outer=size(optima)[3]), :country)
	
	colnames = ["country","Reff_untargeted","dose","Reff_optimized","0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+"]

	names!(res, Symbol.(colnames))
	
	R0 = fill(NaN, length(countries))
	Reff_untargeted = fill(NaN, length(countries))
	
	for i in 1:length(countries)
		p_country = merge(import_params(countries[i], f=f), (;re=re))
		R0[i] = calc_R0(p_country) 
		Reff_untargeted[i] = calc_Reff(fill(1,8), p_country, vac_mode = vac_mode)
	end

	end_points = convert(DataFrame, [R0 zeros(length(countries),1) R0 zeros(length(countries),8); Reff_untargeted ones(length(countries),1) Reff_untargeted ones(length(countries),8)])
	
	insert!(end_points, 1, repeat(countries, outer=2), :country)
	
	names!(end_points, Symbol.(colnames))
	
	return [res; end_points]

end



function optimize_deaths(countries, doses; re=1.0, f=0.5, source="combined", keep="all", scale=0.1, vac_mode="infection", sf_mode="var", nsol=10, iterations=10000, verbose=false)

	optima = fill(NaN, nsol, length(countries), length(doses), 10)
	
	for i in 1:length(countries)
	
		p_country = import_params(countries[i], source=source, sf_mode=sf_mode, f=f)
		
		p_country = merge(p_country, (;re=re))
						
		for j in 1:length(doses)
			
			deaths_country = sum(calc_final_mortality(fill(doses[j], 8), p_country; vac_mode=vac_mode))
			
			optima_country = optimize_vaccination((x,p)->deaths_cost(x,p,doses[j];vac_mode), fill(doses[j], 7), p_country; keep=keep, scale=scale, nsol=nsol, iterations=iterations)
			
			rem_doses = map(1:nsol) do row
							#calc_rem_doses(optima_country[row, 2:end], p_country.S, doses[j])
							calc_rem_doses(optima_country[row, 2:end], p_country.N0, doses[j])
						end
		
			optima[:, i, j, :] = [fill(deaths_country, nsol) optima_country[:,1] rem_doses optima_country[:,2:end] ]	
				
			if verbose
				println("Country: ", countries[i], "; dose: ", doses[j])
				println("Optimum: ", optima[:, i, j, :])
			end
			
		end
		
	end

	return optima
	
end



function format_deaths_optima(optima, countries, doses; vac_mode="infection", re=1.0, f=0.5)

	res = convert(DataFrame, reshape(optima, (size(optima)[1] * size(optima)[2] * size(optima)[3], size(optima)[4])))

	insert!(res, 2, repeat(doses, inner=size(optima)[1] * size(optima)[2]), :dose)

	insert!(res, 1, repeat(countries, inner=size(optima)[1], outer=size(optima)[3]), :country)
	
	colnames = ["country","deaths_untargeted","dose","deaths_optimized","0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+"]

	names!(res, Symbol.(colnames))
	
	deaths_initial = fill(NaN, length(countries))
	deaths_final = fill(NaN, length(countries))
	
	for i in 1:length(countries)
		p_country = merge(import_params(countries[i], f=f), (;re=re))
		deaths_initial[i] = sum(calc_final_mortality(fill(0,8), p_country, vac_mode=vac_mode))
		deaths_final[i] = sum(calc_final_mortality(fill(1,8), p_country, vac_mode=vac_mode))
	end

	end_points = convert(DataFrame, [deaths_initial zeros(length(countries),1) deaths_initial zeros(length(countries),8); deaths_final ones(length(countries),1) deaths_final ones(length(countries),8)])
	
	insert!(end_points, 1, repeat(countries, outer=2), :country)
	
	names!(end_points, Symbol.(colnames))
	
	return [res; end_points]

end



function optimize_hospitalizations(countries, doses; re=1.0, f=f, source="combined", keep="all", scale=0.1, vac_mode="infection", sf_mode="var", nsol=10, iterations=10000, verbose=false)

	optima = fill(NaN, nsol, length(countries), length(doses), 10)
	
	for i in 1:length(countries)
	
		p_country = import_params(countries[i], source=source, sf_mode=sf_mode, f=f)
		
		p_country = merge(p_country, (;re=re))
						
		for j in 1:length(doses)
			
			hospitalizations_country = sum(calc_final_hospitalizations(fill(doses[j], 8), p_country; vac_mode=vac_mode))
			
			optima_country = optimize_vaccination((x,p)->hospitalizations_cost(x,p,doses[j];vac_mode), fill(doses[j], 7), p_country; keep=keep, scale=scale, nsol=nsol, iterations=iterations)
			
			rem_doses = map(1:nsol) do row
							#calc_rem_doses(optima_country[row, 2:end], p_country.S, doses[j])
							calc_rem_doses(optima_country[row, 2:end], p_country.N0, doses[j])
						end
		
			optima[:, i, j, :] = [fill(hospitalizations_country, nsol) optima_country[:,1] rem_doses optima_country[:,2:end] ]	
				
			if verbose
				println("Country: ", countries[i], "; dose: ", doses[j])
				println("Optimum: ", optima[:, i, j, :])
			end
			
		end
		
	end

	return optima
	
end



function format_hospitalizations_optima(optima, countries, doses; vac_mode="infection", re=1.0, f=f)

	res = convert(DataFrame, reshape(optima, (size(optima)[1] * size(optima)[2] * size(optima)[3], size(optima)[4])))

	insert!(res, 2, repeat(doses, inner=size(optima)[1] * size(optima)[2]), :dose)

	insert!(res, 1, repeat(countries, inner=size(optima)[1], outer=size(optima)[3]), :country)
	
	colnames = ["country","hospitalizations_untargeted","dose","hospitalizations_optimized","0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+"]

	names!(res, Symbol.(colnames))
	
	hospitalizations_initial = fill(NaN, length(countries))
	hospitalizations_final = fill(NaN, length(countries))
	
	for i in 1:length(countries)
		p_country = merge(import_params(countries[i], f=f), (;re=re))
		hospitalizations_initial[i] = sum(calc_final_hospitalizations(fill(0,8), p_country, vac_mode=vac_mode))
		hospitalizations_final[i] = sum(calc_final_hospitalizations(fill(1,8), p_country, vac_mode=vac_mode))
	end

	end_points = convert(DataFrame, [hospitalizations_initial zeros(length(countries),1) hospitalizations_initial zeros(length(countries),8); hospitalizations_final ones(length(countries),1) hospitalizations_final ones(length(countries),8)])
	
	insert!(end_points, 1, repeat(countries, outer=2), :country)
	
	names!(end_points, Symbol.(colnames))
	
	return [res; end_points]

end



function summarize_var_optima(optima, countries; trim=true)
	
	out = map([0, 0.025, 0.25, 0.5, 0.75, 0.975, 1]) do q
           quantile3(optima, q; trim)
       end
	
	out = convert(DataFrame, hcat(out...)')
	
	qnames = ["q0", "q025", "q25", "q50", "q75", "q975", "q1"]
	
	colnames = ["R0","HIT","optimum","0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+"]

	names!(out, Symbol.(colnames))
	
	insert!(out, 1, repeat(qnames, inner=size(optima)[2]), :quantile)
	
	insert!(out, 1, repeat(countries, outer=7), :country)
	
	return out[sortperm(out[:, :country]), :]
	
end


function test_median(op)

	med = op[op[:quantile] .== "q50", :]
	
	countries = med[:country]
	
	#println(countries)
	
	x = convert(Array, med[:, 5:12])
	
	#p_countries = map(import_params, med[:country])
	
	Reff = fill(NaN, length(countries))
	
	
	for i in 1:length(countries)
	
		p_country = import_params(countries[i])
				
		Reff[i] = calc_Reff(x[i,:], p_country)

	end
	
	return DataFrame(country = med[:country], Reff = Reff)

end


function quantile3(x, percentile; trim=true)

	q = map(1:size(x)[2]) do col
			
			if trim
				cut_off = findlast(x[:,col,2] .< 1.05 * x[1,col,2])
			else
				cut_off = size(x)[1]
			end
			
			map(1:size(x)[3]) do slice
				quantile!(x[1:cut_off, col, slice], percentile)
			end
		end

	#println(q)

	return hcat(q...)

end

function find_cutoff(x)

	sort!(x, dims=1)
	
	return 

end