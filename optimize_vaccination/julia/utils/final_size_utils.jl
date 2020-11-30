using LinearAlgebra
using NLsolve
include("NGM_utils.jl")
include("..//baseline_params.jl")


# Vectorized final size equation
# x_i = \exp( \sum_j K_{ij} (1 - x_j) )

function fs_res(x; NGM)

	x - exp.(-NGM*(1 .- x))

end


function fs_jac(x; NGM)

	I - diagm(exp.(-NGM*(1 .- x))) * NGM

end


function calc_final_size(NGM)

	f(x) = fs_res(x; NGM=NGM)
	j(x) = fs_jac(x; NGM=NGM)
	
	return nlsolve(f, j, fill(0.5, 16)).zero

end


function calc_final_hospitalizations(x, p; vac_mode="infection")

	# Overall vaccine efficacy is e
	# Relative efficacy for 60+ age groups is re
	x_e = p.e * [x[1:6]' p.re*x[7:8]']

	if vac_mode == "infection"
		S_vac = (1 .- repeat(x_e[:], inner=2) ) .* p.S
		p_vac = merge(p, (;S=S_vac))
	elseif vac_mode == "disease"
		y_vac = (1 .- repeat(x_e[:], inner=2) ) .* p.y
		p_vac = merge(p, (;y=y_vac))
	end

	NGM_FS = calc_NGM_FS(p_vac)
	
	final_size = calc_final_size(NGM_FS)
	
	if vac_mode == "infection"
		return p_vac.S .* (1 .- final_size) .* p.hr' / sum(p_vac.N0) * 1e5
	elseif vac_mode == "disease"
		return p_vac.S .* (1 .- final_size) .* (1 .- repeat(x_e[:], inner=2) ) .* p.hr' / sum(p_vac.N0) * 1e5
	end
	
end


function calc_final_mortality(x, p; vac_mode="infection")

	# Overall vaccine efficacy is e
	# Relative efficacy for 60+ age groups is re
	x_e = p.e * [x[1:6]' p.re*x[7:8]']

	if vac_mode == "infection"
		S_vac = (1 .- repeat(x_e[:], inner=2) ) .* p.S
		p_vac = merge(p, (;S=S_vac))
	elseif vac_mode == "disease"
		y_vac = (1 .- repeat(x_e[:], inner=2) ) .* p.y
		p_vac = merge(p, (;y=y_vac))
	end

	NGM_FS = calc_NGM_FS(p_vac)
	
	final_size = calc_final_size(NGM_FS)
	
	if vac_mode == "infection"
		return p_vac.S .* (1 .- final_size) .* p.mr' / sum(p_vac.N0) * 1e5
	elseif vac_mode == "disease"
		return p_vac.S .* (1 .- final_size) .* (1 .- repeat(x_e[:], inner=2) ) .* p.mr' / sum(p_vac.N0) * 1e5
	end
	
end
