using LinearAlgebra
using UnPack

calc_NGM(p) = calc_NGM!(zeros(length(p.S),length(p.S)), p)

function calc_NGM!(NGM, p)

	@unpack u,c,f,y,te,tc,tp,ts,N0,S,sf = p

	@. begin
		NGM = u * S * c / N0' * (y * (tp + tc) + (1 - y) * f * ts) * sf
	end
	
end


function calc_R0(p; NGM=calc_NGM(p))

	evals = eigvals(NGM)
	
	R0 = maximum(real(evals[isreal.(evals)]))

end


function calc_Reff(x, p; vac_mode="infection")
	
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
		Reff = calc_R0(p_vac)
	
end


function calc_scale_factor(p, R0)

	R0_unscaled = calc_R0(p)

	return R0 / R0_unscaled

end


function calc_HIT(p; R0=calc_R0(p))

	#return (R0 - 1) / (e * R0)
	
	Reff = calc_Reff(fill(1, 8), p)
	
	return (R0 - 1) / (R0 - Reff)

end