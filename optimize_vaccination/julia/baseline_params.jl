# Pre-calibrated susceptibility and clinical fraction

# Consensus fitting
#u = repeat([0.39, 0.38, 0.79, 0.87, 0.80, 0.82, 0.89, 0.74], inner=2)
#y = repeat([0.28, 0.20, 0.26, 0.33, 0.40, 0.49, 0.63, 0.69], inner=2)

# Parametric fitting
u = [0.017781719 0.017781719 0.017781719 0.01826781 0.022270941 0.02828467 0.036056645 0.042069826 0.046073505 0.049071611 0.053326546 0.060093496 0.067362959 0.070109808 0.070595351 0.070595351]
y = fill(0.5, 16)
#u = fill(1, 16)
#y = [0.143576826 0.143576826 0.143576826 0.14861461 0.17884131 0.231738035 0.294710327 0.34256927 0.375314861 0.400503778 0.435768262 0.496221662 0.549118388 0.571788413 0.576826196 0.579345088]

# Empty contact matrix
c = zeros(16,16)

# Relative infectiousness of sub-clinical cases
f = 0.5

# Sojourn times in each infected state
te = 3.0
tp = 2.1
tc = 2.9
ts = 5.0

# Susceptible population
S = ones(16)
N0 = S

# Scaling factor for NGM (i.e. probability of transmission given contact)
sf = 1

# Baseline vaccine efficacy
e = 0.7

# Relative efficacy in 60+ age groups
re = 1.0
#re = 0.5

p_base = (;u,y,c,f,te,tp,tc,ts,S,N0,sf,e,re)

