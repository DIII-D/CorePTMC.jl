"""Hypothesis for temperatures inside the sheath:
    2. Ei = EiSE - e * ΔV;
    3. Ti = TiSE - 2/3 *e *ΔV; true only if distribution is maxwellian inside the sheath. 
    In reality the background ion distribution is skewed and not maxwellian.
    A full kinetic approach would be suitable to estimate thermal forces, this is a placeholder for the moment.
    4. Te = TeSE 
"""

α = 0.71 * Z .^ 2
β = -3 .* (1 - μ - 5 * sqrt(2) * Z ^ 2 * (1.1 * μ ^ 2.5 - 0.35 * μ ^ 1.5) ) ./ (2.6 - 2 * μ + 5.4 * μ ^ 2)

function grad_parallel(x, b_parallel)
    grad_x = grad(x)
    return scalar_product(grad_x, b_parallel)
end

F_thermal = α * grad_parallel(Tₑ) + β * grad_parallel(Tᵦ) # background plasma pushing particles from where it is colder (higher collision frequency) towards volume where it is hotter (lower coll.freq.). 

Δv_thermal = F_thermal / m * dt