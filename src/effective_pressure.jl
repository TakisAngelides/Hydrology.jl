"""
$(TYPEDSIGNATURES)

Update the effective pressure `N` across the grid using a complementary error function 
transition between geometric potential and far-field effective pressure.
"""
function update_N!(HM::HydrologyModel)

    fields, c = HM.fields, HM.c
    
    update_Q!(fields, c) # volumetric water flux [m³ s⁻¹] per conduit in a given grid cell
    update_S_inf!(fields, c) # cross-sectional area of conduits
    update_H!(fields, c) # thickness of conduits 
    update_Po!(fields, c) # ice overburden pressure ρgh
    update_N_inf!(fields, c) # far-field effective pressure valid when we are more upstream than a few kilometers from the grounding line

    @. fields.N.data = (1-erfc((c.sqrtπ*fields.ϕ₀.data)/(2*fields.N_inf.data))) * fields.N_inf.data # effective pressure

    return nothing

end


"""
$(TYPEDSIGNATURES)

Update the hydrostatic ice overburden pressure `Po` based on ice thickness `h`.
"""
function update_Po!(fields, c)
    fields.Po .= c.ρ_i * c.g * fields.h
    return nothing
end


"""
$(TYPEDSIGNATURES)

Update the conduit thickness `H` by calculating separate values for hard and soft beds, 
then interpolating based on the bed heterogeneity indicator `κ`.
"""
function update_H!(fields, c)
    @. fields.H_hard = sqrt(fields.S_inf)
    @. fields.H_soft = c.H_0 + (sqrt(fields.S_inf)/c.F_till - c.H_0) * exp(-fields.Q/c.Q_c)
    @. fields.H = (1 - fields.κ) * fields.H_hard + fields.κ * fields.H_soft
    return nothing
end


"""
$(TYPEDSIGNATURES)

Update the far-field conduit cross-sectional area `S_inf` using the Manning or 
Gauckler-Manning-Strickler flow law.
"""
function update_S_inf!(fields, c)
    @. fields.S_inf = c.K^(-1/c.α) * fields.abs_∇ϕ₀^((1-c.β)/c.α) * fields.Q^(1/c.α) 
    return nothing
end


"""
$(TYPEDSIGNATURES)

Update the far-field effective pressure `N_inf` based on conduit geometry and 
basal velocity, constrained by ice overburden pressure limits.
"""
function update_N_inf!(fields, c)
    # Lower limit of 0.02*Po ≤ N_inf ≤ Po comes from Eq. (20) of Beuler and van Pelt 2015 
    @. fields.N_inf.data = min(max(
        ((fields.H.data/fields.S_inf.data)^2*((c.ρ_i*c.L_w*fields.abs_v_b.data*c.h_b + fields.Q.data*fields.abs_∇ϕ₀.data)/(2.0*c.n^(-c.n)*c.ρ_i*c.L_w*fields.A.data)))^(1.0/c.n),
        0.02*fields.Po.data),
        fields.Po.data)
    return nothing
end


"""
$(TYPEDSIGNATURES)

Update the volumetric water flux per conduit `Q` by scaling the distributed flux `q` 
by the characteristic channel spacing `l_c`.
"""
function update_Q!(fields, c)
    @. fields.Q = fields.q * c.l_c
    return nothing
end