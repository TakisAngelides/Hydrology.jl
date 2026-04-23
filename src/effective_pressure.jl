#################################
# Model: Kazmierczak et al 2024 #
#################################


"""
$(TYPEDSIGNATURES)

Update the effective pressure N across the grid using a complementary error function 
transition between geometric potential and far-field effective pressure.
"""
function update_N!(model::KazmierczakHydroModel, grid::OGRectHydroGrid, state::HydroState)
    update_Q!(model, grid)     # volumetric water flux [m³ s⁻¹] per conduit in a given grid cell
    update_S_inf!(model, grid) # cross-sectional area of conduits
    update_H!(model, grid)     # thickness of conduits 
    update_Po!(model, grid, state)    # ice overburden pressure ρgh
    update_N_inf!(model, grid) # far-field effective pressure valid when we are more upstream than a few kilometers from the grounding line
    @. state.N.data = max(0.0, erf((sqrt(π) * model.ϕ₀.data) / (2 * model.N_inf.data)) * model.N_inf.data) # effective pressure
    fill_halo!(state.N, grid)
    return nothing
end


"""
$(TYPEDSIGNATURES)

Update the hydrostatic ice overburden pressure Po based on ice thickness h.
"""
function update_Po!(model::Union{KazmierczakHydroModel, HABHydroModel}, grid::OGRectHydroGrid, state::HydroState)
    @. model.Po.data = max(model.ρ_i * model.g * state.h.data, 1e5) # minimum limit taken from KORI-ULB model - see KoriInputParams.m at https://github.com/FrankPat/Kori-ULB
    fill_halo!(model.Po, grid)
    return nothing
end


"""
$(TYPEDSIGNATURES)

Update the conduit thickness H by calculating separate values for hard and soft beds, 
then interpolating based on the bed heterogeneity indicator κ.
"""
function update_H!(model::KazmierczakHydroModel, grid::OGRectHydroGrid)
    @. model.H_hard = sqrt(model.S_inf)
    @. model.H_soft = model.H_0 + (sqrt(model.S_inf) / model.F_till - model.H_0) * exp(-model.Q / model.Q_c)
    @. model.H = (1 - model.κ) * model.H_hard + model.κ * model.H_soft
    fill_halo!(model.H, grid)
    return nothing
end


"""
$(TYPEDSIGNATURES)

Update the far-field conduit cross-sectional area S_inf using the Manning or 
Gauckler-Manning-Strickler flow law.
"""
function update_S_inf!(model::KazmierczakHydroModel, grid::OGRectHydroGrid)
    @. model.S_inf = model.K^(-1 / model.α) * model.abs_∇ϕ₀^((1 - model.β) / model.α) * model.Q^(1 / model.α)
    fill_halo!(model.S_inf, grid)
    return nothing
end


"""
$(TYPEDSIGNATURES)

Update the far-field effective pressure N_inf based on conduit geometry and 
basal velocity, constrained by ice overburden pressure limits.
"""
function update_N_inf!(model::KazmierczakHydroModel, grid::OGRectHydroGrid)
    # Lower limit of model.sigmat * Po ≤ N_inf ≤ Po comes from Eq. (20) of Beuler and van Pelt 2015, where e.g. sigmat = 0.02
    @. model.N_inf.data = min(max(
        ((model.H.data/model.S_inf.data)^2*((model.ρ_i*model.L_w*model.abs_v_b.data*model.h_b + model.Q.data*model.abs_∇ϕ₀.data)/(2.0*model.n^(-model.n)*model.ρ_i*model.L_w*model.A_visc.data)))^(1.0/model.n),
        model.sigmat * model.Po.data),
        model.Po.data)
    @. model.N_inf.data[model.S_inf.data .== 0.0] .= model.Po[model.S_inf.data .== 0.0] # if the cross-sectional area of a conduit is zero then the effective pressure is equal to the ice overburden pressure
    fill_halo!(model.N_inf, grid)
    return nothing
end


"""
$(TYPEDSIGNATURES)

Update the volumetric water flux per conduit Q by scaling the distributed flux q 
by the characteristic channel spacing l_c.
"""
function update_Q!(model::KazmierczakHydroModel, grid::OGRectHydroGrid)
    @. model.Q = model.q * model.l_c
    fill_halo!(model.Q, grid)
    return nothing
end


####################################
# Model: Heigh above buyancy (HAB) #
####################################


"""
$(TYPEDSIGNATURES)

Update the effective pressure N across the grid using a complementary error function 
transition between geometric potential and far-field effective pressure.
"""
function update_N!(model::HABHydroModel, grid::OGRectHydroGrid, state::HydroState)

    # Update ice overburden pressure
    update_Po!(model, grid, state)

    # Update water pressure
    update_p_w!(model, grid, state)

    # Update effective pressure
    @. state.N.data = max(model.Po.data - model.p_w.data, model.sigmat * model.Po.data)
    fill_halo!(state.N, grid)
    
    return nothing
end


"""
$(TYPEDSIGNATURES)

Update the water pressure p_w.
"""
function update_p_w!(model::HABHydroModel, grid::OGRectHydroGrid, state::HydroState)

    @. model.p_w.data                          = - model.P_w * model.ρ_sw * model.g * state.b.data 
    @. model.p_w.data[state.mask.data .== 0.0] =   model.P_w * model.ρ_i  * model.g * state.h.data[state.mask.data .== 0.0] # this allocates memory and we can avoid that with a for loop but the current method is faster
    @. model.p_w.data[state.b.data .>= 0.0]    = 0.0
    fill_halo!(model.p_w, grid)

    return nothing

end
