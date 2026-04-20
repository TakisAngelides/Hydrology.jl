"""
$(TYPEDSIGNATURES)

An abstract type for the hydrology model to be simulated. The model can hold revelant constants and model-specific fields. 
"""
abstract type AbstractHydroModel end


"""
$(TYPEDSIGNATURES)

The hydrology model described in Kazmierczak et al 2024 (https://doi.org/10.5194/tc-18-5887-2024). In our implementation here for
the calculations of water flux, we make the assumption that on the grid we have Δx = Δy.
"""
mutable struct KazmierczakHydroModel{T <: AbstractFloat, A} <: AbstractHydroModel

    # Model constants 
    ρ_w    ::T  # Density of water [kg/m³]
    ρ_i    ::T  # Density of ice [kg/m³]
    g      ::T  # Gravitational acceleration [m/s²]
    L_w    ::T  # Latent heat of fusion for ice [J/kg]
    n      ::T  # Glen's flow law exponent (typically 3)
    h_b    ::T  # Typical bed obstacle height [m]
    α      ::T  # Power law exponent for hydraulic transmissivity (m-scale)
    β      ::T  # Power law exponent for hydraulic transmissivity (opening/closing)
    f      ::T  # Darcy-Weisbach friction factor
    F_till ::T  # Till compressibility/yield factor for soft-bed transition
    Q_c    ::T  # Threshold discharge for laminar-to-turbulent transition [m³/s]
    H_0    ::T  # Thickness of canals for soft bed deformation [m]
    l_c    ::T  # Distance between conduits [m]
    K      ::T  # Conductivity coefficient in Darcy–Weisbach relation
    η_w    ::T  # Dynamic viscosity of water [Pa·s]
    Wmin   ::T  # Minimum subglacial water layer thickness [m]
    Wmax   ::T  # Maximum subglacial water layer thickness [m]

    # Geometric potential 
    ϕ₀                   ::A  # Geometric potential [Pa]
    ϕ₀_tmp               ::A  # Temporary storage for potential filling of ϕ₀ to smoothen local minima and avoid stuck water
    minus_∇ϕ₀_x          ::A  # Geometric potential gradient x-component [Pa/m]
    minus_∇ϕ₀_y          ::A  # Geometric potential gradient y-component [Pa/m]
    abs_∇ϕ₀              ::A  # Magnitude of the geometric potential gradient [Pa/m]
    minus_∇ϕ₀_smoothed_x ::A  # Smoothed gradient x-component of the geometric potential [Pa/m]
    minus_∇ϕ₀_smoothed_y ::A  # Smoothed gradient y-component of the geometric potential [Pa/m]
    abs_∇ϕ₀_smoothed     ::A  # Magnitude of the smoothed gradient of the geometric potential [Pa/m]

    # Water flux 
    ψ_out  ::A  # Integrated scalar water flux [m³/s]
    corfac ::A  # Correction factor to go from ψ_out to q
    q      ::A  # Distributed water flux [m²/s]

    # Effective pressure and Bed state 
    Q       ::A  # Volumetric water flux within a conduit [m³/s]
    κ       ::A  # Bed type indicator (0: hard, 1: soft)
    abs_v_b ::A  # Magnitude of basal sliding velocity [m/s]
    A_visc  ::A  # Ice flow law rate factor (Glen's A) [Pa⁻ⁿ s⁻¹]
    S_inf   ::A  # Far-field (away from grounding line) conduit cross-sectional area [m²]
    H_hard  ::A  # Thickness of conduits over a hard bed [m]
    H_soft  ::A  # Thickness of conduits over a soft bed [m]
    H       ::A  # Thickness of conduits [m]
    N_inf   ::A  # Far-field (away from grounding line) effective pressure [Pa]
    Po      ::A  # Ice overburden pressure (ρ_i * g * ice_thickness) [Pa]

end


"""
$(TYPEDSIGNATURES)

The constructor to the Kazmierczak et al 2024 hydrology model. All fields are initialized here to zero, except from 
the viscosity parameter A_visc from Glen's flow, the κ field describing the hardness of the bed, and the absolute value
of the basal velocity of the ice. These three fields are given values from an input file. The user must provide these fields
for the simulation to be able to start.
"""
function KazmierczakHydroModel(
    grid::OGRectHydroGrid,
    κ_in::AbstractArray{<:AbstractFloat},
    abs_v_b_in::AbstractArray{<:AbstractFloat},
    A_visc_in::AbstractArray{<:AbstractFloat},
)

    expected_size = (grid.grid.Nx, grid.grid.Ny)
    
    for (name, arr) in [("κ", κ_in), ("abs_v_b", abs_v_b_in), ("A_visc", A_visc_in)]
        size(arr)[1:2] == expected_size || throw(ArgumentError("$name size $(size(arr)) ≠ grid size $expected_size"))
    end
    
    T  = eltype(grid.grid)
    
    # Physical constants
    ρ_w    = T(1000.0)
    ρ_i    = T(917.0)
    g      = T(9.81)
    L_w    = T(3.34e5)
    n      = T(3.0)
    h_b    = T(0.1)
    α      = T(5//4)
    β      = T(3//2)
    f      = T(0.1)
    F_till = T(1.1)
    Q_c    = T(1.0)
    H_0    = T(0.1)
    l_c    = T(10000.0)
    K      = (T(2)/T(pi))^(T(0.25)) * sqrt((T(pi) + T(2)) / (ρ_w * f))
    η_w = perYear2perSecond(T(1.8e-3)) # viscosity of water (values from KORI-ULB model - see the file KoriInputParams.m at https://github.com/FrankPat/Kori-ULB)
    Wmin = T(1e-8)  # minimum value for water layer thickness W (values from KORI-ULB model)
    Wmax = T(0.015) # maximum value for water layer thickness W (values from KORI-ULB model)
    
    # Geometric potential
    ϕ₀                   = set!(CenterField(grid.grid), 0.0)  # Geometric potential [Pa]
    ϕ₀_tmp               = set!(CenterField(grid.grid), 0.0)  # Geometric potential placeholder for filling
    minus_∇ϕ₀_x          = set!(CenterField(grid.grid), 0.0)  # -∂ϕ₀/∂x [Pa/m]
    minus_∇ϕ₀_y          = set!(CenterField(grid.grid), 0.0)  # -∂ϕ₀/∂y [Pa/m]
    abs_∇ϕ₀              = set!(CenterField(grid.grid), 0.0)  # |∇ϕ₀| [Pa/m]
    minus_∇ϕ₀_smoothed_x = set!(CenterField(grid.grid), 0.0)  # smoothed -∂ϕ₀/∂x [Pa/m]
    minus_∇ϕ₀_smoothed_y = set!(CenterField(grid.grid), 0.0)  # smoothed -∂ϕ₀/∂y [Pa/m]
    abs_∇ϕ₀_smoothed     = set!(CenterField(grid.grid), 0.0)  # |smoothed ∇ϕ₀| [Pa/m]
    
    # Water flux
    ψ_out  = set!(CenterField(grid.grid), 0.0)  # Integrated scalar water flux [m^2/s]
    corfac = set!(CenterField(grid.grid), 0.0)  # Correction factor from ψ_out to q
    q      = set!(CenterField(grid.grid), 0.0)  # Distributed water flux [m^2/s]
    
    # Effective pressure
    Q       = set!(CenterField(grid.grid), 0.0)          # Volumetric water flux per conduit [m^3/s]
    κ       = set!(CenterField(grid.grid), κ_in)         # Bed heterogeneity indicator [0=hard, 1=soft]
    abs_v_b = set!(CenterField(grid.grid), abs_v_b_in)  # Absolute basal velocity [m/s]
    A_visc  = set!(CenterField(grid.grid), A_visc_in)   # Glen's flow law viscosity parameter
    S_inf   = set!(CenterField(grid.grid), 0.0)          # Far-field conduit cross-sectional area [m^2]
    H_hard  = set!(CenterField(grid.grid), 0.0)          # Conduit thickness over hard bed [m]
    H_soft  = set!(CenterField(grid.grid), 0.0)          # Conduit thickness over soft bed [m]
    H       = set!(CenterField(grid.grid), 0.0)          # Conduit thickness [m]
    N_inf   = set!(CenterField(grid.grid), 0.0)          # Far-field effective pressure [Pa]
    Po      = set!(CenterField(grid.grid), 0.0)          # Ice overburden pressure
    
    return KazmierczakHydroModel(
        ρ_w, ρ_i, g, L_w, n, h_b, α, β, f, F_till, Q_c, H_0, l_c, K, η_w, Wmin, Wmax,
        ϕ₀, ϕ₀_tmp, minus_∇ϕ₀_x, minus_∇ϕ₀_y,
        abs_∇ϕ₀, minus_∇ϕ₀_smoothed_x, minus_∇ϕ₀_smoothed_y, abs_∇ϕ₀_smoothed,
        ψ_out, corfac, q,
        Q, κ, abs_v_b, A_visc, S_inf, H_hard, H_soft, H, N_inf, Po
    )

end
