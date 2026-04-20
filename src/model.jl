abstract type AbstractHydroModel end

mutable struct KazmierczakHydroModel{T <: AbstractFloat, A} <: AbstractHydroModel

    # TODO: add meaning to variables with comments

    # Model constants
    ρ_w    ::T
    ρ_i    ::T
    g      ::T
    L_w    ::T
    n      ::T
    h_b    ::T
    α      ::T
    β      ::T
    f      ::T
    F_till ::T
    Q_c    ::T
    H_0    ::T
    l_c    ::T
    K      ::T
    η_w    ::T
    Wmin  ::T
    Wmax  ::T

    # Geometric potential
    ϕ₀                   ::A
    ϕ₀_tmp               ::A
    minus_∇ϕ₀_x          ::A
    minus_∇ϕ₀_y          ::A
    abs_∇ϕ₀              ::A
    minus_∇ϕ₀_smoothed_x ::A
    minus_∇ϕ₀_smoothed_y ::A
    abs_∇ϕ₀_smoothed     ::A

    # Water flux
    ψ_out  ::A
    corfac ::A
    q      ::A

    # Effective pressure
    Q       ::A
    κ       ::A
    abs_v_b ::A
    A_visc  ::A
    S_inf   ::A
    H_hard  ::A
    H_soft  ::A
    H       ::A
    N_inf   ::A
    Po      ::A

end


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
    η_w = perYear2perSecond(T(1.8e-3)) # viscosity of water
    Wmin = T(1e-8)  # minimum value for water layer thickness W
    Wmax = T(0.015) # maximum value for water layer thickness W
    
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
