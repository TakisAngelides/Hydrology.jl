
"""
$(TYPEDSIGNATURES)

Return a struct containing important physical constants.
Comes with default values that can however be changed by the user, for instance by running:

```julia
c = PhysicalConstants(ρ_i = 0.93)   # (kg/m^3)
```

All constants are given in SI units (kilogram, meter, second).
"""
@kwdef struct PhysicalConstants{T <: AbstractFloat}

    ρ_w::T = 1000.0 # density of water [kg/m^3]
    ρ_i::T = 917.0 # density of ice [kg/m^3]
    g::T = 9.81 # graviational acceleration [m/s^2]
    n::T = 3.0 # exponent in Glen's flow law
    L_w::T = 3.34e5 # latent heat of fusion of water [J/kg]
    h_b::T = 0.1 # thickness of bed obstacles [m]
    α::T = 5/4 # exponent in Darcy-Weisbach relation
    β::T = 3/2 # exponent in Darcy-Weisbach relation
    f::T = 0.1 # friction coefficient for a turbulent flow
    F_till::T = 1.1 # factor for the conduit geometry in a till
    Q_c::T = 1.0 # critical water flux in a conduit [m^3/s]
    H_0::T = 0.1 # thickness of canals [m]
    l_c::T = 10000.0 # distance between conduits [m]
    sqrtπ::T = √π
    K::T = (2/π)^(0.25) * sqrt((π + 2) / (ρ_w * f)) # conductivity coefficient in Darcy-Weisbach relation [kg^{1−β} m^{2β−2α+1} s^{2β−3}]

end