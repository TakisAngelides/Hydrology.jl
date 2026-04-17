abstract type AbstractHydroState end

mutable struct HydroState{A} <: AbstractHydroState

    # Inputs
    mask ::A  # grounded ice mask (1 = grounded)
    h    ::A  # ice thickness [m]
    b    ::A  # bedrock elevation [m]
    ṁ    ::A  # basal melt rate [kg/m²/s]

    # Outputs
    N    ::A  # effective pressure [Pa]
    W    ::A  # water layer thickness or storage [m]

end


function HydroState(
    grid::OGRectHydroGrid,
    mask_in::AbstractArray{<:Real},
    h_in::AbstractArray{<:AbstractFloat},
    b_in::AbstractArray{<:AbstractFloat},
    ṁ_in::AbstractArray{<:AbstractFloat},
)

    expected_size = (grid.grid.Nx, grid.grid.Ny)
    
    for (name, arr) in [("mask", mask_in), ("h", h_in), ("b", b_in), ("ṁ", ṁ_in)]
        size(arr)[1:2] == expected_size || throw(ArgumentError("$name size $(size(arr)) ≠ grid size $expected_size"))
    end
        
    # Inputs
    mask = set!(CenterField(grid.grid), mask_in) # grounded ice mask
    h    = set!(CenterField(grid.grid), h_in) # ice thickness
    b    = set!(CenterField(grid.grid), b_in) # bedrock elevation
    ṁ    = set!(CenterField(grid.grid), ṁ_in) # basal melt rate

    # Outputs
    N = set!(CenterField(grid.grid), 0.0) # effective pressure
    W = set!(CenterField(grid.grid), 0.0) # water layer thickness
    
    return HydroState(mask, h, b, ṁ, N, W)

end
