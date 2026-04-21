"""
$(TYPEDSIGNATURES)

Abstract type for a hydrology state to hold fields common to all hydrology models.
"""
abstract type AbstractHydroState end


"""
$(TYPEDSIGNATURES)

The struct for the HydroState which includes physical fields important to hydrology.
"""
mutable struct HydroState{A} <: AbstractHydroState
    # Inputs
    mask ::A  # grounded ice mask (1 = grounded)
    h    ::A  # ice thickness [m]
    b    ::A  # bedrock elevation [m]

    # Outputs
    N    ::A  # effective pressure [Pa]
    W    ::A  # water layer thickness or storage [m]
end


"""
$(TYPEDSIGNATURES)

The constructor to the struct for the HydroState. The user must provide inputs for the ice thickness h, the bedrock elevation b,
the melt rate per unit area in the x-y plane [Kg m⁻² s⁻¹], and the mask which should take the value 1 for grounded ice. The remaining 
fields of the state are initialized to zero, specifically the effective pressure N [Pa] and water layer thickness W [m].

# Arguments

- `grid::OGRectHydroGrid`: grid of the simulation
- `mask_in::AbstractArray{<:Real}`: mask input which should take the value 1 for grounded ice.
- `h_in::AbstractArray{<:AbstractFloat}`: ice thickness input
- `b_in::AbstractArray{<:AbstractFloat}`: bedrock elevation
"""
function HydroState(
    grid::OGRectHydroGrid,
    mask_in::AbstractArray{<:Real},
    h_in::AbstractArray{<:AbstractFloat},
    b_in::AbstractArray{<:AbstractFloat},
)

    # Check that the inputs match the grid size
    expected_size = (grid.grid.Nx, grid.grid.Ny)
    for (name, arr) in [("mask", mask_in), ("h", h_in), ("b", b_in)]
        size(arr)[1:2] == expected_size || throw(ArgumentError("$name size $(size(arr)) ≠ grid size $expected_size"))
    end
        
    # Inputs
    mask = set!(CenterField(grid.grid), mask_in) # grounded ice mask
    h    = set!(CenterField(grid.grid), h_in) # ice thickness
    b    = set!(CenterField(grid.grid), b_in) # bedrock elevation

    # Outputs
    N = set!(CenterField(grid.grid), 0.0) # effective pressure
    W = set!(CenterField(grid.grid), 0.0) # water layer thickness
    
    return HydroState(mask, h, b, N, W)

end
