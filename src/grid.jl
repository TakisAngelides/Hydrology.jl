"""
$(TYPEDSIGNATURES)

Abstract type for the grid of a simulation.
"""
abstract type AbstractHydroGrid end


"""
$(TYPEDSIGNATURES)

A struct for the Oceananigans Rectilinear grid.
"""
struct OGRectHydroGrid <: AbstractHydroGrid
    grid::Oceananigans.RectilinearGrid
end


"""
$(TYPEDSIGNATURES)

The constructor to the struct for the Oceananigans RectilinearGrid.

# Arguments

- `Nx::I`: number of grid cells in the x direction.
- `Ny::I`: 
- `xlims`: tuple specifying the x values of the left-most and right-most edges of the grid in the x direction (e.g. xlims = (0, 1)).
- `ylims`: 

# Keywords

- `T`: type for the physical fields to live on the grid.
    (**Default**: `Float64`)
- `topology`: This specifies the boundary conditions for each of the x, y, z dimensions of the rectilinear grid.
    (**Default**: `(Bounded, Bounded, Flat)`)
- `halo`: tuple specifying the number of halo points in the x, y, z dimensions (e.g. when the z is flat i.e. the fields are not changing in that dimension, then halo = (1, 1) gives one ghost point for x and y to handle boundary conditions)
"""
function OGRectHydroGrid(Nx::I, Ny::I, xlims, ylims; T = Float64, topology = (Bounded, Bounded, Flat), halo = (1, 1)) where {I <: Integer}

    Nx > 0 || throw(ArgumentError("Nx must be positive"))
    Ny > 0 || throw(ArgumentError("Ny must be positive"))

    grid = Oceananigans.RectilinearGrid(T; size = (Nx, Ny), x = xlims, y = ylims, topology = topology, halo = halo)

    return OGRectHydroGrid(grid)
end


"""
$(TYPEDSIGNATURES)

A function to fill the ghost points with the appropriate values for any grid. If this function
is not defined for the specified grid, we inform the user with the error printed below.
"""
function fill_halo!(field, grid::AbstractHydroGrid)
    error("fill_halo! not implemented for $(typeof(grid))")
end


"""
$(TYPEDSIGNATURES)

A function to fill the ghost points of the field input with the appropriate values, specifically for an Oceananigans RectilinearGrid.
"""
function fill_halo!(field, ::OGRectHydroGrid) 
    fill_halo_regions!(field)
end
