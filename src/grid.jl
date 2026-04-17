abstract type AbstractHydroGrid end


struct OGRectHydroGrid <: AbstractHydroGrid
    grid::Oceananigans.RectilinearGrid
end


function OGRectHydroGrid(Nx::I, Ny::I, xlims, ylims; T = Float64, topology = (Bounded, Bounded, Flat), halo = (1, 1)) where {I <: Integer}

    Nx > 0 || throw(ArgumentError("Nx must be positive"))
    Ny > 0 || throw(ArgumentError("Ny must be positive"))

    grid = Oceananigans.RectilinearGrid(T; size = (Nx, Ny), x = xlims, y = ylims, topology = topology, halo = halo)

    return OGRectHydroGrid(grid)
end


function fill_halo!(field, grid::AbstractHydroGrid)
    error("fill_halo! not implemented for $(typeof(grid))")
end


function fill_halo!(field, ::OGRectHydroGrid) 
    fill_halo_regions!(field)
end
