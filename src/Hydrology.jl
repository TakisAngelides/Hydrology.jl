module Hydrology

using MAT
using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using CairoMakie
using BenchmarkTools
using ImageFiltering
using OffsetArrays
using Base.Threads
using SpecialFunctions
using DocStringExtensions

include("constants.jl")
include("utilities.jl")
include("data_loaders.jl")
include("initialize.jl")
include("water_flux.jl")
include("effective_pressure.jl")
include("plotting.jl")

# constants.jl
export PhysicalConstants

# utilities.jl
export compute_lims, perYear2perSecond, Km2m

# data_loaders.jl
export load_kazmierczak2024

# initialize.jl
export HydrologyModel, initialize_grid, initialize_fields, set_initial_fields!, initialize_κ!

# water_flux.jl
export update_q!, update_ϕ₀!, potential_filling!, update_potential_gradients!, update_smoothed_potential_gradients!, accumulate_ψ_out!, update_ψ_out!

# effective_pressure.jl
export update_N!, update_Po!, update_H!, update_S_inf!, update_N_inf!, update_Q!

# plotting.jl
export visualize_grid, visualize_field

end
