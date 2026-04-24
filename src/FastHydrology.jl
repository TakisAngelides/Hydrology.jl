module FastHydrology

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using BenchmarkTools
using ImageFiltering
using OffsetArrays
using Base.Threads
using SpecialFunctions
using DocStringExtensions
using MAT
using NCDatasets
using CairoMakie
using CairoMakie: Reverse

include("grid.jl")
include("model.jl")
include("state.jl")
include("simulation.jl")
include("run.jl")
include("water_flux.jl")
include("effective_pressure.jl")
include("data_loaders.jl")
include("utilities.jl")
include("plotting.jl")

# grid.jl
export AbstractHydroGrid, OGRectHydroGrid, fill_halo!

# model.jl
export AbstractHydroModel, KazmierczakHydroModel, HABHydroModel

# state.jl
export AbstractHydroState, HydroState

# simulations.jl
export AbstractSimulation, TimeSimulation, SteadyStateSimulation

# run.jl
export run!, update_steady_state!

# water_flux.jl
export update_q!, update_ϕ₀!, potential_filling!, update_potential_gradients!, update_smoothed_potential_gradients!, accumulate_ψ_out!, update_ψ_out!, update_W!

# effective_pressure.jl
export update_N!, update_Po!, update_H!, update_S_inf!, update_N_inf!, update_Q!

# data_loaders
export load_Kazmierczak, load_yelmox

# utilities.jl
export compute_lims, perYear2perSecond, perSecond2perYear, Km2m

# plotting.jl
export visualize_grid, visualize_field, mask_field

end
