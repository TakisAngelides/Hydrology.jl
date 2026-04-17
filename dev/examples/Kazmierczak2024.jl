using Hydrology

# Type of the physical fields
T = Float64 

# Path to data; run the file from the repository Hydrology.jl path
path = "$(@__DIR__)/input/Kazmierczak2024/THWAITES2km_m3_HAB_toto.mat"

# Function to load the data
data_loading_function = load_Kazmierczak

# Load the data as they are saved from the repository of Kazmierczak et al 2024 for Fig. S3
Nx, Ny, xlims, ylims, mask, h, b, abs_v_b, A_visc, ṁ_over_ρ_w, κ = load_Kazmierczak(path)

# Build the grid as an Oceananigans rectilinear grid
grid = OGRectHydroGrid(Nx, Ny, xlims, ylims; T = T)

# Build the model of Kazmierczak et al 2024
model = KazmierczakHydroModel(grid, κ, abs_v_b, A_visc)

# Build the hydrology state that is common to all hydrology models
state = HydroState(grid, mask, h, b, ṁ_over_ρ_w*model.ρ_w)

# Create a Simulation struct for Kazmierczak 2024 et al which is a steady state calculation
sim = SteadyStateSimulation(model, grid, state)

# Run the simulation
run!(sim)

visualize_field(model.q, state.mask; plot_title = "q")
visualize_field(state.N, state.mask; plot_title = "N")

println("Finished.")
