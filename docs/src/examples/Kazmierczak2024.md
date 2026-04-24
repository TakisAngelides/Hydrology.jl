```@meta
EditURL = "Kazmierczak2024.jl"
```

# [Kazmierczak et al 2024](@id Kazmierczak2024)
This is an example of how to run the FastHydrology.jl package for the steady state problem of Kazmierczak et al 2024.

````@example Kazmierczak2024
using FastHydrology
using Statistics
using FileIO
using CairoMakie
````

How to load data from a file to get the necessary inputs to initialize the grid, model and state for the Kazmierczak et al 2024 paper.

````@example Kazmierczak2024
T = Float64 # Type of the physical fields
path = "$(@__DIR__)/input/Kazmierczak2024/THWAITES2km_m3_HAB_toto.mat" # Path to data for Fig. S3 of Kazmierczak et al 2024 which focuses on Thwaites glacier.
Nx, Ny, xlims, ylims, mask, h, b, abs_v_b, A_visc, ṁ_over_ρ_w, κ = load_Kazmierczak(path; bed_rheology = :hard)
````

Here we prepare a grid using the Oceananigans rectilinear grid.

````@example Kazmierczak2024
grid = OGRectHydroGrid(Nx, Ny, xlims, ylims; T = T)
````

There is a function to visualize the Oceananigans rectilinear grid.

````@example Kazmierczak2024
fig = visualize_grid(grid)
````

We then build the model using the data from the input file above. The model will hold its model-specific fields, with the rest of the fields common to all models stored in the HydroState below.

````@example Kazmierczak2024
model = KazmierczakHydroModel(grid, κ, abs_v_b, A_visc, ṁ_over_ρ_w);
nothing #hide
````

Further, we build the hydrology state using some of the data from the input file. This state will store the fields common to all hydrology models.

````@example Kazmierczak2024
state = HydroState(grid, mask, h, b);
nothing #hide
````

Finally we can create a simulation struct. Since the Kazmierczak et al 2024 model is a steady state problem we create a SteadyStateSimulation rather than a TimeSimulation.

````@example Kazmierczak2024
sim = SteadyStateSimulation(model, grid, state);
nothing #hide
````

We can now run the simulation which will update the hydrology state, as well as the model-dependent fields defined in the model.

````@example Kazmierczak2024
run!(sim)
````

Applying some rescaling for plotting.

````@example Kazmierczak2024
model.q .*= perSecond2perYear(1.0) / 1e4 # makes q [10⁴ m² a⁻¹]
state.N .*= 1e-6 # makes N [MPa]
model.q .= mask_field(model.q, state.mask, NaN)
state.N .= mask_field(state.N, state.mask, NaN)
state.W .= mask_field(state.W, state.mask, NaN)
display_flag = false
````

Now we can visualize the resulting water flux q [10⁴ m² a⁻¹].

````@example Kazmierczak2024
fig_q = visualize_field(model.q; plot_title = "Distributed water flux q [10⁴ m² a⁻¹]", transpose_data = true, display_flag = display_flag, colorrange = (0, 10))
````

The effective pressure N [MPa].

````@example Kazmierczak2024
fig_N = visualize_field(state.N; plot_title = "Effective pressure N [MPa]", transpose_data = true, display_flag = display_flag, colorrange = (0, 10))
````

The water layer thickness W [m].

````@example Kazmierczak2024
fig_W = visualize_field(state.W; plot_title = "Water thickness W [m]", transpose_data = true, display_flag = display_flag, colorrange = extrema(filter(!isnan, state.W.data)))
````

Using a similar setup to the one above but with model.longcoupwater = 0.0, we can plot results for the whole of Antarctica.

````@example Kazmierczak2024
function plot(f, r)
    img = load("$(@__DIR__)/figures/Kaz24_yelmox_$(f)_$(r)km.png")
    fig = CairoMakie.Figure()
    ax = Axis(fig[1, 1], yreversed = true)
    image!(ax, img')
    ax.aspect = DataAspect()
    hidedecorations!(ax)
    return fig
end
````

Water flux [m²/s] 32km resolution.

````@example Kazmierczak2024
fig = plot("q", "32")
````

Water flux [m²/s] 16km resolution.

````@example Kazmierczak2024
fig = plot("q", "16")
````

Water thickness [m] 32km resolution.

````@example Kazmierczak2024
fig = plot("W", "32")
````

Water thickness [m] 16km resolution.

````@example Kazmierczak2024
fig = plot("W", "16")
````

Effective pressure [Pa] 32km resolution.

````@example Kazmierczak2024
fig = plot("N", "32")
````

Effective pressure [Pa] 16km resolution.

````@example Kazmierczak2024
fig = plot("N", "16")
````

