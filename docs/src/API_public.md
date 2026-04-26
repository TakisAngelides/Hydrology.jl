# Public API

## Grid
```@docs
AbstractHydroGrid
OGRectHydroGrid
```

## Model
```@docs
AbstractHydroModel
KazmierczakHydroModel
HABHydroModel
```

## State
```@docs
AbstractHydroState
HydroState
```

## Simulation
```@docs
AbstractSimulation
TimeSimulation
SteadyStateSimulation
```

## Running Simulations
```@docs
run!
update_steady_state!
```

## Water Flux
```@docs
update_q!
fill_halo!
update_ϕ₀!
potential_filling!
update_potential_gradients!
update_smoothed_potential_gradients!
accumulate_ψ_out!
update_ψ_out!
update_W!
```

## Effective Pressure
```@docs
update_N!
update_Po!
update_H!
update_S_inf!
update_N_inf!
update_Q!
```

## Data Loading
```@docs
load_Kazmierczak
load_yelmox
```

## Utilities
```@docs
compute_lims
perYear2perSecond
perSecond2perYear
Km2m
```

## Plotting
```@docs
visualize_grid
visualize_field
mask_field
```