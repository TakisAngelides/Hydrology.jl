"""
$(TYPEDSIGNATURES)

A function to the run a steady-state simulation.
"""
function run!(sim::SteadyStateSimulation)

    update_steady_state!(sim.model, sim.grid, sim.state)

end


"""
$(TYPEDSIGNATURES)

Update the model and state variables for the Kazmierczak et al 2024 model steady-state calculation. 
We specifically update first the distributed water flux q, then the water layer thickness W, and then the effective pressure N.
"""
function update_steady_state!(model::KazmierczakHydroModel, grid::OGRectHydroGrid, state::HydroState)

    # Update distributed water flux q [m²/s]
    update_q!(model, grid, state)
    
    # Update water layer thickness W [m] stored in the HydroState
    update_W!(model, grid, state)
    
    # Update effective pressure N [Pa]
    update_N!(model, grid, state)

end
