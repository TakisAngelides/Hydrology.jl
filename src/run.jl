"""
$(TYPEDSIGNATURES)

A function to the run a steady-state simulation.
"""
function run!(sim::SteadyStateSimulation)

    update_steady_state!(sim.model, sim.grid, sim.state)

end


#################################
# Model: Kazmierczak et al 2024 #
#################################


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


####################################
# Model: Heigh above buyancy (HAB) #
####################################

"""
$(TYPEDSIGNATURES)

Update the state variable of effective pressure according to Eq. (3) of the paper Kazmierczak et al 2022 (https://doi.org/10.5194/tc-16-4537-2022). 
The model essentially assumes that ocean water inflitrates the ice sheet from the grounding line upwards to grounded ice regions where the bedrock
is below sea level.
"""
function update_steady_state!(model::HABHydroModel, grid::OGRectHydroGrid, state::HydroState)

    # Update effective pressure N [Pa]
    update_N!(model, grid, state)

end
