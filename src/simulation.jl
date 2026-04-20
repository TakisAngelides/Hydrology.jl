"""
$(TYPEDSIGNATURES)

Abstract type for a simulation struct. The simulation should hold important information such as model fields and constants, the grid, and a state with physically relevant fields. 
The method that will perform the computation for the simulation should be implemented for each subtype of AbstractSimulation.
"""
abstract type AbstractSimulation end


"""
$(TYPEDSIGNATURES)

A struct for time evolution simulations. For further details see the corresponding constructor.
"""
struct TimeSimulation <: AbstractSimulation
end


"""
$(TYPEDSIGNATURES)

A struct for steady-state simulations. For further details see the corresponding constructor.
"""
struct SteadyStateSimulation{M <: AbstractHydroModel, G <: AbstractHydroGrid, S <: AbstractHydroState} <: AbstractSimulation
    model::M
    grid::G
    state::S
end


"""

$(TYPEDSIGNATURES)

Constructor to the SteadyStateSimulation struct that holds information for a steady-state simulation.

# Arguments

- `model::M`: This is a subtype of AbstractHydroModel that represents the specific model we want to use for the simulation.
- `grid::G`: The grid on which the simulation should run, which is a subtype of AbstractHydroGrid.
- `state::S`: The state is a subtype of AbstractHydroState and holds the fields common to all hydrology models.
"""
function SteadyStateSimulation(model::M, grid::G, state::S) where 
    {M <: KazmierczakHydroModel, G  <: AbstractHydroGrid, S  <: AbstractHydroState}

    return SteadyStateSimulation{M, G, S}(model, grid, state)
end
