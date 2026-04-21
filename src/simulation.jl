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
