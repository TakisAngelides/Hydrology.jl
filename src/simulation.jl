abstract type AbstractSimulation end

struct TimeSimulation <: AbstractSimulation
    # TODO: implement for models that are not steady state calculations
end

struct SteadyStateSimulation{M <: AbstractHydroModel, G <: AbstractHydroGrid, S <: AbstractHydroState} <: AbstractSimulation

    model::M
    grid::G
    state::S

end

function SteadyStateSimulation(model::M, grid::G, state::S) where 
    {M <: KazmierczakHydroModel, G  <: AbstractHydroGrid, S  <: AbstractHydroState}

    return SteadyStateSimulation{M, G, S}(model, grid, state)

end
