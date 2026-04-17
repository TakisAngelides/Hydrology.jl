function run!(sim::SteadyStateSimulation)

    update_steady_state!(sim.model, sim.grid, sim.state)

end

function update_steady_state!(model::KazmierczakHydroModel, grid::OGRectHydroGrid, state::HydroState)

    update_q!(model, grid, state)
    update_N!(model, grid, state)

end

