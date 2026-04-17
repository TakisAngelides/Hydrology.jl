"""
$(TYPEDSIGNATURES)

Visualize the Oceananigans Rectilinear grid.
"""
function visualize_grid(grid::OGRectHydroGrid)

    grid = grid.grid

    xc = xnodes(grid, Center())
    yc = ynodes(grid, Center())

    dx = xspacings(grid, Center())
    dy = yspacings(grid, Center())

    Nx = length(xc)
    Ny = length(yc)

    quadrants = [
        (1:5, Ny-5:Ny),   # top-left
        (Nx-5:Nx, Ny-5:Ny), # top-right
        (1:5, 1:5),      # bottom-left
        (Nx-5:Nx, 1:5)   # bottom-right
    ]

    titles = ["top left grid corner", "top right grid corner", "bottom left grid corner", "bottom right grid corner"]

    fig = Figure(size=(800,800))
    
    Label(fig[0, 1:2], "Nx = $(Nx), Ny = $(Ny), dx = $(dx[1]), dy = $(dy[1])", halign = :center)
    
    for (idx, (xi, yi)) in enumerate(quadrants)

        ax = Axis(fig[div(idx-1,2)+1, mod(idx-1,2)+1]; xlabel="x", ylabel="y", title=titles[idx], xticklabelsize=10, yticklabelsize=10, xgridvisible = false, ygridvisible = false)

        # scatter centers
        scatter!(ax,repeat(xc[xi], inner=length(yc[yi])), repeat(yc[yi], outer=length(xc[xi])), color=:blue, markersize=4)
    
        # plot rectangles for each cell
        for (i, x) in enumerate(xc[xi])
            for (j, y) in enumerate(yc[yi])
                xs = [x - dx[xi[i]]/2, x + dx[xi[i]]/2, x + dx[xi[i]]/2, x - dx[xi[i]]/2] # left, right, right, left
                ys = [y - dy[yi[j]]/2, y - dy[yi[j]]/2, y + dy[yi[j]]/2, y + dy[yi[j]]/2] # bottom, bottom, top, top
                poly!(ax, xs, ys; color=:transparent, strokewidth=0.5, strokecolor=:red)
            end
        end
    end

    return fig

end


"""
$(TYPEDSIGNATURES)

Visualize a field.
"""
function visualize_field(x, y, data, mask; plot_title = "")

    fig = Figure(size = (900, 700))
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = plot_title, aspect = DataAspect())

    data[mask .!= 1] .= NaN

    if plot_title == "q"
        hm = heatmap!(ax, perSecond2perYear.(data') ./ 1e4; colormap = Reverse(:RdBu), colorrange = (0, 10)) # for q [m² s⁻¹]
    elseif plot_title == "N"
        hm = heatmap!(ax, data' .* 1e-6; colormap = Reverse(:RdBu), colorrange = (0, 10)) # for N [MPa]
    else 
        hm = heatmap!(ax, x, y, data')
    end

    Colorbar(fig[1, 2], hm)

    return fig
end


"""
$(TYPEDSIGNATURES)

Visualize an Oceananigans field.
"""
function visualize_field(field::Oceananigans.Field, mask; kwargs...)

    data = interior(field)[:, :, 1]
    x, y = nodes(field)
    visualize_field(x, y, data, mask; kwargs...)
    
end
