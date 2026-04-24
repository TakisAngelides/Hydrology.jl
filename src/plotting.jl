"""
$(TYPEDSIGNATURES)

Visualize a field.
"""
function visualize_field(x, y, data; 
        plot_title = "", 
        transpose_data = false, 
        colorrange = extrema(filter(!isnan, data)), 
        display_flag = false, 
        colormap = Reverse(:RdBu), 
        colorscale = identity,
        savefig_path = nothing
    )

    fig = Figure(size = (900, 700))
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = plot_title, aspect = DataAspect())

    if transpose_data
        data = data'
    end
    
    hm = heatmap!(ax, x, y, data, colormap = colormap, colorrange = colorrange, colorscale = colorscale)

    Colorbar(fig[1, 2], hm)

    if display_flag
        display(fig)
    end

    if savefig_path !== nothing
        save(savefig_path, fig)
    end

    return fig
end


"""
$(TYPEDSIGNATURES)

Visualize an Oceananigans field.
"""
function visualize_field(field::Oceananigans.Fields.Field; kwargs...)

    data = interior(field)[:, :, 1]
    x, y = nodes(field)
    visualize_field(x, y, data; kwargs...)
    
end


"""
$(TYPEDSIGNATURES)

Set a field to a given input value where the mask is not 1.
"""
function mask_field(field::Oceananigans.Fields.Field, mask, value)

    Nx, Ny = size(mask)
    res = deepcopy(field)
    for j in 1:Ny
        for i in 1:Nx
            if mask[i, j] != 1.0
                res[i, j, 1] = value
            end
        end
    end

    return res
    
end


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
