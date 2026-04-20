using Oceananigans
using MAT
using CairoMakie
using BenchmarkTools
using Oceananigans.BoundaryConditions: fill_halo_regions!
using ImageFiltering
using OffsetArrays
using Base.Threads
using SpecialFunctions

# TODO: add docstrings to all functions
# TODO: convert this file to a package

@kwdef struct PhysicalConstants{T <: AbstractFloat}

    ρ_w::T = 1000.0 # density of water [kg/m^3]
    ρ_i::T = 917.0 # density of ice [kg/m^3]
    g::T = 9.81 # graviational acceleration [m/s^2]
    n::T = 3.0 # exponent in Glen's flow law
    L_w::T = 3.34e5 # latent heat of fusion of water [J/kg]
    h_b::T = 0.1 # thickness of bed obstacles [m]
    α::T = 5/4 # exponent in Darcy-Weisbach relation
    β::T = 3/2 # exponent in Darcy-Weisbach relation
    f::T = 0.1 # friction coefficient for a turbulent flow
    F_till::T = 1.1 # factor for the conduit geometry in a till
    Q_c::T = 1.0 # critical water flux in a conduit [m^3/s]
    H_0::T = 0.1 # thickness of canals [m]
    l_c::T = 10000.0 # distance between conduits [m]
    sqrtπ::T = √π
    K::T = (2/π)^(0.25) * sqrt((π + 2) / (ρ_w * f)) # conductivity coefficient in Darcy-Weisbach relation [kg^{1−β} m^{2β−2α+1} s^{2β−3}]

end

function load_data()

    data_path = "/Users/taange001/Documents/Coding/Code from other people/Kazmierczak et al 2024 Code/results/thwaites_hydro/output/THWAITES2km_m3_HARD_toto.mat"
    data = matread(data_path)

    Nx, Ny = size(data["H"])
    mask = data["MASKo"]
    h = data["H"]
    b = data["B"]
    abs_v_b = data["ub"]
    A = data["A"]
    ṁ_over_ρ_w = data["Bmelt"] ./ 3.154e7 # units should be kg m^-2 s^-1
    
    xc = data["y"].*1000 # multiplication to convert to [m], also the length(data["y"]) = size(data["H"])[1] so I think they just stored them swapped
    yc = data["x"].*1000

    xlims = extrema(xc)
    ylims = extrema(yc)

    dx = xc[2] - xc[1]
    dy = yc[2] - yc[1]

    xlims = (xlims[1] - dx/2, xlims[2] + dx/2)
    ylims = (ylims[1] - dy/2, ylims[2] + dy/2)
    
    return Nx, Ny, xlims, ylims, mask, h, b, abs_v_b, A, ṁ_over_ρ_w
    
end

function initialize_grid(Nx, Ny, xlims, ylims; topology = (Bounded, Bounded, Flat), halo = (1, 1))

    # With bounded topology those dimensions have no-flux boundary conditions for any field that does not specify them
    return RectilinearGrid(topology = topology, size = (Nx, Ny), x = xlims, y = ylims, halo = halo)

end

function initialize_κ(b, Nx, Ny, type = "hard")

    if type == "hard"
        κ = zeros(Nx, Ny)
    elseif type == "soft"
        κ = ones(Nx, Ny)
    else
        κ = zeros(Nx, Ny)
        flag = b .< -1000
        κ[flag] .= 1.0
    end

    return κ

end

function update_ϕ₀!(fields, c)
    @. fields.ϕ₀ = c.ρ_i * c.g * fields.h + c.ρ_w * c.g * fields.b
    fill_halo_regions!(fields.ϕ₀)
end

function update_Po!(fields, c)
    fields.Po .= c.ρ_i * c.g * fields.h 
end

function potential_filling!(fields, c; iterations = 10)

    ϕ₀ = fields.ϕ₀
    fill_halo_regions!(ϕ₀) # TODO: can potentially remove if already called before calling this function
    fields.ϕ₀_tmp .= ϕ₀ # update ϕ₀_tmp to match ϕ₀ before we start the potential filling
    ϕ₀_tmp = fields.ϕ₀_tmp
    Nx, Ny = size(ϕ₀)

    for k in 1:iterations # one iteration might fill a small hole, but many can ensure even deep basins can be filled 
        @inbounds for j in 1:Ny
            for i in 1:Nx

                p = ϕ₀[i, j]

                # Cardinal neighbours of i, j
                # ϕ₀ is an offset array and we use a 1, 1 halo so indexing here is safe and thus we can use @inbounds
                # concretely, the indices of ϕ₀ range from 0 to Nx+1 and from 0 to Ny+1 for the x and y dimensions
                p1, p2 = ϕ₀[i+1, j], ϕ₀[i-1, j]
                p3, p4 = ϕ₀[i, j+1], ϕ₀[i, j-1]
                
                if p < p1 && p < p2 && p < p3 && p < p4 # identify a local minimum where the gradient is effectively pointing inward from all sides
                    ϕ₀_tmp[i, j] = (p1 + p2 + p3 + p4) / 4.0 # Laplacian smoothing step, replacing the pit with the average height of its neighbouring points effectively filling the hole
                end
            end
        end
        ϕ₀ .= ϕ₀_tmp
        fill_halo_regions!(ϕ₀)
    end

    @. fields.h = (fields.ϕ₀ - c.ρ_w * c.g * fields.b)/(c.ρ_i * c.g) # correction to h from potential filling

end

function update_H!(fields, c)
    @. fields.H_hard = sqrt(fields.S_inf)
    @. fields.H_soft = c.H_0 + (sqrt(fields.S_inf)/c.F_till - c.H_0) * exp(-fields.Q/c.Q_c)
    @. fields.H = (1 - fields.κ) * fields.H_hard + fields.κ * fields.H_soft
end

function initialize_fields(mask, h, b, abs_v_b, A, ṁ_over_ρ_w, c, grid)

    fields = (

        # Center fields
        h = CenterField(grid), # ice thickness
        b = CenterField(grid), # bedrock elevation
        abs_v_b = CenterField(grid), # absolute value of basal velocity
        A = CenterField(grid), # viscocity parameter in Glen's flow law
        mask = CenterField(grid), # = 1 when grounded ice where hydrology can exist
        ṁ_over_ρ_w = CenterField(grid), # melting rate divided by density of water
        ϕ₀ = CenterField(grid), # geometric potential
        ϕ₀_tmp = CenterField(grid), # geometric potential place holder to do potential filling 
        ψ_out = CenterField(grid), # integrated scalar water flux
        q = CenterField(grid), # distributed water flux
        Q = CenterField(grid), # volumetric water flux within each conduit
        S_inf = CenterField(grid), # far-field (away from margin) cross-sectional area of conduit
        H_hard = CenterField(grid), # thickness of conduits over hard bed 
        H_soft = CenterField(grid), # thickness of conduits over soft bed
        H = CenterField(grid), # thickness of conduits
        N_inf = CenterField(grid), # far-field effective pressure 
        N = CenterField(grid), # effective pressure
        Po = CenterField(grid), # hydrostatic ice overburden pressure
        abs_∇ϕ₀ = CenterField(grid), # absolute value of the gradient of the geometric potential
        κ = CenterField(grid), # indicator of the heterogenous content of the bed - TODO: should this be updating every iteration of the ice flow model? 
        corfac = CenterField(grid), # correction factor to go from ψ_out to q
        minus_∇ϕ₀_x  = CenterField(grid), # negative of the x-component of the gradient of the geometric potential
        minus_∇ϕ₀_y  = CenterField(grid), # negative of the y-component of the gradient of the geometric potential
        minus_∇ϕ₀_smoothed_x  = CenterField(grid), # negative of the x-component of the gradient of the geometric potential smoothed according to the stress-gradient coupling principle
        minus_∇ϕ₀_smoothed_y  = CenterField(grid), # negative of the y-component of the gradient of the geometric potential smoothed according to the stress-gradient coupling principle
        abs_∇ϕ₀_smoothed = CenterField(grid), # absolute value of the gradient of the geometric potential smoothed according to the stress-gradient coupling principle
        
        )

        fields.h .= h
        fields.b .= b
        fields.abs_v_b .= abs_v_b
        fields.A .= A
        fields.mask .= mask 
        fields.ṁ_over_ρ_w .= ṁ_over_ρ_w 
        
        for name in keys(fields)
            if name ∉ (:h, :b, :abs_v_b, :A, :mask, :ṁ_over_ρ_w)
                fields[name] .= 0.0
            end
        end

        fields.κ .= initialize_κ(b, grid.Nx, grid.Ny, "hard")
        update_Po!(fields, c)
        update_ϕ₀!(fields, c)
        fields.ϕ₀_tmp .= fields.ϕ₀ # does not need halo points update because we will not need to access them in potential filling

    return fields

end

function visualize_grid(grid)

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

    titles = ["top left", "top right", "bottom left", "bottom right"]

    fig = Figure(size=(800,800))
    
    Label(fig[0, 1:2], "Nx = $(Nx), Ny = $(Ny), dx = $(dx[1]), dy = $(dy[1])", halign = :center)
    
    for (idx, (xi, yi)) in enumerate(quadrants)

        ax = Axis(fig[div(idx-1,2)+1, mod(idx-1,2)+1]; xlabel="x", ylabel="y", title=titles[idx])

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

    display(fig)
end

function visualize_field(field, mask; plot_title = "")

    data = interior(field)[:, :, 1]

    x, y = nodes(field)

    fig = Figure(size = (900, 700))
    ax = Axis(fig[1, 1],
              xlabel = "x",
              ylabel = "y",
              title = plot_title,
              aspect = DataAspect())

    for I in findall(mask .!= 1)
        data[I] = NaN
    end

    # hm = heatmap!(ax, data')
    # hm = heatmap!(ax, x, y, data')
    hm = heatmap!(ax, data' .* 3.154e7 ./ 1e4; colormap = Reverse(:RdBu), colorrange = (0, 10)) # for q
    # hm = heatmap!(ax, data' .* 1e-6; colormap = Reverse(:RdBu), colorrange = (0, 10)) # for N

    Colorbar(fig[1, 2], hm)

    display(fig)

    return fig
end

function reflect_index(i, N)
    if i < 1
        return 2 - i # mirror at 1
    elseif i > N
        return 2*N - i # mirror at N
    else
        return i
    end
end

function smooth_potential_gradients!(destination, source, kernel)

    frb = (size(kernel,1)-1) ÷ 2  # integer kernel radius
    Nx, Ny = size(source)

    @threads for j in 1:Ny
        for i in 1:Nx

            tmp_sum = 0.0

            for j_kernel in -frb:frb
                for i_kernel in -frb:frb

                    ii = reflect_index(i + i_kernel, Nx)
                    jj = reflect_index(j + j_kernel, Ny)

                    @inbounds tmp_sum += kernel[i_kernel, j_kernel] * source[ii, jj]

                end
            end

            @inbounds destination[i, j] = tmp_sum

        end
    end

end

"""
    update_potential_gradients!(fields, grid) -> Return type

Compute the x, y gradients of the potential ϕ₀ in `fields`, 
apply a spatial smoothing based on ice thickness horizontal influence (see sec. 8.7.2 and eq. 8.98 of Cuffey & Paterson 2010), 
update halo regions, and store the absolute value of the gradient.

Hence this function updates the fields: minus_∇ϕ₀_x, minus_∇ϕ₀_y, minus_∇ϕ₀_x_tmp, minus_∇ϕ₀_y_tmp, abs_∇ϕ₀ and their halos.

The point of smoothing: A point in the field of the gradient of the geometric potential becomes the weighted sum of itself and its 
neighbours to accommodate the horizontal effect of changing ice thickness and changing bedrock elevation on the water routing, which
is called stress-gradient coupling.

# Arguments

- `fields`: Argument description
- `grid`: Argument description
"""
function update_potential_gradients!(fields, grid)

    # TODO: with regards to the need for smoothing: we say that the water routing follows ∇ϕ, so how can we see that what we are saying 
    # is mathematically inadequate and need to incorporate the horizontal changes in ∇ϕ?

    # Update the potential gradients before smoothing
    fields.minus_∇ϕ₀_x .= -∂x(fields.ϕ₀)
    fields.minus_∇ϕ₀_y .= -∂y(fields.ϕ₀)
    fields.abs_∇ϕ₀ .= sqrt(fields.minus_∇ϕ₀_x^2 + fields.minus_∇ϕ₀_y^2)

    # The water flux at a given point is influenced by variations in ice thickness some distance away. To account for this we perform a convolution of the gradient of the potential
    # such that the influence of nearby points is now incorporated into the value of the gradient of the potential at that point.
    longcoupwater = 5.0
    h_active = @views interior(fields.h, :, :, 1)[fields.mask .== 1] # TODO: this still allocates memory when we do fields.mask .== 1
    h_avg = max(mean(h_active), 10.0) # see Cuffey & Paterson 2010 end of sec. 8.7.2 for the maximum value
    
    # Radius of influence
    scale = h_avg * longcoupwater * 2.0 # represents half the radius of the influence zone of the ice thickness on water
    width = 2.0 * scale # radius of influence, after which the cone ends and the weighting of any points beyond that range is 0
    Δ = grid.Δxᶜᵃᵃ # here we take that dx = dy, assume even spacing, and use the distance between centers # TODO we assume here dx = dy, should we handle it more generally like taking avg?
    if width <= Δ
        scale = Δ / 2.0 + 1.0
    end
    
    # Size of the kernel
    maxlevel = 2 * round(Int, width / Δ - 0.5) + 1 # e.g. 9

    # Filter radius boundary is the radius of the kernel, it calculates how many pixels extend from the center of that kernel window to its edge, for example for 9 x 9 the frb = 4
    frb = Int((maxlevel - 1) / 2) 

    # Filtering kernel
    kernel = zeros(maxlevel, maxlevel)
    for nj in 1:maxlevel, ni in 1:maxlevel
        dist = sqrt((Δ * (ni - frb - 1))^2 + (Δ * (nj - frb - 1))^2) / scale
        kernel[ni, nj] = max(0.0, 1.0 - dist / 2.0)
    end
    kernel ./= sum(kernel) # normalization is important as it ensures that the smoothing doesn't change the total amount of potential but only redistributes it

    # Initialize the views of the arrays for imfilter! # TODO: understand why we cant give the whole data to imfilter! so as to avoid this interior slicing
    @views minus_∇ϕ₀_smoothed_x = interior(fields.minus_∇ϕ₀_smoothed_x)[1:end, :, 1]
    @views minus_∇ϕ₀_smoothed_y = interior(fields.minus_∇ϕ₀_smoothed_y)[:, 1:end, 1]
    @views minus_∇ϕ₀_x = interior(fields.minus_∇ϕ₀_x)[1:end, :, 1]
    @views minus_∇ϕ₀_y = interior(fields.minus_∇ϕ₀_y)[:, 1:end, 1]

    # This slides the cone kernel over every pixel of the gradient and each pixel's new value becomes a weighted average of its neighbors and itself, 
    # at the boundaries we reflect the array when the kernel is outside, imfilter! allocates quite a bit of memory e.g. see discussion here https://discourse.julialang.org/t/how-to-accelerate-the-imfiter-operation/102965/3
    imfilter!(minus_∇ϕ₀_smoothed_x, minus_∇ϕ₀_x, centered(kernel))
    imfilter!(minus_∇ϕ₀_smoothed_y, minus_∇ϕ₀_y, centered(kernel))
    # TODO: Remove if desired: the two below are my custom function for imfilter which is slower but does not allocate (see comment above)
    # smooth_potential_gradients!(minus_∇ϕ₀_smoothed_x, minus_∇ϕ₀_x, centered(kernel))
    # smooth_potential_gradients!(minus_∇ϕ₀_smoothed_y, minus_∇ϕ₀_y, centered(kernel))
    
    # Halo used in the accumulate_ψ_out! so we update them according to the BC of the field
    fill_halo_regions!(fields.minus_∇ϕ₀_smoothed_x)
    fill_halo_regions!(fields.minus_∇ϕ₀_smoothed_y) 

    # Update the absolute value of the gradient of the potential
    @. fields.abs_∇ϕ₀_smoothed = abs(fields.minus_∇ϕ₀_smoothed_x) + abs(fields.minus_∇ϕ₀_smoothed_y) # here we don't need @at (Center, Center, Center) because fields.abs_∇ϕ₀ is initialized already as CenterField
    fill_halo_regions!(fields.abs_∇ϕ₀_smoothed) # used in accumulate_ψ_out!

end

function accumulate_ψ_out!(fields, i, j, grid, w_matrix) # TODO: the w_matrix argument is a test to be removed afterwards

    if fields.ψ_out[i, j] ≥ 0.0 # we initialize all unvisited cells with -1.0 and q is positive semi-definite so if its positive then it has been visited and ψ_out calculated # TODO: why is it positive semi-definite if ṁ can go negative?
        return fields.ψ_out[i, j]
    end

    # TODO: make sure that a negative ṁ case (from refreezing of water) can only reduce the ψ_out to 0 and not make it negative which would be an unphysical disappearance of water into its neighbours
    fields.ψ_out[i, j] = fields.ṁ_over_ρ_w[i, j] * grid.Δxᶜᵃᵃ * grid.Δyᵃᶜᵃ # TODO: this assumes even spacing, if the spacing is not even these will return arrays instead of floats

    for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))

            ni, nj = i+di, j+dj 

            # The normal unit vector (di, dj) as the outgoing normal vector and since we want w to represent incoming flux into i, j, we then take the dot product of minus_∇ϕ₀ with -(dj, dj).
            # In principle, we can have 4 mask == 1 neighbours and only 1 downstream neighbour for which case, in principle, the sum of downstream weights should be 1. The most common case when we have 4 mask == 1
            # neighbours is to have 2 downstream neighbours and again the sum of downstream weights should - in principle - be equal to 1. Given the interpolation of Oceananigans for abs_∇ϕ₀ from minus_∇ϕ₀_x,y and the 
            # way we decided to calculate the weight in general, the sum of downstream weights is close but not equal to 1. Indeed, once an interpolation is present, it is not equal to 1 for any way of the 3 ways below to calculate the weight.
            # When there is a neighbour with mask != 1 but the direction of ∇ϕ₀ would have made it a downstream neighbour, then the sum of downstream weights will be less than 1 and there is loss of water flux into the 
            # sea or to land without grounded ice. Water flux can refreeze and become ice through negative ṁ, but once we lose that water flux this can never happen with that water, so we will have potential ice leakage from the ice flow model 
            # through the hydrology model. # TODO: discuss this with Alex.
            w = -(fields.minus_∇ϕ₀_smoothed_x[ni, nj]*di + fields.minus_∇ϕ₀_smoothed_y[ni, nj]*dj) / fields.abs_∇ϕ₀_smoothed[ni, nj]
        
            if w > 0
                w_matrix[ni, nj] += w # TODO: this is a test to be removed afterwards
                fields.ψ_out[i, j] += accumulate_ψ_out!(fields, ni, nj, grid, w_matrix) * w
            end

    end
    return fields.ψ_out[i, j]
end

function recursive_update_ψ_out!(fields, grid)
    
    for j in 1:grid.Ny
        for i in 1:grid.Nx
            if fields.mask[i, j] == 1
                fields.ψ_out[i, j] = -1.0
            end
        end
    end

    # TODO: This is a test to be removed afterwards
    w_matrix = zeros(grid.Nx, grid.Ny)

    for j in 1:grid.Ny
        for i in 1:grid.Nx
            if fields.mask[i, j] == 1
                accumulate_ψ_out!(fields, i, j, grid, w_matrix) # Argument w_matrix is test to be removed afterwards
            end
        end
    end

    # TODO: This is a test to be removed afterwards
    # for j in 1:grid.Ny
    #     for i in 1:grid.Nx
    #         if fields.mask[i, j] == 1
    #             println(i, " ", j, " ", w_matrix[i, j])
    #         end
    #     end
    # end

end

"""
    update_S_inf!(fields, c) @. fields.S_inf = c.K^(-1/c.α) * fields.abs_∇ϕ₀^((1-c.β)/c.α) * fields.Q^(1/c.α) -> Return type

Update the field representing the cross-sectional area of the conduits. The _inf implies that we used the approximation ∇ϕ ≈ ∇ϕ₀ to derive
the formula for the cross-sectional area of the conduits S_inf. This approximation is valid far from the grounding line as shown from
 Eqs. 5a, 6b of of Kazmierczak et al 2024. 

# Arguments

- `fields`: named tuple of grid fields
- `c`: physical constants
"""
function update_S_inf!(fields, c)
    @. fields.S_inf = c.K^(-1/c.α) * fields.abs_∇ϕ₀^((1-c.β)/c.α) * fields.Q^(1/c.α) 
end


"""
    update_N_inf!(fields, c) -> Return type

Compute the far-field (away from ground line) effective pressure N_inf = Po - Pw.

For the min, max limits used: 

    By definition Neff = Po - Pw ⟹ Neff ≤ Po.
    Enforce N_inf ≥ 0.02*Po to avoid near-flotation (N→0), with 0.02 being heuristic as used in Kori-ULB.

# Arguments

- `fields`: named tuple of grid fields
- `c`: physical constants
"""
function update_N_inf!(fields, c)

    # Lower limit of 0.02*Po ≤ N_inf ≤ Po comes from Eq. (20) of Beuler and van Pelt 2015 where they use a symbol δ = 0.02 # TODO: this lower bound on the far-field effective pressure should be discussed, also why not set the limits on N_eff instead of N_inf?
    @. fields.N_inf.data = min(max(
        ((fields.H.data/fields.S_inf.data)^2*((c.ρ_i*c.L_w*fields.abs_v_b.data*c.h_b + fields.Q.data*fields.abs_∇ϕ₀.data)/(2.0*c.n^(-c.n)*c.ρ_i*c.L_w*fields.A.data)))^(1.0/c.n),
        0.02*fields.Po.data),
        fields.Po.data)

end

"""
    update_N!(fields, c) -> Return type

Update the effective pressure N with the grounding line correction included as described by Eq. 7 of Kazmierczak et al 2024.

erf = 1-erfc.

# Arguments

- `fields`: named tuple of grid fields
- `c`: physical constants
"""
function update_N!(fields, c)
    @. fields.N.data = (1-erfc((c.sqrtπ*fields.ϕ₀.data)/(2*fields.N_inf.data))) * fields.N_inf.data 
end

function main()

    c = PhysicalConstants{Float64}()

    Nx, Ny, xlims, ylims, mask, h, b, abs_v_b, A, ṁ_over_ρ_w = load_data()

    grid = initialize_grid(Nx, Ny, xlims, ylims)

    fields = initialize_fields(mask, h, b, abs_v_b, A, ṁ_over_ρ_w, c, grid)

    # With every ice flow model iteration we need to update the fields: h, b, abs_v_b, A, mask (for grounded ice), ṁ_over_ρ_w and consequently ϕ₀ and Po

    potential_filling!(fields, c) # updates ϕ₀ and consequently also updates h to reflect changes in ϕ₀

    update_potential_gradients!(fields, grid) 

    recursive_update_ψ_out!(fields, grid)
     
    # Budd et al 1996, and Le Brocq et al 2006, 2009 discuss the corfac and conversion. The corfac is to make the flux independent of the grid orientation and then 
    # with that we can think of the width that the flux q is based on - to be defined - as Δx. The correction factor is defined corfac = |cosθ| + |sinθ| where θ is the direction
    # of flow w.r.t. the x-axis. We then have scalar flux magnitude = flux through cell / corfac. Then q = scalar flux magnitude / Δx.
    fields.corfac .= fields.abs_∇ϕ₀_smoothed / sqrt(fields.minus_∇ϕ₀_smoothed_x^2 + fields.minus_∇ϕ₀_smoothed_y^2)
    
    # Limits on q are heuristic and chosen by Frank Pattyn for numerical stability.
    @. fields.q.data = min(max(fields.ψ_out.data / (fields.corfac.data * grid.Δxᶜᵃᵃ), 0), 1e5) # TODO: here we use dx by assuming dx = dy and even grid spacing, should we handle it more generally?

    # # Effective pressure
    
    @. fields.Q = fields.q * c.l_c # volumetric water flux in a single channel
    update_S_inf!(fields, c)
    update_H!(fields, c)
    update_Po!(fields, c)    
    update_N_inf!(fields, c)
    update_N!(fields, c)

    # visualize_field(fields.abs_∇ϕ₀, fields.mask)
    # visualize_field(fields.minus_∇ϕ₀_x, fields.mask)
    # visualize_field(fields.minus_∇ϕ₀_y, fields.mask)
    # visualize_field(fields.ψ_out, fields.mask)
    visualize_field(fields.q, fields.mask)
    # visualize_field(fields.N, fields.mask)

    # Neff = CenterField(grid)
    # data_path = "/Users/takisangelides/Documents/Post-doc/Kazmierczak et al 2024 Code/results/thwaites_hydro/output/THWAITES2km_m3_HARD_toto.mat"
    # data = matread(data_path)
    # Neff .= data["Neff"]
    # visualize_field(Neff, fields.mask)
    
    return

end

main()