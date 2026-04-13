"""
$(TYPEDSIGNATURES)

Update the water flux q via a recursive algorithm.
"""
function update_q!(HM::HydrologyModel)

    fields, c, grid = HM.fields, HM.c, HM.grid

    # Compute ϕ₀ and prepare a temporary place-holder for the potential filling
    update_ϕ₀!(fields, c)
    fields.ϕ₀_tmp .= fields.ϕ₀ # update ϕ₀_tmp to match ϕ₀ before we start the potential filling
    fill_halo_regions!(fields.ϕ₀_tmp)
    
    # Updates ϕ₀ and consequently also updates h to reflect changes in ϕ₀. It fills the local minima of ϕ₀ to avoid water getting stuck in there.
    potential_filling!(fields, c, grid) 

    # Update the gradients of the geometric potential ϕ₀
    update_potential_gradients!(fields)

    # Smoothen the gradients of ϕ₀ to incorporate the concept of stress-gradient coupling.
    update_smoothed_potential_gradients!(fields, grid) 

    # Compute ψ_out via a recursive algorithm that traverses upsteam to compute ψ_out at the highest ϕ₀ and then come back down again, for each unvisited point on the grounded ice grid.
    update_ψ_out!(fields, grid)
     
    # This is the factor to go from ψ_out to q. It can be derived using the definition of ψ_out = integral of q ⋅ n with n the outward normal 
    # and the line integral is over the whole sides of the grid cell that water leaves the cell, which can be 1 side for when q points to any of the 4 cardinal neighbours's center or 2 sides.
    fields.corfac .= fields.abs_∇ϕ₀_smoothed / sqrt(fields.minus_∇ϕ₀_smoothed_x^2 + fields.minus_∇ϕ₀_smoothed_y^2) 
    
    # Limits on q are heuristic and chosen by Frank Pattyn for numerical stability.
    @. fields.q.data = min(max(fields.ψ_out.data / (fields.corfac.data * grid.Δxᶜᵃᵃ), 0), 1e5) # TODO: here we use dx by assuming dx = dy and even grid spacing, should we handle it more generally?

    return nothing

end


"""
$(TYPEDSIGNATURES)

Update the geometric potential ϕ₀ and set the values for the halo points.
"""
function update_ϕ₀!(fields, c)
    
    @. fields.ϕ₀ = c.ρ_i * c.g * fields.h + c.ρ_w * c.g * fields.b
    fill_halo_regions!(fields.ϕ₀)

    return nothing
    
end


"""
$(TYPEDSIGNATURES)

Updates ϕ₀ and consequently also updates h to reflect changes in ϕ₀. It fills the local minima of ϕ₀ to avoid water getting stuck in there.
"""
function potential_filling!(fields, c, grid; iterations = 10)

    ϕ₀ = fields.ϕ₀
    ϕ₀_tmp = fields.ϕ₀_tmp

    for k in 1:iterations # one iteration might fill a small hole, but many can ensure even deep basins can be filled 
        @inbounds for j in 1:grid.Ny
            for i in 1:grid.Nx

                p = ϕ₀[i, j]

                # Cardinal neighbours of i, j. ϕ₀ is an offset array and we use a 1, 1 halo so indexing here is safe and thus we can use @inbounds
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

    return nothing

end


"""
$(TYPEDSIGNATURES)

Compute (the negative of) the gradients of the geometric potential ϕ₀ and its absolute value.
"""
function update_potential_gradients!(fields)

    fields.minus_∇ϕ₀_x .= -∂x(fields.ϕ₀)
    fields.minus_∇ϕ₀_y .= -∂y(fields.ϕ₀)
    fields.abs_∇ϕ₀ .= sqrt(fields.minus_∇ϕ₀_x^2 + fields.minus_∇ϕ₀_y^2)

    return nothing

end


"""
$(TYPEDSIGNATURES)


"""
function update_smoothed_potential_gradients!(fields, grid)

    # TODO: with regards to the need for smoothing: we say that the water routing follows ∇ϕ, so how can we see that what we are saying 
    # is mathematically inadequate and need to incorporate the horizontal changes in ∇ϕ?

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

    return nothing

end

"""
$(TYPEDSIGNATURES)

Helper function to the recursive function to calculate the ψ_out for every grid cell that has grounded ice.
"""
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


"""
$(TYPEDSIGNATURES)

Recursive function to calculate the ψ_out for every grid cell that has grounded ice.
We initialize ψ_out to -1 since it is by definition positive semi-definite and hence we 
know that if a grid point has negative ψ_out, it is still unvisited. 
"""
function update_ψ_out!(fields, grid)
    
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

    return nothing

end
