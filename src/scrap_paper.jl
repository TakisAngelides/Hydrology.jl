using MAT
using LinearAlgebra
using BenchmarkTools
using SpecialFunctions
using WGLMakie
using CairoMakie
using ImageFiltering
using Statistics
using ColorSchemes
import ColorSchemes.Colors
using DSP
using Base.Threads
using OffsetArrays
using Random

data_path = "/Users/taange001/Documents/Coding/Code from other people/Kazmierczak et al 2024 Code/results/thwaites_hydro/output/THWAITES2km_m3_HARD_toto.mat"
data = matread(data_path)

println(keys(data))

# MASK = data_THWAITES2km["MASK"]

# display(heatmap(MASK))

# mask = MASK .== 1

# display(mask)

# display(data["dx"])

# println(unique(MASK))

# display(heatmap(mask))

# b = data_THWAITES2km["B"]

# # println(size(mask))

# # println(size(b))

# display(heatmap(b))

# h = data_THWAITES2km["H"]

# display(heatmap(h))

# # println(size(h))

# v, vx, vy, MASKo = data_THWAITES2km["v"], data_THWAITES2km["vx"], data_THWAITES2km["vy"], data_THWAITES2km["MASKo"]

# display(heatmap(v))
# display(heatmap(vx))
# display(heatmap(vy))
# display(heatmap(MASKo))
# display(heatmap(data["p"]'))
# display(data["p"])

# Nx, Ny = size(data["H"])
# mask = data["MASKo"]
# Neff = data["Neff"]
# for j in 1:Ny
#     for i in 1:Nx
#         if mask[i, j] != 1
#             Neff[i, j] = NaN
#         end
#     end
# end
# display(heatmap(Neff' ./ 1e6; colormap = Reverse(:RdBu), colorrange = (0, 10)))

# Plots.heatmap(data["A"])
# flw = data["flw"]./10^4
# display(Plots.heatmap(flw, c=:redblue))

# -----------------------------------------------------------------------

# Checking fig 4a of https://www.sciencedirect.com/science/article/pii/S0098300406000781

function accumulate_q!(q, i, j, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀) 

    if q[i, j] ≥ 0.0
        return q[i, j]
    end

    q[i, j] = ṁ_over_ρ_w[i, j] * Δ²

    for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))

            ni, nj = i+di, j+dj 

            if (ni < 1) || (ni > size(q, 1)) || (nj < 1) || (nj > size(q, 2)) 
                continue
            end

            if (di > 0) || (dj > 0) # according to d-grid if the neighbour is north or east we need the gradient at i, j otherwise at ni, nj
                w = -(minus_∇ϕ₀_x[i, j]*di + minus_∇ϕ₀_y[i, j]*dj) / abs_∇ϕ₀[i, j]
            else
                w = -(minus_∇ϕ₀_x[ni, nj]*di + minus_∇ϕ₀_y[ni, nj]*dj) / abs_∇ϕ₀[ni, nj]
            end
            w = -(minus_∇ϕ₀_x[ni, nj]*di + minus_∇ϕ₀_y[ni, nj]*dj) / abs_∇ϕ₀[ni, nj]
            w = -(minus_∇ϕ₀_x[i, j]*di + minus_∇ϕ₀_y[i, j]*dj) / abs_∇ϕ₀[i, j]
            if w > 0
                q[i, j] += accumulate_q!(q, ni, nj, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀) * w
            end
    end
    return q[i, j]
end

function recursive_update_q!(q, Nx, Ny, mask, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀)

    q[mask .== 1] .= -1.0 

    for j in 1:Ny
        for i in 1:Nx
            if mask[i, j] == 1 
                accumulate_q!(q, i, j, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀)
            end
        end
    end

end

function potential_filling!(ϕ₀, Nx, Ny, iterations=10)

    pot_next = copy(ϕ₀) # to prevent directional bias, we use the old points for all the new points, hence we need a copy
    for k in 1:iterations # one iteration might fill a small hole, but many can ensure even deep basins can be filled 
        for j in 2:Ny-1
            for i in 2:Nx-1

                p = ϕ₀[i, j]

                # Cardinal neighbours of i, j
                p1, p2 = ϕ₀[i+1, j], ϕ₀[i-1, j]
                p3, p4 = ϕ₀[i, j+1], ϕ₀[i, j-1]
                
                if p < p1 && p < p2 && p < p3 && p < p4 # identify a local minimum where the gradient is effectively pointing inward from all sides
                    pot_next[i, j] = (p1 + p2 + p3 + p4) / 4.0 # Laplacian smoothing step, replacing the pit with the average height of its neighbouring points effectively filling the hole
                end
            end
        end
        ϕ₀ .= pot_next
    end
end

function update_potential_gradients_smoothed!(minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, ϕ₀, h, mask, Nx, Ny, Δ)

    for j in 1:Ny-1, i in 1:Nx-1
        minus_∇ϕ₀_x[i, j] = -(ϕ₀[i+1, j] - ϕ₀[i, j]) / (2 * Δ)
        minus_∇ϕ₀_y[i, j] = -(ϕ₀[i, j+1] - ϕ₀[i, j]) / (2 * Δ)
    end

    longcoupwater = 1.0 
    h_avg = mean(h[mask .== 1])
    
    scale = h_avg * longcoupwater * 2.0 # represents the radius of the influence zone of the ice thickness on water
    width = 2.0 * scale
    if width <= Δ
        scale = Δ / 2.0 + 1.0
    end
    
    maxlevel = 2 * round(Int, width / Δ - 0.5) + 1 # e.g. 9
    # Filter radius boundary is the radius of the kernel, it calculates how many pixels extend from the center of that kernel window to its edge
    # for example for 9 x 9 the frb = 4
    frb = Int((maxlevel - 1) / 2) 

    kernel = zeros(maxlevel, maxlevel)
    for nj in 1:maxlevel, ni in 1:maxlevel
        dist = sqrt((Δ * (ni - frb - 1))^2 + (Δ * (nj - frb - 1))^2) / scale
        kernel[ni, nj] = max(0.0, 1.0 - dist / 2.0)
    end
    kernel ./= sum(kernel) # normalization is important as it ensures that the smoothing doesn't change the total amount of potential but only redistributes it

    # This slides the cone kernel over every pixel of the gradient and each pixel's new value becomes a weighted average of its neighbors
    minus_∇ϕ₀_x .= imfilter(minus_∇ϕ₀_x, centered(kernel), "reflect")
    minus_∇ϕ₀_y .= imfilter(minus_∇ϕ₀_y, centered(kernel), "reflect")

    @. abs_∇ϕ₀ = abs(minus_∇ϕ₀_x) + abs(minus_∇ϕ₀_y)
end

function iterative_update_q!(q, Nx, Ny, mask, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, ϕ₀)

    q .= 0.0
    indices = Tuple{Int, Int}[]
    for j in 1:Ny, i in 1:Nx
        if mask[i, j] == 1
            q[i, j] = ṁ_over_ρ_w[i, j] * Δ²
            push!(indices, (i, j))
        end
    end

    sort!(indices, by = idx -> ϕ₀[idx[1], idx[2]], rev = true)
    
    for (i, j) in indices
    
        for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))

            ni, nj = i + di, j + dj

            if (1 <= ni <= Nx) && (1 <= nj <= Ny) && mask[ni, nj] == 1

                if (di > 0) || (dj > 0) # according to d-grid if the neighbour is north or east we need the gradient at i, j otherwise at ni, nj
                    w = (minus_∇ϕ₀_x[i, j]*di + minus_∇ϕ₀_y[i, j]*dj) / abs_∇ϕ₀[i, j]
                else
                    w = (minus_∇ϕ₀_x[ni, nj]*di + minus_∇ϕ₀_y[ni, nj]*dj) / abs_∇ϕ₀[ni, nj]
                end

                if w > 0
                    q[ni, nj] += q[i, j] * w
                end

            end
        end
    end
end

function update_potential!(ϕ₀, h, b, ρ_w, ρ_i, g)
    @. ϕ₀ = ρ_i * g * h + ρ_w * g * b
end

function main()

    Nx, Ny = 9, 9
    Δ = 1.0
    Δ² = Δ^2
    mask = ones(Nx, Ny)
    ṁ_over_ρ_w = zeros(Nx, Ny)
    ṁ_over_ρ_w[8, 8] = 1.0
    h = zeros(Nx, Ny)
    b = zeros(Nx, Ny)
    deg = 45.0
    for j in 1:Ny, i in 1:Nx
        b[i, j] = i*cos(deg*π/180) + j*sin(deg*π/180)
    end
    
    ρ_w, ρ_i, g = 1000.0, 917.0, 9.81
    
    ϕ₀ = zeros(Nx, Ny)
    minus_∇ϕ₀_x = zeros(Nx, Ny)
    minus_∇ϕ₀_y = zeros(Nx, Ny)
    abs_∇ϕ₀ = zeros(Nx, Ny)
    q = zeros(Nx, Ny)
    Q = zeros(Nx, Ny)

    # Distributed water flux

    update_potential!(ϕ₀, h, b, ρ_w, ρ_i, g)
    
    potential_filling!(ϕ₀, Nx, Ny, 10) 

    @. h = (ϕ₀ - ρ_w*g*b)/(ρ_i*g) # correction to h from potential filling

    gdsmag = zeros(Nx, Ny)
    gx_raw = zeros(Nx, Ny)
    gy_raw = zeros(Nx, Ny)
    for j in 1:Ny-1, i in 1:Nx-1
        gx_raw[i, j] = (ϕ₀[i+1, j] - ϕ₀[i, j]) / (2 * Δ)
        gy_raw[i, j] = (ϕ₀[i, j+1] - ϕ₀[i, j]) / (2 * Δ)
    end
    @. gdsmag = abs(gx_raw) + abs(gy_raw)

    update_potential_gradients_smoothed!(minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, ϕ₀, h, mask, Nx, Ny, Δ)

    recursive_update_q!(q, Nx, Ny, mask, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀)

    corfac = ones(Nx, Ny)
    denom = @. sqrt(minus_∇ϕ₀_x^2 + minus_∇ϕ₀_y^2)
    for i in 1:Nx, j in 1:Ny
        if gdsmag[i, j] > 0 && denom[i, j] > 0
            corfac[i, j] = gdsmag[i, j] / denom[i, j]
        end
    end
    @. q = min(max(q / corfac / Δ, 0.0), 1e5) # distributed water flux
    
    # Plotting

    [x[mask .!= 1] .= NaN for x in (b, h, Q, q)]

    WGLMakie.activate!()
    fig = Figure(size = (1000, 800))

    ax = Axis3(
        fig[1, 1],
        title = "Warner $(deg)°",
        azimuth = -0.6π,
        elevation = 0.2π,
        perspectiveness = 0.5,
        aspect = (1, 1, 1)
    )

    surface!(ax, 1:Nx, 1:Ny, b,
        color = q,
        colormap = :Greys,
        shading = NoShading,
        colorrange = (0, 0.4)
    )

    wireframe!(ax, 1:Nx, 1:Ny, b,
        color = (:black, 0.3),
        linewidth = 1
    )

    zlims!(ax, 0, 20)
    hidedecorations!(ax)

    display(fig)

    CairoMakie.activate!()

    fig = Figure(size = (1000, 800))

    ax1 = Axis(fig[1, 1], aspect = DataAspect(), xgridvisible = false, ygridvisible = false)
    hm1 = heatmap!(ax1, q')
    Colorbar(fig[1, 2], hm1)

    display(fig)

end

# main()

# -----------------------------------------------------------------------

# Convolution

# u = rand(346, 288) 
# v = rand(45, 45)

# out1 = zeros(346, 288)
# out2 = zeros(346+45-1, 288+45-1) # cannot handle boundary conditions
# out3 = zeros(346, 288)

# function reflect_index(i, N)
#     if i < 1
#         return 2 - i # mirror at 1
#     elseif i > N
#         return 2*N - i # mirror at N
#     else
#         return i
#     end
# end

# function smooth_potential_gradients!(source, destination, kernel)

#     frb = (size(kernel,1)-1) ÷ 2  # integer kernel radius
#     Nx, Ny = size(source)
    
#     @threads for j in 1:Ny
#         for i in 1:Nx

#             tmp_sum = 0.0

#             for j_kernel in -frb:frb
#                 for i_kernel in -frb:frb

#                     ii = reflect_index(i + i_kernel, Nx)
#                     jj = reflect_index(j + j_kernel, Ny)

#                     @inbounds tmp_sum += kernel[i_kernel, j_kernel] * source[ii, jj]

#                 end
#             end

#             @inbounds destination[i, j] = tmp_sum

#         end
#     end

# end

# # @btime conv!($out2, $u, $v) # 1.440 ms (31 allocations: 2.05 MiB)

# @btime smooth_potential_gradients!(u, out1, centered(v)) #  with 16 threads

# @btime imfilter!(out3, u, centered(v), "reflect") # 1.948 ms (112 allocations: 11.64 MiB)

# println(isapprox(out1, out3)) # true - test passed

# -----------------------------------------------------------------------

# Taking Oceananigans fields to a non-integer power and applying min and max limits all while using @. syntax

# @kwdef struct PhysicalConstants{T <: AbstractFloat}

#     n::T = 3.0

# end

# c = PhysicalConstants()

# Random.seed!(10)

# Nx, Ny = 500, 500

# grid = RectilinearGrid(topology = (Bounded, Bounded, Flat), size = (Nx, Ny), x = (0, 4), y = (0, 4), halo = (1, 1))

# field = CenterField(grid)
# field_tmp = CenterField(grid)

# a = rand(Nx, Ny)
# field .= a
# Oceananigans.fill_halo_regions!(field)

# b = rand(Nx, Ny)
# field_tmp .= b
# Oceananigans.fill_halo_regions!(field)

# # display(field.data)

# @btime @. field.data = min(max((field.data/field_tmp.data)^(1.0/c.n), 0.2*field_tmp.data), field_tmp.data)

# # display(field.data)

# -----------------------------------------------------------------------

println("Finished.")
