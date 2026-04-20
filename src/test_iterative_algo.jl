using BenchmarkTools
using SpecialFunctions
using CairoMakie
using MAT
using ImageFiltering
using Statistics

function accumulate_Q!(Q, i, j, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, w_matrix) 

    if Q[i, j] ≥ 0.0
        return Q[i, j]
    end

    Q[i, j] = ṁ_over_ρ_w[i, j] * Δ²

    for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))

            ni, nj = i+di, j+dj 

            if (ni < 1) || (ni > size(Q, 1)) || (nj < 1) || (nj > size(Q, 2)) 
                continue
            end

            # if (di > 0) || (dj > 0) 
            #     w = -(minus_∇ϕ₀_x[i, j]*di + minus_∇ϕ₀_y[i, j]*dj) / abs_∇ϕ₀[i, j]
            # else
            #     w = -(minus_∇ϕ₀_x[ni, nj]*di + minus_∇ϕ₀_y[ni, nj]*dj) / abs_∇ϕ₀[ni, nj]
            # end

            w = -(minus_∇ϕ₀_x[ni, nj]*di + minus_∇ϕ₀_y[ni, nj]*dj) / abs_∇ϕ₀[ni, nj]

            # w = -(minus_∇ϕ₀_x[i, j]*di + minus_∇ϕ₀_y[i, j]*dj) / abs_∇ϕ₀[i, j]

            if w > 0
                w_matrix[ni, nj] += w
                Q[i, j] += accumulate_Q!(Q, ni, nj, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, w_matrix) * w
            end
    end
    return Q[i, j]
end

function recursive_update_Q!(Q, Nx, Ny, mask, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀)

    Q[mask .== 1] .= -1.0 

    w_matrix = zeros(Nx, Ny)

    for j in 1:Ny
        for i in 1:Nx
            if mask[i, j] == 1 
                accumulate_Q!(Q, i, j, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, w_matrix)
            end
        end
    end

    # for element in w_matrix
    #     if 0.9 < element < 1.1
    #         println(element)
    #     end
    # end

end

"""
    potential_filling!(ϕ₀, Nx, Ny, iterations=10) -> Return type

Subglacial topography is often "noisy." Digital Elevation Models (DEMs) 
frequently contain sinks (pits)—single grid cells where the potential is 
lower than all four neighbors. In a routing model, water is "greedy": it only moves downhill. 
If it hits a pit, it gets stuck, and the flux Q for the rest of the glacier downstream becomes zero. 
This function "lifts" the bottom of these pits so water can continue to flow.

# Arguments

- `ϕ₀`: Argument description
- `Nx`: Argument description
- `Ny`: Argument description
- `iterations`: Argument description
    (**Default**: `10`)
"""
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

"""
    update_potential_gradients_smoothed!(minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, ϕ₀, h, mask, Nx, Ny, Δ) -> Return type

The distance over which subglacial water "feels" the ice pressure is proportional to the ice thickness h.

The kernel weight square matrix follows the logic that the center has the highest weight, and the weight drops to zero at the edges.
It follows the equation W(d) = max(0, 1 - d/(2*scale)), where d is the distance from the center and scale is the radius of influence of h on water.

# Arguments

- `minus_∇ϕ₀_x`: Argument description
- `minus_∇ϕ₀_y`: Argument description
- `abs_∇ϕ₀`: Argument description
- `ϕ₀`: Argument description
- `h`: Argument description
- `mask`: Argument description
- `Nx`: Argument description
- `Ny`: Argument description
- `Δ`: Argument description
"""
function update_potential_gradients_smoothed!(minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, ϕ₀, h, mask, Nx, Ny, Δ)

    for j in 2:Ny-1, i in 2:Nx-1
        minus_∇ϕ₀_x[i, j] = -(ϕ₀[i+1, j] - ϕ₀[i-1, j]) / (2*Δ)
        minus_∇ϕ₀_y[i, j] = -(ϕ₀[i, j+1] - ϕ₀[i, j-1]) / (2*Δ)
    end

    longcoupwater = 5.0 
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

function iterative_update_Q!(Q, Nx, Ny, mask, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, ϕ₀)

    Q .= 0.0
    indices = Tuple{Int, Int}[]
    for j in 1:Ny, i in 1:Nx
        if mask[i, j] == 1
            Q[i, j] = ṁ_over_ρ_w[i, j] * Δ²
            push!(indices, (i, j))
        end
    end

    sort!(indices, by = idx -> ϕ₀[idx[1], idx[2]], rev = true)
    
    for (i, j) in indices

        weight_sum = 0.0
        num_downstream_neighbours = 0
        downstream_neighbours_list = []

        for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))

            ni, nj = i + di, j + dj

            if (1 <= ni <= Nx) && (1 <= nj <= Ny) && mask[ni, nj] == 1

                if (di < 0) || (dj < 0) # according to d-grid if the neighbour is north or east we need the gradient at i, j otherwise at ni, nj
                    w = (minus_∇ϕ₀_x[i, j]*di + minus_∇ϕ₀_y[i, j]*dj) / abs_∇ϕ₀[i, j]
                else
                    w = (minus_∇ϕ₀_x[ni, nj]*di + minus_∇ϕ₀_y[ni, nj]*dj) / abs_∇ϕ₀[ni, nj]
                end

                w = (minus_∇ϕ₀_x[ni, nj]*di + minus_∇ϕ₀_y[ni, nj]*dj) / abs_∇ϕ₀[ni, nj]

                w = (minus_∇ϕ₀_x[i, j]*di + minus_∇ϕ₀_y[i, j]*dj) / abs_∇ϕ₀[i, j]

                if w > 0
                    Q[ni, nj] += Q[i, j] * w
                    weight_sum += w
                    num_downstream_neighbours += 1
                    push!(downstream_neighbours_list, (di, dj))
                end

            end
        end

        # println(weight_sum, " $num_downstream_neighbours $i $j $downstream_neighbours_list")
        
        # if !isapprox(tmp_sum, 1)
        #     l = [mask[i+1, j], mask[i-1, j], mask[i, j+1], mask[i, j-1]]
        #     if l == [1.0, 1.0, 1.0, 1.0]
        #         println(l)
        #     end
        # end 
        
        # if num_downstream_neighbours != 2
        #     l = [mask[i+1, j], mask[i-1, j], mask[i, j+1], mask[i, j-1]]
        #     if l == [1.0, 1.0, 1.0, 1.0]
        #         println(l)
        #     end
        #     # println(l)
        # end 

        # if num_downstream_neighbours == 2
        #     if !isapprox(weight_sum, 1)
        #         println(weight_sum)
        #     end
        # end

        # l = [mask[i+1, j], mask[i-1, j], mask[i, j+1], mask[i, j-1]]
        # if num_downstream_neighbours == 1 && sum(l) == 4.0
        #     println(weight_sum, " $num_downstream_neighbours $i $j $downstream_neighbours_list")
        # end

    end
end

function update_potential!(ϕ₀, h, b, ρ_w, ρ_i, g)
    @. ϕ₀ = ρ_i * g * h + ρ_w * g * b
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

function update_H!(H, H_hard, H_soft, S_inf, H_0, F_till, Q, Q_c, κ)
    @. H_hard = sqrt(S_inf)
    @. H_soft = H_0 + (sqrt(S_inf)/F_till - H_0) * exp(-Q/Q_c)
    @. H = (1 - κ) * H_hard + κ * H_soft
end

function update_S_inf!(S_inf, K, α, β, abs_∇ϕ₀, Q)
    @. S_inf = K^(-1/α) * abs_∇ϕ₀^((1-β)/α) * Q^(1/α) 
end

function update_N_inf!(N_inf, H, S_inf, ρ_i, L_w, abs_v_b, h_b, Q, abs_∇ϕ₀, n, A)
    @. N_inf = ((H/S_inf)^2*((ρ_i*L_w*abs_v_b*h_b + Q*abs_∇ϕ₀)/(2.0*n^(-n)*ρ_i*L_w*A)))^(1.0/n)
end

function update_N!(N, N_inf, ϕ₀)
    @. N = erf(sqrt(π)*ϕ₀/(2*max(N_inf, 1e-14))) * N_inf
end

function main()

    data_path = "/Users/takisangelides/Documents/Post-doc/Kazmierczak et al 2024 Code/results/thwaites_hydro/output/THWAITES2km_m3_HARD_toto.mat"
    data = matread(data_path)

    Nx, Ny = size(data["H"])
    Δ = 2000.0
    Δ² = Δ^2
    mask = data["MASKo"]
    ṁ_over_ρ_w = data["Bmelt"] ./ 3.154e7 
    h = data["H"]
    b = data["B"]
    abs_v_b = data["ub"]
    
    ρ_w, ρ_i, g = 1000.0, 917.0, 9.81
    n, L_w, h_b = 3.0, 3.35e5, 0.1
    A = data["A"]
    α, β, f = 5/4, 3/2, 0.1
    K = (2/π)^(0.25) * sqrt((π + 2) / (ρ_w * f))
    F_till, Q_c, H_0 = 1.1, 1.0, 0.1

    ϕ₀ = zeros(Nx, Ny)
    minus_∇ϕ₀_x = zeros(Nx, Ny)
    minus_∇ϕ₀_y = zeros(Nx, Ny)
    abs_∇ϕ₀ = zeros(Nx, Ny)
    Q = zeros(Nx, Ny)
    S_inf = zeros(Nx, Ny)
    H_hard = zeros(Nx, Ny)
    H_soft = zeros(Nx, Ny)
    H = zeros(Nx, Ny)
    N_inf = zeros(Nx, Ny)
    N = zeros(Nx, Ny)
    κ = initialize_κ(b, Nx, Ny)

    update_potential!(ϕ₀, h, b, ρ_w, ρ_i, g)
    
    # potential_filling!(ϕ₀, Nx, Ny, 10) 

    update_potential_gradients_smoothed!(minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, ϕ₀, h, mask, Nx, Ny, Δ)

    # iterative_update_Q!(Q, Nx, Ny, mask, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, ϕ₀)
    recursive_update_Q!(Q, Nx, Ny, mask, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀)

    numer = @. abs(minus_∇ϕ₀_x) + abs(minus_∇ϕ₀_y)
    denom = @. sqrt(minus_∇ϕ₀_x^2 + minus_∇ϕ₀_y^2)
    corfac = @. numer / denom

    @. Q = min(max(Q / corfac / Δ, 0.0), 1e5)

    update_S_inf!(S_inf, K, α, β, abs_∇ϕ₀, Q)
    update_H!(H, H_hard, H_soft, S_inf, H_0, F_till, Q, Q_c, κ)
    update_N_inf!(N_inf, H, S_inf, ρ_i, L_w, abs_v_b, h_b, Q, abs_∇ϕ₀, n, A)
    update_N!(N, N_inf, ϕ₀)

    [x[mask .!= 1] .= NaN for x in (b, h, Q, N, κ, abs_v_b)]

    # fig = Figure(size = (1000, 800))

    # ax1 = Axis(fig[1, 1], title="Bed Elevation b", aspect = DataAspect(), xgridvisible = false, ygridvisible = false)
    # hm1 = Makie.heatmap!(ax1, b')
    # Colorbar(fig[1, 2], hm1, label="b [m]")

    # ax2 = Axis(fig[1, 3], title="Ice Thickness h", aspect = DataAspect(), xgridvisible = false, ygridvisible = false)
    # hm2 = Makie.heatmap!(ax2, h')
    # Colorbar(fig[1, 4], hm2, label="h [m]")

    # ax3 = Axis(fig[1, 5], title="Absolute basal velocity", aspect = DataAspect(), xgridvisible = false, ygridvisible = false)
    # hm3 = Makie.heatmap!(ax3, abs_v_b')
    # Colorbar(fig[1, 6], hm3, label="v_b [m/s]")

    # ax4 = Axis(fig[2, 1], title="Bed rheology κ", aspect = DataAspect(), xgridvisible = false, ygridvisible = false)
    # hm4 = Makie.heatmap!(ax4, κ')
    # Colorbar(fig[2, 2], hm4, label="κ")

    # ax5 = Axis(fig[2, 3], title="Effective Pressure N", aspect = DataAspect(), xgridvisible = false, ygridvisible = false)
    # hm5 = Makie.heatmap!(ax5, N'./10^6)
    # Colorbar(fig[2, 4], hm5, label="N [MPa]")

    # ax6 = Axis(fig[2, 5], title="Log scale\nScalar Water Flux Q", aspect = DataAspect(), xgridvisible = false, ygridvisible = false)
    # hm6 = Makie.heatmap!(ax6, Q' .* 3.154e7 ./ 1e4, colormap = Reverse(:RdBu), colorrange = (0, 10))
    # Colorbar(fig[2, 6], hm6, label="Q [m³/a]")

    # colgap!(fig.layout, 15)
    # fig

    display(Makie.heatmap(Q' .* 3.154e7 ./ 1e4, colormap = Reverse(:RdBu), colorrange = (0, 10)))

    return fig
end

main()