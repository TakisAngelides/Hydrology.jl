using BenchmarkTools
using SpecialFunctions
using MAT
using ImageFiltering
using Statistics
using ColorSchemes
import ColorSchemes.Colors
using CairoMakie
# using DSP
# using WGLMakie

function accumulate_ψ_out!(ψ_out, i, j, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, Nx, Ny) 

    if ψ_out[i, j] ≥ 0.0 # we initialize all unvisited cells with -1.0 and q is positive semi-definite so if its positive then it has been visited and ψ_out calculated
        return ψ_out[i, j]
    end

    ψ_out[i, j] = ṁ_over_ρ_w[i, j] * Δ²

    for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))

            ni, nj = i+di, j+dj 

            if (ni < 1) || (ni > Nx) || (nj < 1) || (nj > Ny) 
                continue
            end

            # TODO: Choose the correct calculation of weight
            
            # if (di > 0) || (dj > 0) 
            #     w = -(minus_∇ϕ₀_x[i, j]*di + minus_∇ϕ₀_y[i, j]*dj) / abs_∇ϕ₀[i, j]
            # else
            #     w = -(minus_∇ϕ₀_x[ni, nj]*di + minus_∇ϕ₀_y[ni, nj]*dj) / abs_∇ϕ₀[ni, nj]
            # end
            
            # w = -(minus_∇ϕ₀_x[i, j]*di + minus_∇ϕ₀_y[i, j]*dj) / abs_∇ϕ₀[i, j]
            
            w = -(minus_∇ϕ₀_x[ni, nj]*di + minus_∇ϕ₀_y[ni, nj]*dj) / abs_∇ϕ₀[ni, nj]

            if w > 0
                ψ_out[i, j] += accumulate_ψ_out!(ψ_out, ni, nj, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, Nx, Ny) * w
            end
    end
    return ψ_out[i, j]
end

function recursive_update_ψ_out!(ψ_out, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, active_indices, Nx, Ny)
    
    for I in active_indices
        ψ_out[I] = -1.0
    end

    for I in active_indices
        i = I[1]
        j = I[2]
        accumulate_ψ_out!(ψ_out, i, j, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, Nx, Ny)
    end

end

"""
    potential_filling!(ϕ₀, Nx, Ny, iterations=10) -> Return type

A hollow-filling algorithm must be applied to the topography to ensure 
there are no sinks for the flux. Sinks are likely to be spurious features in the data, 
due to processing and hence can be filled. Hollow filling is carried out iteratively 
by assigning the ‘‘sink’’ cell the average of its four neighbours.

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
function update_potential_gradients_smoothed!(minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, ϕ₀, h, Nx, Ny, Δ, active_indices)

    # TODO: Make efficient once established correct
    # TODO: what boundary conditions should we apply? I think in your code you use periodic but that the end-points are anyway always kept to zero?
    for j in 1:Ny, i in 1:Nx-1
        minus_∇ϕ₀_x[i, j] = -(ϕ₀[i+1, j] - ϕ₀[i, j]) / (Δ)
    end
    minus_∇ϕ₀_x[end, :] = minus_∇ϕ₀_x[end-1, :] # the edge point gets the same gradient as the penultimate point - doesn't change results that much
    
    for j in 1:Ny-1, i in 1:Nx
        minus_∇ϕ₀_y[i, j] = -(ϕ₀[i, j+1] - ϕ₀[i, j]) / (Δ)
    end
    minus_∇ϕ₀_y[:, end] = minus_∇ϕ₀_y[:, end-1] 

    longcoupwater = 5.0
    h_avg = max(mean(view(h, active_indices)), 10.0)
    
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

function update_Po!(Po, ρ_i, g, h)
    @. Po = max(ρ_i * g * h, 1e5)
end

function update_N_inf!(N_inf, H, S_inf, ρ_i, L_w, abs_v_b, h_b, Q, abs_∇ϕ₀, n, A, Po)
    # TODO: why do we set these min and max limits?
    @. N_inf = max(min(((H/S_inf)^2*((ρ_i*L_w*abs_v_b*h_b + Q*abs_∇ϕ₀)/(2.0*n^(-n)*ρ_i*L_w*A)))^(1.0/n), Po), 0.02*Po)
end

"""
    update_N!(N, N_inf, ϕ₀) -> Return type

1−erfc(x)=erf(x) this is how they have it in the Kazmierczak et al 2024 code, probably for numerical stability

# Arguments

- `N`: Argument description
- `N_inf`: Argument description
- `ϕ₀`: Argument description
"""
function update_N!(N, N_inf, ϕ₀)
    @. N = (1 - erfc(sqrt(π)*ϕ₀/(2*N_inf))) * N_inf 
end

function plot_3d_flux(b, q)
    
    fig = Figure(size = (3000, 3000))
    
    ax = Axis3(fig[1, 1], 
        title = "Subglacial Water Flux over Bed Topography",
        aspect = (1, 1, 0.2),
        perspectiveness = 0.5,
        elevation = 0.5,
        azimuth = 1.4
    )

    hidedecorations!(ax)
    hidespines!(ax)
    ax.yzpanelvisible = false
    ax.xzpanelvisible = false
    ax.xypanelvisible = false

    plt = surface!(ax, b', 
        color = q' .* 3.154e7 ./ 1e4, 
        colormap = Reverse(:RdBu), 
        shading = true, 
        nan_color = :transparent,
    )

    Colorbar(fig[1, 2], plt, label = "q_w [m²/a]")
    
    return fig
end

function main()

    data_path = "/Users/takisangelides/Documents/Post-doc/Kazmierczak et al 2024 Code/results/thwaites_hydro/output/THWAITES2km_m3_HARD_toto.mat"
    data = matread(data_path)

    Nx, Ny = size(data["H"])
    Δ = 2000.0
    Δ² = Δ^2
    mask = data["MASKo"]
    active_indices = findall(mask .== 1)
    inactive_indices = findall(mask .!= 1)
    ṁ_over_ρ_w = data["Bmelt"] ./ 3.154e7
    h = data["H"]
    b = data["B"]
    abs_v_b = data["ub"]
    
    ρ_w, ρ_i, g = 1000.0, 917.0, 9.81
    n, L_w, h_b = 3.0, 3.35e5, 0.1
    A = data["A"]
    α, β, f = 5/4, 3/2, 0.1
    K = (2/π)^(0.25) * sqrt((π + 2) / (ρ_w * f))
    F_till, Q_c, H_0, l_c = 1.1, 1.0, 0.1, 1e4

    ϕ₀ = zeros(Nx, Ny)
    minus_∇ϕ₀_x = zeros(Nx, Ny)
    minus_∇ϕ₀_y = zeros(Nx, Ny)
    abs_∇ϕ₀ = zeros(Nx, Ny)
    ψ_out = zeros(Nx, Ny)
    q = zeros(Nx, Ny)
    Q = zeros(Nx, Ny)
    S_inf = zeros(Nx, Ny)
    H_hard = zeros(Nx, Ny)
    H_soft = zeros(Nx, Ny)
    H = zeros(Nx, Ny)
    N_inf = zeros(Nx, Ny)
    N = zeros(Nx, Ny)
    κ = initialize_κ(b, Nx, Ny)
    Po = zeros(Nx, Ny)

    # Distributed water flux

    update_potential!(ϕ₀, h, b, ρ_w, ρ_i, g)
    
    potential_filling!(ϕ₀, Nx, Ny, 10) 

    @. h = (ϕ₀ - ρ_w*g*b)/(ρ_i*g) # correction to h from potential filling

    update_potential_gradients_smoothed!(minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, ϕ₀, h, Nx, Ny, Δ, active_indices)

    recursive_update_ψ_out!(ψ_out, ṁ_over_ρ_w, Δ², minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, active_indices, Nx, Ny)

    corfac = @. abs_∇ϕ₀ / sqrt(minus_∇ϕ₀_x^2 + minus_∇ϕ₀_y^2)
    # Distributed water flux - dividing by Δ is heuristic to go from scalar volumetric flux out of the cell ψ_out
    # to flux density q with units m²/s and its coming from https://doi.org/10.3189/S0260305500013215 W. F. Budd and R. C. Warner 1996
    @. q = min(max(ψ_out / (corfac * Δ), 0), 1e5)

    # Effective pressure

    @. ϕ₀ = ρ_i * g * h + ρ_w * g * b
    for j in 1:Ny, i in 1:Nx-1
        minus_∇ϕ₀_x[i, j] = -(ϕ₀[i+1, j] - ϕ₀[i, j]) / (Δ)
    end
    minus_∇ϕ₀_x[end, :] = minus_∇ϕ₀_x[end-1, :] # the edge point gets the same gradient as the penultimate point - doesn't change results that much
    for j in 1:Ny-1, i in 1:Nx
        minus_∇ϕ₀_y[i, j] = -(ϕ₀[i, j+1] - ϕ₀[i, j]) / (Δ)
    end
    minus_∇ϕ₀_y[:, end] = minus_∇ϕ₀_y[:, end-1] 
    @. abs_∇ϕ₀ = sqrt(minus_∇ϕ₀_x^2 + minus_∇ϕ₀_y^2)

    @. Q = q * l_c # volumetric water flux in a single channel
    update_S_inf!(S_inf, K, α, β, abs_∇ϕ₀, Q)
    update_H!(H, H_hard, H_soft, S_inf, H_0, F_till, Q, Q_c, κ)
    update_Po!(Po, ρ_i, g, h)    
    update_N_inf!(N_inf, H, S_inf, ρ_i, L_w, abs_v_b, h_b, Q, abs_∇ϕ₀, n, A, Po)
    update_N!(N, N_inf, ϕ₀)

    # Plotting

    for I in inactive_indices
        b[I] = h[I] = Q[I] = N[I] = κ[I] = abs_v_b[I] = q[I] = abs_∇ϕ₀[I] = NaN
    end

    parula_custom = ColorScheme([Colors.RGB(62/255, 39/255, 169/255), Colors.RGB(53/255, 121/255, 254/255), Colors.RGB(16/255, 191/255, 187/255), Colors.RGB(200/255, 194/255, 41/255), Colors.RGB(250/255, 252/255, 21/255)])
    light_gray = Colors.RGB(0.95, 0.95, 0.95)

    # WGLMakie.activate!()    
    # fig_3d = plot_3d_flux(b, q)
    # display(fig_3d)

    # CairoMakie.activate!()    
    # fig = Figure(size = (1000, 800))

    # ax1 = Axis(fig[1, 1], title = "Bed Elevation b", aspect = DataAspect(), xgridvisible = false, ygridvisible = false, backgroundcolor = light_gray)
    # hm1 = Makie.heatmap!(ax1, b'; colormap = :oleron, colorrange = (-2000, 2000))
    # Colorbar(fig[1, 2], hm1, label = "b [m]")

    # ax2 = Axis(fig[1, 3], title = "Ice Thickness h", aspect = DataAspect(), xgridvisible = false, ygridvisible = false, backgroundcolor = light_gray)
    # hm2 = Makie.heatmap!(ax2, h'; colormap = parula_custom, colorrange = (0, 3000))
    # Colorbar(fig[1, 4], hm2, label = "h [m]")

    # ax3 = Axis(fig[1, 5], title = "Absolute basal velocity", aspect = DataAspect(), xgridvisible = false, ygridvisible = false, backgroundcolor = light_gray)
    # hm3 = Makie.heatmap!(ax3, log10.(abs_v_b'); colormap = Reverse(:RdBu), colorrange = (0, 3))
    # Colorbar(fig[1, 6], hm3, label = "log₁₀(v_b) [m/s]")

    # ax4 = Axis(fig[2, 1], title = "Bed rheology κ", aspect = DataAspect(), xgridvisible = false, ygridvisible = false, backgroundcolor = light_gray)
    # hm4 = Makie.heatmap!(ax4, κ')
    # Colorbar(fig[2, 2], hm4, label = "κ")

    # ax5 = Axis(fig[2, 3], title = "Effective Pressure N", aspect = DataAspect(), xgridvisible = false, ygridvisible = false, backgroundcolor = light_gray)
    # hm5 = Makie.heatmap!(ax5, N' .* 1e-6; colormap = Reverse(:RdBu), colorrange = (0, 10))
    # Colorbar(fig[2, 4], hm5, label = "N [MPa]")

    # ax6 = Axis(fig[2, 5], title = "Distributed water flux", aspect = DataAspect(), xgridvisible = false, ygridvisible = false, backgroundcolor = light_gray)
    # hm6 = Makie.heatmap!(ax6, q' .* 3.154e7 ./ 1e4; colormap = Reverse(:RdBu)) #, colorrange = (0, 10))
    # Colorbar(fig[2, 6], hm6, label = "q_w [10⁴ m²/a]")

    # colgap!(fig.layout, 15)
    # display(fig)

    # CairoMakie.activate!()    
    # fig = Figure(size = (1000, 800))

    # ax5 = Axis(fig[1, 1], title = "Effective Pressure N", aspect = DataAspect(), xgridvisible = false, ygridvisible = false, backgroundcolor = light_gray)
    # hm5 = Makie.heatmap!(ax5, N' .* 1e-6; colormap = Reverse(:RdBu), colorrange = (0, 10))
    # Colorbar(fig[1, 2], hm5, label = "N [MPa]")

    # ax6 = Axis(fig[1, 3], title = "Distributed water flux", aspect = DataAspect(), xgridvisible = false, ygridvisible = false, backgroundcolor = light_gray)
    # hm6 = Makie.heatmap!(ax6, q' .* 3.154e7 ./ 1e4; colormap = Reverse(:RdBu)) #, colorrange = (0, 10))
    # Colorbar(fig[1, 4], hm6, label = "q_w [10⁴ m²/a]")

    # colgap!(fig.layout, 15)
    # display(fig)

    # CairoMakie.activate!()    
    fig = Figure(size = (900, 700))
    ax1 = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "q old code", aspect = DataAspect())
    hm1 = heatmap!(ax1, q' .* 3.154e7 ./ 1e4; colormap = Reverse(:RdBu), colorrange = (0, 10))
    Colorbar(fig[1, 2], hm1)

    # fig = Figure(size = (900, 700))
    # ax1 = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "N old code", aspect = DataAspect())
    # hm1 = heatmap!(ax1, N' ./ 1e6; colormap = Reverse(:RdBu), colorrange = (0, 10))
    # Colorbar(fig[1, 2], hm1)

    display(fig)

    return nothing
    
end

main()
