"""
    FastHydrology.jl

A numerical model for subglacial water dynamics and effective pressure evolution.
Based on the approach from Kazmierczak et al. 2024.
"""

using BenchmarkTools
using SpecialFunctions
using MAT
using ImageFiltering
using Statistics
using ColorSchemes
import ColorSchemes.Colors
using CairoMakie

# ==============================================================================
# Constants
# ==============================================================================

const WATER_DENSITY = 1000.0  # kg/m³
const ICE_DENSITY = 917.0     # kg/m³
const GRAVITY = 9.81          # m/s²
const MANNING_EXPONENT = 3.0
const LATENT_HEAT_WATER = 3.35e5  # J/kg
const BED_THICKNESS = 0.1     # m
const MANNING_COEFFICIENT_EXPONENT = (5/4)
const BED_FRICTION_EXPONENT = (3/2)
const FRICTION_FACTOR = 0.1
const TILL_FACTOR = 1.1
const CRITICAL_DISCHARGE = 1.0
const INITIAL_CAVITY_HEIGHT = 0.1
const COUPLING_LENGTH = 1e4   # m
const LONG_COUPLING_WATER = 5.0
const SECONDS_PER_YEAR = 3.154e7
const MIN_PRESSURE_FRACTION = 0.02

# ==============================================================================
# Water Flux Calculation
# ==============================================================================

"""
    accumulate_ψ_out!(ψ_out, i, j, ṁ_over_ρ_w, Δ², ∇ϕ, abs_∇ϕ, Nx, Ny)

Recursively accumulate water potential outflow from a cell using melt rate and
pressure gradients. Cells are marked with -1.0 before visiting; positive values
indicate computed cells (memoization).

Arguments:
- `ψ_out`: Output flux array (updated in-place)
- `i, j`: Current cell indices
- `ṁ_over_ρ_w`: Melt rate per water density
- `Δ²`: Grid spacing squared
- `∇ϕ`: Tuple of (minus_∇ϕ₀_x, minus_∇ϕ₀_y)
- `abs_∇ϕ`: Absolute potential gradient magnitude
- `Nx, Ny`: Grid dimensions
"""
function accumulate_ψ_out!(ψ_out, i, j, ṁ_over_ρ_w, Δ², ∇ϕ, abs_∇ϕ, Nx, Ny)
    # Already visited
    ψ_out[i, j] ≥ 0.0 && return ψ_out[i, j]

    # Initialize with melt source
    ψ_out[i, j] = max(0.0, ṁ_over_ρ_w[i, j] * Δ²)

    # Sum contributions from neighboring cells
    for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))
        ni, nj = i + di, j + dj

        # Check bounds
        (ni < 1 || ni > Nx || nj < 1 || nj > Ny) && continue

        # Calculate weight based on pressure gradient direction
        minus_∇ϕ₀_x, minus_∇ϕ₀_y = ∇ϕ
        
        # Safety check for division by zero
        abs_grad = abs_∇ϕ[ni, nj]
        if abs_grad ≤ 1e-12
            continue
        end
        
        w = -(minus_∇ϕ₀_x[ni, nj] * di + minus_∇ϕ₀_y[ni, nj] * dj) / abs_grad

        # Positive weight means flow from neighbor into current cell
        if w > 0
            ψ_out[i, j] += accumulate_ψ_out!(ψ_out, ni, nj, ṁ_over_ρ_w, Δ², ∇ϕ, abs_∇ϕ, Nx, Ny) * w
        end
    end

    return ψ_out[i, j]
end

"""
    recursive_update_ψ_out!(ψ_out, ṁ_over_ρ_w, Δ², ∇ϕ, abs_∇ϕ, mask, Nx, Ny)

Initialize and compute water potential outflow for all active cells.
"""
function recursive_update_ψ_out!(ψ_out, ṁ_over_ρ_w, Δ², ∇ϕ, abs_∇ϕ, mask, Nx, Ny)
    # Mark all active cells as unvisited
    ψ_out[mask .== 1.0] .= -1.0

    # Compute outflow from all active cells
    for j in 1:Ny, i in 1:Nx
        if mask[i, j] == 1.0
            accumulate_ψ_out!(ψ_out, i, j, ṁ_over_ρ_w, Δ², ∇ϕ, abs_∇ϕ, Nx, Ny)
        end
    end
end

# ==============================================================================
# Potential and Topography
# ==============================================================================

"""
    potential_filling!(ϕ₀, Nx, Ny; iterations=10)

Apply iterative hollow-filling to remove spurious sinks in topography.
Uses Laplacian smoothing to raise local minima.

Arguments:
- `ϕ₀`: Potential array (modified in-place)
- `Nx, Ny`: Grid dimensions
- `iterations`: Number of smoothing iterations (default: 10)
"""
function potential_filling!(ϕ₀, Nx, Ny; iterations=10)
    pot_next = copy(ϕ₀)

    for _ in 1:iterations
        for j in 2:Ny-1, i in 2:Nx-1
            p = ϕ₀[i, j]
            p_neighbors = (ϕ₀[i+1, j], ϕ₀[i-1, j], ϕ₀[i, j+1], ϕ₀[i, j-1])

            # Fill local minima with neighbor average
            if all(p .< p_neighbors)
                pot_next[i, j] = mean(p_neighbors)
            end
        end
        ϕ₀ .= pot_next
    end
end

"""
    update_potential!(ϕ₀, h, b, ρ_w, ρ_i, g)

Compute hydraulic potential from ice thickness and bed topography.
ϕ₀ = ρᵢ·g·h + ρw·g·b
"""
function update_potential!(ϕ₀, h, b, ρ_w, ρ_i, g)
    @. ϕ₀ = ρ_i * g * h + ρ_w * g * b
end

"""
    update_potential_gradients_smoothed!(minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀, 
                                        ϕ₀, h, Nx, Ny, Δ, mask)

Compute smoothed pressure gradients using a cone-shaped kernel weighted by
ice thickness. Smoothing radius is proportional to mean ice thickness.
"""
function update_potential_gradients_smoothed!(minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀,
                                             ϕ₀, h, Nx, Ny, Δ, mask)
    # Compute unsmoothed gradients
    for j in 1:Ny, i in 2:Nx-1
        minus_∇ϕ₀_x[i, j] = -(ϕ₀[i+1, j] - ϕ₀[i-1, j]) / (2 * Δ)
    end
    minus_∇ϕ₀_x[1, :] .= minus_∇ϕ₀_x[2, :]
    minus_∇ϕ₀_x[end, :] .= minus_∇ϕ₀_x[end-1, :]

    for j in 2:Ny-1, i in 1:Nx
        minus_∇ϕ₀_y[i, j] = -(ϕ₀[i, j+1] - ϕ₀[i, j-1]) / (2 * Δ)
    end
    minus_∇ϕ₀_y[:, 1] .= minus_∇ϕ₀_y[:, 2]
    minus_∇ϕ₀_y[:, end] .= minus_∇ϕ₀_y[:, end-1]

    # Determine smoothing kernel based on ice thickness
    h_avg = max(mean(h[mask .== 1.0]), 10.0)
    scale = h_avg * LONG_COUPLING_WATER * 2.0
    width = 2.0 * scale
    width ≤ Δ && (scale = Δ / 2.0 + 1.0)

    kernel_size = 2 * round(Int, width / Δ - 0.5) + 1
    frb = Int((kernel_size - 1) / 2)

    # Create cone kernel: W(d) = max(0, 1 - d/(2*scale))
    kernel = zeros(kernel_size, kernel_size)
    for nj in 1:kernel_size, ni in 1:kernel_size
        dist = sqrt((Δ * (ni - frb - 1))^2 + (Δ * (nj - frb - 1))^2) / scale
        kernel[ni, nj] = max(0.0, 1.0 - dist / 2.0)
    end
    kernel ./= sum(kernel)  # Normalize to preserve total potential

    # Apply smoothing filter
    minus_∇ϕ₀_x .= imfilter(minus_∇ϕ₀_x, centered(kernel), "reflect")
    minus_∇ϕ₀_y .= imfilter(minus_∇ϕ₀_y, centered(kernel), "reflect")

    # Magnitude of gradient
    @. abs_∇ϕ₀ = abs(minus_∇ϕ₀_x) + abs(minus_∇ϕ₀_y)
end

# ==============================================================================
# Effective Pressure and Cavity Height
# ==============================================================================

"""
    initialize_κ(b, Nx, Ny; type="hard")

Initialize substrate type indicator κ.
- "hard": κ = 0 everywhere (rigid bed)
- "soft": κ = 1 everywhere (deformable bed)
- "mixed": κ = 1 for deep areas (b < -1000 m), 0 elsewhere
"""
function initialize_κ(b, Nx, Ny; type="hard")
    if type == "hard"
        return zeros(Nx, Ny)
    elseif type == "soft"
        return ones(Nx, Ny)
    else  # "mixed"
        κ = zeros(Nx, Ny)
        κ[b .< -1000] .= 1.0
        return κ
    end
end

"""
    update_Po!(Po, ρ_i, g, h)

Compute overburden pressure from ice thickness.
Po = max(ρᵢ·g·h, 1e5 Pa)
"""
function update_Po!(Po, ρ_i, g, h)
    @. Po = max(ρ_i * g * h, 1e5)
end

"""
    update_S_inf!(S_inf, K, α, β, abs_∇ϕ₀, Q)

Compute cavity opening rate (width * height) at equilibrium.
S∞ = K^(-1/α) · (∇ϕ)^((1-β)/α) · Q^(1/α)
"""
function update_S_inf!(S_inf, K, α, β, abs_∇ϕ₀, Q)
    @. S_inf = K^(-1 / α) * abs_∇ϕ₀^((1 - β) / α) * Q^(1 / α)
    # Ensure S_inf is non-negative and has minimum value for stability
    @. S_inf = max(S_inf, 1e-12)
end

"""
    update_H!(H, H_hard, H_soft, S_inf, H_0, F_till, Q, Q_c, κ)

Compute cavity height as weighted blend of hard-bed and soft-bed solutions.
"""
function update_H!(H, H_hard, H_soft, S_inf, H_0, F_till, Q, Q_c, κ)
    @. H_hard = sqrt(S_inf)
    @. H_soft = H_0 + (sqrt(S_inf) / F_till - H_0) * exp(-Q / Q_c)
    @. H = (1 - κ) * H_hard + κ * H_soft
    # Ensure H is non-negative
    @. H = max(H, 1e-12)
end

"""
    update_N_inf!(N_inf, H, S_inf, ρ_i, L_w, abs_v_b, h_b, Q, abs_∇ϕ₀, n, A, Po)

Compute equilibrium effective pressure using basal velocity and meltwater flux.
Bounded by [0.02·Po, Po].
"""
function update_N_inf!(N_inf, H, S_inf, ρ_i, L_w, abs_v_b, h_b, Q, abs_∇ϕ₀, n, A, Po)
    numerator = ρ_i * L_w * abs_v_b * h_b .+ Q .* abs_∇ϕ₀
    denominator = 2.0 * n^(-n) * ρ_i * L_w * A
    
    # Avoid division by zero: ensure S_inf > 0
    cavity_ratio = @. H / max(S_inf, 1e-12)

    @. N_inf = ((cavity_ratio^2 * numerator / denominator)^(1.0 / n))
    @. N_inf = max(min(N_inf, Po), MIN_PRESSURE_FRACTION * Po)
end

"""
    update_N!(N, N_inf, ϕ₀)

Compute effective pressure from equilibrium pressure and potential.
Uses erf(x) for numerical stability.
N = erf(√π·ϕ₀/(2·N∞)) · N∞
"""
function update_N!(N, N_inf, ϕ₀)
    @. N = max(0.0, erf(sqrt(π) * ϕ₀ / (2 * N_inf)) * N_inf)
end

# ==============================================================================
# Main Simulation
# ==============================================================================

"""
    main(data_path; visualize=true)

Run subglacial hydrology model on Thwaites Glacier data.
Computes water flux and effective pressure.
"""
function main(data_path::String; visualize=true)
    # Load data
    data = matread(data_path)
    Nx, Ny = size(data["H"])
    Δ = 2000.0
    Δ² = Δ^2

    # Extract fields
    mask = data["MASKo"]
    ṁ_over_ρ_w = data["Bmelt"] ./ SECONDS_PER_YEAR
    h = copy(data["H"])
    b = data["B"]
    abs_v_b = data["ub"]
    A = data["A"]

    # Manning-Strickler coefficient
    K = (2 / π)^0.25 * sqrt((π + 2) / (WATER_DENSITY * FRICTION_FACTOR))

    # Initialize fields
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
    Po = zeros(Nx, Ny)
    κ = initialize_κ(b, Nx, Ny)

    # ========== Compute Distributed Water Flux ==========
    update_potential!(ϕ₀, h, b, WATER_DENSITY, ICE_DENSITY, GRAVITY)
    potential_filling!(ϕ₀, Nx, Ny; iterations=10)
    @. h = (ϕ₀ - WATER_DENSITY * GRAVITY * b) / (ICE_DENSITY * GRAVITY)

    update_potential_gradients_smoothed!(minus_∇ϕ₀_x, minus_∇ϕ₀_y, abs_∇ϕ₀,
                                         ϕ₀, h, Nx, Ny, Δ, mask)

    ∇ϕ = (minus_∇ϕ₀_x, minus_∇ϕ₀_y)
    recursive_update_ψ_out!(ψ_out, ṁ_over_ρ_w, Δ², ∇ϕ, abs_∇ϕ₀, mask, Nx, Ny)

    # Correct for coordinate system - avoid division by zero
    corfac = zeros(Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        grad_mag = sqrt(minus_∇ϕ₀_x[i, j]^2 + minus_∇ϕ₀_y[i, j]^2)
        if grad_mag > 1e-12
            corfac[i, j] = abs_∇ϕ₀[i, j] / grad_mag
        else
            corfac[i, j] = 1.0
        end
    end
    @. q = min(max(ψ_out / (corfac * Δ), 0.0), 1e5)

    # ========== Compute Effective Pressure ==========
    @. ϕ₀ = ICE_DENSITY * GRAVITY * h + WATER_DENSITY * GRAVITY * b

    # Recompute gradients without smoothing
    for j in 1:Ny, i in 2:Nx-1
        minus_∇ϕ₀_x[i, j] = -(ϕ₀[i+1, j] - ϕ₀[i-1, j]) / (2 * Δ)
    end
    minus_∇ϕ₀_x[1, :] .= minus_∇ϕ₀_x[2, :]
    minus_∇ϕ₀_x[end, :] .= minus_∇ϕ₀_x[end-1, :]

    for j in 2:Ny-1, i in 1:Nx
        minus_∇ϕ₀_y[i, j] = -(ϕ₀[i, j+1] - ϕ₀[i, j-1]) / (2 * Δ)
    end
    minus_∇ϕ₀_y[:, 1] .= minus_∇ϕ₀_y[:, 2]
    minus_∇ϕ₀_y[:, end] .= minus_∇ϕ₀_y[:, end-1]

    @. abs_∇ϕ₀ = sqrt(minus_∇ϕ₀_x^2 + minus_∇ϕ₀_y^2)

    # Compute effective pressure
    @. Q = q * COUPLING_LENGTH
    update_S_inf!(S_inf, K, MANNING_COEFFICIENT_EXPONENT, BED_FRICTION_EXPONENT, abs_∇ϕ₀, Q)
    update_H!(H, H_hard, H_soft, S_inf, INITIAL_CAVITY_HEIGHT, TILL_FACTOR, Q, CRITICAL_DISCHARGE, κ)
    update_Po!(Po, ICE_DENSITY, GRAVITY, h)
    update_N_inf!(N_inf, H, S_inf, ICE_DENSITY, LATENT_HEAT_WATER, abs_v_b, BED_THICKNESS,
                  Q, abs_∇ϕ₀, MANNING_EXPONENT, A, Po)
    update_N!(N, N_inf, ϕ₀)

    # Convert to annual flux (m²/year)
    @. q *= 60^2 * 24 * 365.25 / 1e4
    q[mask .!= 1.0] .= NaN

    # Visualize
    if visualize
        display(heatmap(q', colormap=Reverse(:RdBu), colorrange=(0, 10)))
    end

    return (; q, N, h, mask)
end

# ==============================================================================
# Entry point
# ==============================================================================

data_path = "/Users/taange001/Documents/Coding/FastHydrology.jl/local_experiments/input/Kazmierczak2024/THWAITES2km_m3_HAB_toto.mat"
main(data_path)