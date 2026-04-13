"""
$(TYPEDSIGNATURES)

A struct containing all the components needed to run the subglacial hydrology model.

# Fields
$(TYPEDFIELDS)
"""
struct HydrologyModel{T<:AbstractFloat, G, F}
    c::PhysicalConstants{T} # physical constants
    grid::G # Oceananigans rectilinear grid
    fields::F # named tuple of Oceananigans fields
end


"""
$(TYPEDSIGNATURES)

Constructor to create a [`HydrologyModel`](@ref) from a data file. The `loader` function
must return `Nx, Ny, xlims, ylims, mask, h, b, abs_v_b, A, ṁ_over_ρ_w`. For example:

```julia
model = HydrologyModel("path/to/file.mat", load_kazmierczak2024)
```
"""
function HydrologyModel(path::String, loader::Function, c::PhysicalConstants{T}) where {T <: AbstractFloat}
    Nx, Ny, xlims, ylims, mask, h, b, abs_v_b, A, ṁ_over_ρ_w, bed_rheology = loader(path)
    grid = initialize_grid(Nx, Ny, xlims, ylims; T = T)
    fields = initialize_fields(grid)
    set_initial_fields!(fields, T.(mask), T.(h), T.(b), T.(abs_v_b), T.(A), T.(ṁ_over_ρ_w), bed_rheology)
    return HydrologyModel(c, grid, fields)
end


"""
$(TYPEDSIGNATURES)

Initialize an Oceananigans rectilinear grid with float type `T`.
"""
function initialize_grid(Nx, Ny, xlims, ylims; T = Float64, topology = (Bounded, Bounded, Flat), halo = (1, 1))
    return RectilinearGrid(T, topology = topology, size = (Nx, Ny), x = xlims, y = ylims, halo = halo)
end


"""
$(TYPEDSIGNATURES)

Initialize all Oceananigans fields on `grid`. All fields are initialized to zero
and then set to their initial values by [`set_initial_fields!`](@ref).
"""
function initialize_fields(grid)
    fields = (
        h                    = CenterField(grid), # ice thickness [m]
        b                    = CenterField(grid), # bedrock elevation [m]
        abs_v_b              = CenterField(grid), # absolute basal velocity [m/s]
        A                    = CenterField(grid), # Glen's flow law viscosity parameter
        mask                 = CenterField(grid), # grounded ice mask (1 = grounded)
        ṁ_over_ρ_w           = CenterField(grid), # basal melt rate / water density [m/s]
        ϕ₀                   = CenterField(grid), # geometric potential [Pa]
        ϕ₀_tmp               = CenterField(grid), # geometric potential placeholder for filling
        ψ_out                = CenterField(grid), # integrated scalar water flux [m^2/s]
        q                    = CenterField(grid), # distributed water flux [m^2/s]
        Q                    = CenterField(grid), # volumetric water flux per conduit [m^3/s]
        S_inf                = CenterField(grid), # far-field conduit cross-sectional area [m^2]
        H_hard               = CenterField(grid), # conduit thickness over hard bed [m]
        H_soft               = CenterField(grid), # conduit thickness over soft bed [m]
        H                    = CenterField(grid), # conduit thickness [m]
        N_inf                = CenterField(grid), # far-field effective pressure [Pa]
        N                    = CenterField(grid), # effective pressure [Pa]
        Po                   = CenterField(grid), # hydrostatic ice overburden pressure [Pa]
        abs_∇ϕ₀              = CenterField(grid), # |∇ϕ₀| [Pa/m]
        κ                    = CenterField(grid), # bed heterogeneity indicator [0=hard, 1=soft]
        corfac               = CenterField(grid), # correction factor from ψ_out to q
        minus_∇ϕ₀_x          = CenterField(grid), # -∂ϕ₀/∂x [Pa/m]
        minus_∇ϕ₀_y          = CenterField(grid), # -∂ϕ₀/∂y [Pa/m]
        minus_∇ϕ₀_smoothed_x = CenterField(grid), # smoothed -∂ϕ₀/∂x [Pa/m]
        minus_∇ϕ₀_smoothed_y = CenterField(grid), # smoothed -∂ϕ₀/∂y [Pa/m]
        abs_∇ϕ₀_smoothed     = CenterField(grid), # |smoothed ∇ϕ₀| [Pa/m]
    )
    return fields
end


"""
$(TYPEDSIGNATURES)

Set the initial values of all fields in `fields` from the input data arrays.
All fields not in the input data are initialized to `zero(T)`. The float type
`T` is inferred from `T = eltype(h)`, where h is the ice thickness array.
"""
function set_initial_fields!(fields, mask, h, b, abs_v_b, A, ṁ_over_ρ_w, bed_rheology)

    T = eltype(h)

    for name in keys(fields)
        if name ∉ (:h, :b, :abs_v_b, :A, :mask, :ṁ_over_ρ_w)
            fields[name] .= zero(T)
        end
    end

    fields.h .= h
    fields.b .= b
    fields.abs_v_b .= abs_v_b
    fields.A .= A
    fields.mask .= mask
    fields.ṁ_over_ρ_w .= ṁ_over_ρ_w
    initialize_κ!(fields, bed_rheology)
    
    return nothing
end


"""
$(TYPEDSIGNATURES)

Initialize the bed hardness indicator `κ` in `fields` based on `bed_rheology`:
- `:hard` → `κ = 0` everywhere
- `:soft` → `κ = 1` everywhere  
- `:mixed` → `κ = 1` where bedrock elevation `b < -1000 m`, `0` elsewhere

Throws an error for unrecognized `type`.
"""
function initialize_κ!(fields, bed_rheology::Symbol = :hard)
    if bed_rheology == :hard
        fields.κ .= zero(eltype(fields.κ))
    elseif bed_rheology == :soft
        fields.κ .= one(eltype(fields.κ))
    elseif bed_rheology == :mixed
        fields.κ .= zero(eltype(fields.κ))
        fields.κ[fields.b .< -1000] .= one(eltype(fields.κ))
    else
        error("Unknown bed type: $type. Choose from :hard, :soft, :mixed.")
    end
    return nothing
end
