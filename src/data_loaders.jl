"""
$(TYPEDSIGNATURES)

Load a `.mat` file from the Kazmierczak et al. (2024) Thwaites hydrology output
and return the processed fields needed for the simulation.

# Arguments
- `path::String`: path to the `.mat` file.

# Returns
- `Nx, Ny`: grid dimensions.
- `xlims, ylims`: grid extent in meters (cell-centered, with half-cell padding to extend to grid edge faces).
- `mask`: grounded ice mask.
- `h`: ice thickness (m).
- `b`: bed elevation (m).
- `abs_v_b`: basal velocity magnitude (m/s).
- `A`: viscocity parameter in Glen's flow law.
- `ṁ_over_ρ_w`: basal melt rate divided by water density (kg m⁻² s⁻¹).
"""
function load_Kazmierczak(path::String)
    
    data = matread(path)
    Nx, Ny = size(data["H"])
    mask = data["MASKo"]
    h = data["H"]
    b = data["B"]
    abs_v_b = data["ub"]
    A_visc = data["A"]
    ṁ_over_ρ_w = perYear2perSecond.(data["Bmelt"]) # They stored this variable in per year units.

    # Note: x and y are swapped in the file, and converted from km to m
    xc = Km2m.(data["y"]) 
    yc = Km2m.(data["x"])
    xlims, ylims = compute_lims(xc, yc)

    bed_rheology = :hard

    function initialize_κ!(Nx, Ny, b, bed_rheology::Symbol = :hard)
        
        T = eltype(b)

        if bed_rheology == :hard 
            κ = zeros(T, Nx, Ny) 
        elseif bed_rheology == :soft 
            κ = one(T, Nx, Ny) 
        elseif bed_rheology == :mixed 
            κ = zeros(T, Nx, Ny) 
            κ[b .< -1000] .= T(1.0)
        end

        return κ

    end

    κ = initialize_κ!(Nx, Ny, b, bed_rheology) 

    return Nx, Ny, xlims, ylims, mask, h, b, abs_v_b, A_visc, ṁ_over_ρ_w, κ
    
end
