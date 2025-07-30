# MOCNTS/src/MOCSolver/SourceUpdate.jl

module SourceUpdate

using ...Types: Problem
using ..SolverState: State
using ...NuclearDataManager: MicroscopicXS
using ...GeometryTracer.FlattenedGeometry: FlatGeometry

export update_source!

"""
Updates the total source in each FSR based on the current scalar flux.
Q_g = (1/k_eff) * χ_g * Σ_{g'} νΣ_f,g' * ϕ_g'  +  Σ_{g' -> g} * ϕ_g'
"""
function update_source!(state::State,
                        problem::Problem,
                        cross_sections::Dict{String, MicroscopicXS},
                        flat_geometry::FlatGeometry)

    num_fsrs = state.num_fsrs
    num_groups = state.num_groups
    
    # Fission spectrum (χ). For now, assume all fission neutrons are born in group 1.
    fission_spectrum = zeros(num_groups)
    fission_spectrum[1] = 1.0

    # Reset the source from the previous iteration
    state.total_source .= 0.0

    # Loop over every FSR to calculate its source
    for i in 1:num_fsrs
        fsr = flat_geometry.fsrs[i]
        material = problem.materials[fsr.material_id]
        
        # --- 1. Calculate the total fission source in this FSR ---
        total_fission_production = 0.0
        for g_prime in 1:num_groups # Loop over source groups g'
            # Calculate macroscopic fission production cross section: νΣ_f
            macro_nu_fission_xs = 0.0
            for (nuclide, atom_density) in material.atom_densities
                # This assumes a single temperature for now
                xs = cross_sections["$(nuclide)_300.0K"]
                macro_nu_fission_xs += atom_density * xs.nu_fission_xs[g_prime]
            end
            
            total_fission_production += macro_nu_fission_xs * state.scalar_flux[i, g_prime]
        end

        # Distribute the fission source over all energy groups using the fission spectrum
        for g in 1:num_groups
            fission_source = (1.0 / state.k_eff) * fission_spectrum[g] * total_fission_production
            state.total_source[i, g] += fission_source
        end

        # --- 2. Calculate the scattering source into each group g ---
        for g in 1:num_groups # Loop over destination groups g
            scattering_source = 0.0
            for g_prime in 1:num_groups # Loop over source groups g'
                # Calculate macroscopic scattering cross section: Σ_{g' -> g}
                macro_scatter_xs = 0.0
                for (nuclide, atom_density) in material.atom_densities
                    xs = cross_sections["$(nuclide)_300.0K"]
                    macro_scatter_xs += atom_density * xs.scatter_matrix[g, g_prime]
                end
                
                scattering_source += macro_scatter_xs * state.scalar_flux[i, g_prime]
            end
            state.total_source[i, g] += scattering_source
        end
    end
end

end # module SourceUpdate