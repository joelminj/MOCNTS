# MOCNTS/src/MOCSolver/KEigenvalueUpdate.jl

module KEigenvalueUpdate

using ...Types: Problem
using ..SolverState: State
using ...NuclearDataManager: MicroscopicXS
using ...GeometryTracer.FlattenedGeometry: FlatGeometry

export update_keff!

"""
Updates the k-eigenvalue based on the current flux distribution.
k_eff = (Total Fission Production) / (Total Absorption)
"""
function update_keff!(state::State,
                      problem::Problem,
                      cross_sections::Dict{String, MicroscopicXS},
                      flat_geometry::FlatGeometry)
    total_fission_production = 0.0
    total_absorption = 0.0
    # Print diagnostics for first 5 iterations
    if state.old_k_eff == 1.0 || isnan(state.k_eff) || state.k_eff == 0.0
        println("[DIAG] update_keff! called. Iteration likely early or unstable.")
    end
    for i in 1:state.num_fsrs
        fsr = flat_geometry.fsrs[i]
        material = problem.materials[fsr.material_id]
        volume = flat_geometry.volumes[i]
        for g in 1:state.num_groups
            macro_nu_fission_xs = 0.0
            macro_absorption_xs = 0.0
            for (nuclide, atom_density) in material.atom_densities
                xs = cross_sections["$(nuclide)_300.0K"]
                macro_nu_fission_xs += atom_density * xs.nu_fission_xs[g]
                macro_absorption_xs += atom_density * xs.absorption_xs[g]
            end
            flux_times_volume = state.scalar_flux[i, g] * volume
            total_fission_production += macro_nu_fission_xs * flux_times_volume
            total_absorption += macro_absorption_xs * flux_times_volume
            # Print diagnostic for first 3 iterations and first FSRs
            if state.old_k_eff == 1.0 && i <= 2 && g <= 2
                println("[DIAG] FSR $i, Group $g: macro_nu_fission_xs=$(macro_nu_fission_xs), macro_absorption_xs=$(macro_absorption_xs), flux_times_volume=$(flux_times_volume)")
            end
        end
    end
    println("[DIAG] total_fission_production=$(total_fission_production), total_absorption=$(total_absorption)")
    if total_absorption == 0.0
        println("[DIAG] WARNING: total_absorption is zero in update_keff! This will cause NaN k_eff.")
    end
    state.k_eff = state.old_k_eff * (total_fission_production / (total_absorption == 0.0 ? 1e-12 : total_absorption))
end

end # module KEigenvalueUpdate