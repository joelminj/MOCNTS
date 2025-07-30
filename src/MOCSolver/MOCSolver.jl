# MOCNTS/src/MOCSolver/MOCSolver.jl

module MOCSolver

using ..ProblemManager
using ..NuclearDataManager: MicroscopicXS
using ..GeometryTracer.FlattenedGeometry
using ..GeometryTracer.RayTracer: Track

export solve

mutable struct SolverState
    num_groups::Int
    num_fsrs::Int
    scalar_flux::Matrix{Float64}
    old_scalar_flux::Matrix{Float64}
    total_source::Matrix{Float64}
    k_eff::Float64
    old_k_eff::Float64
end

function solve(problem::Problem,
               flat_geometry::FlatGeometry,
               cross_sections::Dict{String, MicroscopicXS},
               tracks::Vector{Track})

    println("Initializing MOC solver...")
    num_groups = 2
    num_fsrs = length(flat_geometry.fsrs)
    max_iter = problem.settings.max_iterations
    tolerance = problem.settings.tolerance

    scalar_flux = ones(Float64, num_fsrs, num_groups)
    old_scalar_flux = copy(scalar_flux)
    total_source = zeros(Float64, num_fsrs, num_groups)
    
    state = SolverState(num_groups, num_fsrs, scalar_flux, old_scalar_flux,
                        total_source, 1.0, 1.0)

    println("Starting power iterations...")
    for i in 1:max_iter
        update_source!(state, problem, cross_sections, flat_geometry)
        transport_sweep!(state, problem, cross_sections, tracks, flat_geometry)
        update_keff!(state, problem, cross_sections, flat_geometry)

        flux_err, k_eff_err = calculate_error(state)
        println("Iteration: $i, k_eff = $(round(state.k_eff, digits=6)), Flux Err = $(round(flux_err, digits=8)), k_eff Err = $(round(k_eff_err, digits=8))")
        
        if flux_err < tolerance && k_eff_err < tolerance
            println("Solution converged after $i iterations. ✅")
            break
        end
        if i == max_iter
            println("Warning: Solution did not converge after $max_iter iterations.")
        end

        state.old_scalar_flux .= state.scalar_flux
        state.old_k_eff = state.k_eff
    end

    println("Solver finished.")
    return state
end

function update_source!(state::SolverState, problem::Problem, cross_sections::Dict{String, MicroscopicXS}, flat_geometry::FlatGeometry)
    num_groups = state.num_groups
    fission_spectrum = zeros(num_groups)
    fission_spectrum[1] = 1.0
    state.total_source .= 0.0

    for i in 1:state.num_fsrs
        fsr = flat_geometry.fsrs[i]
        material = problem.materials[fsr.material_id]
        
        total_fission_production = 0.0
        for g_prime in 1:num_groups
            macro_nu_fission_xs = sum(atom_density * cross_sections["$(nuclide)_300.0K"].nu_fission_xs[g_prime] for (nuclide, atom_density) in material.atom_densities)
            total_fission_production += macro_nu_fission_xs * state.scalar_flux[i, g_prime]
        end

        for g in 1:num_groups
            state.total_source[i, g] += (1.0 / state.old_k_eff) * fission_spectrum[g] * total_fission_production
            
            scattering_source = 0.0
            for g_prime in 1:num_groups
                macro_scatter_xs = sum(atom_density * cross_sections["$(nuclide)_300.0K"].scatter_matrix[g, g_prime] for (nuclide, atom_density) in material.atom_densities)
                scattering_source += macro_scatter_xs * state.scalar_flux[i, g_prime]
            end
            state.total_source[i, g] += scattering_source
        end
    end
end

function transport_sweep!(state::SolverState, problem::Problem, cross_sections::Dict{String, MicroscopicXS}, tracks::Vector{Track}, flat_geometry::FlatGeometry)
    num_groups = state.num_groups
    num_azim = problem.settings.num_azim
    angular_weight = π / num_azim
    state.scalar_flux .= 0.0

    for g in 1:num_groups
        for track in tracks
            psi_in = 0.0
            for segment in track.segments
                fsr_id = segment.fsr_id
                L = segment.length
                fsr = flat_geometry.fsrs[fsr_id]
                material = problem.materials[fsr.material_id]
                volume = flat_geometry.volumes[fsr_id]
                Q = state.total_source[fsr_id, g]

                macro_total_xs = sum(atom_density * cross_sections["$(nuclide)_300.0K"].total_xs[g] for (nuclide, atom_density) in material.atom_densities)

                if macro_total_xs <= 0.0 || volume <= 0.0
                    continue
                end

                exp_term = exp(-macro_total_xs * L)
                psi_out = psi_in * exp_term + (Q / macro_total_xs) * (1.0 - exp_term)
                
                # --- CORRECTED FLUX TALLY ---
                avg_psi = (psi_in - psi_out) / (macro_total_xs * L) + (Q / macro_total_xs)
                delta_phi = angular_weight * L * avg_psi / volume
                state.scalar_flux[fsr_id, g] += delta_phi

                psi_in = psi_out
            end
        end
    end
end

function update_keff!(state::SolverState, problem::Problem, cross_sections::Dict{String, MicroscopicXS}, flat_geometry::FlatGeometry)
    total_fission_production = 0.0
    total_absorption = 0.0
    
    for i in 1:state.num_fsrs
        fsr = flat_geometry.fsrs[i]
        material = problem.materials[fsr.material_id]
        volume = flat_geometry.volumes[i]
        for g in 1:state.num_groups
            macro_nu_fission_xs = sum(atom_density * cross_sections["$(nuclide)_300.0K"].nu_fission_xs[g] for (nuclide, atom_density) in material.atom_densities)
            macro_absorption_xs = sum(atom_density * cross_sections["$(nuclide)_300.0K"].absorption_xs[g] for (nuclide, atom_density) in material.atom_densities)
            
            flux_times_volume = state.scalar_flux[i, g] * volume
            total_fission_production += macro_nu_fission_xs * flux_times_volume
            total_absorption += macro_absorption_xs * flux_times_volume
        end
    end

    # --- CORRECTED EIGENVALUE CALCULATION ---
    if total_absorption > 0.0
        state.k_eff = total_fission_production / total_absorption
    end
end

function calculate_error(state::SolverState)
    flux_diff = state.scalar_flux .- state.old_scalar_flux
    flux_norm = sqrt(sum(flux_diff.^2))
    k_eff_err = abs(state.k_eff - state.old_k_eff)
    return flux_norm, k_eff_err
end

end # module MOCSolver