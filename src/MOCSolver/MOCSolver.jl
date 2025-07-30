# MOCNTS/src/MOCSolver/MOCSolver.jl

module MOCSolver

using ..Types: Problem  # Use shared types
using ..NuclearDataManager: MicroscopicXS
using ..GeometryTracer
using ..GeometryTracer.FlattenedGeometry # To get the number of FSRs
using ..GeometryTracer.RayTracer: Track

export solve

# ... SolverState struct remains the same ...
mutable struct SolverState
    num_groups::Int
    num_fsrs::Int
    scalar_flux::Matrix{Float64}
    old_scalar_flux::Matrix{Float64}
    total_source::Matrix{Float64}
    k_eff::Float64
    old_k_eff::Float64
end


"""
The main public function to solve the transport problem.
It initializes the solver state and runs the iterative process.
"""
function solve(problem::Problem,
               flat_geometry::FlatGeometry,
               cross_sections::Dict{String, MicroscopicXS},
               tracks::Vector{Track})

    #======= Initialization =======#
    println("Initializing MOC solver...")

    # For now, assume 2 energy groups. A robust implementation would get this
    # from the cross_sections data.
    num_groups = 2
    num_fsrs = length(flat_geometry.fsrs)
    max_iter = problem.settings.max_iterations
    tolerance = problem.settings.tolerance

    # Initialize all flux and source arrays
    scalar_flux = ones(Float64, num_fsrs, num_groups) # Start with a guess of 1.0
    old_scalar_flux = copy(scalar_flux)
    total_source = zeros(Float64, num_fsrs, num_groups)
    
    state = SolverState(num_groups,
                        num_fsrs,
                        scalar_flux,
                        old_scalar_flux,
                        total_source,
                        1.0, # Initial guess for k_eff
                        1.0)

    println("Starting power iterations...")
    for i in 1:max_iter
        # Update the fission and scattering source
        update_source!(state, problem, cross_sections, flat_geometry)

        # Print diagnostic: check for NaN in total_source
        if any(isnan.(state.total_source))
            println("[DIAG] NaN detected in total_source at iteration $i")
            println(state.total_source)
        end

        # Perform the transport sweep to update the scalar flux
        transport_sweep!(state, problem, cross_sections, tracks, flat_geometry)

        # Print diagnostic: check for NaN in scalar_flux
        if any(isnan.(state.scalar_flux))
            println("[DIAG] NaN detected in scalar_flux at iteration $i")
            println(state.scalar_flux)
        end

        # Update the eigenvalue (k_eff)
        update_keff!(state, problem, cross_sections, flat_geometry)
        # Clamp k_eff to avoid zero/negative/NaN
        if isnan(state.k_eff) || state.k_eff <= 0.0
            println("[DIAG] k_eff clamped from $(state.k_eff) to 1e-6 at iteration $i")
            state.k_eff = 1e-6
        end

        # Normalize scalar flux to prevent runaway growth
        flux_sum = sum(state.scalar_flux)
        if flux_sum != 0.0 && !isnan(flux_sum)
            state.scalar_flux .= state.scalar_flux ./ flux_sum
        end

        # Under-relaxation (damping) of flux update
        alpha = 0.7
        state.scalar_flux .= alpha * state.scalar_flux .+ (1 - alpha) * state.old_scalar_flux

        # Check for convergence
        flux_err, k_eff_err = calculate_error(state)
        is_converged = flux_err < tolerance && k_eff_err < tolerance
        println("Iteration: $i, k_eff = $(round(state.k_eff, digits=6)), Flux Err = $(round(flux_err, digits=8)), k_eff Err = $(round(k_eff_err, digits=8))")
        if is_converged
            println("Solution converged after $i iterations. ✅")
            break
        end
        if i == max_iter
            println("Warning: Solution did not converge after $max_iter iterations.")
        end
        # Prepare for the next iteration
        state.old_scalar_flux .= state.scalar_flux
        state.old_k_eff = state.k_eff
    end

    println("Solver finished.")
    return state
end

#======= Internal Helper Functions (to be implemented) =======#

"""
Updates the total source in each FSR based on the current scalar flux.
Q_g = (1/k_eff) * χ_g * Σ_{g'} νΣ_f,g' * ϕ_g'  +  Σ_{g' -> g} * ϕ_g'
"""
function update_source!(state::SolverState,
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


"""
Performs the transport sweep for all angles and groups. This function updates
the scalar flux in the solver state.
"""
function transport_sweep!(state::SolverState,
                          problem::Problem,
                          cross_sections::Dict{String, MicroscopicXS},
                          tracks::Vector{Track},
                          flat_geometry::FlatGeometry)

    num_groups = state.num_groups
    num_azim = problem.settings.num_azim

    # The weight for each azimuthal angle in a 2D quadrature set
    angular_weight = π / num_azim

    # Reset scalar flux for this iteration before accumulating
    state.scalar_flux .= 0.0
    # Track skipped segments for summary and unique warnings
    skipped_segments = Dict{Tuple{Int,Int},Int}()
    total_skipped = 0

    # Loop over all energy groups
    for g in 1:num_groups
        # Loop over all tracks
        for track_idx in 1:length(tracks)
            track = tracks[track_idx]
            # Assume vacuum boundary conditions: flux entering the geometry is zero
            psi_in = 0.0

            # Sweep along the track, segment by segment
            for seg_idx in 1:length(track.segments)
                segment = track.segments[seg_idx]
                fsr_id = segment.fsr_id
                L = segment.length
                
                # Get material and FSR data
                fsr = flat_geometry.fsrs[fsr_id]
                material = problem.materials[fsr.material_id]
                volume = flat_geometry.volumes[fsr_id]
                
                # Get the total source for this FSR and group
                Q = state.total_source[fsr_id, g]

                # Calculate macroscopic total cross section: Σ_t
                macro_total_xs = 0.0
                for (nuclide, atom_density) in material.atom_densities
                    xs = cross_sections["$(nuclide)_300.0K"]
                    macro_total_xs += atom_density * xs.total_xs[g]
                end

                # Print diagnostic for first 2 tracks and segments
                if track_idx <= 2 && seg_idx <= 2 && g <= 2
                    println("[DIAG] Track $track_idx, Segment $seg_idx, FSR $fsr_id, Group $g: macro_total_xs=$(macro_total_xs), volume=$(volume), Q=$(Q)")
                end

                # Guard against zero cross section or volume
                if macro_total_xs == 0.0 || volume == 0.0
                    key = (fsr_id, g)
                    skipped_segments[key] = get(skipped_segments, key, 0) + 1
                    total_skipped += 1
                    # Only print the first warning for each unique (fsr_id, group)
                    if skipped_segments[key] == 1
                        println("[WARN] Skipping segment: macro_total_xs=$(macro_total_xs), volume=$(volume), fsr_id=$(fsr_id), group=$(g)")
                    end
                    continue
                end

                # Solve the 1D transport equation for this segment
                exp_term = exp(-macro_total_xs * L)
                psi_out = psi_in * exp_term + (Q / macro_total_xs) * (1.0 - exp_term)

                # Tally the contribution to the scalar flux for this FSR
                # Δϕ = w/V * (ψ_in - ψ_out) / Σ_t + w/V * (Q/Σ_t) * L
                delta_phi = angular_weight * ((psi_in - psi_out) / (macro_total_xs * volume) + 
                                              (Q * L) / (macro_total_xs * volume))
                
                state.scalar_flux[fsr_id, g] += delta_phi

                # The outgoing flux of this segment is the incoming flux for the next
                psi_in = psi_out
            end
        end
    end
    # Print summary of skipped segments
    if total_skipped > 0
        println("[SUMMARY] Skipped $total_skipped segments with zero cross section or volume.")
        println("[SUMMARY] Unique (fsr_id, group) skipped: ", collect(keys(skipped_segments)))
    end
end

"""
Updates the k-eigenvalue based on the current flux distribution.
k_eff = (Total Fission Production) / (Total Absorption)
"""
function update_keff!(state::SolverState,
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

"""
Calculates the relative error between the current and previous iteration's
scalar flux and k_eff.
"""
function calculate_error(state::SolverState)
    # Calculate L2 norm of the flux difference
    flux_diff = state.scalar_flux .- state.old_scalar_flux
    flux_norm = sqrt(sum(flux_diff.^2))
    
    # Calculate relative error for k_eff
    k_eff_err = abs(state.k_eff - state.old_k_eff) / abs(state.k_eff)
    
    return flux_norm, k_eff_err
end


function check_convergence(state::SolverState, tolerance::Float64)
    # To be implemented
    return false # For now, never converge
end

end # module MOCSolver