# MOCNTS/src/MOCSolver/TransportSweep.jl

module TransportSweep

using ...Types: Problem
using ..SolverState: State
using ...NuclearDataManager: MicroscopicXS
using ...GeometryTracer.FlattenedGeometry: FlatGeometry
using ...GeometryTracer.RayTracer: Track

export transport_sweep!

"""
Performs the transport sweep for all angles and groups. This function updates
the scalar flux in the solver state.
"""
function transport_sweep!(state::State,
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

end # module TransportSweep