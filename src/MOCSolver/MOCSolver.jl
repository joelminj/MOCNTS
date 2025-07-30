# MOCNTS/src/MOCSolver/MOCSolver.jl

module MOCSolver

using ..Types: Problem  # Use shared types
using ..NuclearDataManager: MicroscopicXS
using ..GeometryTracer.FlattenedGeometry: FlatGeometry
using ..GeometryTracer.RayTracer: Track

# Include sub-modules
include("SolverState.jl")
include("SourceUpdate.jl")
include("TransportSweep.jl")
include("KEigenvalueUpdate.jl")

using .SolverState
using .SolverState: calculate_error, initialize_state
using .SourceUpdate: update_source!
using .TransportSweep: transport_sweep!
using .KEigenvalueUpdate: update_keff!

export solve

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

    state = initialize_state(num_groups, num_fsrs)

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

end # module MOCSolver