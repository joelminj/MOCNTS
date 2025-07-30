# MOCNTS/src/MOCSolver/SolverState.jl

module SolverState

export State, initialize_state

"""
Holds the current state of the MOC solver iteration
"""
mutable struct State
    num_groups::Int
    num_fsrs::Int
    scalar_flux::Matrix{Float64}
    old_scalar_flux::Matrix{Float64}
    total_source::Matrix{Float64}
    k_eff::Float64
    old_k_eff::Float64
end

"""
Initialize the solver state with default values
"""
function initialize_state(num_groups::Int, num_fsrs::Int)
    # Initialize all flux and source arrays
    scalar_flux = ones(Float64, num_fsrs, num_groups) # Start with a guess of 1.0
    old_scalar_flux = copy(scalar_flux)
    total_source = zeros(Float64, num_fsrs, num_groups)
    
    return State(num_groups,
                 num_fsrs,
                 scalar_flux,
                 old_scalar_flux,
                 total_source,
                 1.0, # Initial guess for k_eff
                 1.0)
end

"""
Calculate the relative error between current and previous iteration
"""
function calculate_error(state::State)
    # Calculate L2 norm of the flux difference
    flux_diff = state.scalar_flux .- state.old_scalar_flux
    flux_norm = sqrt(sum(flux_diff.^2))
    
    # Calculate relative error for k_eff
    k_eff_err = abs(state.k_eff - state.old_k_eff) / abs(state.k_eff)
    
    return flux_norm, k_eff_err
end

end # module SolverState