# MOCNTS/src/ProblemManager/ProblemManager.jl

module ProblemManager

# Import shared types instead of defining our own
using ..Types
using ..Types: Material, SolverSettings, Problem, AbstractUniverse

# Re-export for convenience
export Material, SolverSettings, Problem, AbstractUniverse

end # module ProblemManager