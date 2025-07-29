# MOCNTS/src/ProblemManager/ProblemManager.jl

module ProblemManager

using ..GeometryTracer.Hierarchy: Universe, AbstractUniverse

export Material, SolverSettings, Problem

"""
A struct to hold material data. The user provides the composition and mass 
density. The code calculates and stores the atom densities.
"""
mutable struct Material
    id::Int64
    name::String
    mass_density::Float64   # Mass density in [g/cm³]
    
    # User-provided composition mapping nuclide name to weight fraction
    composition::Dict{String, Float64}

    # Calculated atom densities mapping nuclide name to [atoms/barn-cm]
    atom_densities::Dict{String, Float64}

    # Constructor to initialize with empty atom_densities
    function Material(id::Int64, name::String, density::Float64, comp::Dict{String, Float64})
        new(id, name, density, comp, Dict{String, Float64}())
    end
end

struct SolverSettings
    num_azim::Int64
    num_polar::Int64
    track_spacing::Float64
    tolerance::Float64
    max_iterations::Int64
end


struct Problem
    # A dictionary mapping material IDs to Material objects
    materials::Dict{Int64, Material}

    # A dictionary mapping ALL universe IDs to Universe objects
    universes::Dict{Int64, AbstractUniverse}

    # The ID of the top-level Universe that defines the core geometry
    root_universe_id::Int64 

    # Solver settings
    settings::SolverSettings
end

end # module ProblemManager