# MOCNTS/src/MOCNTS.jl

module MOCNTS

"""
MOCNTS main module - Method of Characteristics Neutron Transport Solver
Includes submodules and re-exports key types for user convenience
"""

# --- Shared Types (must come first) ---
include("Types.jl")

# --- Geometry ---
include("GeometryTracer/GeometryTracer.jl")

# --- Physics and Problem Definition ---
include("ProblemManager/ProblemManager.jl")
include("NuclearDataManager/NuclearDataManager.jl")

# --- Solver ---
include("MOCSolver/MOCSolver.jl")

# Make the modules available to code that uses MOCNTS
using .Types
using .GeometryTracer
using .ProblemManager
using .NuclearDataManager
using .MOCSolver

# Re-export commonly used types and functions for user convenience
export Material, SolverSettings, Problem
export ZCylinder, XPlane, YPlane, Halfspace, RegionIntersection, RegionUnion, RegionComplement, Region
export Cell, Universe
export solve, process_nuclear_data!

# Forward key geometry types for convenience - these provide direct access without module prefix
const ZCylinder = GeometryTracer.ZCylinder
const XPlane = GeometryTracer.XPlane
const YPlane = GeometryTracer.YPlane
const Halfspace = GeometryTracer.Halfspace
const RegionIntersection = GeometryTracer.RegionIntersection
const RegionUnion = GeometryTracer.RegionUnion
const RegionComplement = GeometryTracer.RegionComplement
const Region = GeometryTracer.Primitives.Region
const Cell = GeometryTracer.Hierarchy.Cell
const Universe = GeometryTracer.Hierarchy.Universe

end # module MOCNTS