# MOCNTS/src/MOCNTS.jl

module MOCNTS


"""
MOCNTS main module
Includes submodules and re-exports key types
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


# Forward key geometry types for convenience

const ZCylinder = GeometryTracer.ZCylinder
const XPlane = GeometryTracer.XPlane
const YPlane = GeometryTracer.YPlane
const Halfspace = GeometryTracer.Halfspace
const RegionIntersection = GeometryTracer.RegionIntersection
const RegionUnion = GeometryTracer.RegionUnion
const RegionComplement = GeometryTracer.RegionComplement
const Region = GeometryTracer.Primitives.Region

export ZCylinder, XPlane, YPlane, Halfspace, RegionIntersection, RegionUnion, RegionComplement, Region, Cell, Universe

end # module MOCNTS