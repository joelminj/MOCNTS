# MOCNTS/src/GeometryTracer/GeometryTracer.jl

module GeometryTracer




# Include submodules
include("Primitives.jl")
include("Hierarchy.jl")
include("FlattenedGeometry.jl")
include("Flattener.jl")
include("RayTracer.jl")


# Make submodules available
using .Primitives
using .Hierarchy
using .FlattenedGeometry
using .Flattener
using .RayTracer


# Export key geometry types for convenience
export ZCylinder, XPlane, YPlane, Halfspace, RegionIntersection, RegionUnion, RegionComplement
export Primitives, Hierarchy, FlattenedGeometry, Flattener, RayTracer

end # module GeometryTracer
