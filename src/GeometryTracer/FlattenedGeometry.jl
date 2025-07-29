# MOCNTS/src/GeometryTracer/FlattenedGeometry.jl

module FlattenedGeometry

using ..Primitives: Region

export FSR, FlatGeometry

"""
A single Flat-Source Region (FSR).
This is a simple convex region of space filled with a single material.
"""
struct FSR
    id::Int64
    region::Region
    material_id::Int64
end

"""
A flattened representation of the entire problem geometry.
It contains a simple list of all FSRs.
"""
struct FlatGeometry
    fsrs::Vector{FSR}
    volumes::Vector{Float64} # <-- Add this line
    # We can add spatial search acceleration structures here later
end

end # module FlattenedGeometry