# JuliaMOC/src/GeometryTracer/Hierarchy.jl

module Hierarchy

include("Primitives.jl")

using ..Primitives: Region  # Import our previously defined primitives

# Export list
export Cell, Universe, Lattice, AbstractUniverse

abstract type AbstractUniverse end

#======= Hierarchical Geometry Types =======#

"""
A Cell is a spatial region filled with a specific material or another universe.
It connects a geometric `Region` to a `fill`.
The `fill` is identified by an integer ID. If the ID points to a Universe,
it's a nested geometry. If it points to a Material, it's a terminal cell.
"""

struct Cell
    id::Int64
    region::Region
    fill_id::Int64
end

struct Universe <: AbstractUniverse
    id::Int64
    cells::Dict{Int64, Cell}
end

Universe(id::Int64) = Universe(id, Dict{Int64, Cell}())


"""
A Lattice is a special type of Universe that arranges other Universes
in a regular 2D or 3D grid. This is used to model reactor cores and fuel assemblies.
"""

struct Lattice{N} <: AbstractUniverse  # N is the dimension (e.g., 2 for a 2D lattice)
    id::Int64
    pitch::NTuple{N, Float64}
    universes::Array{Int64, N}
    # Location finding logic can be added here later
end


end # module Hierarchy