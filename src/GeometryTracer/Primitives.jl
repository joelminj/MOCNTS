# JuliaMOC/src/GeometryTracer/Primitives.jl

module Primitives

export ZCylinder, XPlane, YPlane, Halfspace, RegionIntersection, RegionUnion, RegionComplement

#======= Abstract Types =======#
"""
An abstract type for all geometric surfaces, defined by the implicit
equation f(x, y, z) = 0.
"""
abstract type Surface end

"""
An abstract type for all spatial regions, defined by boolean operations
on surface half-spaces.
"""
abstract type Region end


#======= Concrete Surface Types =======#
"""
A plane perpendicular to the x-axis, defined by x - x₀ = 0.
"""
struct XPlane <: Surface
    x₀::Float64
end

"""
A plane perpendicular to the y-axis, defined by y - y₀ = 0.
"""
struct YPlane <: Surface
    y₀::Float64
end

"""
A cylinder parallel to the z-axis, defined by (x - x₀)² + (y - y₀)² - r² = 0.
"""
struct ZCylinder <: Surface
    x₀::Float64
    y₀::Float64
    radius::Float64
end


#======= Region Types (CSG Logic) =======#
"""
A region defined by one side of a surface.
`sense = -1` for f(r) < 0
`sense = +1` for f(r) > 0
"""
struct Halfspace <: Region
    surface::Surface
    sense::Int8 # -1 or +1
end

"""
A region defined as the intersection of multiple sub-regions.
(e.g., region1 AND region2 AND ...)
"""
struct RegionIntersection <: Region
    regions::Vector{Region}
end

"""
A region defined as the union of multiple sub-regions.
(e.g., region1 OR region2 OR ...)
"""
struct RegionUnion <: Region
    regions::Vector{Region}
end

"""
A region defined as the logical NOT of a sub-region.
"""
struct RegionComplement <: Region
    region::Region
end

end # module Primitives