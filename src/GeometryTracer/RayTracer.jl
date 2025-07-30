# MOCNTS/src/GeometryTracer/RayTracer.jl

module RayTracer

using LinearAlgebra
using ..Primitives: Region, Surface, XPlane, YPlane, ZCylinder, Halfspace, RegionIntersection, RegionUnion, RegionComplement
using ..FlattenedGeometry
using ...Types: Problem  # Import Problem type from shared Types module

export Track, Segment, generate_tracks

#======= Data Structures =======#
struct Segment
    length::Float64
    fsr_id::Int64
end

struct Track
    p₀::NTuple{2, Float64}
    ϕ::Float64
    segments::Vector{Segment}
end




#======= Helper Functions =======#

"""
Check if a point (x, y) is inside a given Region.
This is a recursive function that evaluates the CSG tree.
"""
function is_inside(p::NTuple{2, Float64}, r::Region)
    if r isa Halfspace
        val = evaluate(p, r.surface)
        return (val * r.sense) > 0
    elseif r isa RegionIntersection
        return all(is_inside(p, sub_r) for sub_r in r.regions)
    elseif r isa RegionUnion
        return any(is_inside(p, sub_r) for sub_r in r.regions)
    elseif r isa RegionComplement
        return !is_inside(p, r.region)
    end
    return false
end

"""
Evaluate the implicit equation f(x, y) for a surface at a given point.
"""
evaluate(p::NTuple{2, Float64}, s::XPlane) = p[1] - s.x₀
evaluate(p::NTuple{2, Float64}, s::YPlane) = p[2] - s.y₀
evaluate(p::NTuple{2, Float64}, s::ZCylinder) = (p[1] - s.x₀)^2 + (p[2] - s.y₀)^2 - s.radius^2

"""
Calculate the distance from a point p along a direction ϕ to a surface.
Returns infinity if the track is parallel or moves away from the surface.
"""
function distance_to_surface(p::NTuple{2, Float64}, ϕ::Float64, s::Surface)
    Ω = (cos(ϕ), sin(ϕ))
    dist = Inf

    if s isa XPlane
        # Distance = (x_plane - x_p) / cos(ϕ)
        if Ω[1] != 0
            d = (s.x₀ - p[1]) / Ω[1]
            dist = d > 1e-10 ? d : Inf
        end
    elseif s isa YPlane
        # Distance = (y_plane - y_p) / sin(ϕ)
        if Ω[2] != 0
            d = (s.y₀ - p[2]) / Ω[2]
            dist = d > 1e-10 ? d : Inf
        end
    elseif s isa ZCylinder
        # Solve quadratic equation for intersection of line and circle
        x₀, y₀, r = s.x₀, s.y₀, s.radius
        px, py = p[1] - x₀, p[2] - y₀
        Ωx, Ωy = Ω[1], Ω[2]
        
        a = Ωx^2 + Ωy^2
        b = 2 * (px * Ωx + py * Ωy)
        c = px^2 + py^2 - r^2
        
        discriminant = b^2 - 4*a*c
        if discriminant >= 0
            sqrt_d = sqrt(discriminant)
            d1 = (-b + sqrt_d) / (2a)
            d2 = (-b - sqrt_d) / (2a)
            # Return the smallest positive distance
            if d1 > 1e-10; dist = min(dist, d1); end
            if d2 > 1e-10; dist = min(dist, d2); end
        end
    end
    return dist
end


"""
Find the FSR that contains a given point p.
This can be slow for large numbers of FSRs.
"""
function find_fsr_at_point(p::NTuple{2, Float64}, flat_geometry::FlatGeometry)
    for fsr in flat_geometry.fsrs
        if is_inside(p, fsr.region)
            return fsr
        end
    end
    return nothing # Point is outside the geometry
end


#======= Main Function =======#

function generate_tracks(problem::Problem, flat_geometry::FlatGeometry, bounding_box::NTuple{4, Float64})
    # Get parameters from solver settings
    num_azim = problem.settings.num_azim
    spacing = problem.settings.track_spacing
    
    # Unpack the bounding box
    min_x, min_y, max_x, max_y = bounding_box

    dϕ = π / num_azim
    all_tracks = Track[]

    for i in 1:num_azim
        ϕ = (i - 0.5) * dϕ
        for x in min_x:spacing:max_x
            for y in min_y:spacing:max_y
                p_start = (x, y)
                # ...existing code...
                trace_single_track!(all_tracks, p_start, ϕ, flat_geometry)
            end
        end
    end
    
    # ...existing code...
    return all_tracks
end


function trace_single_track!(all_tracks::Vector{Track}, 
                             p_start::NTuple{2, Float64}, 
                             ϕ::Float64,
                             flat_geometry::FlatGeometry)
    segments = Segment[]
    p_current = p_start
    Ω = (cos(ϕ), sin(ϕ))

    # Small step for moving across boundaries
    ϵ = 1e-9

    while true
        current_fsr = find_fsr_at_point(p_current, flat_geometry)
        
        # If the point is outside the geometry, stop tracing this track
        if isnothing(current_fsr)
            break
        end

        # Find the minimum distance to any surface of the current FSR
        min_dist = Inf
        for r in current_fsr.region.regions # Assumes RegionIntersection
            if r isa Halfspace
                dist = distance_to_surface(p_current, ϕ, r.surface)
                min_dist = min(min_dist, dist)
            end
        end

        # If no surface is found ahead, the track escapes
        if isinf(min_dist)
            break
        end

        # Add the segment to the list
        push!(segments, Segment(min_dist, current_fsr.id))
        
        # Move the point to the boundary and a tiny bit across it
        p_current = p_current .+ (min_dist + ϵ) .* Ω
    end
    
    # If we found any segments, add the completed track to our list
    if !isempty(segments)
        push!(all_tracks, Track(p_start, ϕ, segments))
    end
end


end # module RayTracer