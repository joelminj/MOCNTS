# MOCNTS/src/GeometryTracer/Flattener.jl

module Flattener

using ..Primitives: Region, XPlane, YPlane, ZCylinder, Halfspace, RegionIntersection, RegionUnion, RegionComplement
using ..Hierarchy
using ..FlattenedGeometry
using ..ProblemManager

export flatten_geometry

# Needed for type-stable FSR vector
using ..FlattenedGeometry: FSR 

#======= Coordinate Transformation Helpers =======#

# Define how to translate different surfaces. We use multiple dispatch.
transform_surface(s::XPlane, translation::NTuple{2, Float64}) = XPlane(s.x₀ + translation[1])
transform_surface(s::YPlane, translation::NTuple{2, Float64}) = YPlane(s.y₀ + translation[2])
transform_surface(s::ZCylinder, translation::NTuple{2, Float64}) = ZCylinder(s.x₀ + translation[1], s.y₀ + translation[2], s.radius)

# Recursively transform a region by applying the transformation to its surfaces
function transform_region(r::Region, translation::NTuple{2, Float64})
    if r isa Halfspace
        return Halfspace(transform_surface(r.surface, translation), r.sense)
    elseif r isa RegionIntersection
        return RegionIntersection([transform_region(sub_r, translation) for sub_r in r.regions])
    elseif r isa RegionUnion
        return RegionUnion([transform_region(sub_r, translation) for sub_r in r.regions])
    elseif r isa RegionComplement
        return RegionComplement(transform_region(r.region, translation))
    end
end


#======= Main Flattener Implementation =======#

"""
The main public function to flatten the geometry.
It initializes the process by calling the recursive helper on the root universe.
"""
function flatten_geometry(problem::Problem)
    fsrs = FSR[]
    fsr_id_counter = 0
    
    # Get the root universe from the dictionary
    root_universe = problem.universes[problem.root_universe_id]

    # Start the recursion at the top-level universe with zero translation
    (fsrs, _) = flatten_recursive!(fsrs, fsr_id_counter, root_universe, (0.0, 0.0), problem)
    
    println("Geometry flattening complete. Found $(length(fsrs)) FSRs.")
    return FlatGeometry(fsrs, []) # Return with empty volumes for now
end

"""
The recursive helper function that traverses the geometry hierarchy.
"""
function flatten_recursive!(fsrs::Vector{FSR}, 
                           fsr_id::Int64,
                           universe::AbstractUniverse, 
                           translation::NTuple{2, Float64}, 
                           problem::Problem)

    # Case 1: The current universe is a Lattice
    if universe isa Lattice
        lattice = universe
        (nx, ny) = size(lattice.universes)
        half_pitch_x, half_pitch_y = lattice.pitch ./ 2.0
        
        for i in 1:nx, j in 1:ny
            # Calculate the center of the lattice cell in the parent's coordinates
            center_x = (i - (nx + 1) / 2.0) * lattice.pitch[1]
            center_y = (j - (ny + 1) / 2.0) * lattice.pitch[2]
            
            # The new translation is the parent's translation plus this cell's offset
            new_translation = (translation[1] + center_x, translation[2] + center_y)
            
            # Get the ID of the universe filling this lattice cell
            fill_universe_id = lattice.universes[i, j]
            
            # Use problem.universes as the dictionary of universes
            fill_universe = problem.universes[fill_universe_id]
            
            # Recurse into the universe filling this lattice cell
            (fsrs, fsr_id) = flatten_recursive!(fsrs, fsr_id, fill_universe, new_translation, problem)
        end

    # Case 2: The current universe is a simple collection of Cells
    else 
        for (cell_id, cell) in universe.cells
            # Check if the cell is filled with a material (base case) or another universe (recursive step)
            if haskey(problem.materials, cell.fill_id)
                # BASE CASE: This cell is filled with a material.
                
                # Apply the current translation to the cell's region to get its global position
                global_region = transform_region(cell.region, translation)
                
                # Create a new FSR
                fsr_id += 1
                new_fsr = FSR(fsr_id, global_region, cell.fill_id)
                push!(fsrs, new_fsr)
            else
                # RECURSIVE STEP: This cell is filled with another universe.
                
                # Find the universe that fills this cell
                fill_universe = problem.universes[cell.fill_id]
                
                # The translation for the next level is the current translation.
                # Here we assume the inner universe's coordinates are relative to the outer one's origin.
                # If cells could also be translated, we would add that here.
                (fsrs, fsr_id) = flatten_recursive!(fsrs, fsr_id, fill_universe, translation, problem)
            end
        end
    end
    
    return fsrs, fsr_id
end


end # module Flattener