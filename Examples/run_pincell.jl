# MOCNTS/examples/run_pincell.jl

# Import the main module to get access to all the functions and types
using MOCNTS
using MOCNTS.GeometryTracer.Hierarchy: Cell, Universe
using MOCNTS: ZCylinder, XPlane, YPlane, Halfspace, RegionIntersection, RegionUnion, RegionComplement
using LinearAlgebra # For norm()
using MOCNTS.NuclearDataManager: MicroscopicXS
using MOCNTS.Flattener: flatten_geometry, FlatGeometry
using MOCNTS.GeometryTracer: generate_tracks
using MOCNTS: solve
println("Starting MOCNTS Pincell Benchmark...")

#======= 1. DEFINE GEOMETRY =======#
# Surfaces for a simple pincell
fuel_radius = 0.41  # cm
pin_pitch = 1.26    # cm
half_pitch = pin_pitch / 2.0

s_fuel_pin = ZCylinder(0.0, 0.0, fuel_radius)
s_box_left = XPlane(-half_pitch)
s_box_right = XPlane(half_pitch)
s_box_bottom = YPlane(-half_pitch)
s_box_top = YPlane(half_pitch)
# Regions defined by these surfaces


# Moderator is just the box
r_moderator = RegionIntersection([
    Halfspace(s_box_left, 1),
    Halfspace(s_box_right, -1),
    Halfspace(s_box_bottom, 1),
    Halfspace(s_box_top, -1)
])
# Fuel is intersection of box and inside cylinder
r_fuel = RegionIntersection([
    Halfspace(s_fuel_pin, -1),
    r_moderator
])

# Cell 1: Moderator (just the box)
c_moderator = Cell(1, r_moderator, 102)
# Cell 2: Fuel (intersection of box and cylinder)
c_fuel = Cell(2, r_fuel, 101)

# The root universe containing these cells
root_universe = Universe(0, Dict(1 => c_moderator, 2 => c_fuel))


#======= 2. DEFINE MATERIALS =======#

# Materials for ProblemManager (solver)
mat_fuel = MOCNTS.ProblemManager.Material(101, "UO2 Fuel", 10.4, Dict("U-235" => 0.04, "U-238" => 0.96, "O-16" => 2.0)) # Note: Composition is simplified
mat_water = MOCNTS.ProblemManager.Material(102, "Water", 1.0, Dict("H-1" => 2.0, "O-16" => 1.0)) # Note: Composition is simplified
materials = Dict(101 => mat_fuel, 102 => mat_water)

# Materials for GeometryTracer (geometry)
mat_fuel_geom = MOCNTS.GeometryTracer.ProblemManager.Material(101, "UO2 Fuel", 10.4, Dict("U-235" => 0.04, "U-238" => 0.96, "O-16" => 2.0))
mat_water_geom = MOCNTS.GeometryTracer.ProblemManager.Material(102, "Water", 1.0, Dict("H-1" => 2.0, "O-16" => 1.0))
materials_geom = Dict(101 => mat_fuel_geom, 102 => mat_water_geom)


#======= 3. DEFINE SOLVER SETTINGS =======#
# Increase the number of azimuthal angles and decrease ray spacing for better coverage
settings = MOCNTS.ProblemManager.SolverSettings(
    128,    # number of azimuthal angles (increase for better coverage)
    8,      # number of polar angles
    0.05,   # ray spacing (cm) - finer for better coverage
    0.01,   # convergence criterion
    100     # max iterations
)
settings_geom = MOCNTS.GeometryTracer.SolverSettings(
    128,    # number of azimuthal angles (increase for better coverage)
    8,      # number of polar angles
    0.05,   # ray spacing (cm) - finer for better coverage
    0.01,   # convergence criterion
    100     # max iterations
)


#======= 4. ASSEMBLE THE PROBLEM =======#
# In a real case, we'd have a dictionary of all universes. For the pincell,
# the root universe is the only one.
universes = Dict(0 => root_universe)




# Create both problem types
problem_geom = MOCNTS.GeometryTracer.ProblemManager.Problem(materials_geom, universes, 0, settings_geom)
problem = MOCNTS.ProblemManager.Problem(materials, universes, 0, settings)


#======= 5. RUN THE SOLVER PIPELINE =======#
# For now, we are creating a dummy cross-section dictionary.
# In a real run, you would uncomment the call to process_nuclear_data!
println("Step 1: Processing nuclear data (using dummy data for now)...")
# Ensure all nuclides and groups have nonzero cross sections
dummy_xs_data = Dict(
    "U-235_300.0K" => MicroscopicXS(
        "U-235", 300.0,
        [1.0, 1.0],           # total_xs: nonzero for both groups
        [0.3, 0.1],           # absorption_xs
        [0.7, 0.2],           # nu_fission_xs
        [0.2 0.1; 0.1 0.2]    # scatter_matrix
    ),
    "U-238_300.0K" => MicroscopicXS(
        "U-238", 300.0,
        [0.8, 0.8],           # total_xs
        [0.6, 0.2],           # absorption_xs
        [0.0, 0.0],           # nu_fission_xs (U-238 is not fissile)
        [0.1 0.05; 0.05 0.1]  # scatter_matrix
    ),
    "O-16_300.0K" => MicroscopicXS(
        "O-16", 300.0,
        [0.5, 0.5],           # total_xs
        [0.01, 0.01],         # absorption_xs
        [0.0, 0.0],           # nu_fission_xs
        [0.4 0.1; 0.1 0.4]    # scatter_matrix
    ),
    "H-1_300.0K" => MicroscopicXS(
        "H-1", 300.0,
        [1.0, 1.0],           # total_xs
        [0.01, 0.01],         # absorption_xs
        [0.0, 0.0],           # nu_fission_xs
        [0.8 0.5; 0.5 0.8]    # scatter_matrix
    )
)
# processed_xs = process_nuclear_data!(problem) # UNCOMMENT FOR REAL NJOY RUN

println("Step 2: Flattening geometry...")
flat_geometry_incomplete = flatten_geometry(problem_geom)

# Manual FSR Area Calculation (since the flattener doesn't do this yet)
mod_area = pin_pitch^2 - π * fuel_radius^2
fuel_area = π * fuel_radius^2
fsr_volumes = [mod_area, fuel_area] # Moderator first, then fuel

flat_geometry = FlatGeometry(flat_geometry_incomplete.fsrs, fsr_volumes)

# Diagnostic printout of FSRs and their volumes
println("FSR count: ", length(flat_geometry.fsrs))
for (i, fsr) in enumerate(flat_geometry.fsrs)
    println("FSR ", i, ": ", fsr, ", Volume: ", fsr_volumes[i])
end

# Calculate atom densities for all materials (required for macroscopic XS)
for mat in values(materials)
    MOCNTS.NuclearDataManager.calculate_atom_densities!(mat)
end

# Explicitly define bounding box for track generation
bounding_box = (-half_pitch-0.01, -half_pitch-0.01, half_pitch+0.01, half_pitch+0.01)

println("Step 3: Generating tracks...")
# Pass both the problem geometry and the bounding box for track generation
println("Using bounding box: ", bounding_box)
tracks = generate_tracks(problem_geom, flat_geometry, bounding_box)

if isnothing(tracks) || isempty(tracks)
    println("Error: No tracks were generated. Check geometry and ray tracer settings.")
else
    println("Step 4: Starting MOC solver...")
    println("Number of tracks generated: ", length(tracks))
    println("Sample track starting points and angles:")
    for i in 1:min(5, length(tracks))
        println("Track ", i, ": p₀ = ", tracks[i].p₀, ", ϕ = ", tracks[i].ϕ)
    end
    final_state = solve(problem, flat_geometry, dummy_xs_data, tracks)

    println("\n--- 🏁 Final Results 🏁 ---")
    println("Converged k_eff: ", round(final_state.k_eff, digits=6))
    println("Group 1 Flux (Fuel): ", final_state.scalar_flux[1, 1])
    println("Group 2 Flux (Fuel): ", final_state.scalar_flux[1, 2])
    println("Group 1 Flux (Mod):  ", final_state.scalar_flux[2, 1])
    println("Group 2 Flux (Mod):  ", final_state.scalar_flux[2, 2])
end