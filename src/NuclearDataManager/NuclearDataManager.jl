# MOCNTS/src/NuclearDataManager/NuclearDataManager.jl

module NuclearDataManager

using ..AtomicData: N_A
using ..ProblemManager
using ..AtomicData
using ..ProblemManager: Problem, Material

export process_nuclear_data!

# ... MicroscopicXS struct  ...
"""
A struct to hold multi-group microscopic cross section data for a single nuclide.
"""

struct MicroscopicXS
    name::String
    temperature::Float64
    total_xs::Vector{Float64}         # Total cross section per group
    absorption_xs::Vector{Float64}    # Absorption cross section per group
    
    # === NEW FIELDS ===
    nu_fission_xs::Vector{Float64}    # ν * σ_fission per group
    scatter_matrix::Matrix{Float64}   # Scattering probability from group g' to g
end

function process_nuclear_data!(problem::Problem)

    # A dictionary to the problem to store the results
    processed_xs = Dict{String, MicroscopicXS}()

    println("Calculating atom densities...")
    for (id, material) in problem.materials
        calculate_atom_densities!(material)
    end

    println("Collecting unique nuclides and temperatures...")
    (nuclides, temps) = collect_nuclides_and_temps(problem)
    
    # Define paths for our workflow
    njoy_input_path = "xs_data/njoy.inp"
    xs_output_dir = "xs_data"
    
    # === ACTION REQUIRED: UPDATE THESE PATHS ===
    # Use your specific path to the NJOY executable
    njoy_exe = "C:\\Users\\poste\\NJOY2016\\build\\njoy.exe" 
    # Use the path to your NEW directory with plain-text .endf files
    endf_dir = "D:\\MOCNTS\\endf_data"
    # ==========================================

    println("Generating NJOY input deck...")
    generate_njoy_input(nuclides, temps, endf_dir, njoy_input_path)

    println("Running NJOY...")
    run_njoy(njoy_exe, njoy_input_path, xs_output_dir, endf_dir)

    println("Parsing NJOY output...")
    for nuclide in nuclides
        for temp in temps
            gendl_tape_path = joinpath(xs_output_dir, "tape24")
            xs_data = parse_njoy_output(gendl_tape_path, nuclide, temp)
            processed_xs["$(nuclide)_$(temp)K"] = xs_data

            # In a real run, NJOY appends to tape24, so we would need to
            # delete it before running for the next nuclide.
            # For simplicity now, we assume one nuclide or manual cleanup.
        end
    end

    println("Nuclear data processing complete!")
    # The 'processed_xs' dictionary now holds all our data.
    # We will later pass this to the solver.
    return processed_xs 
end

# ... collect_nuclides_and_temps is unchanged ...
function collect_nuclides_and_temps(problem::Problem)
    all_nuclides = Set{String}()
    for (id, material) in problem.materials
        for nuclide_name in keys(material.composition)
            push!(all_nuclides, nuclide_name)
        end
    end
    temperatures = Set([300.0])
    return (all_nuclides, temperatures)
end

# ... generate_njoy_input is unchanged ...
function generate_njoy_input(nuclides::Set{String}, temperatures::Set{Float64}, endf_dir::String, output_path::String)
    input_deck = ""
    group_structure = "1.0e-5 0.625 2.0e7"

    for temp in temperatures
        for nuclide in nuclides
            endf_filename = lowercase(replace(nuclide, "-" => "")) * ".endf"
            input_deck *= """
moder /
20 21 /
reconr /
21 22 /
'$nuclide from $endf_dir/$endf_filename' /
9235 1 20 0.001 /
0.0 /
broadr /
21 22 23 /
9235 1 20 1 /
$temp /
0.0 /
groupр /
21 23 0 24 /
9235 1 20 2 6 2 /
'2-group structure' /
$group_structure /
0.001 /
$temp /
/
"""
        end
    end
    input_deck *= "stop\n"
    
    try
        output_dir = dirname(output_path)
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        open(output_path, "w") do f
            write(f, input_deck)
        end
        println("Successfully wrote NJOY input to $output_path")
    catch e
        println("Error writing NJOY input file: $e")
    end
end


"""
Runs the NJOY executable as a subprocess in a specified working directory.
It assumes the ENDF files are linked or present in that directory.
"""
function run_njoy(njoy_exe_path::String, input_path::String, working_dir::String, endf_dir::String)
    # NJOY reads from standard input and writes to standard output.
    # We will redirect our generated input file to stdin.
    inp_file = open(input_path, "r")
    
    # We will also capture the output and error streams to log files.
    out_log = joinpath(working_dir, "njoy.out")
    err_log = joinpath(working_dir, "njoy.err")

    # The command to be executed. We run it within the working_dir.
    # NJOY needs to find the ENDF files, so we must run it from a directory
    # where they are accessible. We will create symlinks for them.
    
    # Create symbolic links to all necessary ENDF files in the working directory
    for file in readdir(endf_dir)
        ln_path = joinpath(working_dir, file)
        # remove old symlink if it exists
        rm(ln_path, force=true) 
        symlink(joinpath(abspath(endf_dir), file), ln_path)
    end
    
    # Create the command object
    cmd = pipeline(ignorestatus(`$njoy_exe_path`); stdin=inp_file, stdout=out_log, stderr=err_log)
    
    try
        # Change to the working directory to run the command
        cd(working_dir) do
            println("Executing NJOY in directory: $(pwd())")
            process = run(cmd)
            
            if process.exitcode == 0
                println("NJOY executed successfully.")
            else
                println("NJOY failed with exit code $(process.exitcode). Check logs in $working_dir")
            end
        end
    catch e
        println("An error occurred while trying to run NJOY: $e")
        println("Please ensure '$njoy_exe_path' is a valid executable and in your system's PATH.")
    finally
        close(inp_file)
    end
end


# ... (calculate_atom_densities! and the placeholder for parse_njoy_output remain) ...
function calculate_atom_densities!(material::Material)
    for (nuclide, weight_fraction) in material.composition
        if !haskey(ATOMIC_MASS, nuclide)
            error("Atomic mass for nuclide $nuclide not found in library.")
        end
        atomic_mass = ATOMIC_MASS[nuclide]
        number_density = (material.mass_density * weight_fraction * N_A) / atomic_mass
        material.atom_densities[nuclide] = number_density * 1e-24
    end
    return nothing
end

"""
Parses the GENDF output file from NJOY (tape24) to extract multi-group cross sections.
This is a simplified parser for a specific GROUPR output format.
"""
function parse_njoy_output(gendl_tape_path::String, nuclide_name::String, temp::Float64)
    total_xs = Float64[]
    absorption_xs = Float64[]
    lines = readlines(gendl_tape_path)

    # This is a brittle parser and depends on the exact NJOY output format.
    # It looks for specific reaction type codes (MT numbers).
    # MT=1 is the total cross section.
    # MT=27 is the total absorption (fission + capture).

    i = 1
    while i <= length(lines)
        # Find the start of a reaction block
        if occursin("MF=3", lines[i]) # MF=3 contains cross section data
            # The next line contains the reaction type (MT)
            mt_number = parse(Int, split(lines[i+1])[4])

            if mt_number == 1 # Total XS
                # The cross section values are on the line after the next
                xs_values_line = lines[i+3]
                total_xs = [parse(Float64, val) for val in split(xs_values_line)]
            elseif mt_number == 27 # Absorption XS
                xs_values_line = lines[i+3]
                absorption_xs = [parse(Float64, val) for val in split(xs_values_line)]
            end
        end
        i += 1
    end

    if isempty(total_xs) || isempty(absorption_xs)
        error("Failed to parse cross sections for $nuclide_name from $gendl_tape_path. Check njoy.out for errors.")
    end

    return MicroscopicXS(nuclide_name, temp, total_xs, absorption_xs)
end

end # module NuclearDataManager