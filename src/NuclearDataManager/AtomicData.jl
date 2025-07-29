# MOCNTS/src/NuclearDataManager/AtomicData.jl

module AtomicData

export ATOMIC_MASS, N_A

# Avogadro's number in [atoms/mol]
const N_A = 6.02214076e23

# A library of atomic masses for common reactor nuclides in [g/mol]
const ATOMIC_MASS = Dict{String, Float64}(
    "H-1"   => 1.007825,
    "O-16"  => 15.994915,
    "B-10"  => 10.012937,
    "Zr-90" => 89.904704,
    "U-235" => 235.043930,
    "U-238" => 238.050788,
    "Pu-239"=> 239.052163
    # We can add many more nuclides here as needed
)

end # module AtomicData