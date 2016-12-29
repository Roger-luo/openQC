using ExpmV, Expokit
export AQCprocessor

abstract AQCprocessor{N}

include("hamiltonian.jl")
include("aqc.jl")
include("adiaprocessing.jl")
include("demoncooling.jl")