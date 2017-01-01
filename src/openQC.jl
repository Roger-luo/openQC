module openQC

#Package Settings
pkg_dir = Pkg.dir("openQC")

#consts 
include("utils/consts.jl")
include("utils/math.jl")
include("utils/qusat.jl")

####################
# type definitions #
####################
# quantum computation abstracts 
include("base.jl")
# state types
include("states/states.jl")

##############
# Simulators #
##############

# AQC model
include("qc/adiabatic/adiabatic.jl")


end # module
