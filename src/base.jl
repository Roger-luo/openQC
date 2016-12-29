# type definitions

export AbstractQC, AbstractQuState, AbstractQuCircuit

# abstracts
abstract AbstractQC{N}
abstract AbstractQuState{T,N}
abstract AbstractQuCircuit{N} <: AbstractQC{N}

