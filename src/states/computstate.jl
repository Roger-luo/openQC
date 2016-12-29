import Base: copy, *, +, -
export ComputState, rawcoeffs, coeff, BellState

type ComputState{A<:AbstractVector,T,N} <: AbstractQuState{T,N}
    coeffs::A
end

# constructors
ComputState{A<:AbstractVector}(state_vec::A) = ComputState{A,eltype(A),length(state_vec)|>log2|>Int}(state_vec)

coefftype{A,T,N}(state::ComputState{A,T,N}) = A
rawcoeffs(state::ComputState) = state.coeffs
coeff(state::ComputState) = state.coeffs

# Base functions overload
copy(state::ComputState) = ComputState(copy(state.coeffs))
similar_type{Q<:ComputState}(::Type{Q}) = ComputState

# Shortcuts
BellState(n::Int) = ComputState(Complex128[1/sqrt(2^n) for i=1:2^n])

# Reload Operators
(*)(a::Number,B::ComputState) = a*coeff(B)
(*)(B::ComputState,a::Number) = a*B
(*)(A::ComputState,B::ComputState) = coeff(A)*coeff(B)
(+)(A::ComputState,B::ComputState) = coeff(A)+coeff(B)
(-)(A::ComputState,B::ComputState) = coeff(A)-coeff(B)
