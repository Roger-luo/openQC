import Base:show
export DemonCooling, t, gamma, spec, steps, tol, hamiltonian, timer, probability, cooling!

immutable DemonCooling{N}<:AQCprocessor{N}
    t::Real
    gamma::Real
    spectrum::AbstractVector
    H::AbstractMatrix
    k::Int # k steps
    tol::Int # tolerance for error correction, at least $tol zeros
end

t(dc::DemonCooling) = dc.t
gamma(dc::DemonCooling) = dc.gamma
spec(dc::DemonCooling) = dc.spectrum
steps(dc::DemonCooling) = dc.k
tol(dc::DemonCooling) = dc.tol
timer(dc::DemonCooling) = dc.t

function DemonCooling(H::AbstractSparseMatrix; k::Int=5, tol::Int=2)
    maxEigen = eigs(H;nev=1,which=:LR)[1][1]|>real
    spec = eigfact(full(H))[:values]
    return DemonCooling{Int(log2(size(H)[1]))}(π/(2*maxEigen),0,spec,H,k,tol)
end

DemonCooling(H::AbstractMatrix;k::Int=5, tol::Int=2) = DemonCooling(sparse(H);k=k,tol=tol)

function show{N}(io::IO,dc::DemonCooling{N})
    println(io,"----------------------------")
    println(io,"$N qubits demon-like cooling")
    println(io,"----------------------------")
    println(io,"t = $(dc.t)")
    println(io,"gamma = $(dc.gamma)")
    println(io,"for $(dc.k) steps")
    println(io,"tolerence: larger than $(dc.tol) zeros")
end

function run!{N}(dc::DemonCooling{N}, aqc::AQC)
    ret=cooling!(aqc.state, dc)
    run_t = 0.0
    if ret==false
        aqc.state = BellState(N)
        run_t = run!(aqc)
        # aqc.process_list = append!(aqc.process_list[1:aqc.status],aqc.process_list)
        # aqc.state = BellState(N)
    end
    return run_t+t(dc)
end

# return a list of probability for each 
# reference: README
function probability{V<:AbstractVector,T,N}(
    state::ComputState{V,T,N},
    dc::DemonCooling{N}
    )
    ret = zeros(Float64,dc.k+1)
    A = dc.spectrum
    map!(x->0.5*sin(x*dc.t-dc.gamma), A)
    for j = 0:dc.k
        for i = 1:2^N
            ret[j+1] += binomial(dc.k,j)*abs2(state.coeffs[i])*(0.5-A[i])^j*(0.5+A[i])^(dc.k-j)
        end
    end
    return ret  
end

# return a cooled state
function cooling!{V<:AbstractVector,T,N}(
    state::ComputState{V,T,N},
    dc::DemonCooling{N}
    )
    p = probability(state, dc)
    k = steps(dc)
    A = dc.spectrum
    map!(x->0.5*sin(x*dc.t-dc.gamma), A)

    dice = rand()
    if dice<sum(p[1:tol(dc)])
        return false
    end
    for j=tol(dc)+1:k
        if dice<sum(p[1:j+1])
            state.coeffs = 0.5^k*binomial(k,j)*(I-im*exp(im*dc.gamma)*expm(-im*full(dc.H)*dc.t))^j*(I+im*exp(im*dc.gamma)*expm(-im*full(dc.H)*dc.t))^(k-j)*state.coeffs
            normalize!(state.coeffs)
            return true            
        end
    end
    return false
end