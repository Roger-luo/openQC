using JuMP, Ipopt
import Base:show,gamma
export DemonCooling, t, gamma, spec, steps, tol, hamiltonian, timer, probability, cooling!
const MAX_COOLING_FAILURE_COUNT = 200
cooling_module_id = 1

type DemonCooling{N}<:AQCprocessor{N}
    t::Real
    gamma::Real
    H::AbstractMatrix
    k::Int # k steps
    tol::Int # tolerance for error correction, at least $tol zeros
    optimize::Int
    failure_count::Int
end

t(dc::DemonCooling) = dc.t
gamma(dc::DemonCooling) = dc.gamma
spec(dc::DemonCooling) = dc.spectrum
steps(dc::DemonCooling) = dc.k
tol(dc::DemonCooling) = dc.tol
timer(dc::DemonCooling) = dc.t

function DemonCooling(H::AbstractSparseMatrix; k::Int=5, tol::Int=2, optimize=:off)
    if optimize == :static
        t,gamma_ = static_opt_parameters(H)
        return DemonCooling{Int(log2(size(H)[1]))}(t,gamma_,H,k,tol,0,0)
    elseif optimize == :off
        maxEigen = eigs(H;nev=1,which=:LR)[1][1]|>real
        t = π/(2*maxEigen)
        gamma_ = 0
        return DemonCooling{Int(log2(size(H)[1]))}(t,gamma_,H,k,tol,1,0)
    elseif optimize == :runtime
        return DemonCooling{Int(log2(size(H)[1]))}(0.0,0.0,H,k,tol,2,0)
    end
end

function static_opt_parameters(H::AbstractSparseMatrix)
    eigens = eigfact(full(H))
    spec = real(eigens[:values])
    val,idx = findmin(spec)
    m = Model()
    @variable(m, t)
    @variable(m, gamma_)
    @NLobjective(m, Min, exp(sum(1-2*sin(e*t-gamma_) for e in spec)/(1-2*sin(val*t-gamma_))))
    solve(m)
    return getvalue(t),getvalue(gamma_)
end 

function multi_cooling_obj(t,gamma_,e,k,j,p)
    return binomial(k,j)*p*(0.5-sin(e*t-gamma_)^j*(0.5+sin(e*t-gamma_))^(k-j))
end

function ∇multi_cooling_obj(g,t,gamma_,e,k,j,p)
    phi = e*t-gamma_
    g[1] = e*(k-j)*cos(phi)*(0.5-sin(phi))^j*(0.5+sin(phi))^(k-j-1)-e*j*cos(phi)*(0.5-sin(phi))^(j-1)*(0.5+sin(phi))^(k-j)

    g[2] = (j-k)*cos(phi)*(0.5-sin(phi))^j*(0.5+sin(phi))^(k-j-1)+j*cos(phi)*(0.5-sin(phi))^(j-1)*(0.5+sin(phi))^(-j+k)
end

function cooling_obj_register(H::AbstractSparseMatrix,k::Int, tol::Int, state::AbstractVector)
    eigens = eigfact(full(H))
    spec = real(eigens[:values])
    val,idx = findmin(spec)
    vecs = eigens[:vectors]
    p = abs2(state.'*eigens[:vectors])
    ground = vecs[:,idx]
    obj(x,y) = sum(sum(
            multi_cooling_obj(x,y,spec[i],k,j,p[i])*multi_cooling_obj(x,y,spec[idx],k,j,p[idx]) 
            for i in eachindex(spec)) for j =tol:k)
    ∇obj(g,x,y) = sum(sum(
                        ∇multi_cooling_obj(g,x,y,spec[i],k,j,p[i])*multi_cooling_obj(x,y,spec[idx],k,j,p[idx])+∇multi_cooling_obj(g,x,y,spec[idx],k,j,p[idx])*multi_cooling_obj(x,y,spec[i],k,j,p[i])
                        for i in eachindex(spec)) for j =tol:k)
    cooling_obj_name = Symbol("cool_obj$(cooling_module_id)")
    global cooling_module_id += 1
    JuMP.register(cooling_obj_name,2,obj,∇obj)
    return cooling_obj_name
end

function runtime_opt_parameters(dc::DemonCooling,aqc::AQC)
    cooling_obj_name = cooling_obj_register(dc.H,dc.k,dc.tol,coeff(aqc.state))
    m = Model()
    spec = real(eigvals(full(dc.H)))
    minEigen, maxEigen = minimum(spec), maximum(spec)
    @variable(m, t_)
    @variable(m, gamma_)
    @constraint(m, -pi<=minEigen*t_-gamma_<=pi)
    @constraint(m, -pi<=maxEigen*t_-gamma_<=pi)
    @NLobjective(m, Max, $(cooling_obj_name)(t_,gamma_))
    solve(m)
    return getvalue(t_),getvalue(gamma_)
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
    if dc.optimize == 2
        dc.t,dc.gamma = runtime_opt_parameters(dc,aqc)
    end
    ret=cooling!(aqc.state, dc)
    run_t = 0.0
    if ret==false
        dc.failure_count += 1
        aqc.state = BellState(N)
        run_t = run!(aqc)
        # aqc.process_list = append!(aqc.process_list[1:aqc.status],aqc.process_list)
        # aqc.state = BellState(N)
        return run_t+t(dc)
    elseif ret
        return t(dc)
    else
    end
end

# return a list of probability for each 
# reference: README
function probability{V<:AbstractVector,T,N}(
    state::ComputState{V,T,N},
    dc::DemonCooling{N}
    )
    ret = zeros(Float64,dc.k+1)
    eigens = eigfact(full(dc.H))
    spec = real(eigens[:values])
    A = spec
    map!(x->0.5*sin(x*dc.t-dc.gamma), A)
    p = coeff(state).'*eigens[:vectors]
    for j = 0:dc.k
        for i = 1:2^N
            ret[j+1] += binomial(dc.k,j)*abs2(p[i])*(0.5-A[i])^j*(0.5+A[i])^(dc.k-j)
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
    eigens = eigfact(full(dc.H))
    spec = real(eigens[:values])
    A = spec
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