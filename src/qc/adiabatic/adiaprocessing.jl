import Base:show
export AbstractPropMethod, AdiabaticProcessing, QuExpmV, QuExpokit, timer, hamiltonian

abstract AbstractPropMethod

immutable QuExpokit <: AbstractPropMethod
    options::Dict{Symbol, Any}
end

QuExpokit() = QuExpokit(Dict())

function (prob::QuExpokit)(scalar::Real, A::AbstractMatrix, v::Vector)
    return Expokit.expmv(scalar,A,v,m=get(prob.options, :m, 30), tol=get(prob.options, :tol, 1e-7))
end

type QuExpmV <: AbstractPropMethod
    options::Dict{Symbol, Any}
end

QuExpmV() = QuExpmV(Dict())

function (prob::QuExpmV)(scalar::Real, A::AbstractMatrix, v::Vector)
    return ExpmV.expmv(scalar,A,v,M = get(prob.options, :M, []), prec = get(prob.options, :prec, "double"), shift = get(prob.options, :shift, false), full_term = get(prob.options, :full_term, false))
end

immutable AdiabaticProcessing{H<:AbstractMatrix, S<:AbstractPropMethod, N}<:AQCprocessor{N}
    HB::H
    HP::H
    T::AbstractFloat
    prop::S
    dt::Real
    tlist::Tuple

    function AdiabaticProcessing(
        HP::H,
        method::S,
        maxtime::Real, 
        dt::Real,
        tlist::Tuple
        )
        HB = bHamiltonian(N)

        # Similar Matrix
        P = eigenvector(0,N)
        for i=1:2^N-1
            P = [P eigenvector(i,N)]
        end
        invP = inv(P)
        HB = P*HB*invP        
        new(HB,HP,maxtime,method,dt,tlist)
    end
end

AdiabaticProcessing{H<:AbstractMatrix,S<:AbstractPropMethod}(HP::H,method::S;dt=1e-2,tlist=(0.0,1.0),maxtime=1)=AdiabaticProcessing{H,S,Int(log2(size(HP)[1]))}(HP,method,maxtime,dt,tlist)

AdiabaticProcessing(HP::AbstractMatrix) = AdiabaticProcessing(HP,QuExpmV())

function show{H,S,N}(io::IO,ap::AdiabaticProcessing{H,S,N})
    println(io,"------------------------------")
    println(io,"$N qubits Adiabatic Processing")
end

# function print{H,S,N}(ap::AdiabaticProcessing{H,S,N})
#     println("------------------------------")
#     println("$N qubits Adiabatic Processing")
# end

function hamiltonian(ap::AdiabaticProcessing, s::Real)
    return (1-s)*ap.HB+s*ap.HP
end

# return processor output states and time cost
function run!(ap::AdiabaticProcessing, aqc::AQC)
    for t in ap.tlist[1]*ap.T:ap.dt:ap.tlist[2]*ap.T
        s = t/ap.T
        H = hamiltonian(ap, s)
        aqc.state.coeffs = ap.prop(ap.dt, -im*H, coeff(aqc.state))
    end
    aqc.sys_h = hamiltonian(ap, ap.tlist[2])
    return (ap.tlist[2]-ap.tlist[1])*ap.T
end

function timer(ap::AdiabaticProcessing)
    return ap.T
end