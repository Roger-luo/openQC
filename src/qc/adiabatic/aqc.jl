import Base:run, |>, show, start, next, done
export AQC, state, run!, |>, clear_process!, copy_process, timer, print, reset!,success_probability

type AQC{A<:AbstractVector,T,N}<:AbstractQC{N}
    state::ComputState{A,T,N}
    status::Int # current processor index
    sys_h::AbstractMatrix # system hamiltonian
    process_list::AbstractArray
end

AQC(n::Int) = AQC{Array{Complex128,1},Complex128,n}(BellState(n), 0, spzeros(2^n,2^n), Array(AQCprocessor,0))
state(aqc::AQC) = aqc.state
hamiltonian(aqc::AQC) = aqc.sys_h


# TODO: aqc processing type for iterations
function start(aqc::AQC)
    return 1
end

function next(aqc::AQC, state)
    return aqc.process_list[state], state+1
end

function done(aqc::AQC, state)
    return state > length(aqc.process_list)
end

function run!(aqc::AQC)
    run_timer = 0.0
    for each_processor in aqc
        aqc.status += 1
        # @show each_processor
        run_timer += run!(each_processor, aqc)
    end
    return run_timer
end

function show(io::IO,aqc::AQC)
    for each_processor in aqc.process_list
        show(io::IO,each_processor)
    end
end

function (|>)(aqc::AQC,aqcprocessor::AQCprocessor)
    push!(aqc.process_list, aqcprocessor)
    return aqc
end

function clear_process!(aqc::AQC)
    aqc.process_list = Array(AQCprocessor,0)
end

function copy_process(aqc::AQC)
    return deepcopy(aqc.process_list)
end

function timer(aqc::AQC)
    ret = 0.0
    for i=1:aqc.status
        ret += timer(aqc.process_list[i])
    end
    return ret 
end

function reset!{A,T,N}(aqc::AQC{A,T,N})
    aqc.state = BellState(Int(N))
end

function success_probability(aqc::AQC,tg_ground_state::AbstractVector)
    return abs2(dot(coeff(aqc.state),tg_ground_state))
end

function success_probability(aqc::AQC,pH::AbstractMatrix)
    eigen = eigfact(full(pH))
    val,idx = findmin(eigen[:values])
    return abs2(dot(coeff(aqc.state),eigen[:vectors][:,idx]))
end