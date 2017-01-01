###################################################################
#
# This file generates SAT's Instance in classical method for tests
# in quantum computing.
# Ref: arxiv:quant-ph/0007071v1
#
###################################################################


import Base: show, call
import LightGraphs.save
export AbstractBits, AbstractClause, Bits,
        Clause, ECClause, TruthTable, save,
        generate, Instance, AssignNum, readins


abstract AbstractBits
abstract AbstractClause{N}

immutable Bits <: AbstractBits
    value::UInt32

    Bits(n::Integer) = new(n)
end

function (b::Bits)(n::Integer)
    @assert n>0
    return (b.value>>(n-1))&1
end

# function convert(Bits,Integer)
#     return Bits()

"""
TruthTable Example:

|bits value | truth value|
|-----------|------------|
| 000       | 0          |
| 001       | 1          |
| 010       | 0          |
| 011       | 1          |
| 100       | 0          |
| 101       | 1          |
| 110       | 0          |
| 111       | 0          |

TruthTable(0b00101010)
"""
immutable TruthTable
    value::UInt32
end

function (truthtable::TruthTable)(index::Integer)
    return (truthtable.value>>index)&1
end

immutable Clause{N} <: AbstractClause{N}
    table::TruthTable
    ids::AbstractVector{Integer}

    function Clause(table::TruthTable,ids::AbstractVector{Integer})
        @assert length(ids)==N
        sort!(ids)
        new(table,ids)
    end
end

Clause(table::TruthTable,ids::Integer...) = Clause{length(ids)}(table,[ids...])
Clause(num::Integer,table::TruthTable,ids::Integer...) = Clause{num}(table,[ids...])
Clause(num::Integer,table::Integer,ids::Integer...) = Clause{num}(TruthTable(table),[ids...])

function (clause::Clause{N}){N}(assign::Integer)
    return clause.table(assign)
end

################################################################################
#
# Clauses in exact cover problem
#
################################################################################

type ECClause{N} <: AbstractClause{N}
    ids::Vector{Int}

    function ECClause(ids::Vector{Int})
        @assert length(ids) == N
        sort!(ids)
        new(ids)
    end
end

ECClause(ids::Vector{Int}) = ECClause{length(ids)}(ids)
ECClause(ids::Integer...) = ECClause{length(ids)}([ids...])

function (c::ECClause{N}){N}(assign::Integer)
    res = 0
    for i = 1:N
        res += assign&1
        assign = assign>>1
    end
    return Int(res==1)
end

function save{N}(io::IO,clause::ECClause{N})
    write(io,"$(clause.ids[1])")
    for id in clause.ids[2:end]
        write(io,"\t")
        write(io,"$(id)")
    end
    write(io,"\n")
end

"""
Instance is a collection of M clauses:

~~~TeX
C_1 Λ C_2 Λ... C_M
~~~

construct an instance:

- `num` is the number of bits
- `clauses` is the collection of clause
"""
type Instance{M,N}
    c::AbstractVector{AbstractClause{N}}
    ans::AbstractVector{Int}
end

Instance{N}(clauses::AbstractVector{ECClause{N}},ans::AbstractVector{Int}) = Instance{length(clauses),N}(clauses,ans)
Instance{N}(clauses::AbstractVector{ECClause{N}},ans::Int) = Instance{length(clauses),N}(clauses,[ans])
Instance{N}(num::Integer,clauses::AbstractVector{ECClause{N}}) = Instance{num,N}(clauses,Array(Int,0))
Instance{N}(num::Integer,clause::ECClause{N},clauses::ECClause{N}...) = Instance(num,[clause,clauses...])
Instance{N}(num::Integer,clauses::AbstractVector{Clause{N}}) = Instance{num,N}(clauses,Array(Int,0))
Instance{N}(num::Integer,clause::Clause{N},clauses::Clause{N}...) = Instance(num,[clause,clauses...])
Instance{N}(clauses::AbstractVector{ECClause{N}}) = Instance{length(clauses),N}(clauses)
Instance{N}(clause::Clause{N},clauses::Clause{N}...) = Instance([clause,clauses...])

function (clauses::Instance{M,N}){M,N}(assign::Bits)
    res = 1
    for clause in clauses.c
        assignment = 0
        digit = 0
        for id in clause.ids
            assignment += assign(id)<<digit
            digit += 1
        end
        res &= clause(assignment)
    end
    return res
end

function show{M,N}(io::IO,ins::Instance{M,N})
    for clause in ins.c[1:end-1]
        print(io,clause.ids)
        print(io,"Λ")
    end
    print(io,ins.c[end].ids)
    if !isempty(ins.ans)
        print(io,"\tanswer is :")
        print(io, ins.ans)
    end
    print(io,"\n")
end

show{M,N}(ins::Instance{M,N}) = show(STDOUT,ins)

function save{M,N}(io::IO,ins::Instance{M,N})
    write(io,"$(ins.ans[1])")
    for ans in ins.ans[2:end]
        write(io,"$(ans)")
        write(io,"\t")
    end
    write(io,"\n")
    for clause in ins.c
        save(io,clause)
    end
end

function readins(io::IO)
    data = convert(Array{Int,2},readdlm(io))
    clauses = ECClause[ECClause(vec(data[i,:])) for i=1:size(data)[1]]
    @show typeof(clauses[1])
    return Instance(size(data)[2],clauses)
end

function AssignNum{M,N}(instance::Instance{M,N})
    count = 0
    ans = Int[]
    for i = 0:2^M-1
        if instance(Bits(i)) == 1
            count+=1
            push!(ans,i)
        end
    end
    return count,ans
end

import Base: push!

function push!{M,N}(instance::Instance{M,N},clause::ECClause{N})
    push!(instance.c,clause)
end

################################################################################
#
# generate instance with only one assignment satisfied
#
################################################################################

function engine(n::Integer,N::Integer,maxtry=100)
    list = collect(1:n)
    ids = (list|>shuffle)[1:N]

    clauses = [ECClause(ids)]
    ins = Instance(n,clauses)

    for i = 1:maxtry
        pre_ids = ids
        ids = (list|>shuffle)[1:N]
        if sort!(ids) != sort!(pre_ids)
            push!(ins,ECClause(ids))
        end

        asign = AssignNum(ins)
        if asign[1]==1
            append!(ins.ans,asign[2])
            return ins
        end
    end

    return nothing
end

function generate_ins(n::Integer, maxiter::Int, maxbits::Int, N::Int)
    for i=1:maxiter
        ins = engine(n,N)
        t === nothing || return ins
    end
    warn("may not have an assignment\n")
    return nothing
end

function generate(n::Integer;maxiter=1000,maxbits=16,N=3)
    ins_dir = pkg_dir*"/src/utils/data/"*"$(N)_sat_ins"
    if isdir(ins_dir)
        try
            open(ins_dir*"/"*string(n)*".ins","r") do io
                ans = vec(readdlm(IOBuffer(readline(io)),Int))
                data = readdlm(io,Int)
                clauses = ECClause[ECClause(vec(data[i,:])) for i=1:size(data)[1]]
                return Instance{length(clauses),N}(clauses,ans)
            end
        catch
            return generate(ins_dir*"/"*string(n)*".ins", n;maxiter=maxiter,maxbits=maxbits,N=N)
        end
    else
        mkdir(ins_dir)
        return generate(ins_dir*"/"*string(n)*".ins",n,maxiter=maxiter, maxbits=maxbits, N=N)
    end
end

function generate(file::AbstractString, n::Integer; maxiter=1000, maxbits=16, N=3)
    open(file,"w") do io
        ins = generate_ins(n,maxiter,maxbits,N)
        save(io, ins)
        return ins
    end
end

function generate(rng::Range; maxiter=1000, maxbits=16, N=3)
    ins_dir = pkg_dir*"/src/utils/data/"*"$(N)_sat_ins"
    if !isdir(ins_dir)
        mkdir(ins_dir)
        for i in rng
            generate(ins_dir*"/"*string(i)*".ins",i,maxiter=maxiter, maxbits=maxbits, N=N)
        end
    else
        ins_list = map(x->parse(split(x,".")[1]),readdir(ins_dir))
        for i in rng
            if !(i in ins_list)
                generate(ins_dir*"/"*string(i)*".ins",i,maxiter=maxiter, maxbits=maxbits, N=N)
            end
        end
    end
end