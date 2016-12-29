# Graph States

using LightGraphs
import Base: convert
export GraphState, stabilizers, isLCeq

type GraphState
    graph::Graph
end

# operators

function stabilizers(gs::GraphState)
    res = [sigmax for i in vertices(gs.graph)]
    for node in vertices(gs.graph)
        for itr in vertices(gs.graph)
            if itr in neighbors(gs.graph,node)
                res[node] = kron(res[node],sigmaz)
            else
                res[node] = kron(res[node],sigmao)
            end
        end
    end
    return res
end

function isLCeq(gsa::GraphState,gsb::GraphState)
    @assert nv(gsa) == nv(gsb) "Number of nodes does not match"
    n = nv(gsa)
    generator_matrix_left = [speye(n) adjacency_matrix(gsa)]
    generator_matrix_right = [speye(n) adjacency_matrix(gsb)]
    P = [spzeros(n,n) speye(n);speye(n) spzeros(n,n)]
    flag = generator_matrix_left*P*generator_matrix_right.'
    return isempty(nonzeros(flag))
end

function convert(::Type{ComputState},gs::GraphState)
    return eigs(gs,nev=1,which=:SR)[2][:,1]|>ComputState
end
