function eigenvector(index::Integer,n::Integer)
    res = 1
    for i in bin(index,n)
        if i=='0'
            res = kron(res,1/sqrt(2)*[1,1])
        elseif i=='1'
            res = kron(res,1/sqrt(2)*[1,-1])
        end
    end
    return res
end