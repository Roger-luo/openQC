n=4; T = 1

# Cooling Algorithm Test
P = rand(2^n,2^n)
test_cooling_matrix = P*diagm(collect(1:2^n))*inv(P)
eigen = eigfact(test_cooling_matrix)
val,idx = findmin(eigen[:values])
maxEigen = maximum(eigen[:values])
cooling_state = BellState(n)

while !cooling!(cooling_state, DemonCooling(test_cooling_matrix; k=1, tol=0))
    cooling_state = BellState(n)
end

test_cooling_state = BellState(n)

for i=1:1
test_cooling_state.coeffs = normalize((I-im*expm(-im*test_cooling_matrix*π/(2*maxEigen)))*test_cooling_state.coeffs)
end
@test round(abs(dot(coeff(test_cooling_state),eigen[:vectors][:,idx])),5) == round(abs(dot(coeff(cooling_state),eigen[:vectors][:,idx])),5)



# @show round(norm(coeff(aqc.state)[ans[1]+1])^2,3)
# @show abs(dot(coeff(aqc.state),eigen[:vectors][:,idx]))
# @show abs(dot(coeff(cooled_aqc.state),eigen[:vectors][:,idx]))

# show(cooled_aqc)
# print(cooled_aqc)




# test_cooling_state = BellState(n)

# for i=1:200
# test_cooling_state.coeffs = normalize(test_cooling_state.coeffs-im*expm(-im*test_cooling_matrix*0.522697)*test_cooling_state.coeffs)
# end
# @show abs(dot(coeff(test_cooling_state),eigen[:vectors][:,idx]))


# test_cooling_state = BellState(n)

# for i=1:200
# test_cooling_state.coeffs = normalize(test_cooling_state.coeffs-im*exp(im*1.5368)*expm(-im*test_cooling_matrix*2.06334)*test_cooling_state.coeffs)
# end
# @show abs(dot(coeff(test_cooling_state),eigen[:vectors][:,idx]))