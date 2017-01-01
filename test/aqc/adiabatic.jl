
# adiabatic processing test
n=4; T = 3
# ins,ans = generate(n)
# pH = pHamiltonian(ins,n)
# aqc = AQC(n)
# aqc|>AdiabaticProcessing(pH,QuExpmV();dt=1e-2,maxtime=T)
# @test run!(aqc) == T
# @test round(norm(coeff(aqc.state)[ans[1]+1])^2,3) == 0.130

# cooling-assisted processing test
ins = generate(n)
pH = pHamiltonian(ins,n)
eigen = eigfact(full(pH))
val,idx = findmin(eigen[:values])
maxEigen = maximum(eigen[:values])

aqc = AQC(n)
cooled_aqc = AQC(n)

ap1 = AdiabaticProcessing(pH,QuExpmV();dt=1e-2,maxtime=T,tlist=(0.0,1/3))
ap2 = AdiabaticProcessing(pH,QuExpmV();dt=1e-2,maxtime=T,tlist=(1/3,2/3))
ap3 = AdiabaticProcessing(pH,QuExpmV();dt=1e-2,maxtime=T,tlist=(2/3,1.0))
cooled_aqc|>ap1|>DemonCooling(hamiltonian(ap1,ap1.tlist[2]); k=10, tol=3, optimize=:off)
cooled_aqc|>ap2|>DemonCooling(hamiltonian(ap2,ap2.tlist[2]); k=10, tol=3, optimize=:off)
cooled_aqc|>ap3

aqc|>ap1|>ap2|>ap3
@show run!(cooled_aqc)
@show abs(dot(coeff(cooled_aqc.state),eigen[:vectors][:,idx]))
# @show run!(aqc)
# data = @parallel (+) for i=1:1
#   cooled_aqc.state = BellState(n)
#   run_time = run!(cooled_aqc)
#   suc_prob = abs(dot(coeff(cooled_aqc.state),eigen[:vectors][:,idx]))
#   Float64[run_time,suc_prob]
# end

# @show data