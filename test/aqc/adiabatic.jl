using QuSAT

# adiabatic processing test
n=4; T = 3
ins,ans = generate(n)
pH = pHamiltonian(ins,n)
aqc = AQC(n)
aqc|>AdiabaticProcessing(pH,QuExpmV();dt=1e-2,maxtime=T)
@test run!(aqc) == T
@test round(norm(coeff(aqc.state)[ans[1]+1])^2,3) == 0.130

