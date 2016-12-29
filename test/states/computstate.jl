t = ComputState([1,2,3,4])
@test rawcoeffs(t) == [1,2,3,4]
@test coeff(t) == [1,2,3,4]