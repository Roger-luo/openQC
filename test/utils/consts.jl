check_list = [
    :h,
    :hbar,
    :ħ,
    :sigmao,
    :sigmax,
    :sigmay,
    :sigmaz,
    :σ₀,
    :σ₁,
    :σ₂,
    :σ₃
]

for each_symbol in check_list
    @test isdefined(each_symbol)
end