using JSON
export consts

const consts = Dict()

#import NIST data set
const NIST_CODATA_FUNDAMENTAL_PHYSICS_CONSTS = "srd121_allascii_2014.json" 
const NIST_consts = JSON.parsefile("src/utils/data/"*NIST_CODATA_FUNDAMENTAL_PHYSICS_CONSTS)["constant"]
for each_const in NIST_consts
    each_const["Value"] = replace(each_const["Value"],r"( )[^e]","_")
    each_const["Value"] = replace(each_const["Value"]," ","")
    each_const["Value"] = replace(each_const["Value"],".","")
    each_const["Value"] = AbstractFloat(parse(each_const["Value"]))
    if each_const["Uncertainty"]!="(exact)"
        each_const["Uncertainty"] = replace(each_const["Uncertainty"],r"( )[^e]","_")
        each_const["Uncertainty"] = replace(each_const["Uncertainty"]," ","")
        each_const["Uncertainty"] = AbstractFloat(parse(each_const["Uncertainty"]))
    else
        each_const["Uncertainty"] = :exact
    end
end


#####################
# Consts Definition #
#####################
immutable PhyConst{T<:Number}
    value::T 
    uncertainty::Union{T,Symbol} 
    unit::AbstractString
end
# constructors
PhyConst{T<:Number}(value::T, uncertainty::T, unit::AbstractString) = PhyConst{T}(value, uncertainty, unit)

# utils
value(pc::PhyConst) = pc.value
uncertainty(pc::PhyConst) = pc.uncertainty
unit(pc::PhyConst) = pc.unit


# choose dataset
function dataset(name::AbstractString)
    if name=="NIST_CODATA"
        for each_const_elem in NIST_consts
            consts[each_const_elem["Quantity "]] = PhyConst(
                each_const_elem["Value"],
                each_const_elem["Uncertainty"],
                each_const_elem["Unit"]
                )
        end
    else
        throw(ArgumentError("Input dataset: $(name) does exist"))
    end
end

dataset("NIST_CODATA")

# Planck consts
const h = value(consts["Planck constant"])
const hbar = value(consts["Planck constant over 2 pi"])
const ħ = hbar
export h, hbar, ħ

# Pauli Groups
const sigmax = sparse(Complex[0.0 1.0;1.0 0.0])
const sigmay = sparse(Complex[0.0 -1.0im;1.0im 0.0])
const sigmaz = sparse(Complex[1.0 0.0; 0.0 -1.0])
const sigmao = sparse(Complex[1.0 0.0;0.0 1.0])
export sigmax, sigmay,sigmaz, sigmao

const σ₀ = sigmao
const σ₁ = sigmax
const σ₂ = sigmay
const σ₃ = sigmaz
export σ₀,σ₁,σ₂,σ₃