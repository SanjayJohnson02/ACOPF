using Pkg
using PowerModels

my_data = PowerModels.parse_file("pglib_opf_case14_ieee.m")

#Pkg.add(url = "https://github.com/exanauts/ExaModelsPower.jl")

using ExaModelsPower
em, emvars = ExaModelsPower.opf_model("pglib_opf_case14_ieee.m")

#Pkg.add("ExaModels")
#Pkg.add("NLPModelsIpopt")
using ExaModels
using NLPModelsIpopt
result2 = ipopt(em)

data, dicts = ExaModelsPower.parse_ac_power_data("pglib_opf_case14_ieee.m")
println(data[4])