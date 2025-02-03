using Pkg
Pkg.add(["JuMP", "PowerModels", "Interpolations", "ExaModels", "Ipopt", "NLPModelsIpopt", "LinearAlgebra", "CUDA", "MadNLPGPU", "MadNLPHSL"])
Pkg.develop(path = "MadNCL.jl-main")
Pkg.develop(PackageSpec(path="HSL_jll.jl.v2025.1.19/"))
using JuMP, PowerModels, Interpolations, ExaModels, Ipopt, NLPModelsIpopt, LinearAlgebra, CUDA, MadNLPGPU, HSL_jll, MadNLPHSL
