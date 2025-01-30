using Pkg
Pkg.add(["JuMP", "PowerModels", "Interpolations", "ExaModels", "Ipopt", "NLPModelsIpopt", "LinearAlgebra", "CUDA", "MadNLPGPU", "MadNLPHSL"])
Pkg.develop(PackageSpec(path="acopf_scratch/HSL_jll.jl.v2025.1.19/"))
using JuMP, PowerModels, Interpolations, ExaModels, Ipopt, NLPModelsIpopt, LinearAlgebra, CUDA, MadNLPGPU, HSL_jll, MadNLPHSL
