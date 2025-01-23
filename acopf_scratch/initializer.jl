using Pkg
Pkg.add(["JuMP", "PowerModels", "Interpolations", "ExaModels", "Ipopt", "NLPModelsIpopt", "LinearAlgebra", "CUDA", "MadNLPGPU"])
using JuMP, PowerModels, Interpolations, ExaModels, Ipopt, NLPModelsIpopt, LinearAlgebra, CUDA, MadNLPGPU
import HSL_jll