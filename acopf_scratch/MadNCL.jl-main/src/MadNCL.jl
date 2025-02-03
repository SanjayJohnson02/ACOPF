module MadNCL

import LinearAlgebra: norm
import LinearAlgebra: dot, axpy!, mul!
import LinearAlgebra: Symmetric
import SparseArrays: sparsevec, nonzeros
import SparseArrays: SparseMatrixCSC
import Printf: @printf
import NLPModels
import MadNLP
import MadNLP: AbstractExecutionStats, getStatus
import Polynomials: Polynomial, roots

# CUDA related deps
import Atomix
import MadNLPGPU
import CUDA: CUSPARSE, CuVector, CUDABackend
import KernelAbstractions: @kernel, @index, synchronize

include("utils.jl")
include("Models/ncl.jl")
include("Models/implicit_ncl.jl")
include("Models/scaled.jl")

include("Algorithms/ncl.jl")

include("KKT/k1s.jl")
include("KKT/k2r.jl")

include("cuda_wrapper.jl")

end
