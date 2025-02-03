@kernel function _transfer_to_map!(dest, to_map, src)
    k = @index(Global, Linear)
    @inbounds begin
        Atomix.@atomic dest[to_map[k]] += src[k]
    end
end

function MadNLP.transfer!(
    dest::CUSPARSE.CuSparseMatrixCSC{Tv},
    src::MadNLP.SparseMatrixCOO{Tv},
    map::CuVector{Int},
) where {Tv}
    fill!(nonzeros(dest), zero(Tv))
    _transfer_to_map!(CUDABackend())(nonzeros(dest), map, src.V; ndrange=length(map))
    synchronize(CUDABackend())
    return
end

@kernel function _remove_diagonal_kernel!(y, Ap, Ai, Az, x, alpha)
    j = @index(Global, Linear)
    @inbounds for k in Ap[j]:Ap[j+1]-1
        if Ai[k] == j
            y[j] -= alpha * Az[k] * x[j]
            break
        end
    end
end

function MadNCL.symul!(
    y::CuVector{T},
    A::CUSPARSE.CuSparseMatrixCSC{T},
    x::CuVector{T},
    alpha::Number,
    beta::Number,
) where {T}
    m, n = size(A)
    mul!(y, A , x, alpha, beta)
    mul!(y, A', x, alpha, one(T))
    _remove_diagonal_kernel!(CUDABackend())(y, A.colPtr, A.rowVal, A.nzVal, x, alpha; ndrange=n)
    synchronize(CUDABackend())
    return y
end

function MadNLP.build_condensed_aug_coord!(
    kkt::MadNCL.K1sAuglagKKTSystem{T,VT,MT},
) where {T,VT,MT<:CUSPARSE.CuSparseMatrixCSC{T}}
    fill!(kkt.aug_com.nzVal, zero(T))
    if length(kkt.hptr) > 0
        MadNLPGPU._transfer_hessian_kernel!(CUDABackend())(
            kkt.aug_com.nzVal,
            kkt.hptr,
            kkt.hess_com.nzVal;
            ndrange = length(kkt.hptr),
        )
        synchronize(CUDABackend())
    end
    if length(kkt.dptr) > 0
        MadNLPGPU._transfer_hessian_kernel!(CUDABackend())(
            kkt.aug_com.nzVal,
            kkt.dptr,
            kkt.pr_diag_reduced;
            ndrange = length(kkt.dptr),
        )
        synchronize(CUDABackend())
    end
    if length(kkt.ext.jptrptr) > 1 # otherwise error is thrown
        MadNLPGPU._transfer_jtsj_kernel!(CUDABackend())(
            kkt.aug_com.nzVal,
            kkt.jptr,
            kkt.ext.jptrptr,
            kkt.jt_csc.nzVal,
            kkt.diag_buffer;
            ndrange = length(kkt.ext.jptrptr) - 1,
        )
        synchronize(CUDABackend())
    end
    return
end

