#=
    ImplicitNCLModel.

    Implement magic steps for NCL augmented Lagrangian.

=#

struct ImplicitNCLModel{T, VT, VI, M} <: NLPModels.AbstractNLPModel{T, VT}
    nlp::M
    n::Int
    m::Int
    meta::NLPModels.NLPModelMeta{T, VT}
    counters::NLPModels.Counters
    yk::VT   # size [m]
    rk::VT   # size [m]
    ck::VT   # size [m]
    σk::VT   # size [m]
    hashx::Ref{UInt64}
    ρk::Ref{T}
    μk::Ref{T}
    ind_eq::VI
    ind_lb::VI
    ind_ub::VI
    ind_rng::VI
    # (Transposed) Jacobian J [CSR format]
    JTp::VI  # size [n]
    JTj::VI  # size [nnzj]
    JTx::VT  # size [nnzj]
    map_jac::VI
    buffer_jac::VT  # size [nnzj]
    # L = JᵀDJ  [CSR format]
    Lp::VI   # size [n]
    Lj::VI   # size [nnzl]
    Lx::VT   # size [nnzl]
    Dx::VT   # size [m]
end

function ImplicitNCLModel(
    nlp::NLPModels.AbstractNLPModel{T, VT};
    resid::T=zero(T),
    ρ::T=T(10),
    μ::T=T(0.1),
) where {T, VT}
    nlp.meta.minimize || error("only minimization problems are currently supported")

    n = NLPModels.get_nvar(nlp)
    m = NLPModels.get_ncon(nlp)

    y = fill!(similar(VT, m), one(T))
    r = fill!(similar(VT, m), one(T))
    c = fill!(similar(VT, m), one(T))
    σ = fill!(similar(VT, m), one(T))

    lvarx = NLPModels.get_lvar(nlp)
    uvarx = NLPModels.get_uvar(nlp)

    x0 = NLPModels.get_x0(nlp)
    hashx = UInt64(0)

    ind_eq = findall(nlp.meta.lcon .== nlp.meta.ucon)
    ind_low = findall((nlp.meta.lcon .> -Inf) .&& (nlp.meta.ucon .== Inf))
    ind_upp = findall((nlp.meta.lcon .== -Inf) .&& (nlp.meta.ucon .< Inf))
    ind_rng = findall(-Inf .< nlp.meta.lcon .< nlp.meta.ucon .< Inf)

    # Build condensed Jacobian
    nnzj = NLPModels.get_nnzj(nlp)
    Ji = zeros(Int, nnzj)
    Jj = zeros(Int, nnzj)
    NLPModels.jac_structure!(nlp, Ji, Jj)
    Jx = zeros(nnzj)
    # Build mapping
    Jx .= 1:nnzj

    buffer_jac = zeros(nnzj)

    # Build CSR representation of transposed Jacobian
    Jtp, Jtj, Jtx = coo_to_csr(n, m, Jj, Ji, Jx)
    map_jac = Jtx

    # Build CSR representation of Jᵀ * J
    Cp, Cj = build_symbolic_analysis_jtj(n, m, Jtp, Jtj)
    nnz_jtj = Cp[end] - 1
    Cx = zeros(nnz_jtj)
    Dx = zeros(m)

    meta = NLPModels.NLPModelMeta{T, VT}(
        n;
        lvar=lvarx,
        uvar=uvarx,
        x0=x0,
        nnzj=0,
        nnzh=NLPModels.get_nnzh(nlp) + nnz_jtj,
        ncon=0,
        minimize=true,
    )
    cnt = NLPModels.Counters()

    return ImplicitNCLModel{T, VT, Vector{Int}, typeof(nlp)}(
        nlp, n, m, meta, cnt, y, r, c, σ, hashx, ρ, μ,
        ind_eq, ind_low, ind_upp, ind_rng,
        Jtp, Jtj, Jtx, map_jac, buffer_jac, Cp, Cj, Cx, Dx,
    )
end

function set_barrier!(ncl::ImplicitNCLModel, μ)
    ncl.hashx[] = UInt64(0)
    ncl.μk[] = μ
    return
end

function _update!(ncl::ImplicitNCLModel, x)
    idx = hash(x)
    if idx == ncl.hashx[]
        return
    end
    ncl.hashx[] = idx

    lb = NLPModels.get_lcon(ncl.nlp)
    ub = NLPModels.get_ucon(ncl.nlp)
    ρ = ncl.ρk[]
    μ = ncl.μk[]
    y = ncl.yk
    rk = ncl.rk
    ck = ncl.ck
    NLPModels.cons!(ncl.nlp, x, ck)

    _reset_primal_slack_ncl!(rk, ck, y, ρ, μ, lb, ub)
    ncl.σk .= y .- ρ .* ncl.rk
    return
end

function NLPModels.obj(ncl::ImplicitNCLModel{T, VT}, x::VT) where {T, VT <: AbstractVector{T}}
    _update!(ncl, x)
    μ = ncl.μk[]
    ρ = ncl.ρk[]
    lb = NLPModels.get_lcon(ncl.nlp)
    ub = NLPModels.get_ucon(ncl.nlp)

    # Objective
    obj_val = NLPModels.obj(ncl.nlp, x)
    # Augmented Lagrangian term
    obj_aug = -dot(ncl.yk, ncl.rk) + ρ * dot(ncl.rk, ncl.rk) / T(2)
    # Barrier term
    obj_bar = zero(T)
    # lower-bound
    @inbounds for i in ncl.ind_lb
        obj_bar -= μ * log(ncl.rk[i] + ncl.ck[i] - lb[i])
    end
    # upper-bound
    @inbounds for i in ncl.ind_ub
        obj_bar -= μ * log(ub[i] - ncl.rk[i] - ncl.ck[i])
    end
    # range
    @inbounds for i in ncl.ind_rng
        obj_bar -= μ * (log(ub[i] - ncl.rk[i] - ncl.ck[i]) + log(ncl.rk[i] + ncl.ck[i] - lb[i]))
    end
    return obj_val + obj_aug + obj_bar
end

function NLPModels.cons!(ncl::ImplicitNCLModel{T, VT}, x::VT, c::VT) where {T, VT <: AbstractVector{T}}
    _update!(ncl, x)
    return c
end

function NLPModels.grad!(ncl::ImplicitNCLModel{T, VT}, x::VT, g::VT) where {T, VT <: AbstractVector{T}}
    _update!(ncl, x)
    jtv = similar(g)
    NLPModels.grad!(ncl.nlp, x, g)
    NLPModels.jtprod!(ncl.nlp, x, ncl.σk, jtv)
    g .+= jtv
    return g
end

function NLPModels.jprod!(ncl::ImplicitNCLModel{T, VT}, x::VT, v::VT, Jv::VT) where {T, VT <: AbstractVector{T}}
    _update!(ncl, x)
    return Jv
end

function NLPModels.jtprod!(ncl::ImplicitNCLModel{T, VT}, x::VT, v::VT, Jtv::VT) where {T, VT <: AbstractVector{T}}
    _update!(ncl, x)
    return Jtv
end

function NLPModels.jac_structure!(ncl::ImplicitNCLModel, jrows::VI, jcols::VI) where {VI <: AbstractVector{Int}}
    return Int[], Int[]
end

function NLPModels.jac_coord!(ncl::ImplicitNCLModel{T, VT}, x::VT, jac::VT) where {T, VT <: AbstractVector{T}}
    _update!(ncl, x)
    return jac
end

function NLPModels.hess_structure!(ncl::ImplicitNCLModel, hrows::VI, hcols::VI) where {VI <: AbstractVector{Int}}
    n = NLPModels.get_nvar(ncl)
    nnzhx = NLPModels.get_nnzh(ncl.nlp)
    hrowsx = view(hrows, 1:nnzhx)
    hcolsx = view(hcols, 1:nnzhx)
    NLPModels.hess_structure!(ncl.nlp, hrowsx, hcolsx)
    # Unpack condensed Jacobian Jᵀ J from CSR to COO
    cnt = nnzhx
    @inbounds for i in 1:n
        for c in ncl.Lp[i]:ncl.Lp[i+1]-1
            cnt += 1
            hrows[cnt] = i
            hcols[cnt] = ncl.Lj[c]
        end
    end
    return (hrows, hcols)
end

function NLPModels.hess_coord!(
    ncl::ImplicitNCLModel{T, VT},
    x::AbstractVector{T},
    y::AbstractVector{T},
    hess::AbstractVector{T};
    obj_weight::T=one(T),
) where {T, VT <: AbstractVector{T}}
    _update!(ncl, x)
    n, m = NLPModels.get_nvar(ncl.nlp), NLPModels.get_ncon(ncl.nlp)
    ind_eq = ncl.ind_eq
    ind_ineq = setdiff(1:m, ind_eq)
    # Evaluate Hessian of original model
    nnzhx = NLPModels.get_nnzh(ncl.nlp)
    hessx = view(hess, 1:nnzhx)
    NLPModels.hess_coord!(ncl.nlp, x, ncl.σk, hessx; obj_weight=obj_weight)

    # Evaluate Jacobian of original model
    nnzj = NLPModels.get_nnzj(ncl.nlp)
    jac = ncl.buffer_jac
    NLPModels.jac_coord!(ncl.nlp, x, jac)

    # Transfer data to (transposed) Jacobian
    ncl.JTx .= jac[ncl.map_jac]

    # Assemble diagonal
    c, r = ncl.ck, ncl.rk
    c♭, c♯ = NLPModels.get_lcon(ncl.nlp), NLPModels.get_ucon(ncl.nlp)
    ρ, μ = ncl.ρk[], ncl.μk[]
    # Assemble diagonal matrix for NCL
    @inbounds for i in ind_ineq
        Σ = μ / (c[i] + r[i] - c♭[i])^2 + μ / (c[i] + r[i] - c♯[i])^2
        ncl.Dx[i] = ρ * Σ / (ρ + Σ)
    end
    @inbounds for i in ind_eq
        ncl.Dx[i] = ρ
    end
    # Assemble condensed matrix Jᵀ D J
    assemble_condensed_jacobian!(
        n, m,
        ncl.JTp, ncl.JTj, ncl.JTx,
        ncl.Lp, ncl.Lj, ncl.Lx, ncl.Dx,
    )
    hess[nnzhx+1:end] .= ncl.Lx

    return hess
end

