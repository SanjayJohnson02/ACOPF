
function symul!(y::AbstractVector, A::AbstractMatrix, x::AbstractVector, alpha::Number, beta::Number)
    return mul!(y, Symmetric(A, :L), x, alpha, beta)
end

function coo_to_csr(
    n_rows,
    n_cols,
    Ai::AbstractVector{Ti},
    Aj::AbstractVector{Ti},
    Ax::AbstractVector{Tv},
) where {Tv, Ti}
    @assert length(Ai) == length(Aj) == length(Ax)
    nnz = length(Ai)
    Bp = zeros(Ti, n_rows+1)
    Bj = zeros(Ti, nnz)
    Bx = zeros(Tv, nnz)

    nnz = length(Ai)
    @inbounds for n in 1:nnz
        Bp[Ai[n]] += 1
    end

    # cumsum the nnz per row to get Bp
    cumsum = 1
    @inbounds for i in 1:n_rows
        tmp = Bp[i]
        Bp[i] = cumsum
        cumsum += tmp
    end
    Bp[n_rows+1] = nnz + 1

    @inbounds for n in 1:nnz
        i = Ai[n]
        dest = Bp[i]
        Bj[dest] = Aj[n]
        Bx[dest] = Ax[n]
        Bp[i] += 1
    end

    last = 1
    @inbounds for i in 1:n_rows+1
        tmp = Bp[i]
        Bp[i] = last
        last = tmp
    end

    return (Bp, Bj, Bx)
end

function build_symbolic_analysis_jtj(
    n_rows,
    n_cols,
    Jtp::AbstractVector{Ti},
    Jtj::AbstractVector{Ti},
) where {Ti}
    Cp = zeros(Ti, n_rows + 1)
    xb = zeros(UInt8, n_cols)

    # Count nonzeros per rows
    nnz = 0
    @inbounds for i in 1:n_rows
        for c in Jtp[i]:Jtp[i+1]-1
            j = Jtj[c]
            xb[j] = UInt8(1)
        end
        # JᵀJ is symmetric, store only lower triangular part
        for j in 1:i
            for c in Jtp[j]:Jtp[j+1]-1
                k = Jtj[c]
                if xb[k] == 1
                    nnz += 1
                    Cp[i] += 1
                    break
                end
            end
        end
        # Reset to 0
        for c in Jtp[i]:Jtp[i+1]-1
            xb[Jtj[c]] = UInt8(0)
        end
    end
    # cumsum the nnz per row to get Bp
    cumsum = 1
    @inbounds for i in 1:n_rows
        tmp = Cp[i]
        Cp[i] = cumsum
        cumsum += tmp
    end
    Cp[n_rows+1] = nnz + 1

    Cj = zeros(Ti, nnz)
    cnt = 0
    @inbounds for i in 1:n_rows
        for c in Jtp[i]:Jtp[i+1]-1
            j = Jtj[c]
            xb[j] = UInt8(1)
        end
        # JᵀJ is symmetric, store only lower triangular part
        for j in 1:i
            for c in Jtp[j]:Jtp[j+1]-1
                k = Jtj[c]
                if xb[k] == 1
                    cnt += 1
                    Cj[cnt] = j
                    break
                end
            end
        end
        for c in Jtp[i]:Jtp[i+1]-1
            xb[Jtj[c]] = UInt8(0)
        end
    end

    return (Cp, Cj)
end

function assemble_condensed_jacobian!(
    n_rows,
    n_cols,
    Jtp::AbstractVector{Ti},
    Jtj::AbstractVector{Ti},
    Jtx::AbstractVector{Tv},
    Cp::AbstractVector{Ti},
    Cj::AbstractVector{Ti},
    Cx::AbstractVector{Tv},
    Dx::AbstractVector{Tv},
) where {Ti, Tv}

    buffer = zeros(Tv, n_cols)

    @inbounds for i in 1:n_rows
        # Read row i
        for c in Jtp[i]:Jtp[i+1]-1
            j = Jtj[c]
            buffer[j] = Jtx[c] * Dx[j]
        end
        for c in Cp[i]:Cp[i+1]-1
            j = Cj[c]
            Cx[c] = Tv(0)
            for d in Jtp[j]:Jtp[j+1]-1
                k = Jtj[d]
                Cx[c] += buffer[k] * Jtx[d]
            end
        end
        # Reset buffer
        for c in Jtp[i]:Jtp[i+1]-1
            j = Jtj[c]
            buffer[j] = Tv(0)
        end
    end
end

# Magic steps for NCL
function _reset_primal_slack_ncl!(
    rk, ck,
    y, ρ, μ,
    lb, ub,
)
    @assert length(lb) == length(ub)
    m = length(lb)

    has_succeeded = true
    ki = 0
    ## equality
    @inbounds for i in 1:m
        if lb[i] == ub[i]
            # EQ
            rk[i] = -ck[i] + lb[i]
        elseif isfinite(lb[i]) && isinf(ub[i])
            # LB
            ki +=1
            a0 = -μ - y[i] * (ck[i] - lb[i])
            a1 = ρ * (ck[i] - lb[i]) - y[i]
            a2 = ρ
            Δ = a1^2 - 4*a0*a2
            rk[i] = (-a1 + sqrt(Δ)) / (2*a2)
            has_succeeded &= (ck[i] + rk[i] >= lb[i])
        elseif isinf(lb[i]) && isfinite(ub[i])
            # UB
            ki +=1
            a0 = -μ - y[i] * (ck[i] - ub[i])
            a1 = ρ * (ck[i] - ub[i]) - y[i]
            a2 = ρ
            Δ = a1^2 - 4*a0*a2
            rk[i] = (-a1 - sqrt(Δ)) / (2*a2)
            has_succeeded &= (ck[i] + rk[i] <= ub[i])
        elseif isfinite(lb[i]) && isfinite(ub[i])
            # RNG
            ki +=1
            a0 = -μ * (2*ck[i] - lb[i] - ub[i]) - y[i] * (ck[i] - lb[i]) * (ck[i] - ub[i])
            a1 = ρ * (ck[i] - lb[i]) * (ck[i] - ub[i]) - y[i] * (2*ck[i] - lb[i] - ub[i]) - 2*μ
            a2 = ρ * (2*ck[i] - lb[i] - ub[i]) - y[i]
            a3 = ρ
            # N.B.: solve cubic polynomial to find the solution of the magic step
            #       for range constraints. The second root is the only one returning
            #       a feasible slack variable.
            rk[i] = roots(Polynomial([a0, a1, a2, a3]))[2]
            has_succeeded &= (lb[i] <= ck[i] + rk[i] <= ub[i])
        end
    end
    @assert has_succeeded
    return has_succeeded
end

function _reset_dual_slack_ncl!(
    λ, zl, zu,
    rk, sk, y, ρ, μ,
    lb, ub,
)
    @assert length(λ) == length(rk)
    @assert length(zl) == length(zu) == length(sk)
    ind_ineq = findall(lb .< ub)
    ki = 0
    m = length(lb)
    # λ .= y .- ρ .* rk
    zl .= μ ./ (sk .- lb[ind_ineq])
    zu .= μ ./ (ub[ind_ineq] .- sk)
    return
end

