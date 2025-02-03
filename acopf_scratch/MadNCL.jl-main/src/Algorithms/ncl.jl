
#=
    NCLOptions
=#

@kwdef struct NCLOptions{T}
    verbose::Bool = true
    slack_reset::Bool = false
    opt_tol::T = 1e-6
    feas_tol::T = 1e-6
    rho0::T = 1e2
    rhomax::T = 1e12
    max_auglag_iter::Int = 20
end

#=
    NCLSolver
=#

struct NCLSolver{T, VT, M}
    ncl::NCLModel{T, VT, M}
    ipm::MadNLP.MadNLPSolver{T, VT}
    options::NCLOptions{T}
    n::Int
    m::Int
end

function NCLSolver(nlp::NLPModels.AbstractNLPModel{T, VT}; scaling=true, ncl_options=NCLOptions(), ipm_options...) where {T, VT}
    n, m = NLPModels.get_nvar(nlp), NLPModels.get_ncon(nlp)
    ncl = if scaling
        NCLModel(ScaledModel(nlp))
    else
        NCLModel(nlp)
    end
    solver = MadNLP.MadNLPSolver(
        ncl;
        nlp_scaling=false,
        ipm_options...,
    )
    return NCLSolver{T, VT, typeof(ncl.nlp)}(ncl, solver, ncl_options, n, m)
end

#=
    NCLStats
=#

mutable struct NCLStats{T, VT} <: AbstractExecutionStats
    status::MadNLP.Status
    solution::VT
    regularization::VT
    objective::T
    dual_feas::T
    primal_feas::T
    multipliers::VT
    multipliers_L::VT
    multipliers_U::VT
    iter::Int
    counters::MadNLP.MadNLPCounters
end

function NCLStats(solver::NCLSolver{T, VT, M}, status) where {T, VT, M<:NLPModels.AbstractNLPModel}
    n, m = solver.n, solver.m
    ncl = solver.ncl
    x_final = MadNLP.primal(solver.ipm.x)[1:n]
    r_final = MadNLP.primal(solver.ipm.x)[1+n:m+n]
    zl = MadNLP.primal(solver.ipm.zl)[1:n]
    zu = MadNLP.primal(solver.ipm.zu)[1:n]
    return NCLStats(
        status,
        x_final,
        r_final,
        NLPModels.obj(ncl.nlp, x_final),
        solver.ipm.inf_du,
        norm(r_final, Inf),
        copy(solver.ncl.yk),
        zl,
        zu,
        solver.ipm.cnt.k,
        solver.ipm.cnt,
    )
end

function NCLStats(solver::NCLSolver{T, VT, M}, status) where {T, VT, M<:ScaledModel}
    n, m = solver.n, solver.m
    ncl = solver.ncl
    obj_scale, con_scale = ncl.nlp.scaling_obj, ncl.nlp.scaling_cons
    # Unscale solution
    x_final = MadNLP.primal(solver.ipm.x)[1:n]
    r_final = MadNLP.primal(solver.ipm.x)[1+n:m+n] ./ con_scale
    zl = MadNLP.primal(solver.ipm.zl)[1:n] ./ obj_scale
    zu = MadNLP.primal(solver.ipm.zu)[1:n] ./ obj_scale
    y = copy(solver.ipm.y) .* con_scale ./ obj_scale
    return NCLStats(
        status,
        x_final,
        r_final,
        NLPModels.obj(ncl.nlp, x_final) ./ obj_scale,
        solver.ipm.inf_du,
        norm(r_final, Inf),
        y,
        zl,
        zu,
        solver.ipm.cnt.k,
        solver.ipm.cnt,
    )
end

function getStatus(result::NCLStats)
    if result.status == MadNLP.SOLVE_SUCCEEDED
        println("Optimal solution found.")
    elseif result.status == MadNLP.INFEASIBLE_PROBLEM_DETECTED
        println("Convergence to an infeasible point.")
    elseif result.status == MadNLP.MAXIMUM_ITERATIONS_EXCEEDED
        println("Maximum number of iterations reached.")
    else
        println("Unknown return status.")
    end
end



#=
    NCL Algorithm
=#

function _introduce(nx, nr)
    println("MadNCL algorithm\n")

    println("Total number of variables............................:      ", nx)
    println("Total number of constraints..........................:      ", nr)
    println()
end

function _log_header()
    @printf(
        "outer  inner     objective    inf_pr   inf_du    η        μ       ρ \n"
    )
end

function _log_iter(nit, flag, n_inner, obj, inf_pr, inf_du, alpha, mu, rho)
    @printf(
        "%5s%1s %5i %+13.7e %6.2e %6.2e %6.2e %6.1e %6.2e\n",
        nit, flag, n_inner, obj, inf_pr, inf_du, alpha, mu, rho,
    )
end

function get_inf_du(solver::NCLSolver)
    return solver.ipm.inf_du
end

# Slack-reset, aka magic step
function project!(solver::MadNLP.MadNLPSolver; dual=false)
    ncl = solver.nlp
    n = solver.n
    nx, m = ncl.nx, ncl.nr
    yk = ncl.yk
    ρk = ncl.ρk[]
    μk = solver.mu

    lb, ub = NLPModels.get_lcon(ncl.nlp), NLPModels.get_ucon(ncl.nlp)

    ck = solver.c
    xk = view(MadNLP.primal(solver.x), 1:nx)
    rk = view(MadNLP.primal(solver.x), nx+1:nx+m)
    sk = MadNLP.slack(solver.x)

    # update constraint
    NLPModels.cons!(ncl.nlp, xk, ck)

    _reset_primal_slack_ncl!(rk, ck, yk, ρk, μk, lb, ub)
    for (k, i) in enumerate(solver.ind_ineq)
        sk[k] = ck[i] + rk[i]
    end

    if dual
        n_ineq = length(solver.ind_ineq)
        λ = solver.y
        zl = MadNLP.slack(solver.zl)
        zu = MadNLP.slack(solver.zu)
        _reset_dual_slack_ncl!(λ, zl, zu, rk, sk, yk, ρk, μk, lb, ub)
    end
    return
end

function setup!(solver; μ=1e-1, tol=1e-8, slack_reset=true)
    # Update options
    solver.opt.mu_init = μ
    solver.opt.tol = tol
    solver.mu = solver.opt.mu_init
    # Ensure the barrier parameter is fixed
    solver.opt.mu_min = solver.opt.mu_init

    # Slack reset
    if slack_reset
        project!(solver)
    end

    # Refresh values
    solver.obj_val = MadNLP.eval_f_wrapper(solver, solver.x)
    MadNLP.eval_grad_f_wrapper!(solver, solver.f, solver.x)
    MadNLP.eval_cons_wrapper!(solver, solver.c, solver.x)

    # Update filter
    theta = MadNLP.get_theta(solver.c)
    solver.theta_max = 1e4*max(1,theta)
    solver.theta_min = 1e-4*max(1,theta)
    solver.tau = max(solver.opt.tau_min,1-solver.opt.mu_init)
    empty!(solver.filter)
    push!(solver.filter, (solver.theta_max,-Inf))

    return MadNLP.REGULAR
end

function solve!(solver::NCLSolver)
    n, m = solver.n, solver.m
    ncl = solver.ncl
    options = solver.options

    options.verbose && _introduce(n, m)

    # Parameters
    ### Penalty ρ
    ncl.ρk[] = options.rho0
    ρ_max = options.rhomax
    tau_ρ = 10.0
    ### Barrier parameters
    μ = 1e-1
    μ_min = 1e-8
    γ = 0.05
    τ = 0.6 * 2.0 / (1 + γ)
    μ_fac = 0.2
    ### Forcing parameters
    eps_c = 10*μ
    eps_d = 100*μ^(1.0+γ)
    η = 1.0e-1  # initial primal feasibility tolerance

    start_time = time()

    ncl_status = MadNLP.INITIAL

    # Unless the parameter `dual_initialized` is set to `true`,
    # MadNLP computes the initial multipliers y0 in `initialize!`.
    MadNLP.initialize!(solver.ipm)
    # Update initial multiplier using multiplier computed by MadNLP in `initialize!`.
    ncl.yk .= solver.ipm.y

    solver.ipm.status = MadNLP.REGULAR
    cnt_it_inner = 0

    x = view(MadNLP.primal(solver.ipm.x), 1:n)
    r = view(MadNLP.primal(solver.ipm.x), 1+n:n+m)

    primal_feas = solver.ipm.inf_pr
    dual_feas = get_inf_du(solver)

    options.verbose && _log_header()
    flag = " "
    options.verbose && _log_iter(0, flag, 0, solver.ipm.obj_val, primal_feas, dual_feas, η, μ, ncl.ρk[])

    iter = 1
    while iter <= options.max_auglag_iter
        setup!(
            solver.ipm;
            tol=eps_d,
            μ=μ,
            slack_reset=options.slack_reset,
        )

        status = MadNLP.regular!(solver.ipm)
        # Check return status
        has_converged = status == MadNLP.SOLVE_SUCCEEDED
        flag = has_converged ? " " : "r"

        # If solve has not succeeded, we call a feasibility restoration phase.
        if !has_converged
            MadNLP.robust!(solver.ipm)
        end

        # Unpack solution
        primal_feas = norm(r, Inf)
        dual_feas = get_inf_du(solver)

        if primal_feas <= max(η, options.feas_tol)
            # Update multiplier.
            axpy!(-ncl.ρk[], r, ncl.yk)
            # Decrease barrier parameter
            μ = min(μ^τ, μ_fac * μ)
            # Update forcing parameter
            eps_d = 100.0 * μ^(1+γ)
            eps_c = 10.0 * μ
            η = min(μ^1.1, 0.1 * μ)
        else
            # Increase penalty
            ncl.ρk[] = min(ncl.ρk[] * tau_ρ, ρ_max)
        end

        # Log evolution
        ipm_iter = solver.ipm.cnt.k
        obj_val = NLPModels.obj(ncl.nlp, x)
        options.verbose && _log_iter(iter, flag, ipm_iter, obj_val, primal_feas, dual_feas, η, μ, ncl.ρk[])

        if primal_feas <= options.feas_tol && dual_feas <= options.opt_tol
            # Check convergence
            ncl_status = MadNLP.SOLVE_SUCCEEDED
            break
        elseif (ncl.ρk[] >= options.rhomax) && (primal_feas > options.feas_tol)
            # Check infeasibility
            ncl_status = MadNLP.INFEASIBLE_PROBLEM_DETECTED
            break
        end
        iter += 1
    end

    if (iter >= options.max_auglag_iter) || (solver.ipm.status == MadNLP.MAXIMUM_ITERATIONS_EXCEEDED)
        ncl_status = MadNLP.MAXIMUM_ITERATIONS_EXCEEDED
    end

    solver.ipm.cnt.total_time = time() - start_time

    return NCLStats(solver, ncl_status)
end

function madncl(
    nlp::NLPModels.AbstractNLPModel;
    options...
)
    solver = NCLSolver(nlp; options...)
    return solve!(solver)
end
