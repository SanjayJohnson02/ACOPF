#=

    Script to benchmark nonlinear solvers.

    Four classes of problems are used:
    - Tax instances (from Ampl files)
    - ACOPF problems (from PGLib_opf)
    - COPS instances (from COPSBenchmark.jl)
    - Degenerate COPS instances (from COPSBenchmark.jl)

=#

using DelimitedFiles
using LinearAlgebra
using NLPModels
using JuMP
using MadNLPHSL
using MadNCL
using HybridKKT

mkl = false
if mkl
    using MKL
end
lbt = mkl ? "mkl" : "openblas"

function time_linear_solver(solver::MadNLP.MadNLPSolver)
    t_build     = @elapsed MadNLP.build_kkt!(solver.kkt)
    t_factorize = @elapsed MadNLP.factorize!(solver.kkt.linear_solver)
    t_backsolve = @elapsed MadNLP.solve!(solver.kkt, solver.d)
    return (
        time_build_kkt = t_build,
        time_factorization = t_factorize,
        time_backsolve = t_backsolve,
    )
end

time_linear_solver(solver::MadNCL.NCLSolver) = time_linear_solver(solver.ipm)

function get_primal_residual(model::AbstractNLPModel, stats)
    x = stats.solution
    c = NLPModels.cons(model, x)
    cl = NLPModels.get_lcon(model)
    cu = NLPModels.get_ucon(model)
    err = 0.0
    for i in eachindex(c)
        if cl[i] == cu[i]
            err = max(err, abs(c[i] - cl[i]))
        else
            if isfinite(cl[i])
                inf_pr = max(0.0, cl[i] - c[i])
                err = max(err, abs(inf_pr))
            end
            if isfinite(cu[i])
                inf_pr = max(0.0, c[i] - cu[i])
                err = max(err, abs(inf_pr))
            end
        end
    end
    return err
end

# TODO: does not work for max problem or with fixed variables
function get_dual_residual(model::AbstractNLPModel, stats)
    x = stats.solution
    y = stats.multipliers
    g = NLPModels.grad(model, x)
    jtv = NLPModels.jtprod(model, x, y)
    zl = stats.multipliers_L
    zu = stats.multipliers_U
    err = norm(g .+ jtv .- zl .+ zu, Inf)
    return err
end

function get_primal_complementarity_residual(model::AbstractNLPModel, stats)
    x = stats.solution
    xl = NLPModels.get_lvar(model)
    xu = NLPModels.get_uvar(model)
    zl = stats.multipliers_L
    zu = stats.multipliers_U
    err = 0.0
    for i in eachindex(x)
        if isfinite(xl[i])
            inf_compl = (x[i] - xl[i]) * zl[i]
            err = max(err, abs(inf_compl))
        end
        if isfinite(xu[i])
            inf_compl = (xu[i] - x[i]) * zu[i]
            err = max(err, abs(inf_compl))
        end
    end
    return err
end

function get_dual_complementarity_residual(model::AbstractNLPModel, stats)
    x = stats.solution
    y = stats.multipliers
    c = NLPModels.cons(model, x)
    cl = NLPModels.get_lcon(model)
    cu = NLPModels.get_ucon(model)
    err = 0.0
    for i in eachindex(c)
        if cl[i] < cu[i]
            if isfinite(cl[i])
                err = max(err, abs((c[i] - cl[i]) * y[i]))
            end
            if isfinite(cu[i])
                err = max(err, abs((cu[i] - c[i]) * y[i]))
            end
        end
    end
    return err
end

function print_accuracy(nlp::AbstractNLPModel, stats)
    println("* Primal feasibility: ", get_primal_residual(nlp, stats))
    println("* Dual feasibility:   ", get_dual_residual(nlp, stats))
    println("* Primal comp:        ", get_primal_complementarity_residual(nlp, stats))
    println("* Dual comp:          ", get_dual_complementarity_residual(nlp, stats))
end

get_optimizer(model::JuMP.Model) = model.moi_backend.optimizer.model
get_backend(model::JuMP.Model) = model.moi_backend.optimizer.model.solver

function build_madnlp(model::AbstractNLPModel; kkt_system=MadNLP.SparseKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    print_level = verbose ? MadNLP.INFO : MadNLP.ERROR
    solver = MadNLP.MadNLPSolver(
        model;
        linear_solver=linear_solver,
        kkt_system=kkt_system,
        max_iter=max_iter,
        print_level=print_level,
        tol=tol,
    )
    return solver
end

function solve_madnlp(model::AbstractNLPModel; kkt_system=MadNLP.SparseKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    solver = build_madnlp(model; kkt_system, linear_solver, tol, max_iter, verbose)
    stats = MadNLP.solve!(solver)
    # Issue with MadNLP: the solver does not rescale the multipliers
    stats.multipliers .*= solver.cb.con_scale
    return stats
end

function build_madnlp(model::JuMP.Model; kkt_system=MadNLP.SparseKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    nlp = MathOptNLPModel(model)
    build_madnlp(nlp; kkt_system, linear_solver, tol, max_iter, verbose)
end

function solve_madnlp(model::JuMP.Model; kkt_system=MadNLP.SparseKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    nlp = MathOptNLPModel(model)
    solve_madnlp(nlp; kkt_system, linear_solver, tol, max_iter, verbose)
end

function build_madncl(model::AbstractNLPModel; kkt_system=MadNLP.SparseKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    print_level = verbose ? MadNLP.INFO : MadNLP.ERROR
    ncl_options = MadNCL.NCLOptions(
        max_auglag_iter=40,
        feas_tol=tol,
        opt_tol=tol,
        slack_reset=false,
        verbose=verbose,
    )
    solver = MadNCL.NCLSolver(
        model;
        ncl_options=ncl_options,
        linear_solver=linear_solver,
        kkt_system=kkt_system,
        max_iter=max_iter,
        print_level=print_level,
        nlp_scaling=false,
    )
    return solver
end

function solve_madncl(model::AbstractNLPModel; kkt_system=MadNLP.SparseKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    solver = build_madncl(model; kkt_system, linear_solver, tol, max_iter, verbose)
    MadNCL.solve!(solver)
end

function build_madncl(model::JuMP.Model; kkt_system=MadNLP.SparseKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    nlp = MathOptNLPModel(model)
    build_madncl(nlp; kkt_system, linear_solver, tol, max_iter, verbose)
end

function solve_madncl(model::JuMP.Model; kkt_system=MadNLP.SparseKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    nlp = MathOptNLPModel(model)
    solve_madncl(nlp; kkt_system, linear_solver, tol, max_iter, verbose)
end

function build_likkt(model::AbstractNLPModel; kkt_system=MadNLP.SparseCondensedKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    print_level = verbose ? MadNLP.INFO : MadNLP.ERROR
    solver = MadNLP.MadNLPSolver(
        model;
        linear_solver=linear_solver,
        kkt_system=kkt_system,
        equality_treatment=MadNLP.RelaxEquality,
        fixed_variable_treatment=MadNLP.RelaxBound,
        dual_initialized=true,
        max_iter=max_iter,
        print_level=print_level,
        tol=tol,
    )
    return solver
end

function build_likkt(model::JuMP.Model; kkt_system=MadNLP.SparseCondensedKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    nlp = MathOptNLPModel(model)
    build_likkt(nlp; kkt_system, linear_solver, tol, max_iter, verbose)
end

function solve_likkt(model::AbstractNLPModel; kkt_system=MadNLP.SparseCondensedKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    solver = build_likkt(model; kkt_system, linear_solver, tol, max_iter, verbose)
    stats = MadNLP.solve!(solver)
    # Issue with MadNLP: the solver does not rescale the multipliers
    stats.multipliers .*= solver.cb.con_scale
    return stats
end

function solve_likkt(model::JuMP.Model; kkt_system=MadNLP.SparseCondensedKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    nlp = MathOptNLPModel(model)
    solve_likkt(nlp; kkt_system, linear_solver, tol, max_iter, verbose)
end

function build_hykkt(model::AbstractNLPModel; gamma::Float64=1e7, kkt_system=HybridKKT.HybridCondensedKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    print_level = verbose ? MadNLP.INFO : MadNLP.ERROR
    solver = MadNLP.MadNLPSolver(
        model;
        linear_solver=linear_solver,
        kkt_system=kkt_system,
        equality_treatment=MadNLP.EnforceEquality,
        fixed_variable_treatment=MadNLP.MakeParameter,
        max_iter=max_iter,
        print_level=print_level,
        tol=tol,
    )
    solver.kkt.gamma[] = gamma
    return solver
end

function build_hykkt(model::JuMP.Model; gamma::Float64=1e7, kkt_system=MadNLP.SparseCondensedKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    nlp = MathOptNLPModel(model)
    build_hykkt(nlp; gamma, kkt_system, linear_solver, tol, max_iter, verbose)
end

function solve_hykkt(model::AbstractNLPModel; gamma::Float64=1e7, kkt_system=HybridKKT.HybridCondensedKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    solver = build_hykkt(model; gamma, kkt_system, linear_solver, tol, max_iter, verbose)
    stats = MadNLP.solve!(solver)
    # Issue with MadNLP: the solver does not rescale the multipliers
    stats.multipliers .*= solver.cb.con_scale
    return stats
end

function solve_hykkt(model::JuMP.Model; gamma::Float64=1e7, kkt_system=MadNLP.SparseCondensedKKTSystem, linear_solver=Ma27Solver, tol::Float64=1e-8, max_iter::Int=1000, verbose::Bool=false)
    nlp = MathOptNLPModel(model)
    solve_hykkt(nlp; gamma, kkt_system, linear_solver, tol, max_iter, verbose)
end

function run_benchmark(solver, instances; options...)
    n_instances = length(instances)
    results = zeros(n_instances, 12)
    k = 1
    for (instance, args, instance_name) in instances
        @info instance_name
        model = instance(args...)
        tol = haskey(options, :tol) ? options[:tol] : 1e-8
        stats = solver(model; options...)
        results[k, 1]  = Int(stats.status)
        results[k, 2]  = stats.objective
        results[k, 3]  = stats.iter
        results[k, 4]  = stats.counters.total_time
        results[k, 5]  = stats.counters.linear_solver_time
        results[k, 6]  = stats.counters.eval_function_time
        results[k, 7]  = tol
        results[k, 8]  = get_primal_residual(model, stats)
        results[k, 9]  = get_dual_residual(model, stats)
        results[k, 10] = get_primal_complementarity_residual(model, stats)
        results[k, 11] = get_dual_complementarity_residual(model, stats)
        results[k, 12] = norm(results[k, 8:11], Inf)
        k += 1
    end
    names = [n for (f, a, n) in instances]
    cols = [
        "instance",
        "status",
        "obj",
        "it",
        "total",
        "linsolve",
        "AD",
        "tolerance",
        "primal residual",
        "dual residual",
        "primal complementarity residual",
        "dual complementarity residual",
        "maximum residual",
    ]
    return [reshape(cols, 1, 13) ; names results]
end

function run_benchmark(name::String, solver, instances; options...)
    found = false
    pos = 0
    for (k, instance) in enumerate(instances)
        if instance[end] == name
            found = true
            pos = k
        end
    end
    !found && error("The problem $name is not in the collection.")
    results = zeros(1, 12)
    (instance, args, instance_name) = instances[pos]
    @info instance_name
    model = instance(args...)
    tol = haskey(options, :tol) ? options[:tol] : 1e-8
    stats = solver(model; options...)
    results[1, 1]  = Int(stats.status)
    results[1, 2]  = stats.objective
    results[1, 3]  = stats.iter
    results[1, 4]  = stats.counters.total_time
    results[1, 5]  = stats.counters.linear_solver_time
    results[1, 6]  = stats.counters.eval_function_time
    results[1, 7]  = tol
    results[1, 8]  = get_primal_residual(model, stats)
    results[1, 9]  = get_dual_residual(model, stats)
    results[1, 10] = get_primal_complementarity_residual(model, stats)
    results[1, 11] = get_dual_complementarity_residual(model, stats)
    results[1, 12] = norm(results[1, 8:11], Inf)
    names = [name]
    cols = [
        "instance",
        "status",
        "obj",
        "it",
        "total",
        "linsolve",
        "AD",
        "tolerance",
        "primal residual",
        "dual residual",
        "primal complementarity residual",
        "dual complementarity residual",
        "maximum residual",
    ]
    return [reshape(cols, 1, 13) ; names results]
end


function run_kkt(build_solver, instances; options...)
    n_instances = length(instances)
    results = zeros(n_instances, 3)
    k = 1
    for (instance, args, instance_name) in instances
        model = instance(args...)
        solver = build_solver(model; options...)
        solver isa MadNCL.NCLSolver ? MadNCL.solve!(solver) : MadNLP.solve!(solver)

        timings = time_linear_solver(solver)
        results[k, 1] = timings[1]
        results[k, 2] = timings[2]
        results[k, 3] = timings[3]
        k += 1
    end
    names = [n for (f, a, n) in instances]
    cols = [
        "instance",
        "timer_build",
        "timer_factorization",
        "timer_backsolve",
    ]

    return [reshape(cols, 1, 4) ; names results]
end

function run_kkt(name::String, build_solver, instances; options...)
    found = false
    pos = 0
    for (k, instance) in enumerate(instances)
        if instance[end] == name
            found = true
            pos = k
        end
    end
    !found && error("The problem $name is not in the collection.")

    results = zeros(1, 3)
    (instance, args, instance_name) = instances[pos]
    @info instance_name
    model = instance(args...)
    solver = build_solver(model; options...)
    solver isa MadNCL.NCLSolver ? MadNCL.solve!(solver) : MadNLP.solve!(solver)

    timings = time_linear_solver(solver)
    results[1, 1] = timings[1]
    results[1, 2] = timings[2]
    results[1, 3] = timings[3]
    names = [name]
    cols = [
        "instance",
        "timer_build",
        "timer_factorization",
        "timer_backsolve",
    ]

    return [reshape(cols, 1, 4) ; names results]
end
