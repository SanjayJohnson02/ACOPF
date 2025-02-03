
using Pkg.Artifacts
using PowerModels
using ExaModels

PowerModels.silence()
pglib_path = joinpath(artifact"PGLib_opf", "pglib-opf-23.07")

include(joinpath(@__DIR__, "..", "common.jl"))

OPF_INSTANCES = [
    "pglib_opf_case1354_pegase",
    "pglib_opf_case1888_rte",
    "pglib_opf_case1951_rte",
    "pglib_opf_case2000_goc",
    "pglib_opf_case2312_goc",
    "pglib_opf_case2742_goc",
    "pglib_opf_case2848_rte",
    "pglib_opf_case2868_rte",
    "pglib_opf_case2869_pegase",
    "pglib_opf_case3022_goc",
    "pglib_opf_case3970_goc",
    "pglib_opf_case4020_goc",
    "pglib_opf_case4601_goc",
    "pglib_opf_case4619_goc",
    "pglib_opf_case4661_sdet",
    "pglib_opf_case4837_goc",
    "pglib_opf_case4917_goc",
    "pglib_opf_case5658_epigrids",
    "pglib_opf_case6468_rte",
    "pglib_opf_case9241_pegase",
]

function maxload_model(file_name)
    data = PowerModels.parse_file(file_name)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    ngen = length(ref[:gen])
    alpha = ones(ngen)

    model = JuMP.Model()

    JuMP.@variable(model, va[i in keys(ref[:bus])], start=0.0)
    JuMP.@variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)
    JuMP.@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    JuMP.@variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    JuMP.@variable(model, 0 <= λ, start=1.0)

    # Feasibility problem.
    JuMP.@objective(model, Min, -λ)

    # Voltage angle at ref node
    for (i,bus) in ref[:ref_buses]
        JuMP.@constraint(model, va[i] == 0)
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        JuMP.@constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            λ * sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2
        )

        JuMP.@constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) -
            λ * sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2
        )
    end

    # Branch power flow physics and limit constraints
    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        vm_fr = vm[branch["f_bus"]]
        vm_to = vm[branch["t_bus"]]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        ttm = tr^2 + ti^2
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]

        # From side of the branch flow
        JuMP.@constraint(model, p_fr ==  (g+g_fr)/ttm*vm_fr^2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        JuMP.@constraint(model, q_fr == -(b+b_fr)/ttm*vm_fr^2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        JuMP.@constraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        JuMP.@constraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Voltage angle difference limit
        JuMP.@constraint(model, branch["angmin"] <= va_fr - va_to <= branch["angmax"])

        # Apparent power limit, from side and to side
        JuMP.@constraint(model, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
        JuMP.@constraint(model, p_to^2 + q_to^2 <= branch["rate_a"]^2)
    end

    return model
end

MAXOPF_INSTANCES = [
    (() -> ExaModels.ExaModel(maxload_model(joinpath(pglib_path, "$(instance).m"))), (), instance) for instance in OPF_INSTANCES
]

function main(
    name::String="maxopf";
    tol::Float64=1e-8,
    max_iter::Int=1000,
    kkt_benchmark::Bool=false,
    madnlp::Bool=true,
    madncl::Bool=true,
    verbose::Bool=false,
)

    results_dir = joinpath(@__DIR__, "..",  "results") |> normpath
    !isdir(results_dir) && mkpath(results_dir)

    prefix = kkt_benchmark ? "kkt_" : ""
    _madnlp = kkt_benchmark ? build_madnlp : solve_madnlp
    _madncl = kkt_benchmark ? build_madncl : solve_madncl
    _benchmark = kkt_benchmark ? run_kkt : run_benchmark

    all_solvers = []

    if madnlp
        append!(all_solvers, [
            (_madnlp, "madnlp", "cpu", Ma27Solver, "ma27", MadNLP.SparseKKTSystem, "K2"),
            (_madnlp, "madnlp", "cpu", Ma57Solver, "ma57", MadNLP.SparseKKTSystem, "K2"),
        ])
    end

    if madncl
        append!(all_solvers, [
            (_madncl, "madncl", "cpu", Ma27Solver, "ma27", MadNLP.SparseKKTSystem, "K2"),
            (_madncl, "madncl", "cpu", Ma27Solver, "ma27", MadNCL.K2rAuglagKKTSystem, "K2r"),
            (_madncl, "madncl", "cpu", Ma27Solver, "ma27", MadNCL.K1sAuglagKKTSystem, "K1s"),
        ])
    end

    @info "Warm-start -- MAXOPF"
    for (builder, algo, solver, solver_name, kkt, kkt_name) in all_solvers
        results = _benchmark("pglib_opf_case1354_pegase", builder, MAXOPF_INSTANCES; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=false)
    end

    for (builder, algo, solver, solver_name, kkt, kkt_name)  in all_solvers
        @info "Benchmark $(algo)-$(solver_name)-$(kkt_name)"
        if name == "maxopf"
            results = _benchmark(builder, MAXOPF_INSTANCES; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=verbose)
        else
            results = _benchmark(name, builder, MAXOPF_INSTANCES; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=verbose)
        end
        path_results = joinpath(results_dir, "$(prefix)$(name)_$(algo)_$(kkt_name)_$(solver_name).txt")
        open(path_results, "w") do io
            writedlm(io, results)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

