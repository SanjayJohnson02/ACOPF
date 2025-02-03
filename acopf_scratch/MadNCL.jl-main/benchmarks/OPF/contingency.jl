
#=
    Script for contingency screening.

    Check if a given contingency is feasible w.r.t. the base case solution.
    The constraints model the PV/PQ switches and the droop control with
    complementarity constraints.

=#

using Printf
using LinearAlgebra
using DelimitedFiles
using JuMP
using PowerModels
using Ipopt
using NLPModelsJuMP
using MadNLP, MadNLPHSL
using MadNCL
using HSL_jll
using NLPModels
using ExaModels

PowerModels.silence()

include("opf.jl")

function extract_solution(model)
    return (
        va=JuMP.value.(model[:va]),
        vm=JuMP.value.(model[:vm]),
        pg=JuMP.value.(model[:pg]),
        qg=JuMP.value.(model[:qg]),
        p=JuMP.value.(model[:p]),
        q=JuMP.value.(model[:q]),
    )
end

function contingency_problem(file_name, base_case, contingency; load_factor=1.0)
    data = PowerModels.parse_file(file_name)
    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    ngen = length(ref[:gen])
    alpha = ones(ngen)

    model = JuMP.Model()

    JuMP.@variable(model, va[i in keys(ref[:bus])], start=base_case.va[i])
    JuMP.@variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=base_case.vm[i])
    JuMP.@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"], start=base_case.pg[i])
    JuMP.@variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"], start=base_case.qg[i])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"], start=base_case.p[(l, i, j)])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"], start=base_case.q[(l, i, j)])

    JuMP.@variable(model, 0.0 <= ρp[i in keys(ref[:gen])])
    JuMP.@variable(model, 0.0 <= ρn[i in keys(ref[:gen])])
    JuMP.@variable(model, 0.0 <= vp[i in keys(ref[:gen])])
    JuMP.@variable(model, 0.0 <= vn[i in keys(ref[:gen])])

    JuMP.@variable(model, Δ, start=0.0)

    # Feasibility problem.
    JuMP.@objective(model, Min, 0.0)

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
            sum(load_factor * load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2
        )

        JuMP.@constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(load_factor * load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2
        )
    end

    # Set flux to 0 at contingency
    for (l, i, j) in ref[:arcs]
        if l == contingency
            @constraint(model, p[(l, i, j)] == 0.0)
            @constraint(model, q[(l, i, j)] == 0.0)
        end
    end

    # Droop control
    for (j, i) in enumerate(keys(ref[:gen]))
        pmin, pmax = ref[:gen][i]["pmin"], ref[:gen][i]["pmax"]
        @constraint(model, ρp[i] - ρn[i] == pg[i] - base_case.pg[i] - alpha[j] * Δ)
        @constraint(model, ρn[i] * (pmax - pg[i]) <= 0.0)
        @constraint(model, ρp[i] * (pg[i] - pmin) <= 0.0)
    end

    # Voltage magnitude are not adjusted at PV buses
    for g in keys(ref[:gen])
        b = ref[:gen][g]["gen_bus"]
        qmin, qmax = ref[:gen][g]["qmin"], ref[:gen][g]["qmax"]
        @constraint(model, vp[g] - vn[g] == vm[b] - base_case.vm[b])
        if isfinite(qmax)
            @constraint(model, vn[g] * (qmax - qg[g]) <= 0.0)
        end
        if isfinite(qmin)
            @constraint(model, vp[g] * (qg[g] - qmin) <= 0.0)
        end
    end

    # Branch power flow physics and limit constraints
    for (i,branch) in ref[:branch]
        if i == contingency
            continue
        end
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


function screen_contingency(case; tol=1e-6, max_cont=Inf)
    # Base case
    model  = opf_model(case)
    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_attribute(model, "print_level", 0)
    JuMP.set_attribute(model, "max_iter", 500)
    JuMP.set_attribute(model, "hsllib", HSL_jll.libhsl_path)
    JuMP.set_attribute(model, "linear_solver", "ma27")
    JuMP.optimize!(model)
    base_case = extract_solution(model)

    # Contingency
    nl = div(length(model[:q]), 2)
    nK = isfinite(max_cont) ? min(nl, max_cont) : nl
    screen = zeros(nK, 4)
    @info "Screen contingency..."
    for k in 1:nK
        @printf "\rScreen contingency: %d / %d " k nK
        cnt_model = contingency_problem(case, base_case, k)

        nlp = ExaModels.ExaModel(cnt_model)
        snlp = MadNCL.ScaledModel(nlp)
        ncl = MadNCL.NCLModel(snlp)
        ncl.yk .= 0
        ncl.ρk[] = 1e6
        res = madnlp(
            ncl;
            print_level=MadNLP.ERROR,
            # linear_solver=MadNLPGPU.CUDSSSolver,
            linear_solver=Ma57Solver,
            # cholmod_algorithm=MadNLP.LDL,
            max_iter=200,
            tol=1e-8,
            kkt_system=MadNCL.K2rAuglagKKTSystem,
        )

        nvar = NLPModels.get_nvar(nlp)
        ncon = NLPModels.get_ncon(nlp)

        r = res.solution[nvar+1:nvar+ncon]
        screen[k, 1] = Int(res.status)
        screen[k, 2] = res.objective
        screen[k, 3] = norm(r, Inf)
        screen[k, 4] = length(findall(abs.(r) .> tol))
    end
    return screen
end

function demo()
    case = "../matpower/data/case9.m"

    model  = opf_model(case)
    JuMP.set_optimizer(model, MadNLP.Optimizer)
    JuMP.set_attribute(model, "linear_solver", Ma27Solver)
    JuMP.optimize!(model)
    base_case = extract_solution(model)

    index_contingency = 9
    cnt_model = contingency_problem(case, base_case, index_contingency)

    nlp = ExaModels.ExaModel(cnt_model)
    snlp = MadNCL.ScaledModel(nlp)
    ncl = MadNCL.NCLModel(snlp)
    ncl.yk .= 0
    ncl.ρk[] = 1e6
    res = madnlp(
        ncl;
        print_level=MadNLP.INFO,
        linear_solver=Ma27Solver,
        max_iter=200,
        tol=1e-8,
        kkt_system=MadNCL.K2rAuglagKKTSystem,
    )
end

