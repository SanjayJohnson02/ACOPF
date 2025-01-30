function opf_exa_rect(filename, backend, process, tol)
    my_data = PowerModels.parse_file(filename)
    PowerModels.standardize_cost_terms!(my_data, order = 2)
    PowerModels.calc_thermal_limits!(my_data)
    ref = PowerModels.build_ref(my_data)[:it][:pm][:nw][0]

    arcdict = Dict(a => k for (k, a) in enumerate(ref[:arcs]))
    busdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:bus]))
    gendict = Dict(k => i for (i, (k, v)) in enumerate(ref[:gen]))
    branchdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:branch]))


    buses = [
            begin
                bus_loads = [ref[:load][l] for l in ref[:bus_loads][k]]
                bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][k]]
                pd = sum(load["pd"] for load in bus_loads; init = 0.0)
                gs = sum(shunt["gs"] for shunt in bus_shunts; init = 0.0)
                qd = sum(load["qd"] for load in bus_loads; init = 0.0)
                bs = sum(shunt["bs"] for shunt in bus_shunts; init = 0.0)
                (i = busdict[k], pd = pd, gs = gs, qd = qd, bs = bs)
            end for (k, v) in ref[:bus]
        ]

    gens = [
            (
                i = gendict[k],
                cost1 = v["cost"][1],
                cost2 = v["cost"][2],
                cost3 = v["cost"][3],
                bus = busdict[v["gen_bus"]],
            ) for (k, v) in ref[:gen]
        ]

    arcs = [
            (i = k, rate_a = ref[:branch][l]["rate_a"], bus = busdict[i]) for
            (k, (l, i, j)) in enumerate(ref[:arcs])
        ]

    branches = [
            begin
                
                g, b = PowerModels.calc_branch_y(branch)
                tr, ti = PowerModels.calc_branch_t(branch)
                (f_idx = arcdict[i, branch["f_bus"], branch["t_bus"]],
                t_idx = arcdict[i, branch["t_bus"], branch["f_bus"]],
                g = g,
                b = b,
                tr = tr,
                ti = ti,
                ttm = tr^2 + ti^2,
                g_fr = branch["g_fr"],
                b_fr = branch["b_fr"],
                g_to = branch["g_to"],
                b_to = branch["b_to"],
                i = branchdict[i],
                f_bus = busdict[branch["f_bus"]],
                t_bus = busdict[branch["t_bus"]],
                rate_a_sq = branch["rate_a"]^2,
                )
            end for (i, branch) in ref[:branch]
        ]

    ref_buses = [busdict[i] for (i, k) in ref[:ref_buses]]
    vmax = [v["vmax"] for (k, v) in ref[:bus]]
    vmin = [v["vmin"] for (k, v) in ref[:bus]]
    pmax = [v["pmax"] for (k, v) in ref[:gen]]
    pmin = [v["pmin"] for (k, v) in ref[:gen]]
    qmax = [v["qmax"] for (k, v) in ref[:gen]]
    qmin = [v["qmin"] for (k, v) in ref[:gen]]
    rate_a = [ref[:branch][l]["rate_a"] for (k, (l, i, j)) in enumerate(ref[:arcs])]
    angmax = [b["angmax"] for (i, b) in ref[:branch]]
    angmin = [b["angmin"] for (i, b) in ref[:branch]]

    model = ExaCore(Float64; backend = nothing)

    #not sure if these should be initialized at nonzero
    vr = variable(model, length(buses), start = (ones(length(buses))))
    vim = variable(model, length(buses))

    pg = variable(model, length(gens), lvar = (pmin), uvar = (pmax))
    qg = variable(model, length(gens), lvar = (qmin), uvar = (qmax))

    p = variable(model, length(arcs), lvar = (-rate_a), uvar = (rate_a))
    q = variable(model, length(arcs), lvar = (-rate_a), uvar = (rate_a))

    obj = objective(model, pg[g.i]^2*g.cost1+pg[g.i]*g.cost2+g.cost3 for g in gens)



    c0 = constraint(model,  atan(vim[i]/vr[i]) for i in ref_buses)



    c1 = constraint(model, p[branch.f_idx] -(  (branch.g + branch.g_fr)/branch.ttm*(vr[branch.f_bus]^2+vim[branch.f_bus]^2) + 
                    (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vr[branch.f_bus]*vr[branch.t_bus] + vim[branch.f_bus]*vim[branch.t_bus])+
                    (-branch.b*branch.tr-branch.g*branch.ti)/branch.ttm*(vim[branch.f_bus]*vr[branch.t_bus] - vr[branch.f_bus]*vim[branch.t_bus])) for branch in branches)


    c2 = constraint(model, q[branch.f_idx] -( -(branch.b + branch.b_fr)/branch.ttm*(vr[branch.f_bus]^2+vim[branch.f_bus]^2) -
                    (-branch.b*branch.tr - branch.g*branch.ti)/branch.ttm*(vr[branch.f_bus]*vr[branch.t_bus] + vim[branch.f_bus]*vim[branch.t_bus])+
                    (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vim[branch.f_bus]*vr[branch.t_bus] - vr[branch.f_bus]*vim[branch.t_bus])) for branch in branches)


    c3 = constraint(model, p[branch.t_idx] -( (branch.g+branch.g_to)*(vr[branch.t_bus]^2+vim[branch.t_bus]^2) +
                    (-branch.g*branch.tr-branch.b*branch.ti)/branch.ttm*(vr[branch.f_bus]*vr[branch.t_bus] + vim[branch.f_bus]*vim[branch.t_bus]) + 
                    (-branch.b*branch.tr+branch.g*branch.ti)/branch.ttm*(vim[branch.t_bus]*vr[branch.f_bus] - vr[branch.t_bus]*vim[branch.f_bus])) for branch in branches)


    c4 = constraint(model, q[branch.t_idx] -( -(branch.b+branch.b_to)*(vr[branch.t_bus]^2+vim[branch.t_bus]^2) - 
                    (-branch.b*branch.tr + branch.g*branch.ti)/branch.ttm* (vr[branch.f_bus]*vr[branch.t_bus] + vim[branch.f_bus]*vim[branch.t_bus]) +
                    (-branch.g*branch.tr - branch.b*branch.ti)/branch.ttm * (vim[branch.t_bus]*vr[branch.f_bus] - vr[branch.t_bus]*vim[branch.f_bus])) for branch in branches)


    c5 = constraint(model, atan(vim[branch.f_bus]/vr[branch.f_bus]) - atan(vim[branch.t_bus]/vr[branch.t_bus]) for branch in branches; lcon = angmin, ucon = angmax)

    c6 =  constraint(model, p[branch.f_idx]^2 + q[branch.f_idx]^2 - branch.rate_a_sq for branch in branches; lcon = fill!(similar(branches, Float64, length(branches)), -Inf))

    c7 =  constraint(model, p[branch.t_idx]^2 + q[branch.t_idx]^2 - branch.rate_a_sq for branch in branches; lcon = fill!(similar(branches, Float64, length(branches)), -Inf))

    c8 = constraint(model, bus.pd + bus.gs*(vr[bus.i]^2 + vim[bus.i]^2) for bus in buses)
    c8a = constraint!(model, c8, arc.bus => p[arc.i] for arc in arcs)
    c8b = constraint!(model, c8, gen.bus => -pg[gen.i] for gen in gens)


    c9 = constraint(model, bus.qd - bus.bs*(vr[bus.i]^2 + vim[bus.i]^2) for bus in buses)
    c9a = constraint!(model, c9, arc.bus => q[arc.i] for arc in arcs)
    c9b = constraint!(model, c9, gen.bus => -qg[gen.i] for gen in gens)

    c10 = constraint(model, vr[bus.i]^2+vim[bus.i]^2 for bus in buses; lcon = vmin.^2, ucon = vmax.^2) 
    if process == "gpu"
        result = madnlp(ExaModel(model), tol = tol, linear_solver = Ma57Solver)
    elseif process == "cpu"
        result = ipopt(ExaModel(model), tol = tol, linear_solver = "ma57")
    end
    return [solution(result, vr), result.objective, solution(result, pg)]
end

filename = "acopf_scratch/pglib_opf_case1354_pegase.m"
cpu_sol = opf_exa_rect(filename, nothing, "cpu", 1e-6)
cpu_sol2 = opf_exa_rect(filename, nothing, "cpu", 1e-8)
gpu_sol = opf_exa_rect(filename, CUDABackend(), "gpu", 1e-6)
gpu_sol2 = opf_exa_rect(filename, CUDABackend(), "gpu", 1e-8)


println(filename)

println("Error in gpu solution between 1e-4 and 1e-8 tol")
println(norm(Array(gpu_sol[3]) - Array(gpu_sol2[3]))/norm(Array(gpu_sol[3])))
println("Error in cpu solution between 1e-4 and 1e-8 tol")
println(norm(Array(cpu_sol[3]) - Array(cpu_sol2[3]))/norm(Array(cpu_sol[3])))


println("tol for gpu tol = 1e-4")
println("vr tol ", norm(cpu_sol2[1] - Array(gpu_sol[1]))/norm(cpu_sol2[1]))
println("obj tol ", norm(cpu_sol2[2] - (gpu_sol[2]))/cpu_sol2[2])
println("pg tol ", norm(cpu_sol2[3] - Array(gpu_sol[3]))/norm(cpu_sol2[3]))

println("tol for gpu tol = 1e-8")
println("vr tol ", norm(cpu_sol2[1] - Array(gpu_sol2[1]))/norm(cpu_sol2[1]))
println("obj tol ", norm(cpu_sol2[2] - (gpu_sol2[2]))/cpu_sol2[2])
println("pg tol ", norm(cpu_sol2[3] - Array(gpu_sol2[3]))/norm(cpu_sol2[3]))

