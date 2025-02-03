function opf_exa_rect_mp(filename, backend, process, tol)
    my_data = PowerModels.parse_file(filename)
    PowerModels.standardize_cost_terms!(my_data, order = 2)
    PowerModels.calc_thermal_limits!(my_data)
    ref = PowerModels.build_ref(my_data)[:it][:pm][:nw][0]

    interval_split = my_data["time_elapsed"]
    total_hrs = 24

    # from RTS 96 paper
    summer_wkdy_hour_scalar = [
        .64, .60, .58, .56, .56, .58, .64, .76, .87, .95, .99, 1.0, .99, 1.0, 1.0,
        .97, .96, .96, .93, .92, .92, .93, .87, .72, .64]

    get_profile=scale(interpolate(summer_wkdy_hour_scalar,BSpline(Linear())),0:total_hrs)


    summer_wkdy_qrtr_scalar = zeros(Int(total_hrs/interval_split))
    for i in 0:length(summer_wkdy_qrtr_scalar) - 1
        summer_wkdy_qrtr_scalar[i+1] = get_profile(i*interval_split)
    end


    arcdict = Dict(a => k for (k, a) in enumerate(ref[:arcs]))
    busdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:bus]))
    gendict = Dict(k => i for (i, (k, v)) in enumerate(ref[:gen]))
    branchdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:branch]))


    buses = [
            begin
                bus_loads = [ref[:load][l] for l in ref[:bus_loads][k]]
                bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][k]]
                pd = sum(load["pd"] for load in bus_loads; init = 0.0)*summer_wkdy_qrtr_scalar[t]
                gs = sum(shunt["gs"] for shunt in bus_shunts; init = 0.0)
                qd = sum(load["qd"] for load in bus_loads; init = 0.0)*summer_wkdy_qrtr_scalar[t]
                bs = sum(shunt["bs"] for shunt in bus_shunts; init = 0.0)
                (i = busdict[k], pd = pd, gs = gs, qd = qd, bs = bs, t=t)
            end for  t in eachindex(summer_wkdy_qrtr_scalar), (k, v) in ref[:bus]]
        
    println(typeof(buses), size(buses))
    gens = [
            (
                i = gendict[k],
                cost1 = v["cost"][1],
                cost2 = v["cost"][2],
                cost3 = v["cost"][3],
                bus = busdict[v["gen_bus"]],
                t = t
            ) for  t in eachindex(summer_wkdy_qrtr_scalar), (k, v) in ref[:gen]
        ]

    arcs = [
            (i = k, rate_a = ref[:branch][l]["rate_a"], bus = busdict[i], t=t) for
             t in eachindex(summer_wkdy_qrtr_scalar), (k, (l, i, j)) in enumerate(ref[:arcs])
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
                t = t
                )
            end  for t in eachindex(summer_wkdy_qrtr_scalar), (i, branch) in ref[:branch]
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

    model = ExaCore(Float64; backend = backend)

    vr = variable(model, length(busdict), length(summer_wkdy_qrtr_scalar), start = ones(size(buses)))
    vim = variable(model, length(busdict), length(summer_wkdy_qrtr_scalar))

    pg = variable(model, length(gendict), length(summer_wkdy_qrtr_scalar), lvar = repeat(pmin, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(pmax, length(summer_wkdy_qrtr_scalar), 1))
    qg = variable(model, length(gendict),length(summer_wkdy_qrtr_scalar), lvar = repeat(qmin, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(qmax, length(summer_wkdy_qrtr_scalar), 1))

    p = variable(model, length(arcdict), length(summer_wkdy_qrtr_scalar), lvar = -repeat(rate_a, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(rate_a, length(summer_wkdy_qrtr_scalar), 1))
    q = variable(model, length(arcdict), length(summer_wkdy_qrtr_scalar), lvar = -repeat(rate_a, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(rate_a, length(summer_wkdy_qrtr_scalar), 1))

    obj = objective(model, pg[g.i, g.t]^2*g.cost1+pg[g.i, g.t]*g.cost2+g.cost3 for g in gens)

    for i in ref_buses
        c0 = constraint(model,  atan(vim[i, t]/vr[i, t]) for t in eachindex(summer_wkdy_qrtr_scalar))
    end



    c1 = constraint(model, p[branch.f_idx, branch.t] -(  (branch.g + branch.g_fr)/branch.ttm*(vr[branch.f_bus, branch.t]^2+vim[branch.f_bus, branch.t]^2) + 
                    (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vr[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] + vim[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t])+
                    (-branch.b*branch.tr-branch.g*branch.ti)/branch.ttm*(vim[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] - vr[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t])) for branch in branches)


    c2 = constraint(model, q[branch.f_idx, branch.t] -( -(branch.b + branch.b_fr)/branch.ttm*(vr[branch.f_bus, branch.t]^2+vim[branch.f_bus, branch.t]^2) -
                    (-branch.b*branch.tr - branch.g*branch.ti)/branch.ttm*(vr[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] + vim[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t])+
                    (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vim[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] - vr[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t])) for branch in branches)


    c3 = constraint(model, p[branch.t_idx, branch.t] -( (branch.g+branch.g_to)*(vr[branch.t_bus, branch.t]^2+vim[branch.t_bus, branch.t]^2) +
                    (-branch.g*branch.tr-branch.b*branch.ti)/branch.ttm*(vr[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] + vim[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t]) + 
                    (-branch.b*branch.tr+branch.g*branch.ti)/branch.ttm*(vim[branch.t_bus, branch.t]*vr[branch.f_bus, branch.t] - vr[branch.t_bus, branch.t]*vim[branch.f_bus, branch.t])) for branch in branches)


    c4 = constraint(model, q[branch.t_idx, branch.t] -( -(branch.b+branch.b_to)*(vr[branch.t_bus, branch.t]^2+vim[branch.t_bus, branch.t]^2) - 
                    (-branch.b*branch.tr + branch.g*branch.ti)/branch.ttm* (vr[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] + vim[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t]) +
                    (-branch.g*branch.tr - branch.b*branch.ti)/branch.ttm * (vim[branch.t_bus, branch.t]*vr[branch.f_bus, branch.t] - vr[branch.t_bus, branch.t]*vim[branch.f_bus, branch.t])) for branch in branches)



    c5 = constraint(model, atan(vim[branch.f_bus, branch.t]/vr[branch.f_bus, branch.t]) - atan(vim[branch.t_bus, branch.t]/vr[branch.t_bus, branch.t]) for branch in branches; lcon = repeat(angmin, length(summer_wkdy_qrtr_scalar), 1), ucon = repeat(angmax, length(summer_wkdy_qrtr_scalar), 1))

    c6 =  constraint(model, p[branch.f_idx, branch.t]^2 + q[branch.f_idx, branch.t]^2 - branch.rate_a_sq for branch in branches; lcon = fill!(similar(branches, Float64, size(branches)), -Inf))

    c7 =  constraint(model, p[branch.t_idx, branch.t]^2 + q[branch.t_idx, branch.t]^2 - branch.rate_a_sq for branch in branches; lcon = fill!(similar(branches, Float64, size(branches)), -Inf))

    c8 = constraint(model, bus.pd + bus.gs*(vr[bus.i, bus.t]^2 + vim[bus.i, bus.t]^2) for bus in buses)
    c8a = constraint!(model, c8, (arc.t, arc.bus) => p[arc.i, arc.t] for arc in arcs)
    c8b = constraint!(model, c8, (gen.t, gen.bus) => -pg[gen.i, gen.t] for gen in gens)


    c9 = constraint(model, bus.qd - bus.bs*(vr[bus.i, bus.t]^2 + vim[bus.i, bus.t]^2) for bus in buses)
    c9a = constraint!(model, c9, (arc.t, arc.bus) => q[arc.i, arc.t] for arc in arcs)
    c9b = constraint!(model, c9, (gen.t, gen.bus) => -qg[gen.i, gen.t] for gen in gens)


    c10 = constraint(model, vr[bus.i, bus.t]^2+vim[bus.i, bus.t]^2 for bus in buses; lcon = repeat(vmin, length(summer_wkdy_qrtr_scalar), 1).^2, ucon = repeat(vmax, length(summer_wkdy_qrtr_scalar), 1).^2) 


    if process == "gpu"
        result = madnlp(ExaModel(model), tol = tol)
    elseif process == "cpu"
        result = ipopt(ExaModel(model), tol = tol, linear_solver = "ma27")
    end
    return [solution(result, vr), result.objective, solution(result, pg)] 
end