#=
using ExaModels
using NLPModelsIpopt
using PowerModels
using Interpolations
using Ipopt 
using LinearAlgebra 
=# 

function opf_exa_rect_storage(filename)
    my_data = PowerModels.parse_file(filename)
    n_bus = length(my_data["bus"])
    n_storage = length(my_data["storage"])

    base_MVA = my_data["baseMVA"]

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

    # Create a list of dictionaries for each bus
    bus_list = [
        Dict("i" => 0, "pd" => zeros(length(summer_wkdy_qrtr_scalar)), "gs" => 0.0, "qd" => zeros(length(summer_wkdy_qrtr_scalar)), "bs" => 0, "bus_type" => 0, "arcs" => [], "gen_idx" => []) for _ in 1:n_bus
    ]


    vmaxs = zeros(n_bus)
    vmins = zeros(n_bus)


    for key in keys(my_data["bus"])
        i = parse(Int, key)
        bus_list[i]["i"] = i
        bus_list[i]["bus_type"] = my_data["bus"][key]["bus_type"]
        vmaxs[i] = my_data["bus"][key]["vmax"]
        vmins[i] = my_data["bus"][key]["vmin"]
        bus_list[i]["arcs"] = []
    end


    n_gen = length(my_data["gen"])

    #bus # aligns with bus list index
    gen_names = ["i", "cost1", "cost2", "cost3", "bus"]

    # Create a list of dictionaries for each bus
    gen_list = [
        Dict(name => 0.0 for name in gen_names) for _ in 1:n_gen
    ]

    pmins = zeros(n_gen)
    pmaxs = zeros(n_gen)
    qmins = zeros(n_gen)
    qmaxs = zeros(n_gen)


    for key in keys(my_data["gen"])
        i = parse(Int, key)
        gen_list[i]["i"] = i
        cost = my_data["gen"][key]["cost"]
        if length(cost) == 3
            gen_list[i]["cost1"] = cost[1]
            gen_list[i]["cost2"] = cost[2]
            gen_list[i]["cost3"] = cost[3]
        elseif length(cost) == 2
            gen_list[i]["cost2"] = cost[1]
            gen_list[i]["cost3"] = cost[2]
        elseif length(cost) == 1
            gen_list[i]["cost3"] = cost[1]
        end

        #hard coded for time interval

        push!(bus_list[my_data["gen"][key]["gen_bus"]]["gen_idx"], i)
        gen_list[i]["bus"] = my_data["gen"][key]["gen_bus"]
        pmins[i] = my_data["gen"][key]["pmin"]
        pmaxs[i] = my_data["gen"][key]["pmax"]
        qmins[i] = my_data["gen"][key]["qmin"]
        qmaxs[i] = my_data["gen"][key]["qmax"]
    end


    arc_list = []
    branch_list = []
    rate_as = Vector{Float64}()
    angmaxs = Vector{Float64}()
    angmins = Vector{Float64}()

    dict_names = ["i", "rate_a", "from_bus", "to_bus"]

    #need to consider arc, include a to and from dict, add a list to each bus indicating what arc indexes it uses

    #hard coded adjustment to rate a
    for key in keys(my_data["branch"])
        branch_dict = Dict{Any, Any}()
        branch_dict["f_bus"] = my_data["branch"][key]["f_bus"]
        branch_dict["t_bus"] = my_data["branch"][key]["t_bus"]
        branch_dict["a_rate_sq"] = (my_data["branch"][key]["rate_a"])^2
        y = pinv( my_data["branch"][key]["br_r"] + im* my_data["branch"][key]["br_x"])
        tr =  my_data["branch"][key]["tap"]*cos( my_data["branch"][key]["shift"])
        ti =  my_data["branch"][key]["tap"]*sin( my_data["branch"][key]["shift"])

        #for conciseness, could choose to add computations of these variables rather than output them raw
        branch_dict["g"] = real(y)
        branch_dict["b"] = imag(y)
        branch_dict["tr"] = tr
        branch_dict["ti"] = ti
        branch_dict["ttm"] = tr^2+ti^2
        branch_dict["g_fr"] = my_data["branch"][key]["g_fr"]
        branch_dict["g_to"] = my_data["branch"][key]["g_to"]
        branch_dict["b_fr"] = my_data["branch"][key]["b_fr"]
        branch_dict["b_to"] = my_data["branch"][key]["b_to"]

        arc_dict = Dict(name => 0.0 for name in dict_names) 
        arc_dict["rate_a"] = my_data["branch"][key]["rate_a"]
        arc_dict["from_bus"] = my_data["branch"][key]["f_bus"]
        arc_dict["to_bus"] = my_data["branch"][key]["t_bus"]
        arc_dict["i"] = length(arc_list) + 1
        push!(arc_list, arc_dict)
        push!(bus_list[my_data["branch"][key]["f_bus"]]["arcs"], length(arc_list))
        branch_dict["f_idx"] = length(arc_list)

        #done twice to account for both ways of the arc
        push!(rate_as, my_data["branch"][key]["rate_a"])
        push!(rate_as, my_data["branch"][key]["rate_a"])

        arc_dict = Dict(name => 0.0 for name in dict_names) 
        arc_dict["rate_a"] = my_data["branch"][key]["rate_a"]
        arc_dict["from_bus"] = my_data["branch"][key]["t_bus"]
        arc_dict["to_bus"] = my_data["branch"][key]["f_bus"]
        arc_dict["i"] = length(arc_list) + 1
        push!(arc_list, arc_dict)
        push!(bus_list[my_data["branch"][key]["t_bus"]]["arcs"], length(arc_list))
        branch_dict["t_idx"] = length(arc_list)

        push!(branch_list, branch_dict)
        push!(angmaxs, my_data["branch"][key]["angmax"])
        push!(angmins, my_data["branch"][key]["angmin"])

        
    end

    for key in keys(my_data["shunt"])
        bus_list[my_data["shunt"][key]["shunt_bus"]]["gs"] = my_data["shunt"][key]["gs"]
        bus_list[my_data["shunt"][key]["shunt_bus"]]["bs"] = my_data["shunt"][key]["bs"]
    end




    for key in keys(my_data["load"])
        pd_list = zeros(length(summer_wkdy_qrtr_scalar))
        qd_list = zeros(length(summer_wkdy_qrtr_scalar))
        for t in eachindex(summer_wkdy_qrtr_scalar)
            pd_list[t] = my_data["load"][key]["pd"]*summer_wkdy_qrtr_scalar[t]
            qd_list[t] = my_data["load"][key]["qd"]*summer_wkdy_qrtr_scalar[t]
        end
        bus_list[my_data["load"][key]["load_bus"]]["pd"] = pd_list
        bus_list[my_data["load"][key]["load_bus"]]["qd"] = qd_list
    end

    storage_list = [
        Dict("c" => 0, "Einit" => 0, "Emax" => 0.0, "Pcmax" => 0, "Pdmax" => 0, "etac" => 0, "etad" => 0, "Zr" => 0, "Zim" => 0, "Srating" => 0, "Pexts" => 0, "Qexts" => 0, "bus" => 0) for _ in 1:n_storage
    ]


    for key in eachindex(my_data["storage"])
        c = parse(Int, key)
        storage_list[c]["c"] = c
        storage_list[c]["Einit"] = my_data["storage"][key]["energy"]
        storage_list[c]["Emax"] = my_data["storage"][key]["energy_rating"]
        storage_list[c]["Pcmax"] = my_data["storage"][key]["charge_rating"]

        #this is slightly different from what paper states
        storage_list[c]["Pdmax"] = my_data["storage"][key]["discharge_rating"]
        #storage_list[i]["Pdmax"] = 0.75

        #these are swapped form data and paper?
        storage_list[c]["etac"] = my_data["storage"][key]["charge_efficiency"]
        #storage_list[i]["etac"] = .85
        storage_list[c]["etad"] = my_data["storage"][key]["discharge_efficiency"]
        #storage_list[i]["etad"] = 0.9

        storage_list[c]["Zr"] = my_data["storage"][key]["r"]
        storage_list[c]["Zim"] = my_data["storage"][key]["x"]
        storage_list[c]["Srating"] = my_data["storage"][key]["thermal_rating"]
        storage_list[c]["Pexts"] = my_data["storage"][key]["ps"]
        storage_list[c]["Qexts"] = my_data["storage"][key]["qs"]
        storage_list[c]["bus"] = my_data["storage"][key]["storage_bus"]

    end


    gen_tuple_array = Array{NamedTuple{(:i, :t, :cost1, :cost2, :cost3, :bus), Tuple{Int, Int, Float64, Float64, Float64, Int}}}(undef, length(summer_wkdy_qrtr_scalar), length(gen_list))
    for t in eachindex(summer_wkdy_qrtr_scalar)
        for dict in gen_list
            gen_tuple = (i = dict["i"], t = t, cost1 = dict["cost1"], cost2 = dict["cost2"], cost3 = dict["cost3"], bus = dict["bus"])
            gen_tuple_array[t, Int(dict["i"])] = gen_tuple
        end
    end

    bus_tuple_array = Array{NamedTuple{(:i, :t, :pd, :gs, :qd, :bs, :bus_type), Tuple{Int, Int, Float64, Float64, Float64, Float64, Float64}}}(undef, length(summer_wkdy_qrtr_scalar), length(bus_list))
    for t in eachindex(summer_wkdy_qrtr_scalar)
        for dict in bus_list
            bus_tuple = (i = Int(dict["i"]), t = t, pd = dict["pd"][t], gs = dict["gs"], qd = dict["qd"][t], bs = dict["bs"], bus_type = dict["bus_type"])
            bus_tuple_array[t, Int(dict["i"])] = bus_tuple
        end
    end

    branch_tuple_array = Array{NamedTuple{(:i, :t, :f_idx, :t_idx, :f_bus, :t_bus, :g, :b, :tr, :ti, :ttm, :g_fr, :g_to, :b_fr, :b_to, :a_rate_sq), Tuple{Int, Int, Int, Int, Int, Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}}(undef, length(summer_wkdy_qrtr_scalar), length(branch_list))
    for t in eachindex(summer_wkdy_qrtr_scalar)
        for i in eachindex(branch_list)
            branch_tuple = (i = i, t=t, f_idx = branch_list[i]["f_idx"], t_idx = branch_list[i]["t_idx"], f_bus = branch_list[i]["f_bus"], t_bus = branch_list[i]["t_bus"], g = branch_list[i]["g"], b = branch_list[i]["b"], tr = branch_list[i]["tr"], ti = branch_list[i]["ti"],
                            ttm = branch_list[i]["ttm"], g_fr = branch_list[i]["g_fr"], g_to = branch_list[i]["g_to"], b_fr = branch_list[i]["b_fr"], b_to = branch_list[i]["b_to"], a_rate_sq = branch_list[i]["a_rate_sq"])
            branch_tuple_array[t, i] = branch_tuple
        end
    end

    arc_tuple_array = Array{NamedTuple{(:i, :t, :from_bus, :to_bus, :rate_a), Tuple{Int,Int,  Int, Int, Float64}}}(undef, length(summer_wkdy_qrtr_scalar), length(arc_list))
    for t in eachindex(summer_wkdy_qrtr_scalar)
        for dict in arc_list
            arc_tuple = (i = Int(dict["i"]), t = t, from_bus = Int(dict["from_bus"]), to_bus = Int(dict["to_bus"]), rate_a = dict["rate_a"])
            arc_tuple_array[t, Int(dict["i"])] = arc_tuple
        end
    end

    stor_tuple_array = Array{NamedTuple{(:c, :t, :Einit, :Emax, :Pcmax, :Pdmax, :etac, :etad, :Zr, :Zim, :Srating, :Pexts, :Qexts, :bus), Tuple{Int,Int,  Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Int}}}(undef, length(summer_wkdy_qrtr_scalar), length(storage_list))
    for t in eachindex(summer_wkdy_qrtr_scalar)
        for dict in storage_list
            stor_tuple = (c = Int(dict["c"]), t = t, Einit = dict["Einit"], Emax = dict["Emax"], Pcmax = dict["Pcmax"], Pdmax = dict["Pdmax"], etac = dict["etac"], etad = dict["etad"], Zr = dict["Zr"], Zim = dict["Zim"], Srating = dict["Srating"], Pexts = dict["Pexts"], Qexts = dict["Qexts"], bus = Int(dict["bus"]))
            stor_tuple_array[t, Int(dict["c"])] = stor_tuple
        end
    end

    model = ExaCore(; backend = CUDABackend())

    #not sure if these should be initialized at nonzero
    vr = variable(model, length(bus_list), length(summer_wkdy_qrtr_scalar), start = ones(size(bus_tuple_array)))
    vim = variable(model, length(bus_list), length(summer_wkdy_qrtr_scalar))

    pg = variable(model, length(gen_list), length(summer_wkdy_qrtr_scalar), lvar = repeat(pmins, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(pmaxs, length(summer_wkdy_qrtr_scalar), 1))
    qg = variable(model, length(gen_list),length(summer_wkdy_qrtr_scalar), lvar = repeat(qmins, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(qmaxs, length(summer_wkdy_qrtr_scalar), 1))

    p = variable(model, length(arc_list), length(summer_wkdy_qrtr_scalar), lvar = -repeat(rate_as, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(rate_as, length(summer_wkdy_qrtr_scalar), 1))
    q = variable(model, length(arc_list), length(summer_wkdy_qrtr_scalar), lvar = -repeat(rate_as, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(rate_as, length(summer_wkdy_qrtr_scalar), 1))

    pstc = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar), lvar = zeros(size(stor_tuple_array)), uvar = [stor.Pcmax for stor in stor_tuple_array])

    
    pstd = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar), lvar = zeros(size(stor_tuple_array)), uvar = [stor.Pdmax for stor in stor_tuple_array])#[[stor.Pdmax for stor in stors] for stors in stor_tuple_array])


    pst = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar))

    qst = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar))
    I2 = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar), lvar = zeros(size(stor_tuple_array)))
    qint = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar), lvar = -[stor.Srating for stor in stor_tuple_array], uvar = [stor.Srating for stor  in stor_tuple_array])
    E = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar), lvar = zeros(size(stor_tuple_array)), uvar = [stor.Emax for stor in stor_tuple_array])
    obj = objective(model, pg[g.i, g.t]^2*g.cost1+pg[g.i, g.t]*g.cost2+g.cost3 for g in gen_tuple_array)



    #this is hard coded, in this model it doesn't matter but in the future it may?
    c0 = constraint(model,  atan(vim[4, t]/vr[4, t]) for t in eachindex(summer_wkdy_qrtr_scalar))



    c1 = constraint(model, p[branch.f_idx, branch.t] -(  (branch.g + branch.g_fr)/branch.ttm*(vr[branch.f_bus, branch.t]^2+vim[branch.f_bus, branch.t]^2) + 
                    (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vr[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] + vim[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t])+
                    (-branch.b*branch.tr-branch.g*branch.ti)/branch.ttm*(vim[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] - vr[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t])) for branch in branch_tuple_array)


    c2 = constraint(model, q[branch.f_idx, branch.t] -( -(branch.b + branch.b_fr)/branch.ttm*(vr[branch.f_bus, branch.t]^2+vim[branch.f_bus, branch.t]^2) -
                    (-branch.b*branch.tr - branch.g*branch.ti)/branch.ttm*(vr[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] + vim[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t])+
                    (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vim[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] - vr[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t])) for branch in branch_tuple_array)


    c3 = constraint(model, p[branch.t_idx, branch.t] -( (branch.g+branch.g_to)*(vr[branch.t_bus, branch.t]^2+vim[branch.t_bus, branch.t]^2) +
                    (-branch.g*branch.tr-branch.b*branch.ti)/branch.ttm*(vr[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] + vim[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t]) + 
                    (-branch.b*branch.tr+branch.g*branch.ti)/branch.ttm*(vim[branch.t_bus, branch.t]*vr[branch.f_bus, branch.t] - vr[branch.t_bus, branch.t]*vim[branch.f_bus, branch.t])) for branch in branch_tuple_array)


    c4 = constraint(model, q[branch.t_idx, branch.t] -( -(branch.b+branch.b_to)*(vr[branch.t_bus, branch.t]^2+vim[branch.t_bus, branch.t]^2) - 
                    (-branch.b*branch.tr + branch.g*branch.ti)/branch.ttm* (vr[branch.f_bus, branch.t]*vr[branch.t_bus, branch.t] + vim[branch.f_bus, branch.t]*vim[branch.t_bus, branch.t]) +
                    (-branch.g*branch.tr - branch.b*branch.ti)/branch.ttm * (vim[branch.t_bus, branch.t]*vr[branch.f_bus, branch.t] - vr[branch.t_bus, branch.t]*vim[branch.f_bus, branch.t])) for branch in branch_tuple_array)



    c5 = constraint(model, atan(vim[branch.f_bus, branch.t]/vr[branch.f_bus, branch.t]) - atan(vim[branch.t_bus, branch.t]/vr[branch.t_bus, branch.t]) for branch in branch_tuple_array; lcon = repeat(angmins, length(summer_wkdy_qrtr_scalar), 1), ucon = repeat(angmaxs, length(summer_wkdy_qrtr_scalar), 1))

    c6 =  constraint(model, p[branch.f_idx, branch.t]^2 + q[branch.f_idx, branch.t]^2 - branch.a_rate_sq for branch in branch_tuple_array; lcon = fill!(similar(branch_tuple_array, Float64, size(branch_tuple_array)), -Inf))

    c7 =  constraint(model, p[branch.t_idx, branch.t]^2 + q[branch.t_idx, branch.t]^2 - branch.a_rate_sq for branch in branch_tuple_array; lcon = fill!(similar(branch_tuple_array, Float64, size(branch_tuple_array)), -Inf))

    c8 = constraint(model, bus.pd + bus.gs*(vr[bus.i, bus.t]^2 + vim[bus.i, bus.t]^2) for bus in bus_tuple_array)
    c8a = constraint!(model, c8, (arc.t, arc.from_bus) => p[arc.i, arc.t] for arc in arc_tuple_array)
    c8b = constraint!(model, c8, (gen.t, gen.bus) => -pg[gen.i, gen.t] for gen in gen_tuple_array)
    c8c = constraint!(model, c8, (stor.t, stor.bus) => pst[stor.c, stor.t] for stor in stor_tuple_array)


    c9 = constraint(model, bus.qd - bus.bs*(vr[bus.i, bus.t]^2 + vim[bus.i, bus.t]^2) for bus in bus_tuple_array)
    c9a = constraint!(model, c9, (arc.t, arc.from_bus) => q[arc.i, arc.t] for arc in arc_tuple_array)
    c9b = constraint!(model, c9, (gen.t, gen.bus) => -qg[gen.i, gen.t] for gen in gen_tuple_array)
    c9c = constraint!(model, c9, (stor.t, stor.bus) => qst[stor.c, stor.t] for stor in stor_tuple_array)


    c10 = constraint(model, vr[bus.i, bus.t]^2+vim[bus.i, bus.t]^2 for bus in bus_tuple_array; lcon = repeat(vmins, length(summer_wkdy_qrtr_scalar), 1).^2, ucon = repeat(vmaxs, length(summer_wkdy_qrtr_scalar), 1).^2) 

    c11 = constraint(model, pst[stor.c, stor.t] + pstd[stor.c, stor.t] - pstc[stor.c, stor.t] - stor.Pexts - stor.Zr*I2[stor.c, stor.t] for stor in stor_tuple_array)

    c12 = constraint(model, qst[stor.c, stor.t] - qint[stor.c, stor.t] - stor.Qexts - stor.Zim*I2[stor.c, stor.t] for stor in stor_tuple_array)

    c13 = constraint(model, pst[stor.c, stor.t]^2 + qst[stor.c, stor.t]^2 - (vr[stor.bus, stor.t]^2 + vim[stor.bus, stor.t]^2)*I2[stor.c, stor.t] for stor in stor_tuple_array)

    c14 = constraint(model, E[stor.c, stor.t] - E[stor.c, stor.t - 1] - interval_split*(stor.etac*pstc[stor.c, stor.t] - pstd[stor.c, stor.t]/stor.etad) for stor in stor_tuple_array[2:length(summer_wkdy_qrtr_scalar), :])

    c15 = constraint(model, E[stor.c, stor.t] - stor.Einit - interval_split*(stor.etac*pstc[stor.c, stor.t] - pstd[stor.c, stor.t]/stor.etad) for stor in stor_tuple_array[1, :])

    c16 = constraint(model, pst[stor.c, stor.t]^2 + qst[stor.c, stor.t]^2 - stor.Srating^2 for stor in stor_tuple_array; lcon = fill!(similar(stor_tuple_array, Float64, size(stor_tuple_array)), -Inf))

    c17 = constraint(model, pstd[stor.c, stor.t] - pstc[stor.c, stor.t] for stor in stor_tuple_array; lcon = -[stor.Srating for stor in stor_tuple_array], ucon = [stor.Srating for stor in stor_tuple_array])

    c18 = constraint(model, pstc[stor.c, stor.t]*pstd[stor.c, stor.t] for stor in stor_tuple_array)
    return ExaModel(model)
end

model = opf_exa_rect_storage("pglib_opf_case14_ieee_mod.m")
madnlp(model)