#coords options: polar, rect
#model_type options: static, multiperiod, storage
#solve_method options: jump, exa

function opf_combined(filename, coords, model_type, solve_method)

    my_data = PowerModels.parse_file(filename)
    base_MVA = my_data["baseMVA"]
    n_bus = length(my_data["bus"])

    #hardcoded
    refbus_list = []

    if model_type == "multiperiod" || model_type == "storage"
        multiperiod = true
    else
        multiperiod = false
    end

    if multiperiod
        n_storage = length(my_data["storage"])
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
        Dict("i" => 0, "pd" => zeros(length(summer_wkdy_qrtr_scalar)), "gs" => 0.0, "qd" => zeros(length(summer_wkdy_qrtr_scalar)), "bs" => 0, "bus_type" => 0, "arcs" => [], "gen_idx" => [], "storage_idx" => []) for _ in 1:n_bus
    ]
    else
        bus_list = [
            Dict("i" => 0, "pd" => 0, "gs" => 0.0, "qd" => 0, "bs" => 0, "bus_type" => 0, "arcs" => [], "gen_idx" => []) for _ in 1:n_bus
        ]
    end

    vmaxs = zeros(n_bus)
    vmins = zeros(n_bus)


    for key in keys(my_data["bus"])
        i = parse(Int, key)
        bus_list[i]["i"] = i
        bus_list[i]["bus_type"] = my_data["bus"][key]["bus_type"]
        vmaxs[i] = my_data["bus"][key]["vmax"]
        vmins[i] = my_data["bus"][key]["vmin"]
        bus_list[i]["arcs"] = []

        if bus_list[i]["bus_type"] == 3
            push!(refbus_list, i)
        end
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

        
        gen_list[i]["pstart"] = my_data["gen"][key]["pg"]
        gen_list[i]["qstart"] = my_data["gen"][key]["qg"]

        push!(bus_list[my_data["gen"][key]["gen_bus"]]["gen_idx"], i)
        gen_list[i]["bus"] = my_data["gen"][key]["gen_bus"]
        pmins[i] = my_data["gen"][key]["pmin"]
        pmaxs[i] = my_data["gen"][key]["pmax"]
        qmins[i] = my_data["gen"][key]["qmin"]
        qmaxs[i] = my_data["gen"][key]["qmax"]
    end


    arc_list = []
    branch_list = []
    rate_as = []
    angmaxs = []
    angmins = []

    dict_names = ["i", "rate_a", "from_bus", "to_bus"]

    #need to consider arc, include a to and from dict, add a list to each bus indicating what arc indexes it uses

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



    if multiperiod
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
    else
        for key in keys(my_data["load"])
            bus_list[my_data["load"][key]["load_bus"]]["pd"] = my_data["load"][key]["pd"]
            bus_list[my_data["load"][key]["load_bus"]]["qd"] = my_data["load"][key]["qd"]
        end
    end

    if model_type == "storage"
        storage_list = [
            Dict("c" => 0, "Einit" => 0, "Emax" => 0.0, "Pcmax" => 0, "Pdmax" => 0, "etac" => 0, "etad" => 0, "Zr" => 0, "Zim" => 0, "Srating" => 0, "Pexts" => 0, "Qexts" => 0, "bus" => 0) for _ in 1:n_storage
        ]


        for key in eachindex(my_data["storage"])
            c = parse(Int, key)
            storage_list[c]["c"] = c
            storage_list[c]["Einit"] = my_data["storage"][key]["energy"]
            storage_list[c]["Emax"] = my_data["storage"][key]["energy_rating"]
            storage_list[c]["Pcmax"] = my_data["storage"][key]["charge_rating"]

            storage_list[c]["Pdmax"] = my_data["storage"][key]["discharge_rating"]

            storage_list[c]["etac"] = my_data["storage"][key]["charge_efficiency"]
            storage_list[c]["etad"] = my_data["storage"][key]["discharge_efficiency"]

            storage_list[c]["Zr"] = my_data["storage"][key]["r"]
            storage_list[c]["Zim"] = my_data["storage"][key]["x"]
            storage_list[c]["Srating"] = my_data["storage"][key]["thermal_rating"]
            storage_list[c]["Pexts"] = my_data["storage"][key]["ps"]
            storage_list[c]["Qexts"] = my_data["storage"][key]["qs"]
            storage_list[c]["bus"] = my_data["storage"][key]["storage_bus"]

            push!(bus_list[my_data["storage"][key]["storage_bus"]]["storage_idx"], c)
        end
    end

    if solve_method == "exa"

        model = ExaCore()

        #not sure if these should be initialized at nonzero
        if multiperiod
            gen_tuple_array = Array{NamedTuple{(:i, :t, :cost1, :cost2, :cost3, :bus, :pstart, :qstart), Tuple{Int, Int, Float64, Float64, Float64, Int, Float64, Float64}}}(undef, length(summer_wkdy_qrtr_scalar), length(gen_list))
            for t in eachindex(summer_wkdy_qrtr_scalar)
                for dict in gen_list
                    gen_tuple = (i = dict["i"], t = t, cost1 = dict["cost1"], cost2 = dict["cost2"], cost3 = dict["cost3"], bus = dict["bus"], pstart = dict["pstart"], qstart = dict["qstart"])
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

            refbus_tuple_array = Array{NamedTuple{(:i, :t, :bus), Tuple{Int,Int,  Int}}}(undef, length(summer_wkdy_qrtr_scalar), length(refbus_list))
            for t in eachindex(summer_wkdy_qrtr_scalar)
                for refbus_i in eachindex(refbus_list)
                    refbus_tuple = (i = refbus_i, t = t, bus = refbus_list[refbus_i])
                    refbus_tuple_array[t, refbus_i] = refbus_tuple
                end
            end

            if model_type == "storage"
                stor_tuple_array = Array{NamedTuple{(:c, :t, :Einit, :Emax, :Pcmax, :Pdmax, :etac, :etad, :Zr, :Zim, :Srating, :Pexts, :Qexts, :bus), Tuple{Int,Int,  Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Int}}}(undef, length(summer_wkdy_qrtr_scalar), length(storage_list))
                for t in eachindex(summer_wkdy_qrtr_scalar)
                    for dict in storage_list
                        stor_tuple = (c = Int(dict["c"]), t = t, Einit = dict["Einit"], Emax = dict["Emax"], Pcmax = dict["Pcmax"], Pdmax = dict["Pdmax"], etac = dict["etac"], etad = dict["etad"], Zr = dict["Zr"], Zim = dict["Zim"], Srating = dict["Srating"], Pexts = dict["Pexts"], Qexts = dict["Qexts"], bus = Int(dict["bus"]))
                        stor_tuple_array[t, Int(dict["c"])] = stor_tuple
                    end
                end
            end
            if coords == "rect"
                vr = variable(model, length(bus_list), length(summer_wkdy_qrtr_scalar), start = ones(size(bus_tuple_array)))
                vim = variable(model, length(bus_list), length(summer_wkdy_qrtr_scalar))
            elseif coords == "polar"
                vm = variable(model, length(bus_list), length(summer_wkdy_qrtr_scalar), start = ones(size(bus_tuple_array)), lvar = repeat(vmins, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(vmaxs, length(summer_wkdy_qrtr_scalar), 1))
                va = variable(model, length(bus_list), length(summer_wkdy_qrtr_scalar))      
            end          

            pg = variable(model, length(gen_list), length(summer_wkdy_qrtr_scalar), start = [gen.pstart for gen in gen_tuple_array], lvar = repeat(pmins, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(pmaxs, length(summer_wkdy_qrtr_scalar), 1))
            qg = variable(model, length(gen_list),length(summer_wkdy_qrtr_scalar), start = [gen.qstart for gen in gen_tuple_array], lvar = repeat(qmins, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(qmaxs, length(summer_wkdy_qrtr_scalar), 1))

            p = variable(model, length(arc_list), length(summer_wkdy_qrtr_scalar), lvar = -repeat(rate_as, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(rate_as, length(summer_wkdy_qrtr_scalar), 1))
            q = variable(model, length(arc_list), length(summer_wkdy_qrtr_scalar), lvar = -repeat(rate_as, length(summer_wkdy_qrtr_scalar), 1), uvar = repeat(rate_as, length(summer_wkdy_qrtr_scalar), 1))

            obj = objective(model, pg[g.i, g.t]^2*g.cost1+pg[g.i, g.t]*g.cost2+g.cost3 for g in gen_tuple_array)


            if coords == "rect"
                c0 = constraint(model,  atan(vim[refbus.bus, refbus.t]/vr[refbus.bus, refbus.t]) for refbus in refbus_tuple_array)



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


                c8 = constraint(model, bus.pd + bus.gs*(vr[bus.i, bus.t]^2 + vim[bus.i, bus.t]^2) for bus in bus_tuple_array)
                


                c9 = constraint(model, bus.qd - bus.bs*(vr[bus.i, bus.t]^2 + vim[bus.i, bus.t]^2) for bus in bus_tuple_array)
                


                c10 = constraint(model, vr[bus.i, bus.t]^2+vim[bus.i, bus.t]^2 for bus in bus_tuple_array; lcon = repeat(vmins, length(summer_wkdy_qrtr_scalar), 1).^2, ucon = repeat(vmaxs, length(summer_wkdy_qrtr_scalar), 1).^2) 
            
            elseif coords == "polar"
                c0 = constraint(model,  va[refbus.bus, refbus.t] for refbus in refbus_tuple_array)

                c1 = constraint(model, p[branch.f_idx, branch.t] -(  (branch.g + branch.g_fr)/branch.ttm*(vm[branch.f_bus, branch.t]^2) + 
                                (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vm[branch.f_bus, branch.t]*vm[branch.t_bus, branch.t]*cos(va[branch.f_bus, branch.t] - va[branch.t_bus, branch.t]))+
                                (-branch.b*branch.tr-branch.g*branch.ti)/branch.ttm*(vm[branch.f_bus, branch.t]*vm[branch.t_bus, branch.t]*sin(va[branch.f_bus, branch.t] - va[branch.t_bus, branch.t]))) for branch in branch_tuple_array)

                c2 = constraint(model, q[branch.f_idx, branch.t] -( -(branch.b + branch.b_fr)/branch.ttm*(vm[branch.f_bus, branch.t]^2) -
                                (-branch.b*branch.tr - branch.g*branch.ti)/branch.ttm*(vm[branch.f_bus, branch.t]*vm[branch.t_bus, branch.t]*cos(va[branch.f_bus, branch.t] - va[branch.t_bus, branch.t]))+
                                (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vm[branch.f_bus, branch.t]*vm[branch.t_bus, branch.t]*sin(va[branch.f_bus, branch.t] - va[branch.t_bus, branch.t]))) for branch in branch_tuple_array)
                                
                c3 = constraint(model, p[branch.t_idx, branch.t] -( (branch.g+branch.g_to)*(vm[branch.t_bus, branch.t]^2) +
                                (-branch.g*branch.tr-branch.b*branch.ti)/branch.ttm*(vm[branch.f_bus, branch.t]*vm[branch.t_bus, branch.t]*cos(va[branch.t_bus, branch.t] - va[branch.f_bus, branch.t])) + 
                                (-branch.b*branch.tr+branch.g*branch.ti)/branch.ttm*(vm[branch.f_bus, branch.t]*vm[branch.t_bus, branch.t]*sin(va[branch.t_bus, branch.t] - va[branch.f_bus, branch.t]))) for branch in branch_tuple_array)

                c4 = constraint(model, q[branch.t_idx, branch.t] -( -(branch.b+branch.b_to)*(vm[branch.t_bus, branch.t]^2) - 
                                (-branch.b*branch.tr + branch.g*branch.ti)/branch.ttm* (vm[branch.f_bus, branch.t]*vm[branch.t_bus, branch.t]*cos(va[branch.t_bus, branch.t] - va[branch.f_bus, branch.t])) +
                                (-branch.g*branch.tr - branch.b*branch.ti)/branch.ttm * (vm[branch.f_bus, branch.t]*vm[branch.t_bus, branch.t]*sin(va[branch.t_bus, branch.t] - va[branch.f_bus, branch.t]))) for branch in branch_tuple_array)

                c5 = constraint(model, va[branch.f_bus, branch.t] - va[branch.t_bus, branch.t] for branch in branch_tuple_array; lcon = repeat(angmins, length(summer_wkdy_qrtr_scalar), 1), ucon = repeat(angmaxs, length(summer_wkdy_qrtr_scalar), 1))

                c8 = constraint(model, bus.pd + bus.gs*(vm[bus.i, bus.t]^2) for bus in bus_tuple_array)

                c9 = constraint(model, bus.qd - bus.bs*(vm[bus.i, bus.t]^2) for bus in bus_tuple_array)
            end

            c6 =  constraint(model, p[branch.f_idx, branch.t]^2 + q[branch.f_idx, branch.t]^2 - branch.a_rate_sq for branch in branch_tuple_array; lcon = fill!(similar(branch_tuple_array, Float64, size(branch_tuple_array)), -Inf))

            c7 =  constraint(model, p[branch.t_idx, branch.t]^2 + q[branch.t_idx, branch.t]^2 - branch.a_rate_sq for branch in branch_tuple_array; lcon = fill!(similar(branch_tuple_array, Float64, size(branch_tuple_array)), -Inf))
            
            c8a = constraint!(model, c8, (arc.t, arc.from_bus) => p[arc.i, arc.t] for arc in arc_tuple_array)
            c8b = constraint!(model, c8, (gen.t, gen.bus) => -pg[gen.i, gen.t] for gen in gen_tuple_array)

            c9a = constraint!(model, c9, (arc.t, arc.from_bus) => q[arc.i, arc.t] for arc in arc_tuple_array)
            c9b = constraint!(model, c9, (gen.t, gen.bus) => -qg[gen.i, gen.t] for gen in gen_tuple_array)

            if model_type == "storage"
                pstc = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar), lvar = zeros(size(stor_tuple_array)), uvar = [stor.Pcmax for stor in stor_tuple_array])
    
                
                pstd = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar), lvar = zeros(size(stor_tuple_array)), uvar = [stor.Pdmax for stor in stor_tuple_array])#[[stor.Pdmax for stor in stors] for stors in stor_tuple_array])
    
    
                pst = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar))
    
                qst = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar))
                I2 = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar), lvar = zeros(size(stor_tuple_array)))
                qint = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar), lvar = -[stor.Srating for stor in stor_tuple_array], uvar = [stor.Srating for stor  in stor_tuple_array])
                E = variable(model, length(storage_list), length(summer_wkdy_qrtr_scalar), lvar = zeros(size(stor_tuple_array)), uvar = [stor.Emax for stor in stor_tuple_array])
                
                c8c = constraint!(model, c8, (stor.t, stor.bus) => pst[stor.c, stor.t] for stor in stor_tuple_array)
                c9c = constraint!(model, c9, (stor.t, stor.bus) => qst[stor.c, stor.t] for stor in stor_tuple_array)

                c11 = constraint(model, pst[stor.c, stor.t] + pstd[stor.c, stor.t] - pstc[stor.c, stor.t] - stor.Pexts - stor.Zr*I2[stor.c, stor.t] for stor in stor_tuple_array)

                c12 = constraint(model, qst[stor.c, stor.t] - qint[stor.c, stor.t] - stor.Qexts - stor.Zim*I2[stor.c, stor.t] for stor in stor_tuple_array)

                if coords == "rect"
                    c13 = constraint(model, pst[stor.c, stor.t]^2 + qst[stor.c, stor.t]^2 - (vr[stor.bus, stor.t]^2 + vim[stor.bus, stor.t]^2)*I2[stor.c, stor.t] for stor in stor_tuple_array)
                elseif coords == "polar"
                    c13 = constraint(model, pst[stor.c, stor.t]^2 + qst[stor.c, stor.t]^2 - (vm[stor.bus, stor.t]^2)*I2[stor.c, stor.t] for stor in stor_tuple_array)
                end

                c14 = constraint(model, E[stor.c, stor.t] - E[stor.c, stor.t - 1] - interval_split*(stor.etac*pstc[stor.c, stor.t] - pstd[stor.c, stor.t]/stor.etad) for stor in stor_tuple_array[2:length(summer_wkdy_qrtr_scalar), :])

                c15 = constraint(model, E[stor.c, stor.t] - stor.Einit - interval_split*(stor.etac*pstc[stor.c, stor.t] - pstd[stor.c, stor.t]/stor.etad) for stor in stor_tuple_array[1, :])

                c16 = constraint(model, pst[stor.c, stor.t]^2 + qst[stor.c, stor.t]^2 - stor.Srating^2 for stor in stor_tuple_array; lcon = fill!(similar(stor_tuple_array, Float64, size(stor_tuple_array)), -Inf))

                c17 = constraint(model, pstd[stor.c, stor.t] - pstc[stor.c, stor.t] for stor in stor_tuple_array; lcon = -[stor.Srating for stor in stor_tuple_array], ucon = [stor.Srating for stor in stor_tuple_array])

                c18 = constraint(model, pstc[stor.c, stor.t]*pstd[stor.c, stor.t] for stor in stor_tuple_array)
            end
        #if static
        else            
            gen_tuple_list = Vector{NamedTuple{(:i, :cost1, :cost2, :cost3, :bus, :pstart, :qstart), Tuple{Int, Float64, Float64, Float64, Int, Float64, Float64}}}()
            for dict in gen_list
                gen_tuple = (i = dict["i"], cost1 = dict["cost1"], cost2 = dict["cost2"], cost3 = dict["cost3"], bus = dict["bus"], pstart = dict["pstart"], qstart = dict["qstart"])
                push!(gen_tuple_list, gen_tuple)
            end


            bus_tuple_list = Vector{NamedTuple{(:i, :pd, :gs, :qd, :bs, :bus_type), Tuple{Int, Float64, Float64, Float64, Float64, Float64}}}()
            for dict in bus_list
                bus_tuple = (i = Int(dict["i"]), pd = dict["pd"], gs = dict["gs"], qd = dict["qd"], bs = dict["bs"], bus_type = dict["bus_type"])
                push!(bus_tuple_list, bus_tuple)
            end

            branch_tuple_list = Vector{NamedTuple{(:i, :f_idx, :t_idx, :f_bus, :t_bus, :g, :b, :tr, :ti, :ttm, :g_fr, :g_to, :b_fr, :b_to, :a_rate_sq), Tuple{Int, Int, Int, Int, Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}}()
            for i in eachindex(branch_list)
                branch_tuple = (i = i, f_idx = branch_list[i]["f_idx"], t_idx = branch_list[i]["t_idx"], f_bus = branch_list[i]["f_bus"], t_bus = branch_list[i]["t_bus"], g = branch_list[i]["g"], b = branch_list[i]["b"], tr = branch_list[i]["tr"], ti = branch_list[i]["ti"],
                                ttm = branch_list[i]["ttm"], g_fr = branch_list[i]["g_fr"], g_to = branch_list[i]["g_to"], b_fr = branch_list[i]["b_fr"], b_to = branch_list[i]["b_to"], a_rate_sq = branch_list[i]["a_rate_sq"])
                push!(branch_tuple_list, branch_tuple)
            end

            arc_tuple_list = Vector{NamedTuple{(:i, :from_bus, :to_bus, :rate_a), Tuple{Int, Int, Int, Float64}}}()
            for dict in arc_list
                arc_tuple = (i = Int(dict["i"]), from_bus = Int(dict["from_bus"]), to_bus = Int(dict["to_bus"]), rate_a = dict["rate_a"])
                push!(arc_tuple_list, arc_tuple)
            end

            if coords == "rect"
                vr = variable(model, length(bus_list), start = ones(length(bus_list)))
                vim = variable(model, length(bus_list))
            elseif coords == "polar"
                vm = variable(model, length(bus_list), start = ones(length(bus_list)), lvar = vmins, uvar = vmaxs)
                va = variable(model, length(bus_list))
            end

            pg = variable(model, length(gen_list), start = [gen.pstart for gen in gen_tuple_list], lvar = pmins, uvar = pmaxs)
            qg = variable(model, length(gen_list), start = [gen.qstart for gen in gen_tuple_list], lvar = qmins, uvar = qmaxs)

            p = variable(model, length(arc_list), lvar = -rate_as, uvar = rate_as)
            q = variable(model, length(arc_list), lvar = -rate_as, uvar = rate_as)

            obj = objective(model, pg[g.i]^2*g.cost1+pg[g.i]*g.cost2+g.cost3 for g in gen_tuple_list)

            if coords == "rect"

                c0 = constraint(model,  atan(vim[refbus]/vr[refbus]) for refbus in refbus_list)

                c1 = constraint(model, p[branch.f_idx] -(  (branch.g + branch.g_fr)/branch.ttm*(vr[branch.f_bus]^2+vim[branch.f_bus]^2) + 
                                (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vr[branch.f_bus]*vr[branch.t_bus] + vim[branch.f_bus]*vim[branch.t_bus])+
                                (-branch.b*branch.tr-branch.g*branch.ti)/branch.ttm*(vim[branch.f_bus]*vr[branch.t_bus] - vr[branch.f_bus]*vim[branch.t_bus])) for branch in branch_tuple_list)


                c2 = constraint(model, q[branch.f_idx] -( -(branch.b + branch.b_fr)/branch.ttm*(vr[branch.f_bus]^2+vim[branch.f_bus]^2) -
                                (-branch.b*branch.tr - branch.g*branch.ti)/branch.ttm*(vr[branch.f_bus]*vr[branch.t_bus] + vim[branch.f_bus]*vim[branch.t_bus])+
                                (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vim[branch.f_bus]*vr[branch.t_bus] - vr[branch.f_bus]*vim[branch.t_bus])) for branch in branch_tuple_list)


                c3 = constraint(model, p[branch.t_idx] -( (branch.g+branch.g_to)*(vr[branch.t_bus]^2+vim[branch.t_bus]^2) +
                                (-branch.g*branch.tr-branch.b*branch.ti)/branch.ttm*(vr[branch.f_bus]*vr[branch.t_bus] + vim[branch.f_bus]*vim[branch.t_bus]) + 
                                (-branch.b*branch.tr+branch.g*branch.ti)/branch.ttm*(vim[branch.t_bus]*vr[branch.f_bus] - vr[branch.t_bus]*vim[branch.f_bus])) for branch in branch_tuple_list)


                c4 = constraint(model, q[branch.t_idx] -( -(branch.b+branch.b_to)*(vr[branch.t_bus]^2+vim[branch.t_bus]^2) - 
                                (-branch.b*branch.tr + branch.g*branch.ti)/branch.ttm* (vr[branch.f_bus]*vr[branch.t_bus] + vim[branch.f_bus]*vim[branch.t_bus]) +
                                (-branch.g*branch.tr - branch.b*branch.ti)/branch.ttm * (vim[branch.t_bus]*vr[branch.f_bus] - vr[branch.t_bus]*vim[branch.f_bus])) for branch in branch_tuple_list)




                c5 = constraint(model, atan(vim[branch.f_bus]/vr[branch.f_bus]) - atan(vim[branch.t_bus]/vr[branch.t_bus]) for branch in branch_tuple_list; lcon = angmins, ucon = angmaxs)

                
                c8 = constraint(model, bus.pd + bus.gs*(vr[bus.i]^2 + vim[bus.i]^2) for bus in bus_tuple_list)


                c9 = constraint(model, bus.qd - bus.bs*(vr[bus.i]^2 + vim[bus.i]^2) for bus in bus_tuple_list)

                c10 = constraint(model, vr[bus.i]^2+vim[bus.i]^2 for bus in bus_tuple_list; lcon = vmins.^2, ucon = vmaxs.^2) 
            
            elseif coords == "polar"
                c0 = constraint(model,  va[refbus] for refbus in refbus_list)

                c1 = constraint(model, p[branch.f_idx] -(  (branch.g + branch.g_fr)/branch.ttm*(vm[branch.f_bus]^2) + 
                                (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vm[branch.f_bus]*vm[branch.t_bus]*cos(va[branch.f_bus] - va[branch.t_bus]))+
                                (-branch.b*branch.tr-branch.g*branch.ti)/branch.ttm*(vm[branch.f_bus]*vm[branch.t_bus]*sin(va[branch.f_bus] - va[branch.t_bus]))) for branch in branch_tuple_list)

                c2 = constraint(model, q[branch.f_idx] -( -(branch.b + branch.b_fr)/branch.ttm*(vm[branch.f_bus]^2)-
                                (-branch.b*branch.tr - branch.g*branch.ti)/branch.ttm*(vm[branch.f_bus]*vm[branch.t_bus]*cos(va[branch.f_bus] - va[branch.t_bus]))+
                                (-branch.g*branch.tr + branch.b*branch.ti)/branch.ttm*(vm[branch.f_bus]*vm[branch.t_bus]*sin(va[branch.f_bus] - va[branch.t_bus]))) for branch in branch_tuple_list)

                c3 = constraint(model, p[branch.t_idx] -( (branch.g+branch.g_to)*(vm[branch.t_bus]^2) +
                                (-branch.g*branch.tr-branch.b*branch.ti)/branch.ttm*(vm[branch.f_bus]*vm[branch.t_bus]*cos(va[branch.t_bus] - va[branch.f_bus])) + 
                                (-branch.b*branch.tr+branch.g*branch.ti)/branch.ttm*(vm[branch.f_bus]*vm[branch.t_bus]*sin(va[branch.t_bus] - va[branch.f_bus]))) for branch in branch_tuple_list)

                c4 = constraint(model, q[branch.t_idx] -( -(branch.b+branch.b_to)*(vm[branch.t_bus]^2) - 
                                (-branch.b*branch.tr + branch.g*branch.ti)/branch.ttm* (vm[branch.f_bus]*vm[branch.t_bus]*cos(va[branch.t_bus] - va[branch.f_bus])) +
                                (-branch.g*branch.tr - branch.b*branch.ti)/branch.ttm * (vm[branch.f_bus]*vm[branch.t_bus]*sin(va[branch.t_bus] - va[branch.f_bus])) ) for branch in branch_tuple_list)

                c5 = constraint(model, va[branch.f_bus] - va[branch.t_bus] for branch in branch_tuple_list; lcon = angmins, ucon = angmaxs)

                c8 = constraint(model, bus.pd + bus.gs*(vm[bus.i]^2) for bus in bus_tuple_list)
                c9 = constraint(model, bus.qd - bus.bs*(vm[bus.i]^2) for bus in bus_tuple_list)

            end

            c6 =  constraint(model, p[branch.f_idx]^2 + q[branch.f_idx]^2 - branch.a_rate_sq for branch in branch_tuple_list; lcon = fill!(similar(branch_tuple_list, Float64, length(branch_tuple_list)), -Inf))

            c7 =  constraint(model, p[branch.t_idx]^2 + q[branch.t_idx]^2 - branch.a_rate_sq for branch in branch_tuple_list; lcon = fill!(similar(branch_tuple_list, Float64, length(branch_tuple_list)), -Inf))

            
            c8a = constraint!(model, c8, arc.from_bus => p[arc.i] for arc in arc_tuple_list)
            c8b = constraint!(model, c8, gen.bus => -pg[gen.i] for gen in gen_tuple_list)
    
            
            c9a = constraint!(model, c9, arc.from_bus => q[arc.i] for arc in arc_tuple_list)
            c9b = constraint!(model, c9, gen.bus => -qg[gen.i] for gen in gen_tuple_list)

        end
        return ExaModel(model)
    elseif solve_method == "jump"
        model = Model(Ipopt.Optimizer)

        if multiperiod
            if coords == "rect"
                vr = @variable(model, vr[i = 1:length(bus_list), t = 1:length(summer_wkdy_qrtr_scalar)], start = 1)
                vim = @variable(model, vim[i=1:length(bus_list), t = 1:length(summer_wkdy_qrtr_scalar)])
            elseif coords == "polar"
                vm = @variable(model, vmins[i] <= vm[i = 1:length(bus_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= vmaxs[i], start = 1)
                va = @variable(model, va[i=1:length(bus_list), t = 1:length(summer_wkdy_qrtr_scalar)])
            end

            pg = @variable(model, pmins[i] <= pg[i = 1:length(gen_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= pmaxs[i], start = gen_list[i]["pstart"])

            qg = @variable(model, qmins[i] <= qg[i =1:length(gen_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= qmaxs[i], start = gen_list[i]["qstart"])

            p = @variable(model, -rate_as[i] <= p[i =1:length(arc_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= rate_as[i])

            q = @variable(model, -rate_as[i] <= q[i=1:length(arc_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= rate_as[i])


            obj = @objective(model, Min, sum(pg[i,t]^2*gen_list[i]["cost1"]+pg[i,t]*gen_list[i]["cost2"]+gen_list[i]["cost3"] for i in eachindex(gen_list), t in eachindex(summer_wkdy_qrtr_scalar)))

            if coords == "rect"
                c0 = @constraint(model,  [i = eachindex(refbus_list), t = eachindex(summer_wkdy_qrtr_scalar)], atan(vim[i, t]/vr[i, t]) == 0)

                c1 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], p[branch_list[i]["f_idx"], t] -(  (branch_list[i]["g"] + branch_list[i]["g_fr"])/branch_list[i]["ttm"]*(vr[branch_list[i]["f_bus"], t]^2+vim[branch_list[i]["f_bus"], t]^2) + 
                                (-branch_list[i]["g"]*branch_list[i]["tr"] + branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vr[branch_list[i]["f_bus"], t]*vr[branch_list[i]["t_bus"], t] + vim[branch_list[i]["f_bus"], t]*vim[branch_list[i]["t_bus"], t])+
                                (-branch_list[i]["b"]*branch_list[i]["tr"]-branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vim[branch_list[i]["f_bus"], t]*vr[branch_list[i]["t_bus"], t] - vr[branch_list[i]["f_bus"], t]*vim[branch_list[i]["t_bus"], t])) == 0)


                c2 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], q[branch_list[i]["f_idx"], t] -( -(branch_list[i]["b"] + branch_list[i]["b_fr"])/branch_list[i]["ttm"]*(vr[branch_list[i]["f_bus"], t]^2+vim[branch_list[i]["f_bus"], t]^2) -
                                (-branch_list[i]["b"]*branch_list[i]["tr"] - branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vr[branch_list[i]["f_bus"], t]*vr[branch_list[i]["t_bus"], t] + vim[branch_list[i]["f_bus"], t]*vim[branch_list[i]["t_bus"], t])+
                                (-branch_list[i]["g"]*branch_list[i]["tr"] + branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vim[branch_list[i]["f_bus"], t]*vr[branch_list[i]["t_bus"], t] - vr[branch_list[i]["f_bus"], t]*vim[branch_list[i]["t_bus"], t])) == 0)


                c3 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], p[branch_list[i]["t_idx"], t] -( (branch_list[i]["g"]+branch_list[i]["g_to"])*(vr[branch_list[i]["t_bus"], t]^2+vim[branch_list[i]["t_bus"], t]^2) +
                                (-branch_list[i]["g"]*branch_list[i]["tr"]-branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vr[branch_list[i]["f_bus"], t]*vr[branch_list[i]["t_bus"], t] + vim[branch_list[i]["f_bus"], t]*vim[branch_list[i]["t_bus"], t]) + 
                                (-branch_list[i]["b"]*branch_list[i]["tr"]+branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vim[branch_list[i]["t_bus"], t]*vr[branch_list[i]["f_bus"], t] - vr[branch_list[i]["t_bus"], t]*vim[branch_list[i]["f_bus"], t])) == 0)


                c4 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], q[branch_list[i]["t_idx"], t] -( -(branch_list[i]["b"]+branch_list[i]["b_to"])*(vr[branch_list[i]["t_bus"], t]^2+vim[branch_list[i]["t_bus"], t]^2) - 
                                (-branch_list[i]["b"]*branch_list[i]["tr"] + branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]* (vr[branch_list[i]["f_bus"], t]*vr[branch_list[i]["t_bus"], t] + vim[branch_list[i]["f_bus"], t]*vim[branch_list[i]["t_bus"], t]) +
                                (-branch_list[i]["g"]*branch_list[i]["tr"] - branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"] * (vim[branch_list[i]["t_bus"], t]*vr[branch_list[i]["f_bus"], t] - vr[branch_list[i]["t_bus"], t]*vim[branch_list[i]["f_bus"], t])) == 0)

                c5 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], angmins[i] <= atan(vim[branch_list[i]["f_bus"], t]/vr[branch_list[i]["f_bus"], t]) - atan(vim[branch_list[i]["t_bus"], t]/vr[branch_list[i]["t_bus"], t]) <= angmaxs[i])
                

                c10 = @constraint(model, [i = eachindex(bus_list), t = eachindex(summer_wkdy_qrtr_scalar)], vmins[i]^2 <= vr[i, t]^2+vim[i, t]^2 <= vmaxs[i]^2)
            elseif coords == "polar"
                c0 = @constraint(model,  [i = eachindex(refbus_list), t = eachindex(summer_wkdy_qrtr_scalar)],va[i, t] == 0)

                c1 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], p[branch_list[i]["f_idx"], t] -(  (branch_list[i]["g"] + branch_list[i]["g_fr"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]^2) + 
                                (-branch_list[i]["g"]*branch_list[i]["tr"] + branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*cos(va[branch_list[i]["f_bus"], t] - va[branch_list[i]["t_bus"], t]))+
                                (-branch_list[i]["b"]*branch_list[i]["tr"]-branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*sin(va[branch_list[i]["f_bus"], t] - va[branch_list[i]["t_bus"], t]))) == 0)

                c2 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], q[branch_list[i]["f_idx"], t] -( -(branch_list[i]["b"] + branch_list[i]["b_fr"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]^2) -
                                (-branch_list[i]["b"]*branch_list[i]["tr"] - branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*cos(va[branch_list[i]["f_bus"], t] - va[branch_list[i]["t_bus"], t]))+
                                (-branch_list[i]["g"]*branch_list[i]["tr"] + branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*sin(va[branch_list[i]["f_bus"], t] - va[branch_list[i]["t_bus"], t]))) == 0)

                c3 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], p[branch_list[i]["t_idx"], t] -( (branch_list[i]["g"]+branch_list[i]["g_to"])*(vm[branch_list[i]["t_bus"], t]^2) +
                                (-branch_list[i]["g"]*branch_list[i]["tr"]-branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*cos(va[branch_list[i]["t_bus"], t] - va[branch_list[i]["f_bus"], t])) + 
                                (-branch_list[i]["b"]*branch_list[i]["tr"]+branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*sin(va[branch_list[i]["t_bus"], t] - va[branch_list[i]["f_bus"], t]))) == 0)
                
                c4 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], q[branch_list[i]["t_idx"], t] -( -(branch_list[i]["b"]+branch_list[i]["b_to"])*(vm[branch_list[i]["t_bus"], t]^2) - 
                                (-branch_list[i]["b"]*branch_list[i]["tr"] + branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*cos(va[branch_list[i]["t_bus"], t] - va[branch_list[i]["f_bus"], t])) +
                                (-branch_list[i]["g"]*branch_list[i]["tr"] - branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"] *(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*sin(va[branch_list[i]["t_bus"], t] - va[branch_list[i]["f_bus"], t]))) == 0)

                c5 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], angmins[i] <= va[branch_list[i]["f_bus"], t] - va[branch_list[i]["f_bus"], t] <= angmaxs[i])


            end 

            c6 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], p[branch_list[i]["f_idx"], t]^2 + q[branch_list[i]["f_idx"], t]^2 - branch_list[i]["a_rate_sq"] <= 0)

            c7 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], p[branch_list[i]["t_idx"], t]^2 + q[branch_list[i]["t_idx"], t]^2 - branch_list[i]["a_rate_sq"] <= 0)

            

            if model_type == "storage"

                pstc = @variable(model, 0 <= pstc[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= storage_list[i]["Pcmax"])
    
                pstd = @variable(model, 0 <= pstd[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= storage_list[i]["Pdmax"])
    
                pst = @variable(model, pst[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)])
    
                #this is maybe supposed to have bounds??
                qst = @variable(model, qst[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)])
    
                I2 = @variable(model, 0 <= I2[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)])
    
                qint = @variable(model, -storage_list[i]["Srating"] <= qint[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= storage_list[i]["Srating"])
    
                E = @variable(model, 0 <= E[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= storage_list[i]["Emax"])
    
                #only for use in MINLP version
                #zc = @variable(model, zc[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)], Bin)

                if coords == "rect"
                    c8 = @constraint(model, [i = eachindex(bus_list), t = eachindex(summer_wkdy_qrtr_scalar)], bus_list[i]["pd"][t] + bus_list[i]["gs"]*(vr[i, t]^2+vim[i, t]^2) + sum(p[j, t] for j in bus_list[i]["arcs"]) - 
                                sum(pg[k, t] for k in bus_list[i]["gen_idx"]) + sum(pst[l, t] for l in bus_list[i]["storage_idx"]) == 0)

                    c9 = @constraint(model, [i = eachindex(bus_list), t = eachindex(summer_wkdy_qrtr_scalar)], bus_list[i]["qd"][t] - bus_list[i]["bs"]*(vr[i, t]^2+vim[i, t]^2)  + sum(q[j, t] for j in bus_list[i]["arcs"]) - 
                                sum(qg[k, t] for k in bus_list[i]["gen_idx"]) + sum(qst[l, t] for l in bus_list[i]["storage_idx"])== 0)
                    
                    c13 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pst[c, t]^2 + qst[c, t]^2 == (vr[storage_list[c]["bus"], t]^2 + vim[storage_list[c]["bus"], t]^2)*I2[c, t])
                elseif coords == "polar"
                    c8 = @constraint(model, [i = eachindex(bus_list), t = eachindex(summer_wkdy_qrtr_scalar)], bus_list[i]["pd"][t] + bus_list[i]["gs"]*(vm[i, t]^2) + sum(p[j, t] for j in bus_list[i]["arcs"]) - 
                                sum(pg[k, t] for k in bus_list[i]["gen_idx"]) + sum(pst[l, t] for l in bus_list[i]["storage_idx"]) == 0)

                    c9 = @constraint(model, [i = eachindex(bus_list), t = eachindex(summer_wkdy_qrtr_scalar)], bus_list[i]["qd"][t] - bus_list[i]["bs"]*(vm[i, t]^2)  + sum(q[j, t] for j in bus_list[i]["arcs"]) - 
                                sum(qg[k, t] for k in bus_list[i]["gen_idx"]) + sum(qst[l, t] for l in bus_list[i]["storage_idx"])== 0)
                    
                    c13 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pst[c, t]^2 + qst[c, t]^2 == (vm[storage_list[c]["bus"], t]^2)*I2[c, t])
                end
                c11 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pst[c, t] + pstd[c, t] - pstc[c, t] == storage_list[c]["Pexts"] + storage_list[c]["Zr"]*I2[c, t])

                c12 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], qst[c, t] == qint[c, t] + storage_list[c]["Qexts"] + storage_list[c]["Zim"]*I2[c, t])


                c14 = @constraint(model, [c = eachindex(storage_list), t = 2:length(summer_wkdy_qrtr_scalar)], E[c, t] - E[c, t-1] == interval_split*(storage_list[c]["etac"]*pstc[c, t]- pstd[c, t]/storage_list[c]["etad"]))

                c15 = @constraint(model, [c = eachindex(storage_list)], E[c, 1] - storage_list[c]["Einit"] == interval_split*(storage_list[c]["etac"]*pstc[c, 1]- pstd[c, 1]/storage_list[c]["etad"]))

                c16 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pst[c, t]^2 + qst[c, t]^2 <= storage_list[c]["Srating"]^2)

                c17 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], -storage_list[c]["Srating"] <= pstd[c, t] - pstc[c, t] <= storage_list[c]["Srating"])

                c18 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pstc[c, t]*pstd[c, t] == 0)
            else
                if coords == "rect"
                    c8 = @constraint(model, [i = eachindex(bus_list), t = eachindex(summer_wkdy_qrtr_scalar)], bus_list[i]["pd"][t] + bus_list[i]["gs"]*(vr[i, t]^2+vim[i, t]^2) + sum(p[j, t] for j in bus_list[i]["arcs"]) - 
                                sum(pg[k, t] for k in bus_list[i]["gen_idx"]) == 0)

                    c9 = @constraint(model, [i = eachindex(bus_list), t = eachindex(summer_wkdy_qrtr_scalar)], bus_list[i]["qd"][t] - bus_list[i]["bs"]*(vr[i, t]^2+vim[i, t]^2)  + sum(q[j, t] for j in bus_list[i]["arcs"]) - 
                                sum(qg[k, t] for k in bus_list[i]["gen_idx"])== 0)
                elseif coords == "polar"
                    c8 = @constraint(model, [i = eachindex(bus_list), t = eachindex(summer_wkdy_qrtr_scalar)], bus_list[i]["pd"][t] + bus_list[i]["gs"]*(vm[i, t]^2) + sum(p[j, t] for j in bus_list[i]["arcs"]) - 
                                sum(pg[k, t] for k in bus_list[i]["gen_idx"]) == 0)
                    c9 = @constraint(model, [i = eachindex(bus_list), t = eachindex(summer_wkdy_qrtr_scalar)], bus_list[i]["qd"][t] - bus_list[i]["bs"]*(vm[i, t]^2)  + sum(q[j, t] for j in bus_list[i]["arcs"]) - 
                                sum(qg[k, t] for k in bus_list[i]["gen_idx"])== 0)
                end
            end
        elseif model_type == "static"

            if coords == "rect"
                vr = @variable(model, vr[i = 1:length(bus_list)], start = 1)
                vim = @variable(model, vim[i=1:length(bus_list)])
            elseif coords == "polar"
                vm = @variable(model, vmins[i] <= vm[i = 1:length(bus_list)] <= vmaxs[i], start = 1)
                va = @variable(model, va[i=1:length(bus_list)])
            end

            pg = @variable(model, pmins[i] <= pg[i = 1:length(gen_list)] <= pmaxs[i], start = gen_list[i]["pstart"])
        
            qg = @variable(model, qmins[i] <= qg[i =1:length(gen_list)] <= qmaxs[i], start = gen_list[i]["qstart"])
        
            p = @variable(model, -rate_as[i] <= p[i =1:length(arc_list)] <= rate_as[i])
        
            q = @variable(model, -rate_as[i] <= q[i=1:length(arc_list)] <= rate_as[i])
        
            obj = @objective(model, Min, sum(pg[i]^2*gen_list[i]["cost1"]+pg[i]*gen_list[i]["cost2"]+gen_list[i]["cost3"] for i in eachindex(pg)))
        
        
            if coords == "rect"
                c0 = @constraint(model,  [i = eachindex(refbus_list)], atan(vim[i]/vr[i]) == 0)
            
            
            
                c1 = @constraint(model, [i = eachindex(branch_list)], p[branch_list[i]["f_idx"]] -(  (branch_list[i]["g"] + branch_list[i]["g_fr"])/branch_list[i]["ttm"]*(vr[branch_list[i]["f_bus"]]^2+vim[branch_list[i]["f_bus"]]^2) + 
                                (-branch_list[i]["g"]*branch_list[i]["tr"] + branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vr[branch_list[i]["f_bus"]]*vr[branch_list[i]["t_bus"]] + vim[branch_list[i]["f_bus"]]*vim[branch_list[i]["t_bus"]])+
                                (-branch_list[i]["b"]*branch_list[i]["tr"]-branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vim[branch_list[i]["f_bus"]]*vr[branch_list[i]["t_bus"]] - vr[branch_list[i]["f_bus"]]*vim[branch_list[i]["t_bus"]])) == 0)
            
            
                c2 = @constraint(model, [i = eachindex(branch_list)], q[branch_list[i]["f_idx"]] -( -(branch_list[i]["b"] + branch_list[i]["b_fr"])/branch_list[i]["ttm"]*(vr[branch_list[i]["f_bus"]]^2+vim[branch_list[i]["f_bus"]]^2) -
                                (-branch_list[i]["b"]*branch_list[i]["tr"] - branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vr[branch_list[i]["f_bus"]]*vr[branch_list[i]["t_bus"]] + vim[branch_list[i]["f_bus"]]*vim[branch_list[i]["t_bus"]])+
                                (-branch_list[i]["g"]*branch_list[i]["tr"] + branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vim[branch_list[i]["f_bus"]]*vr[branch_list[i]["t_bus"]] - vr[branch_list[i]["f_bus"]]*vim[branch_list[i]["t_bus"]])) == 0)
            
            
                c3 = @constraint(model, [i = eachindex(branch_list)], p[branch_list[i]["t_idx"]] -( (branch_list[i]["g"]+branch_list[i]["g_to"])*(vr[branch_list[i]["t_bus"]]^2+vim[branch_list[i]["t_bus"]]^2) +
                                (-branch_list[i]["g"]*branch_list[i]["tr"]-branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vr[branch_list[i]["f_bus"]]*vr[branch_list[i]["t_bus"]] + vim[branch_list[i]["f_bus"]]*vim[branch_list[i]["t_bus"]]) + 
                                (-branch_list[i]["b"]*branch_list[i]["tr"]+branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vim[branch_list[i]["t_bus"]]*vr[branch_list[i]["f_bus"]] - vr[branch_list[i]["t_bus"]]*vim[branch_list[i]["f_bus"]])) == 0)
            
            
                c4 = @constraint(model, [i = eachindex(branch_list)], q[branch_list[i]["t_idx"]] -( -(branch_list[i]["b"]+branch_list[i]["b_to"])*(vr[branch_list[i]["t_bus"]]^2+vim[branch_list[i]["t_bus"]]^2) - 
                                (-branch_list[i]["b"]*branch_list[i]["tr"] + branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]* (vr[branch_list[i]["f_bus"]]*vr[branch_list[i]["t_bus"]] + vim[branch_list[i]["f_bus"]]*vim[branch_list[i]["t_bus"]]) +
                                (-branch_list[i]["g"]*branch_list[i]["tr"] - branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"] * (vim[branch_list[i]["t_bus"]]*vr[branch_list[i]["f_bus"]] - vr[branch_list[i]["t_bus"]]*vim[branch_list[i]["f_bus"]])) == 0)
            
                c5 = @constraint(model, [i = eachindex(branch_list)], angmins[i] <= atan(vim[branch_list[i]["f_bus"]]/vr[branch_list[i]["f_bus"]]) - atan(vim[branch_list[i]["t_bus"]]/vr[branch_list[i]["t_bus"]]) <= angmaxs[i])
        
                c8 = @constraint(model, [i = eachindex(bus_list)], bus_list[i]["pd"] + bus_list[i]["gs"]*(vr[i]^2+vim[i]^2) + sum(p[j] for j in bus_list[i]["arcs"]) - sum(pg[k] for k in bus_list[i]["gen_idx"]) == 0)
        
                c9 = @constraint(model, [i = eachindex(bus_list)], bus_list[i]["qd"] - bus_list[i]["bs"]*(vr[i]^2+vim[i]^2)  + sum(q[j] for j in bus_list[i]["arcs"]) - sum(qg[k] for k in bus_list[i]["gen_idx"]) == 0)
            
                c10 = @constraint(model, [i = eachindex(bus_list)], vmins[i]^2 <= vr[i]^2+vim[i]^2 <= vmaxs[i]^2)
            elseif coords == "polar"
                c0 = @constraint(model,  [i = eachindex(refbus_list)], va[i] == 0)

                c1 = @constraint(model, [i = eachindex(branch_list)], p[branch_list[i]["f_idx"]] -(  (branch_list[i]["g"] + branch_list[i]["g_fr"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"]]^2) + 
                                (-branch_list[i]["g"]*branch_list[i]["tr"] + branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"]]*vm[branch_list[i]["t_bus"]]*cos(va[branch_list[i]["f_bus"]] - va[branch_list[i]["t_bus"]]))+
                                (-branch_list[i]["b"]*branch_list[i]["tr"]-branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"]]*vm[branch_list[i]["t_bus"]]*sin(va[branch_list[i]["f_bus"]] - va[branch_list[i]["t_bus"]]))) == 0)

                c2 = @constraint(model, [i = eachindex(branch_list)], q[branch_list[i]["f_idx"]] -( -(branch_list[i]["b"] + branch_list[i]["b_fr"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"]]^2) -
                                (-branch_list[i]["b"]*branch_list[i]["tr"] - branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"]]*vm[branch_list[i]["t_bus"]]*cos(va[branch_list[i]["f_bus"]] - va[branch_list[i]["t_bus"]]))+
                                (-branch_list[i]["g"]*branch_list[i]["tr"] + branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"]]*vm[branch_list[i]["t_bus"]]*sin(va[branch_list[i]["f_bus"]] - va[branch_list[i]["t_bus"]]))) == 0)

                c3 = @constraint(model, [i = eachindex(branch_list)], p[branch_list[i]["t_idx"]] -( (branch_list[i]["g"]+branch_list[i]["g_to"])*(vm[branch_list[i]["t_bus"]]^2) +
                                (-branch_list[i]["g"]*branch_list[i]["tr"]-branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"]]*vm[branch_list[i]["t_bus"]]*cos(va[branch_list[i]["t_bus"]] - va[branch_list[i]["f_bus"]])) + 
                                (-branch_list[i]["b"]*branch_list[i]["tr"]+branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"]]*vm[branch_list[i]["t_bus"]]*sin(va[branch_list[i]["t_bus"]] - va[branch_list[i]["f_bus"]]))) == 0)
            
                c4 = @constraint(model, [i = eachindex(branch_list)], q[branch_list[i]["t_idx"]] -( -(branch_list[i]["b"]+branch_list[i]["b_to"])*(vm[branch_list[i]["t_bus"]]^2) - 
                                (-branch_list[i]["b"]*branch_list[i]["tr"] + branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]* (vm[branch_list[i]["f_bus"]]*vm[branch_list[i]["t_bus"]]*cos(va[branch_list[i]["t_bus"]] - va[branch_list[i]["f_bus"]])) +
                                (-branch_list[i]["g"]*branch_list[i]["tr"] - branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"] *(vm[branch_list[i]["f_bus"]]*vm[branch_list[i]["t_bus"]]*sin(va[branch_list[i]["t_bus"]] - va[branch_list[i]["f_bus"]]))) == 0)
            
                c5 = @constraint(model, [i = eachindex(branch_list)], angmins[i] <= va[branch_list[i]["f_bus"]] - va[branch_list[i]["t_bus"]] <= angmaxs[i])

                c8 = @constraint(model, [i = eachindex(bus_list)], bus_list[i]["pd"] + bus_list[i]["gs"]*(vm[i]^2) + sum(p[j] for j in bus_list[i]["arcs"]) - sum(pg[k] for k in bus_list[i]["gen_idx"]) == 0)
        
                c9 = @constraint(model, [i = eachindex(bus_list)], bus_list[i]["qd"] - bus_list[i]["bs"]*(vm[i]^2)  + sum(q[j] for j in bus_list[i]["arcs"]) - sum(qg[k] for k in bus_list[i]["gen_idx"]) == 0)


            end
            c6 = @constraint(model, [i = eachindex(branch_list)], 0 >= p[branch_list[i]["f_idx"]]^2 + q[branch_list[i]["f_idx"]]^2 - branch_list[i]["a_rate_sq"])
        
            c7 = @constraint(model, [i = eachindex(branch_list)], 0 >= p[branch_list[i]["t_idx"]]^2 + q[branch_list[i]["t_idx"]]^2 - branch_list[i]["a_rate_sq"])
        
           
        
        
        end
        
        return model
    end
end

model = opf_combined("pglib_opf_case14_ieee.m", "polar", "static", "exa")
ipopt(model)
#optimize!(model)

