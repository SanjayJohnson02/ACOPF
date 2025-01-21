#=
using Pkg
Pkg.add("PowerModels")
Pkg.add("Interpolations")
Pkg.add("JuMP")
Pkg.add("Ipopt")
Pkg.add(["Juniper", "NLopt"])
using Juniper

using Interpolations
using PowerModels
using JuMP
using Ipopt 
using LinearAlgebra 
=#

#build nested dictionary
function opf_jump_polar_storage(filename)
    my_data = PowerModels.parse_file(filename)
    n_bus = length(my_data["bus"])
    n_storage = length(my_data["storage"])

    base_MVA = my_data["baseMVA"]

    interval_split = my_data["time_elapsed"]
    #interval_split = 1
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
    gen_names = ["i", "cost1", "cost2", "cost3", "bus", "pstart", "qstart"]

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
        Dict("i" => 0, "Einit" => 0, "Emax" => 0.0, "Pcmax" => 0, "Pdmax" => 0, "etac" => 0, "etad" => 0, "Zr" => 0, "Zim" => 0, "Srating" => 0, "Pexts" => 0, "Qexts" => 0, "bus" => 0) for _ in 1:n_storage
    ]
    for key in eachindex(my_data["storage"])
        i = parse(Int, key)
        storage_list[i]["i"] = i
        storage_list[i]["Einit"] = my_data["storage"][key]["energy"]
        storage_list[i]["Emax"] = my_data["storage"][key]["energy_rating"]
        storage_list[i]["Pcmax"] = my_data["storage"][key]["charge_rating"]

        #this is slightly different from what paper states
        storage_list[i]["Pdmax"] = my_data["storage"][key]["discharge_rating"]
        #storage_list[i]["Pdmax"] = 0.75

        #these are swapped form data and paper?
        storage_list[i]["etac"] = my_data["storage"][key]["charge_efficiency"]
        #storage_list[i]["etac"] = .85
        storage_list[i]["etad"] = my_data["storage"][key]["discharge_efficiency"]
        #storage_list[i]["etad"] = 0.9

        storage_list[i]["Zr"] = my_data["storage"][key]["r"]
        storage_list[i]["Zim"] = my_data["storage"][key]["x"]
        storage_list[i]["Srating"] = my_data["storage"][key]["thermal_rating"]
        storage_list[i]["Pexts"] = my_data["storage"][key]["p_loss"]
        storage_list[i]["Qexts"] = my_data["storage"][key]["q_loss"]
        storage_list[i]["bus"] = my_data["storage"][key]["storage_bus"]

        push!(bus_list[my_data["storage"][key]["storage_bus"]]["storage_idx"], i)

    end




    #optimizer = Juniper.Optimizer
    #nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

    # Create a JuMP model using Juniper
    #model = Model(optimizer_with_attributes(optimizer, "nl_solver"=>nl_solver))

    model = Model(Ipopt.Optimizer)

    set_optimizer_attribute(model, "max_iter", 500)

    #not sure if these should be initialized at nonzero
    va = @variable(model, va[i = 1:length(bus_list), t = 1:length(summer_wkdy_qrtr_scalar)])
    vm = @variable(model, vmins[i] <= vm[i=1:length(bus_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= vmaxs[i], start = 1)

    pg = @variable(model, pmins[i] <= pg[i = 1:length(gen_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= pmaxs[i], start = gen_list[i]["pstart"])

    qg = @variable(model, qmins[i] <= qg[i =1:length(gen_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= qmaxs[i], start = gen_list[i]["qstart"])

    p = @variable(model, -rate_as[i] <= p[i =1:length(arc_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= rate_as[i])

    q = @variable(model, -rate_as[i] <= q[i=1:length(arc_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= rate_as[i])

    pstc = @variable(model, 0 <= pstc[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= storage_list[i]["Pcmax"])

    pstd = @variable(model, 0 <= pstd[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= storage_list[i]["Pdmax"])

    pst = @variable(model, -storage_list[i]["Srating"] <= pst[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= storage_list[i]["Srating"])

    #this is maybe supposed to have bounds??
    qst = @variable(model, -storage_list[i]["Srating"] <= qst[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= storage_list[i]["Srating"])

    #equivalent to I^2
    I2 = @variable(model, 0 <= I2[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)])

    qint = @variable(model, -storage_list[i]["Srating"] <= qint[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= storage_list[i]["Srating"])

    E = @variable(model, 0 <= E[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)] <= storage_list[i]["Emax"])

    #only for use in MINLP version
    #zc = @variable(model, zc[i = 1:length(storage_list), t = 1:length(summer_wkdy_qrtr_scalar)], Bin)


    obj = @objective(model, Min, sum(pg[i,t]^2*gen_list[i]["cost1"]+pg[i,t]*gen_list[i]["cost2"]+gen_list[i]["cost3"] for i in eachindex(gen_list), t in eachindex(summer_wkdy_qrtr_scalar)))

    #this is hard coded, in this model it doesn't matter but in the future it may?
    c0 = @constraint(model,  [t = eachindex(summer_wkdy_qrtr_scalar)], va[4, t] == 0)



    c1 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], p[branch_list[i]["f_idx"], t] -(  (branch_list[i]["g"] + branch_list[i]["g_fr"])/branch_list[i]["ttm"]*vm[branch_list[i]["f_bus"], t]^2 + 
                    (-branch_list[i]["g"]*branch_list[i]["tr"] + branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*cos(va[branch_list[i]["f_bus"], t] - va[branch_list[i]["t_bus"], t]))+
                    (-branch_list[i]["b"]*branch_list[i]["tr"]-branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*sin(va[branch_list[i]["f_bus"], t] - va[branch_list[i]["t_bus"], t]))) == 0)


    c2 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], q[branch_list[i]["f_idx"], t] -( -(branch_list[i]["b"] + branch_list[i]["b_fr"])/branch_list[i]["ttm"]*vm[branch_list[i]["f_bus"], t]^2 -
                    (-branch_list[i]["b"]*branch_list[i]["tr"] - branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*cos(va[branch_list[i]["f_bus"], t] - va[branch_list[i]["t_bus"], t]))+
                    (-branch_list[i]["g"]*branch_list[i]["tr"] + branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*sin(va[branch_list[i]["f_bus"], t] - va[branch_list[i]["t_bus"], t]))) == 0)


    c3 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], p[branch_list[i]["t_idx"], t] -( (branch_list[i]["g"]+branch_list[i]["g_to"])*vm[branch_list[i]["t_bus"], t]^2 +
                    (-branch_list[i]["g"]*branch_list[i]["tr"]-branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*cos(va[branch_list[i]["t_bus"], t] - va[branch_list[i]["f_bus"], t])) + 
                    (-branch_list[i]["b"]*branch_list[i]["tr"]+branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]*(vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*sin(va[branch_list[i]["t_bus"], t] - va[branch_list[i]["f_bus"], t]))) == 0)


    c4 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], q[branch_list[i]["t_idx"], t] -( -(branch_list[i]["b"]+branch_list[i]["b_to"])*vm[branch_list[i]["t_bus"], t]^2 - 
                    (-branch_list[i]["b"]*branch_list[i]["tr"] + branch_list[i]["g"]*branch_list[i]["ti"])/branch_list[i]["ttm"]* (vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*cos(va[branch_list[i]["t_bus"], t] - va[branch_list[i]["f_bus"], t])) +
                    (-branch_list[i]["g"]*branch_list[i]["tr"] - branch_list[i]["b"]*branch_list[i]["ti"])/branch_list[i]["ttm"] * (vm[branch_list[i]["f_bus"], t]*vm[branch_list[i]["t_bus"], t]*sin(va[branch_list[i]["t_bus"], t] - va[branch_list[i]["f_bus"], t]))) == 0)

    c5 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], angmins[i] <= va[branch_list[i]["f_bus"], t] -va[branch_list[i]["t_bus"], t]  <= angmaxs[i])


    c6 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], p[branch_list[i]["f_idx"], t]^2 + q[branch_list[i]["f_idx"], t]^2 - branch_list[i]["a_rate_sq"] <= 0)

    c7 = @constraint(model, [i = eachindex(branch_list), t = eachindex(summer_wkdy_qrtr_scalar)], p[branch_list[i]["t_idx"], t]^2 + q[branch_list[i]["t_idx"], t]^2 - branch_list[i]["a_rate_sq"] <= 0)

    c8 = @constraint(model, [i = eachindex(bus_list), t = eachindex(summer_wkdy_qrtr_scalar)], bus_list[i]["pd"][t] + bus_list[i]["gs"]*(vm[i, t]^2) + sum(p[j, t] for j in bus_list[i]["arcs"]) - 
                    sum(pg[k, t] for k in bus_list[i]["gen_idx"]) + sum(pst[l, t] for l in bus_list[i]["storage_idx"]) == 0)

    c9 = @constraint(model, [i = eachindex(bus_list), t = eachindex(summer_wkdy_qrtr_scalar)], bus_list[i]["qd"][t] - bus_list[i]["bs"]*(vm[i, t]^2)  + sum(q[j, t] for j in bus_list[i]["arcs"]) - 
                    sum(qg[k, t] for k in bus_list[i]["gen_idx"]) + sum(qst[l, t] for l in bus_list[i]["storage_idx"])== 0)

                    
    c11 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pst[c, t] + pstd[c, t] - pstc[c, t] == storage_list[c]["Pexts"] + storage_list[c]["Zr"]*I2[c, t])

    c12 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], qst[c, t] == qint[c, t] + storage_list[c]["Qexts"] + storage_list[c]["Zim"]*I2[c, t])

    c13 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pst[c, t]^2 + qst[c, t]^2 == vm[storage_list[c]["bus"], t]^2*I2[c, t])

    c14 = @constraint(model, [c = eachindex(storage_list), t = 2:length(summer_wkdy_qrtr_scalar)], E[c, t] - E[c, t-1] == interval_split*(storage_list[c]["etac"]*pstc[c, t]- pstd[c, t]/storage_list[c]["etad"]))

    c15 = @constraint(model, [c = eachindex(storage_list)], E[c, 1] - storage_list[c]["Einit"] == interval_split*(storage_list[c]["etac"]*pstc[c, 1]- pstd[c, 1]/storage_list[c]["etad"]))

    c16 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pst[c, t]^2 + qst[c, t]^2 <= storage_list[c]["Srating"]^2)

    c17 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], -storage_list[c]["Srating"] <= pstd[c, t] - pstc[c, t] <= storage_list[c]["Srating"])

    #adding variable Pstor
    #c17b = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pstor[c, t] == pstd[c, t] - pstc[c, t])

    c18 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pstc[c, t]*pstd[c, t] == 0)

    #for use in MINLP
    #c18 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pstc[c, t] <= storage_list[c]["Pcmax"]*zc[c, t])
    #c18b = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pstd[c, t] <= storage_list[c]["Pdmax"]*(1 - zc[c, t]))


    #=
    c19 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pstc[c, t] == 0)
    c20 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pstd[c, t] == 0)
    c21 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], pst[c, t] == 0)
    c22 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], qst[c, t] == 0)
    c23 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], qint[c, t] == 0)
    c24 = @constraint(model, [c = eachindex(storage_list), t = eachindex(summer_wkdy_qrtr_scalar)], I[c, t] == 0)
    =#
    return model
end
#optimize!(model)

