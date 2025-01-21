#=
using Pkg
Pkg.add("PowerModels")

using PowerModels
using JuMP
using Ipopt=#

#build nested dictionary
function opf_jump_rect(filename)

    my_data = PowerModels.parse_file(filename)

    n_bus = length(my_data["bus"])

    base_MVA = my_data["baseMVA"]


    # Create a list of dictionaries for each bus
    bus_list = [
        Dict("i" => 0, "pd" => 0.0, "gs" => 0.0, "qd" => 0, "bs" => 0, "bus_type" => 0, "arcs" => [], "gen_idx" => []) for _ in 1:n_bus
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
        branch_dict["a_rate_sq"] = my_data["branch"][key]["rate_a"]^2
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
        bus_list[my_data["load"][key]["load_bus"]]["pd"] = my_data["load"][key]["pd"]
        bus_list[my_data["load"][key]["load_bus"]]["qd"] = my_data["load"][key]["qd"]
    end



    model = Model(Ipopt.Optimizer)

    #not sure if these should be initialized at nonzero
    vr = @variable(model, vr[i = 1:length(bus_list)], start = 1)
    vim = @variable(model, vim[i=1:length(bus_list)])

    pg = @variable(model, pmins[i] <= pg[i = 1:length(gen_list)] <= pmaxs[i])

    qg = @variable(model, qmins[i] <= qg[i =1:length(gen_list)] <= qmaxs[i])

    p = @variable(model, -rate_as[i] <= p[i =1:length(arc_list)] <= rate_as[i])

    q = @variable(model, -rate_as[i] <= q[i=1:length(arc_list)] <= rate_as[i])

    obj = @objective(model, Min, sum(pg[i]^2*gen_list[i]["cost1"]+pg[i]*gen_list[i]["cost2"]+gen_list[i]["cost3"] for i in eachindex(pg)))


    #this is hard coded, in this model it doesn't matter but in the future it may?
    c0 = @constraint(model,  atan(vim[4]/vr[4]) == 0)



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

    c6 = @constraint(model, [i = eachindex(branch_list)], 0 >= p[branch_list[i]["f_idx"]]^2 + q[branch_list[i]["f_idx"]]^2 - branch_list[i]["a_rate_sq"])

    c7 = @constraint(model, [i = eachindex(branch_list)], 0 >= p[branch_list[i]["t_idx"]]^2 + q[branch_list[i]["t_idx"]]^2 - branch_list[i]["a_rate_sq"])

    c8 = @constraint(model, [i = eachindex(bus_list)], bus_list[i]["pd"] + bus_list[i]["gs"]*(vr[i]^2+vim[i]^2) + sum(p[j] for j in bus_list[i]["arcs"]) - sum(pg[k] for k in bus_list[i]["gen_idx"]) == 0)

    c9 = @constraint(model, [i = eachindex(bus_list)], bus_list[i]["qd"] - bus_list[i]["bs"]*(vr[i]^2+vim[i]^2)  + sum(q[j] for j in bus_list[i]["arcs"]) - sum(qg[k] for k in bus_list[i]["gen_idx"]) == 0)

    c10 = @constraint(model, [i = eachindex(bus_list)], vmins[i]^2 <= vr[i]^2+vim[i]^2 <= vmaxs[i]^2)


    return model
end
