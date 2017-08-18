export
    post_ac_opf

"""
Given a JuMP model and a PowerModels network data structure, 
Builds an AC-OPF formulation of the given and returns the JuMP model
"""
function post_ac_opf(data::Dict{String,Any}, model=Model())
    ref = PMs.build_ref(data)

    @variable(model, t[i in keys(ref[:bus])])
    @variable(model, ref[:bus][i]["vmin"] <= v[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)

    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    @variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])


    pmin = Dict([(a, 0.0) for a in ref[:arcs_dc]])
    pmax = Dict([(a, 0.0) for a in ref[:arcs_dc]])
    qmin = Dict([(a, 0.0) for a in ref[:arcs_dc]])
    qmax = Dict([(a, 0.0) for a in ref[:arcs_dc]])
    for (l,i,j) in ref[:arcs_from_dc]
        qmin[(l,i,j)] =  ref[:dcline][l]["qminf"]
        qmax[(l,i,j)] =  ref[:dcline][l]["qmaxf"]
        qmin[(l,j,i)] =  ref[:dcline][l]["qmint"]
        qmax[(l,j,i)] =  ref[:dcline][l]["qmaxt"]
        pmin[(l,i,j)] =  ref[:dcline][l]["pminf"]
        pmax[(l,i,j)] =  ref[:dcline][l]["pmaxf"]
        pmin[(l,j,i)] =  ref[:dcline][l]["pmint"]
        pmax[(l,j,i)] =  ref[:dcline][l]["pmaxt"]
    end

    # pmin, qmin computations need to be moved to ref in power models
    @variable(model, pmin[(l,i,j)] <= p_dc[(l,i,j) in ref[:arcs_dc]] <= pmax[(l,i,j)])
    @variable(model, qmin[(l,i,j)] <= q_dc[(l,i,j) in ref[:arcs_dc]] <= qmax[(l,i,j)])
    
    
    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    @objective(model, Min, 
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
        sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    )

    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @constraint(model, t[i] == 0)
    end

    
    for (i,bus) in ref[:bus]
        # Bus KCL
        @constraint(model, sum(p[a] for a in ref[:bus_arcs][i]) + sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == sum(pg[g] for g in ref[:bus_gens][i]) - bus["pd"] - bus["gs"]*v[i]^2)
        @constraint(model, sum(q[a] for a in ref[:bus_arcs][i]) + sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == sum(qg[g] for g in ref[:bus_gens][i]) - bus["qd"] + bus["bs"]*v[i]^2)
    end

    for (i,branch) in ref[:branch]
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]

        t_fr = t[f_bus]
        t_to = t[t_bus]

        f_idx = (i, f_bus, t_bus)
        t_idx = (i, t_bus, f_bus)

        # Line Flow
        g, b = PMs.calc_branch_y(branch)
        tr, ti = PMs.calc_branch_t(branch)
        c = branch["br_b"]
        tm = branch["tap"]^2

        #constraint_ohms_yt_from(pm, branch)
        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        v_fr = v[f_bus]
        v_to = v[t_bus]
        t_fr = t[f_bus]
        t_to = t[t_bus]

        @NLconstraint(model, p_fr == g/tm*v_fr^2 + (-g*tr+b*ti)/tm*(v_fr*v_to*cos(t_fr-t_to)) + (-b*tr-g*ti)/tm*(v_fr*v_to*sin(t_fr-t_to)) )
        @NLconstraint(model, q_fr == -(b+c/2)/tm*v_fr^2 - (-b*tr-g*ti)/tm*(v_fr*v_to*cos(t_fr-t_to)) + (-g*tr+b*ti)/tm*(v_fr*v_to*sin(t_fr-t_to)) )

        #constraint_ohms_yt_to(pm, branch)
        @NLconstraint(model, p_to == g*v_to^2 + (-g*tr-b*ti)/tm*(v_to*v_fr*cos(t_to-t_fr)) + (-b*tr+g*ti)/tm*(v_to*v_fr*sin(t_to-t_fr)) )
        @NLconstraint(model, q_to == -(b+c/2)*v_to^2 - (-b*tr+g*ti)/tm*(v_to*v_fr*cos(t_fr-t_to)) + (-g*tr-b*ti)/tm*(v_to*v_fr*sin(t_to-t_fr)) )


        # Phase Angle Difference Limit
        #constraint_phase_angle_difference(pm, branch)
        @constraint(model, t_fr - t_to <= branch["angmax"])
        @constraint(model, t_fr - t_to >= branch["angmin"])

        # Apparent Power Limit, From and To
        @constraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
        @constraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_bus = dcline["f_bus"]
        t_bus = dcline["t_bus"]
        f_idx = (i, f_bus, t_bus)
        t_idx = (i, t_bus, f_bus)
        loss0 = dcline["loss0"]
        loss1 = dcline["loss1"]

        p_fr = p_dc[f_idx]
        p_to = p_dc[t_idx]

        @constraint(model, (1-loss1) * p_fr + (p_to - loss0) == 0)
    end

    return model
end



