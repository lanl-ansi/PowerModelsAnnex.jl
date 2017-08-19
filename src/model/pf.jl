export post_ac_pf, post_dc_pf

"""
Given a JuMP model and a PowerModels network data structure, 
Builds an AC-PF formulation of the given data and returns the JuMP model
"""
function post_ac_pf(data::Dict{String,Any}, model=Model())
    ref = PMs.build_ref(data)
    epsilon = 0.00001

    @variable(model, t[i in keys(ref[:bus])])
    @variable(model, v[i in keys(ref[:bus])] >= 0, start=1.0)

    @variable(model, pg[i in keys(ref[:gen])])
    @variable(model, qg[i in keys(ref[:gen])])

    @variable(model, p[(l,i,j) in ref[:arcs]])
    @variable(model, q[(l,i,j) in ref[:arcs]])

    @variable(model, p_dc[(l,i,j) in ref[:arcs_dc]])
    @variable(model, q_dc[(l,i,j) in ref[:arcs_dc]])

    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @constraint(model, t[i] == 0)
        @constraint(model, v[i] == bus["vm"])
    end

    for (i,bus) in ref[:bus]
        # Bus KCL
        @constraint(model, 
            sum(p[a] for a in ref[:bus_arcs][i]) + 
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == 
            sum(pg[g] for g in ref[:bus_gens][i]) - bus["pd"] - bus["gs"]*v[i]^2
        )
        @constraint(model, 
            sum(q[a] for a in ref[:bus_arcs][i]) + 
            sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == 
            sum(qg[g] for g in ref[:bus_gens][i]) - bus["qd"] + bus["bs"]*v[i]^2
        )

        # PV Bus Constraints
        if length(ref[:bus_gens][i]) > 0 && !(i in keys(ref[:ref_buses]))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            # @constraint(model, v[i] == bus["vm"])
            # soft equality needed becouse vm in file may not be precice enough to ensure feasiblity
            @constraint(model, v[i] <= bus["vm"] + epsilon)
            @constraint(model, v[i] >= bus["vm"] - epsilon)

            for j in ref[:bus_gens][i]
                @constraint(model, pg[j] == ref[:gen][j]["pg"])
            end
        end
    end

    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        v_fr = v[branch["f_bus"]]
        v_to = v[branch["t_bus"]]
        t_fr = t[branch["f_bus"]]
        t_to = t[branch["t_bus"]]

        # Line Flow
        g, b = PMs.calc_branch_y(branch)
        tr, ti = PMs.calc_branch_t(branch)
        c = branch["br_b"]
        tm = branch["tap"]^2

        @NLconstraint(model, p_fr == g/tm*v_fr^2 + (-g*tr+b*ti)/tm*(v_fr*v_to*cos(t_fr-t_to)) + (-b*tr-g*ti)/tm*(v_fr*v_to*sin(t_fr-t_to)) )
        @NLconstraint(model, q_fr == -(b+c/2)/tm*v_fr^2 - (-b*tr-g*ti)/tm*(v_fr*v_to*cos(t_fr-t_to)) + (-g*tr+b*ti)/tm*(v_fr*v_to*sin(t_fr-t_to)) )

        @NLconstraint(model, p_to == g*v_to^2 + (-g*tr-b*ti)/tm*(v_to*v_fr*cos(t_to-t_fr)) + (-b*tr+g*ti)/tm*(v_to*v_fr*sin(t_to-t_fr)) )
        @NLconstraint(model, q_to == -(b+c/2)*v_to^2 - (-b*tr+g*ti)/tm*(v_to*v_fr*cos(t_fr-t_to)) + (-g*tr-b*ti)/tm*(v_to*v_fr*sin(t_to-t_fr)) )
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        @constraint(model, p_dc[f_idx] == dcline["pf"])
        @constraint(model, p_dc[t_idx] == dcline["pt"])

        # @constraint(model, v[dcline["f_bus"]] == dcline["vf"])
        # soft equality needed becouse vm in file may not be precice enough to ensure feasiblity
        @constraint(model, v[dcline["f_bus"]] <= dcline["vf"] + epsilon)
        @constraint(model, v[dcline["f_bus"]] >= dcline["vf"] - epsilon)

        # @constraint(model, v[dcline["t_bus"]] == dcline["vt"])
        # soft equality needed becouse vm in file may not be precice enough to ensure feasiblity
        @constraint(model, v[dcline["t_bus"]] <= dcline["vt"] + epsilon)
        @constraint(model, v[dcline["t_bus"]] >= dcline["vt"] - epsilon)
    end

    return model
end


"""
Given a JuMP model and a PowerModels network data structure, 
Builds an DC-PF formulation of the given data and returns the JuMP model
"""
function post_dc_pf(data::Dict{String,Any}, model=Model())
    ref = PMs.build_ref(data)

    @variable(model, t[i in keys(ref[:bus])])

    @variable(model, pg[i in keys(ref[:gen])])

    @variable(model, p[(l,i,j) in ref[:arcs_from]])

    p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
    p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))

    @variable(model, p_dc[(l,i,j) in ref[:arcs_dc]])

    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @constraint(model, t[i] == 0)
    end

    for (i,bus) in ref[:bus]
        # Bus KCL
        @constraint(model, 
            sum(p_expr[a] for a in ref[:bus_arcs][i]) + 
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == 
            sum(pg[g] for g in ref[:bus_gens][i]) - bus["pd"] - bus["gs"]*1.0^2
        )

        # PV Bus Constraints
        if length(ref[:bus_gens][i]) > 0 && !(i in keys(ref[:ref_buses]))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            for j in ref[:bus_gens][i]
                @constraint(model, pg[j] == ref[:gen][j]["pg"])
            end
        end
    end

    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        t_fr = t[branch["f_bus"]]
        t_to = t[branch["t_bus"]]

        # Line Flow
        g, b = PMs.calc_branch_y(branch)

        @constraint(model, p_fr == -b*(t_fr - t_to))
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        @constraint(model, p_dc[f_idx] == dcline["pf"])
        @constraint(model, p_dc[t_idx] == dcline["pt"])
    end

    return model
end
