export post_ac_pf, post_soc_pf, post_dc_pf

"""
Given a JuMP model and a PowerModels network data structure, 
Builds an AC-PF formulation of the given data and returns the JuMP model
"""
function post_ac_pf(data::Dict{String,Any}, model=Model())
    @assert !InfrastructureModels.ismultinetwork(data)
    ref = PMs.build_ref(data)[:nw][0]

    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, vm[i in keys(ref[:bus])] >= 0, start=1.0)

    @variable(model, pg[i in keys(ref[:gen])])
    @variable(model, qg[i in keys(ref[:gen])])

    @variable(model, p[(l,i,j) in ref[:arcs]])
    @variable(model, q[(l,i,j) in ref[:arcs]])

    @variable(model, p_dc[(l,i,j) in ref[:arcs_dc]])
    @variable(model, q_dc[(l,i,j) in ref[:arcs_dc]])

    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @assert bus["bus_type"] == 3
        @constraint(model, va[i] == 0)
        @constraint(model, vm[i] == bus["vm"])
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Bus KCL
        @constraint(model, 
            sum(p[a] for a in ref[:bus_arcs][i]) + 
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == 
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2
        )
        @constraint(model, 
            sum(q[a] for a in ref[:bus_arcs][i]) + 
            sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == 
            sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2
        )

        # PV Bus Constraints
        if length(ref[:bus_gens][i]) > 0 && !(i in keys(ref[:ref_buses]))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            @constraint(model, vm[i] == bus["vm"])

            for j in ref[:bus_gens][i]
                @constraint(model, pg[j] == ref[:gen][j]["pg"])
            end
        end
    end

    for (i,branch) in ref[:branch]
        # AC Line Flow Constraint
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

        # Line Flow
        g, b = PMs.calc_branch_y(branch)
        tr, ti = PMs.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2

        @NLconstraint(model, p_fr ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        @NLconstraint(model, q_fr == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        @NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        @NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        @constraint(model, p_dc[f_idx] == dcline["pf"])
        @constraint(model, p_dc[t_idx] == dcline["pt"])

        f_bus = ref[:bus][dcline["f_bus"]]
        if f_bus["bus_type"] == 1
            @constraint(model, vm[dcline["f_bus"]] == f_bus["vm"])
        end

        t_bus = ref[:bus][dcline["t_bus"]]
        if t_bus["bus_type"] == 1
            @constraint(model, vm[dcline["t_bus"]] == t_bus["vm"])
        end
    end

    return model
end


"""
Given a JuMP model and a PowerModels network data structure, 
Builds an SOC-PF formulation of the given data and returns the JuMP model
"""
function post_soc_pf(data::Dict{String,Any}, model=Model())
    @assert !InfrastructureModels.ismultinetwork(data)
    ref = PMs.build_ref(data)[:nw][0]

    @variable(model, w[i in keys(ref[:bus])] >= 0, start=1.001)
    @variable(model, wr[bp in keys(ref[:buspairs])], start=1.0)
    @variable(model, wi[bp in keys(ref[:buspairs])])

    @variable(model, pg[i in keys(ref[:gen])])
    @variable(model, qg[i in keys(ref[:gen])])

    @variable(model, p[(l,i,j) in ref[:arcs]])
    @variable(model, q[(l,i,j) in ref[:arcs]])

    @variable(model, p_dc[(l,i,j) in ref[:arcs_dc]])
    @variable(model, q_dc[(l,i,j) in ref[:arcs_dc]])

    for (i,j) in keys(ref[:buspairs])
        # Voltage Product Relaxation
        @constraint(model, wr[(i,j)]^2 + wi[(i,j)]^2 <= w[i]*w[j])
    end

    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @assert bus["bus_type"] == 3
        @constraint(model, w[i] == bus["vm"]^2)
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Bus KCL
        @constraint(model, 
            sum(p[a] for a in ref[:bus_arcs][i]) + 
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == 
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*w[i]
        )
        @constraint(model, 
            sum(q[a] for a in ref[:bus_arcs][i]) + 
            sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == 
            sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts)*w[i]
        )

        # PV Bus Constraints
        if length(ref[:bus_gens][i]) > 0 && !(i in keys(ref[:ref_buses]))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            @constraint(model, w[i] == bus["vm"]^2)
            for j in ref[:bus_gens][i]
                @constraint(model, pg[j] == ref[:gen][j]["pg"])
            end
        end
    end

    for (i,branch) in ref[:branch]
        # AC Line Flow Constraint
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        bp_idx = (branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        w_fr = w[branch["f_bus"]]
        w_to = w[branch["t_bus"]]
        wr_br = wr[bp_idx]
        wi_br = wi[bp_idx]

        g, b = PMs.calc_branch_y(branch)
        tr, ti = PMs.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2

        @constraint(model, p_fr ==  (g+g_fr)/tm*w_fr + (-g*tr+b*ti)/tm*(wr_br) + (-b*tr-g*ti)/tm*(wi_br) )
        @constraint(model, q_fr == -(b+b_fr)/tm*w_fr - (-b*tr-g*ti)/tm*(wr_br) + (-g*tr+b*ti)/tm*(wi_br) )

        @constraint(model, p_to ==  (g+g_to)*w_to + (-g*tr-b*ti)/tm*(wr_br) + (-b*tr+g*ti)/tm*(-wi_br) )
        @constraint(model, q_to == -(b+b_to)*w_to - (-b*tr+g*ti)/tm*(wr_br) + (-g*tr-b*ti)/tm*(-wi_br) )
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        @constraint(model, p_dc[f_idx] == dcline["pf"])
        @constraint(model, p_dc[t_idx] == dcline["pt"])

        f_bus = ref[:bus][dcline["f_bus"]]
        if f_bus["bus_type"] == 1
            @constraint(model, w[dcline["f_bus"]] == f_bus["vm"]^2)
        end

        t_bus = ref[:bus][dcline["t_bus"]]
        if t_bus["bus_type"] == 1
            @constraint(model, w[dcline["t_bus"]] == t_bus["vm"]^2)
        end
    end

    return model
end


"""
Given a JuMP model and a PowerModels network data structure, 
Builds an DC-PF formulation of the given data and returns the JuMP model
"""
function post_dc_pf(data::Dict{String,Any}, model=Model())
    @assert !InfrastructureModels.ismultinetwork(data)
    ref = PMs.build_ref(data)[:nw][0]

    @variable(model, va[i in keys(ref[:bus])])

    @variable(model, pg[i in keys(ref[:gen])])

    @variable(model, p[(l,i,j) in ref[:arcs_from]])

    p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
    p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))

    @variable(model, p_dc[(l,i,j) in ref[:arcs_dc]])

    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @constraint(model, va[i] == 0)
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Bus KCL
        @constraint(model, 
            sum(p_expr[a] for a in ref[:bus_arcs][i]) + 
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) == 
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*1.0^2
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
        # AC Line Flow Constraint
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        # Line Flow
        g, b = PMs.calc_branch_y(branch)

        @constraint(model, p_fr == -b*(va_fr - va_to))
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
