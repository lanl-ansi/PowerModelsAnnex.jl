function net2pmc_gen(net::Network; baseMVA=100)
    gen_dict = Dict{String, Any}()
    to_warn = Int[]
    for r in eachrow(gen(net))
        old = Dict()
        if in("gen", keys(pmc(net))) && string(r[:element_id]) in keys(pmc(net)["gen"])
            old = pmc(net)["gen"][string(r[:element_id])]
        end
        coeffs = [0.0]
        there = !isempty(old)
        if (!there && ismissing(r[:cost]))
            push!(to_warn, r[:element_id])
        end
        gen_dict[string(r[:element_id])] = Dict(
            "index" => (!there || !ismissing(r[:element_id])) ? r[:element_id] : old["index"],
            "gen_bus" => (!there || !ismissing(r[:bus])) ? r[:bus] : old["gen_bus"],
            "pg" => (!there || !ismissing(r[:gen_p])) ? r[:gen_p] : old["pg"],
            "pmax" => (!there || !ismissing(r[:p_max])) ? r[:p_max] : old["pmax"],
            "pmin" => (!there || !ismissing(r[:p_min])) ? r[:p_min] : old["pmin"],
            "qg" => (!there || !ismissing(r[:gen_q])) ? r[:gen_q] : old["qg"],
            "qmax" => (!there || !ismissing(r[:q_max])) ? r[:q_max] : old["qmax"],
            "qmin" => (!there || !ismissing(r[:q_min])) ? r[:q_min] : old["qmin"],
            "startup" => (!there || !ismissing(r[:startup_cost])) ? r[:startup_cost] : old["startup"],
            "gen_status" => (!there || !ismissing(r[:status])) ? r[:status] : old["gen_status"],
            "ramp_10" => (!there || !ismissing(r[:ramp])) ? r[:ramp] : old["ramp_10"],
            "model" => !there  ? -1 : old["model"], # Default or old value
            "ncost" => !there  ? 0 : old["ncost"], # Default or old value
            "cost" =>  !there ? coeffs : old["cost"], # Default or old value
            "qc1max" => 0.0,
            "qc2max" => 0.0,
            "mbase" => baseMVA,
            "pc2" => 0.0,
            "pc1" => 0.0,
            "ramp_q" => 0.0,
            "ramp_30" => 0.0,
            "apf" => 0.0,
            "shutdown" => 0.0,
            "ramp_agc" => 0.0,
            "vg" => 1.01,
            "qc1min" => 0.0,
            "qc2min" => 0.0,
        )
        if !ismissing(r[:cost]) && !there
            if  isa(cost_gen(net)[r[:cost], :coeffs], PolynomialCost)
                gen_dict[string(r[:element_id])]["cost"] = costcurve2pmc(cost_gen(net)[r[:cost], :coeffs])
                gen_dict[string(r[:element_id])]["model"] = 2
                gen_dict[string(r[:element_id])]["ncost"] = n_cost(cost_gen(net)[r[:cost], :coeffs])
            elseif isa(cost_gen(net)[r[:cost], :coeffs], PWLCost)
                gen_dict[string(r[:element_id])]["cost"] = costcurve2pmc(cost_gen(net)[r[:cost], :coeffs])
                gen_dict[string(r[:element_id])]["model"] = 1
                gen_dict[string(r[:element_id])]["ncost"] = Int(length(gen_dict[string(r[:element_id])]["cost"])/2)
            end
        end
    end
    if !isempty(to_warn)
        warn(LOGGER, "There is no cost data for generator $to_warn). Default cost will be assigned.")
    end
    return gen_dict
end

function net2pmc_ps_load(net::Network; baseMVA=100)
    ps_load_dict = Dict{String, Any}()
    # Here we also need to add the price sensitive loads as generators
    id_offset = maximum(gen(net)[!, :element_id])
    to_warn = Int[]
    for r in eachrow(ps_load(net))
        old = Dict()
        if string(r[:element_id] + id_offset) in keys(pmc(net)["gen"])
            old = pmc(net)["gen"][string(r[:element_id])]
        end
        coeffs = [0.0]
        there = !isempty(old)
        if (!there && ismissing(r[:cost]))
            push!(to_warn, r[:element_id])
        end
        ps_load_dict[string(r[:element_id] + id_offset)] = Dict(
            "index" => (!there || !ismissing(r[:element_id])) ? r[:element_id] + id_offset : old["index"],
            "gen_bus" => (!there || !ismissing(r[:bus])) ? r[:bus] : old["gen_bus"],
            "pmax" => (!there || !ismissing(r[:load_max])) ? r[:load_max] : old["pmax"],
            "cost" => !there ? coeffs : old["cost"], # Default or old value
            "model" => !there  ? -1 : old["model"], # Default or old value
            "ncost" => !there  ? 0 : old["ncost"], # Default or old value
            "gen_status" => (!there || !ismissing(r[:status])) ? r[:status] : old["gen_status"],
            "ramp_10" => 0.0,
            "startup" => 0.0,
            "qg" => 0.0,
            "qmax" => 0.0,
            "qmin" => 0.0,
            "pmin" => 0.0,
            "qc1max" => 0.0,
            "model" => 2,
            "qc2max" => 0.0,
            "mbase" => baseMVA,
            "pc2" => 0.0,
            "pc1" => 0.0,
            "ramp_q" => 0.0,
            "ramp_30" => 0.0,
            "apf" => 0.0,
            "shutdown" => 0.0,
            "ramp_agc" => 0.0,
            "vg" => 1.01,
            "qc1min" => 0.0,
            "qc2min" => 0.0,
        )
        # TODO: Check the signs of ps load costs
        if !ismissing(r[:cost]) && !there
            if isa(cost_load(net)[r[:cost], :coeffs], PolynomialCost)
                # Load is negative generation. By serving load, the total cost receives a negative contribution
                # due to the amount payed by the utility.
                # The cost curve has to be flipped along the x axis.
                ps_load_dict[string(r[:element_id] + id_offset)]["cost"] = -costcurve2pmc(cost_load(net)[r[:cost], :coeffs])
                aux = ps_load_dict[string(r[:element_id] + id_offset)]["pmin"]
                ps_load_dict[string(r[:element_id] + id_offset)]["pmin"] = -ps_load_dict[string(r[:element_id] + id_offset)]["pmax"] # The cost curve has to be flipped along the y axis.
                ps_load_dict[string(r[:element_id] + id_offset)]["pmax"] = -aux
                ps_load_dict[string(r[:element_id] + id_offset)]["model"] = 2
                ps_load_dict[string(r[:element_id] + id_offset)]["ncost"] = n_cost(cost_load(net)[r[:cost], :coeffs])
            elseif isa(cost_load(net)[r[:cost], :coeffs], PWLCost)
                ps_load_dict[string(r[:element_id] + id_offset)]["cost"] = -costcurve2pmc(cost_load(net)[r[:cost], :coeffs])
                aux = ps_load_dict[string(r[:element_id] + id_offset)]["pmin"]
                ps_load_dict[string(r[:element_id] + id_offset)]["pmin"] = -ps_load_dict[string(r[:element_id] + id_offset)]["pmax"] # The cost curve has to be flipped along the y axis.
                ps_load_dict[string(r[:element_id] + id_offset)]["pmax"] = -aux
                ps_load_dict[string(r[:element_id] + id_offset)]["model"] = 1
                ps_load_dict[string(r[:element_id] + id_offset)]["ncost"] = Int(length(ps_load_dict[string(r[:element_id] + id_offset)]["cost"])/2)
            end
        end
    end
    if !isempty(to_warn)
        warn(LOGGER, "There is no cost data for ps_load $to_warn). Default cost will be assigned.")
    end
    return ps_load_dict
end
