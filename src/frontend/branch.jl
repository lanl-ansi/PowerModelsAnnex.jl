function net2pmc_branch(net::Network)
    branch_dict = Dict{String, Any}()
    for r in eachrow(line(net))
        old = Dict()
        if in("branch", keys(pmc(net))) && string(r[:element_id]) in keys(pmc(net)["branch"])
            old = pmc(net)["branch"][string(r[:element_id])]
        end
        there = !isempty(old)
        branch_dict[string(r[:element_id])] = Dict(
            "index" => (!there || !ismissing(r[:element_id])) ? r[:element_id] : old["index"],
            "br_r" => (!there || !ismissing(r[:element_id])) ? r[:resistance] : old["br_r"],
            "rate_a" => (!there || !ismissing(r[:element_id])) ? r[:rate_a] : old["rate_a"],
            "rate_b" => (!there || !ismissing(r[:element_id])) ? r[:rate_b] : old["rate_b"],
            "rate_c" => (!there || !ismissing(r[:element_id])) ? r[:rate_c] : old["rate_c"],
            "transformer" => (!there || !ismissing(r[:element_id])) ? r[:transformer] : old["transformer"],
            "f_bus" => (!there || !ismissing(r[:element_id])) ? r[:from_bus] : old["f_bus"],
            "t_bus" => (!there || !ismissing(r[:element_id])) ? r[:to_bus] : old["t_bus"],
            "br_status" => (!there || !ismissing(r[:element_id])) ? r[:status] : old["br_status"],
            "br_x" => (!there || !ismissing(r[:element_id])) ? r[:reactance] : old["br_x"],
            "br_b" => (!there || !ismissing(r[:element_id])) ? r[:susceptance] : old["br_b"],
            "angmin" => (!there || !ismissing(r[:element_id])) ? r[:ang_min] : old["angmin"],
            "angmax" => (!there || !ismissing(r[:element_id])) ? r[:ang_max] : old["angmax"],
            "b_to" => (!there || !ismissing(r[:element_id])) ? r[:b_to] : old["b_to"],
            "b_fr" => (!there || !ismissing(r[:element_id])) ? r[:b_fr] : old["b_fr"],
            "g_to" => (!there || !ismissing(r[:element_id])) ? r[:g_to] : old["g_to"],
            "g_fr" => (!there || !ismissing(r[:element_id])) ? r[:g_fr] : old["g_fr"],
            "shift" => 0.0,
            "tap" => 1.0,
        )
    end
    return branch_dict
end
