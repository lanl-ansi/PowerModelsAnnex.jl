include("branch.jl")
include("bus.jl")
include("generators.jl")
include("load.jl")

"""
    network2pmc(net::Network)

Return a PowerModels network model (dictionary) based on the information contained within
a Network `net`.
"""
function network2pmc(
    net::Network;
    per_unit::Bool=true,
    name::AbstractString="pmc",
    baseMVA::Int=100,
    multinetwork::Bool=false,
    version::Int=2,
)
    outnet = Dict(
        "per_unit" => per_unit, # We may want to check if it is ok to use a deafult here.
        "name" => name,
        "baseMVA" => baseMVA, # We likely don't want a default value here.
        "multinetwork" => multinetwork,
        "version" => version,
        "dcline" => Dict{String, Any}() # We don't use this.
    )
    stripunits!(net)

    gen_dict = net2pmc_gen(net; baseMVA=baseMVA)
    ps_load_dict = net2pmc_ps_load(net; baseMVA=baseMVA)
    outnet["gen"] = merge(gen_dict, ps_load_dict)
    outnet["branch"] = net2pmc_branch(net)
    outnet["bus"] = net2pmc_bus(net)
    outnet["load"] = net2pmc_load(net)

    # Add extra keys that `Network` does not have. These are simply placeholders
    outnet["shunt"] = Dict{String,Any}()
    outnet["storage"] = Dict{String,Any}()
    return outnet
end

"""
    build_pmc!(net::Network)

Read information contained in `net` and uses it to populate `net.pmc`. This will overwrite
any information previously contained in there.
"""
function build_pmc!(net::Network)
    updated = network2pmc(net)
    for k in keys(updated)
        net.pmc[k] = updated[k]
    end
    #PowerModels.check_network_data(updated) # Small problem, it's the cause of
    # extreme verbosity.
end
