using DataFrames
using Memento
using Missings
using Unitful
include("units.jl")
using .Units

export
    Network,
    add_bus!,
    add_gen!,
    add_load!,
    add_pi_load!,
    add_ps_load!,
    add_line!,
    add_cost_gen!,
    add_cost_load!,
    network2pmc,
    build_pmc!,
    converged,
    PolynomialCost,
    PWLCost,
    CostCurve

# These dictionaries store which quantities are relevant to us when building a network model.
# The keys are the titles of the columns we are going to use, while the values are the corresponding
# keys in a PowerModels network.
# Whenever we don't want to import a column from the PowerModels network, the dictionary establishes
# a default value to be used. Default value is `missing` for most quantites, except for status,
# for which we use `1`, meaning the element is on, and for `element_id`, for which we use -1

"""
    bus_columns

:element_id -> id number of the bus
:name ->  name of the bus
:voltage -> bus operating voltage
:volt_max -> maximum bus voltage
:volt_min -> minimum bus voltage
:bus_type -> Bus type
:base_kv -> Base voltage
:coords -> Bus geocoordinates in latitude and longitude
:zone -> Bus zone
:area -> Bus Area
"""
const bus_columns = Dict(
    :element_id => :index,
    :name => :bus_name,
    :voltage => :vm,
    :volt_max => :vmax,
    :volt_min => :vmin,
    :bus_type => :bus_type,
    :base_kv => :base_kv,
    :coords => missing, # The cases don't come with location data.
    :zone => :zone,
    :area => :area,
)

"""
    gen_columns

:element_id -> generator ID number
:bus -> ID of the bus where the generator is placed
:gen_p -> active power generated
:p_max -> maximum active power
:p_min -> minimum active power
:gen_q -> reactive power generated
:q_max -> maximum reactive power
:q_min -> minimum reactive power
:cost -> ID number of the associated cost
:startup_cost -> startup cost
:status -> 0=IDLE, 1=ACTIVE
:ramp -> 10-minute ramping constraint
"""
const gen_columns = Dict(
    :element_id => :index,
    :bus => :gen_bus,
    :gen_p => :pg, # Active power
    :p_max => :pmax,
    :p_min => :pmin,
    :gen_q => :qg, # Reactive power
    :q_max => :qmax,
    :q_min => :qmin,
    :cost => missing, # We don't want to store the cost function here, just its id.
    :startup_cost => :startup,
    :status => :gen_status,
    :ramp => :ramp_10,
)

"""
    pi_load_columns

:element_id -> ID number of the load
:bus -> ID number of the bus where the load is placed
:load_p -> active load value
:load_q -> reactive load value
:status -> 0=IDLE, 1=ACTIVE
"""
const pi_load_columns = Dict(
    :element_id => :index,
    :bus => :load_bus,
    :load_p => :pd,
    :load_q => :qd,
    :status => :status,
)

"""
    ps_load_columns

:element_id -> ID number of the load
:bus -> ID number of the bus where the load is placed
:load -> active load value
:load_max -> maximum active load
:cost -> ID number of the associated cost
:status -> 0=IDLE, 1=ACTIVE
"""
const ps_load_columns = Dict( # This is meant for our use, not for importing powermodels cases.
    :element_id => -1,
    :bus => missing,
    :load => missing, # Active power
    :load_max => missing,
    :cost => missing,
    :status => 1,
)

"""
    line_columns

:element_id -> ID number of the line
:rate_a -> maximum current under rate A
:rate_b -> maximum current under rate B
:rate_c -> maximum current under rate C
:transformer -> is transformer?
:from_bus -> first end-bus ID
:to_bus -> second end-bus ID
:status -> 0=IDLE, 1=ACTIVE
:resistance -> line electric resistance
:reactance -> line reactance
:susceptance -> line susceptance
:ang_min -> minimum voltage angle difference
:ang_max -> maximum voltage angle difference
:q_to -> Reactive power reaching `to_bus`
:q_from -> Reactive power leaving `from_bus`
:p_to -> Real power reaching `to_bus`
:p_from -> Real power leaving `from_bus`
"""
const line_columns = Dict(
    :element_id => :index,
    :rate_a => :rate_a,
    :rate_b => :rate_b,
    :rate_c => :rate_c,
    :transformer => :transformer,
    :from_bus => :f_bus,
    :to_bus => :t_bus,
    :status => :br_status,
    :resistance => :br_r,
    :reactance => :br_x,
    :susceptance => :br_b,
    :ang_min => :angmin,
    :ang_max => :angmax,
    :q_to => :q_to,
    :q_from => :q_fr,
    :p_to => :p_to,
    :p_from => :p_from,
)

"""
    cost_columns

:element_id -> ID number of the cost (defaults to -1)
:coeffs -> coefficients of the cost function (defaults to `:cost`)
:model -> Piecewise linear (1) or polynomial (2), defaults to polynomial (2)
:ncost -> Number of coefficients in the polynomial cost or number of pairs (MWh, USD) in piecewise linear costs.
"""
const cost_columns = Dict(
    :element_id => -1,
    :coeffs => :cost,
    :model => :model,
    :ncost => :ncost,
)

"""
    Network

Basic type for storing a network as a group of DataFrames plus a PowerModels network case.

- Constructors:

    Network(; bus, gen, pi_load, ps_load, line, cost_gen, cost_load, pmc)

Create a `Network` containing the specified dataframes. Specifying a PowerModel case `pmc`
will solely store the case inside the `Network`, without using it to populate the
dataframes. If that is desired, another constructor should be used.

    Network(pmc::Dict{String, Any})

Create a `Network` from a PowerModel network.

    Network(path::AbstractString)

Create a `Network` from a Matpower case file.

    Network(case::Symbol)

Create a `Network` from a case. All cases are stored inside the `case_library` folder. In
order to list all valid case names, use `PowerModelsAPI.list_cases()`.
"""
struct Network
    bus::DataFrame
    gen::DataFrame
    pi_load::DataFrame
    ps_load::DataFrame
    line::DataFrame
    cost_gen::DataFrame
    cost_load::DataFrame

    pmc::Dict{String,Any} # PowerModels case
    results::Dict{String,Any} # Results from OPF. Always start empty.

    function Network(;
            bus::DataFrame = DataFrame(
                [[] for i in 1:length(bus_columns)],
                sort(collect(keys(bus_columns)))
            ),
            gen::DataFrame = DataFrame(
                [[] for i in 1:length(gen_columns)],
                sort(collect(keys(gen_columns)))
            ),
            pi_load::DataFrame = DataFrame(
                [[] for i in 1:length(pi_load_columns)],
                sort(collect(keys(pi_load_columns)))
            ),
            ps_load::DataFrame = DataFrame(
                [[] for i in 1:length(ps_load_columns)],
                sort(collect(keys(ps_load_columns)))
            ),
            line::DataFrame = DataFrame(
                [[] for i in 1:length(line_columns)],
                sort(collect(keys(line_columns)))
            ),
            cost_gen::DataFrame = DataFrame(
                [[] for i in 1:length(cost_columns)],
                sort(collect(keys(cost_columns)))
            ),
            cost_load::DataFrame = DataFrame(
                [[] for i in 1:length(cost_columns)],
                sort(collect(keys(cost_columns)))
            ),
            pmc::Dict{String, Any} = Dict{String,Any}(),
            results::Dict{String, Any} = Dict{String,Any}(),
        )
        return new(bus, gen, pi_load, ps_load, line, cost_gen, cost_load, pmc, results)
    end
end


"""
    build_df_from_pmc(columns::Dict{Symbol, <: Any}, block::Dict{String,Any})

Return a DataFrame as built from a PowerModels case dictionary. `columns` specifies which
quantites to grab and from where. `block` corresponds to a dictionary from a PowerModels
case.
"""
function build_df_from_pmc(columns::Dict{Symbol, <: Any}, block::Dict{String,Any})
    entries = Dict()
    for k in keys(columns)
        if isa(columns[k], Symbol)
            s = string(columns[k])
            entries[k] = Vector{Any}( # This is necessary for to avoid having columns
            # initialized as `NAType` and thus refusing other types. Not ideal, though.
                    [(s in keys(block[i]) ? block[i][s] : missing) for i in keys(block)]
                )
        else
            entries[k] = fill(columns[k], length(block))
        end
    end
    return allowmissing!(DataFrame(entries))
end

function Network(pmc::Dict)
    gen_df = build_df_from_pmc(gen_columns, pmc["gen"])
    gen_df[:cost] = collect(1:size(gen_df)[1]) # Since the only place where costs are stored
    # in the pmc is in pmc["gen"], the order should automatically be preserved.
    allowmissing!(gen_df)
    # When importing data from matpower we need to convert the generator costs to the data
    # type appropriate for the cost curve
    aux = deepcopy(pmc["gen"])
    to_warn = []
    for k in keys(aux)
        row = aux[k]["cost"]
        if aux[k]["model"] == 1
            if length(row) % 2 != 0
                error(LOGGER, "The type selected (PWL) and the cost data are incompatible")
            end
            tmp_ncost = length(row) / 2
            idx = range(1, tmp_ncost)
            mws = row[2 * idx - 1] #* u"MW*hr"
            costs = row[2 * idx] #* u"USD"
            aux[k]["cost"] = PWLCost(mw=mws, cost=costs)
            aux[k]["ncost"] = tmp_ncost
        elseif aux[k]["model"] == 2
            aux[k]["cost"] = PolynomialCost(row)
            aux[k]["ncost"] = length(row)
        else
            push!(to_warn, k)
        end
    end
    if !isempty(to_warn)
        warn(LOGGER, "Unable to set the costs for keys $to_warn. Please check.")
    end
    cost_gen_df = build_df_from_pmc(cost_columns, aux)
    cost_gen_df[:element_id] = collect(1:size(cost_gen_df)[1])
    allowmissing!(cost_gen_df)
    pi_load_df = build_df_from_pmc(pi_load_columns, pmc["load"])
    allowmissing!(pi_load_df)

    return Network(
        bus = build_df_from_pmc(bus_columns, pmc["bus"]),
        gen = gen_df,
        pi_load = pi_load_df,
        ps_load = DataFrames.DataFrame(
            [[] for i in 1:length(ps_load_columns)],
            sort(collect(keys(ps_load_columns)))
        ),
        line = build_df_from_pmc(line_columns, pmc["branch"]),
        cost_gen = cost_gen_df,
        cost_load = DataFrames.DataFrame(
            [[] for i in 1:length(cost_columns)],
            sort(collect(keys(cost_columns)))
        ),
        pmc = pmc,
        results = Dict{String,Any}(),
    )
    applyunits!(net)

    return net
end


"""
    applyunits!(net::Network)

This function annotates the Network structure with physical unit annotations.
For example, for the load column, it is originally in Float64. This function
will convert it to a type representing the units "u"MWh"). See Units.jl
for more information. We are assuming energy units as opposed to power units.
"""
function applyunits!(net::Network)
    net.pi_load[:load_p] = Array{UnitfulMissing}(net.pi_load[:load_p]*u"MWh")
    net.pi_load[:load_q] = Array{UnitfulMissing}(net.pi_load[:load_q]*u"MVARh")
    net.ps_load[:load] = Array{UnitfulMissing}(net.ps_load[:load]*u"MWh")
    net.ps_load[:load_max] = Array{UnitfulMissing}(net.ps_load[:load_max]*u"MWh")
    net.gen[:gen_p] = Array{UnitfulMissing}(net.gen[:gen_p]*u"MWh")
    net.gen[:p_max] = Array{UnitfulMissing}(net.gen[:p_max]*u"MWh")
    net.gen[:p_min] = Array{UnitfulMissing}(net.gen[:p_min]*u"MWh")
    net.gen[:gen_q] = Array{UnitfulMissing}(net.gen[:gen_q]*u"MVARh")
    net.gen[:q_max] = Array{UnitfulMissing}(net.gen[:q_max]*u"MVARh")
    net.gen[:q_min] = Array{UnitfulMissing}(net.gen[:q_min]*u"MVARh")
    net.gen[:startup_cost] = Array{UnitfulMissing}(net.gen[:startup_cost]*u"USD")
    net.gen[:ramp] = Array{UnitfulMissing}(net.gen[:ramp]*u"MWh/minute")
    net.bus[:base_kv] = Array{UnitfulMissing}(net.bus[:base_kv]*u"kV")
    #TODO: add cost units
    net.line[:p_to] = Array{UnitfulMissing}(net.line[:p_to]*u"MWh")
    net.line[:p_from] = Array{UnitfulMissing}(net.line[:rate_a]*u"MWh")
    net.line[:q_to] = Array{UnitfulMissing}(net.line[:rate_a]*u"MVARh")
    net.line[:q_from] = Array{UnitfulMissing}(net.line[:rate_a]*u"MVARh")
    net.line[:rate_a] = Array{UnitfulMissing}(net.line[:rate_a]*u"A")
    net.line[:rate_b] = Array{UnitfulMissing}(net.line[:rate_b]*u"A")
    net.line[:rate_c] = Array{UnitfulMissing}(net.line[:rate_c]*u"A")
    net.line[:resistance] = Array{UnitfulMissing}(net.line[:resistance]*u"Ω")
    net.line[:reactance] = Array{UnitfulMissing}(net.line[:reactance]*u"Ω")
    net.line[:susceptance] = Array{UnitfulMissing}(net.line[:susceptance]*u"S")
end

function stripunits!(net::Network)
    # Incomplete method to be extended as we annotate units
    net.pi_load[:load_p] = Array{Union{Missings.Missing,Float64}}(
        fustrip(net.pi_load[:load_p])
    )
    net.pi_load[:load_q] = Array{Union{Missings.Missing,Float64}}(
        fustrip(net.pi_load[:load_q])
    )
    net.ps_load[:load] = Array{Union{Missings.Missing,Float64}}(fustrip(net.ps_load[:load]))
    net.ps_load[:load_max] = Array{Union{Missings.Missing,Float64}}(
        fustrip(net.ps_load[:load_max])
    )
    net.gen[:gen_p] = Array{Union{Missings.Missing,Float64}}(fustrip(net.gen[:gen_p]))
    net.gen[:p_max] = Array{Union{Missings.Missing,Float64}}(fustrip(net.gen[:p_max]))
    net.gen[:p_min] = Array{Union{Missings.Missing,Float64}}(fustrip(net.gen[:p_min]))
    net.gen[:gen_q] = Array{Union{Missings.Missing,Float64}}(fustrip(net.gen[:gen_q]))
    net.gen[:q_max] = Array{Union{Missings.Missing,Float64}}(fustrip(net.gen[:q_max]))
    net.gen[:q_min] = Array{Union{Missings.Missing,Float64}}(fustrip(net.gen[:q_min]))
    net.gen[:startup_cost] = Array{Union{Missings.Missing,Float64}}(
        fustrip(net.gen[:startup_cost])
    )
    net.gen[:ramp] = Array{Union{Missings.Missing,Float64}}(fustrip(net.gen[:ramp]))
    net.bus[:base_kv] = Array{Union{Missings.Missing,Float64}}(fustrip(net.bus[:base_kv]))
    #TODO: strip cost units
    net.line[:p_to] = Array{Union{Missings.Missing,Float64}}(fustrip(net.line[:p_to]))
    net.line[:p_fr] = Array{Union{Missings.Missing,Float64}}(fustrip(net.line[:p_from]))
    net.line[:q_to] = Array{Union{Missings.Missing,Float64}}(fustrip(net.line[:q_to]))
    net.line[:q_from] = Array{Union{Missings.Missing,Float64}}(fustrip(net.line[:q_from]))
    net.line[:rate_a] = Array{Union{Missings.Missing,Float64}}(fustrip(net.line[:rate_a]))
    net.line[:rate_b] = Array{Union{Missings.Missing,Float64}}(fustrip(net.line[:rate_b]))
    net.line[:rate_c] = Array{Union{Missings.Missing,Float64}}(fustrip(net.line[:rate_c]))
    net.line[:resistance] = Array{Union{Missings.Missing,Float64}}(
        fustrip(net.line[:resistance])
    )
    net.line[:reactance] = Array{Union{Missings.Missing,Float64}}(
        fustrip(net.line[:reactance])
    )
    net.line[:susceptance] = Array{Union{Missings.Missing,Float64}}(
        fustrip(net.line[:susceptance])
    )
end

function matpower2pmc(path::AbstractString)
    return setlevel!(getlogger(PowerModels), "error") do
        PowerModels.parse_file(path)
    end
end

function Network(casepath::AbstractString)
    return Network(matpower2pmc(casepath))
end

"""
    add_bus!(
        net::Network;
        bus_type::Int=0,
        base_kv::Union{Missings.Missing,<:Number}=missing,
        name::AbstractString="none",
        element_id::Int=-1,
        voltage::Union{Missings.Missing,<:Number}=missing,
        volt_max::Union{Missings.Missing,<:Number}=missing,
        volt_min::Union{Missings.Missing,<:Number}=missing,
        coords::Union{Missings.Missing,Tuple{Float64,Float64}}=missing,
    )

Add a new bus to Network `net`. If `name` and `element_id` are not provided, reasonable
values will be automatically chosen.
"""
function add_bus!(
    net::Network;
    bus_type::Int=0,
    base_kv::Union{Missings.Missing,<:Number}=missing,
    name::AbstractString="none",
    element_id::Int=-1,
    voltage::Union{Missings.Missing,<:Number}=missing,
    volt_max::Union{Missings.Missing,<:Number}=missing,
    volt_min::Union{Missings.Missing,<:Number}=missing,
    coords::Union{Missings.Missing,Tuple{Float64,Float64}}=missing,
    area::Union{Missings.Missing,Int}=missing,
    zone::Union{Missings.Missing,Int}=missing,
)
    ids = net.bus[:element_id]
    if element_id == -1 || element_id in ids
        control = false
        if element_id in ids
            w = "Desired bus id, $element_id, is already in use. "
            control = true
        end
        element_id = isempty(ids) ? 1 : maximum(ids) + 1
        if control
            warn(LOGGER, w * "Using id = $element_id instead.")
        end
    end
    name == "none" ? name = "bus_" * string(element_id) : nothing
    push!(
        net.bus,
        Any[
            area,
            base_kv,
            bus_type,
            coords,
            element_id,
            name,
            volt_max,
            volt_min,
            voltage,
            zone
        ]
    )
end

"""
    add_gen!(
        net::Network;
        element_id::Int=-1,
        bus::Union{Missings.Missing,Int}=missing,
        gen_p::Union{Missings.Missing,<:Number}=missing,
        p_max::Union{Missings.Missing,<:Number}=missing,
        p_min::Union{Missings.Missing,<:Number}=missing,
        gen_q::Union{Missings.Missing,<:Number}=missing,
        q_max::Union{Missings.Missing,<:Number}=missing,
        q_min::Union{Missings.Missing,<:Number}=missing,
        cost::Union{Missings.Missing,Int}=missing,
        startup_cost::Union{Missings.Missing,<:Number}=missing,
        status::Int= 1,
        ramp::Union{Missings.Missing,<:Number}=missing,
    )

Add a new generator to a Network `net`. In case `element_id` and `status` are not provided,
a reasonable value will be chosen for `element_id`, while `status` wil be assumed as `1`
(the element is active).
"""
function add_gen!(
    net::Network;
    element_id::Int=-1,
    bus::Union{Missings.Missing,Int}=missing,
    gen_p::Union{Missings.Missing,<:Number}=missing,
    p_max::Union{Missings.Missing,<:Number}=missing,
    p_min::Union{Missings.Missing,<:Number}=missing,
    gen_q::Union{Missings.Missing,<:Number}=missing,
    q_max::Union{Missings.Missing,<:Number}=missing,
    q_min::Union{Missings.Missing,<:Number}=missing,
    cost::Union{Missings.Missing,Int}=missing,
    startup_cost::Union{Missings.Missing,<:Number}=missing,
    status::Int= 1,
    ramp::Union{Missings.Missing,<:Number}=missing,
)
    ids = net.gen[:element_id]
    if element_id == -1 || element_id in ids
        control = false
        if element_id in ids
            w = "Desired bus id, $element_id, is already in use. "
            control = true
        end
        element_id = isempty(ids) ? 1 : maximum(ids) + 1
        if control
            warn(LOGGER, w * "Using id = $element_id instead.")
        end
    end
    push!(net.gen, (Any[
        bus,
        cost,
        element_id,
        gen_p,
        gen_q,
        p_max,
        p_min,
        q_max,
        q_min,
        ramp,
        startup_cost,
        status,
    ]))
end

"""
    add_pi_load!(
        net::Network;
        element_id::Int=-1,
        bus::Union{Missings.Missing,Int}=missing,
        load_p::Union{Missings.Missing,<:Number}=missing,
        load_q::Union{Missings.Missing,<:Number}=missing,
        status::Int=1,
    )

Add a new price insensitive load to a Network `net`. In case `element_id` and `status` are
not provided, a reasonable value will be chosen for `element_id`, while `status` wil be
assumed as `1` (the element is active).
"""
function add_pi_load!(
    net::Network;
    element_id::Int=-1,
    bus::Union{Missings.Missing,Int}=missing,
    load_p::Union{Missings.Missing,<:Number}=missing,
    load_q::Union{Missings.Missing,<:Number}=missing,
    status::Int=1,
)
    ids = net.pi_load[:element_id]
    if element_id == -1 || element_id in ids
        control = false
        if element_id in ids && element_id != -1
            w = "Desired bus id, $element_id, is already in use. "
            control = true
        end
        element_id = isempty(ids) ? 1 : maximum(ids) + 1
        if control
            warn(LOGGER, w * "Using id = $element_id instead.")
        end
    end
    push!(net.pi_load, Any[bus, element_id, load_p, load_q, status])
end

"""
    add_load!(
        net::Network;
        element_id::Int=-1,
        bus::Union{Missings.Missing,Int}=missing,
        load_p::Union{Missings.Missing,<:Number}=missing,
        load_q::Union{Missings.Missing,<:Number}=missing,
        status::Int=1,
    )

Add a new price insensitive load to a Network `net`. In case `element_id` and `status` are
not provided, a reasonable value will be chosen for `element_id`, while `status` wil be
assumed as `1` (the element is active).
"""
function add_load!(
    net::Network;
    element_id::Int=-1,
    bus::Union{Missings.Missing,Int}=missing,
    load_p::Union{Missings.Missing,<:Number}=missing,
    load_q::Union{Missings.Missing,<:Number}=missing,
    status::Int=1,
)
    add_pi_load!(
        net,
        element_id=element_id,
        bus=bus,
        load_p=load_p,
        load_q=load_q,
        status=status
    )
end


"""
    add_ps_load!(
        net::Network;
        element_id::Int=-1,
        bus::Union{Missings.Missing,Int}=missing,
        load::Union{Missings.Missing,<:Number}=missing,
        load_max::Union{Missings.Missing,<:Number}=missing,
        cost::Union{Missings.Missing,Int}=missing,
        status::Int=1,
    )

Add a new price sensitive load to a Network `net`. In case `element_id` and `status` are
not provided, a reasonable value will be chosen for `element_id`, while `status` wil be
assumed as `1` (the element is active).
"""
function add_ps_load!(
    net::Network;
    element_id::Int=-1,
    bus::Union{Missings.Missing,Int}=missing,
    load::Union{Missings.Missing,<:Number}=missing,
    load_max::Union{Missings.Missing,<:Number}=missing,
    cost::Union{Missings.Missing,Int}=missing,
    status::Int=1,
)
    ids = net.ps_load[:element_id]
    if element_id == -1 || element_id in ids
        control = false
        if element_id in ids && element_id != -1
            w = "Desired bus id, $element_id, is already in use. "
            control = true
        end
        element_id = isempty(ids) ? 1 : maximum(ids) + 1
        if control
            warn(LOGGER, w * "Using id = $element_id instead.")
        end
    end
    push!(net.ps_load, Any[bus, cost, element_id, load, load_max, status])
end

"""
    add_line!(
        net::Network;
        rate_a::Union{Missings.Missing,<:Number}= missing,
        rate_b::Union{Missings.Missing,<:Number}= missing,
        rate_c::Union{Missings.Missing,<:Number}= missing,
        transformer::Bool= false,
        from_bus::Union{Missings.Missing,Int}= missing,
        to_bus::Union{Missings.Missing,Int}= missing,
        status::Int=1,
        element_id::Int=-1,
        resistance::Union{Missings.Missing,<:Number}= missing,
        reactance::Union{Missings.Missing,<:Number}= missing,
        susceptance::Union{Missings.Missing,<:Number}= missing,
        ang_min::Float64=-1.0,
        ang_max::Float64=1.0,
        p_to::Float64=0.0,
        p_from::Float64=0.0,
        q_to::Float64=0.0,
        q_from::Float64=0.0,
    )

Add a new line to a Network `net`. In case `element_id`, `transformer` and `status` are
not provided, a reasonable value will be chosen for `element_id`, while `status` wil be
assumed as `1` (the element is active) and `transformer` will be assumed as `false` (the
line is not a transformer).
"""
function add_line!(
    net::Network;
    rate_a::Union{Missings.Missing,<:Number}= missing,
    rate_b::Union{Missings.Missing,<:Number}= missing,
    rate_c::Union{Missings.Missing,<:Number}= missing,
    transformer::Bool= false,
    from_bus::Union{Missings.Missing,Int}= missing,
    to_bus::Union{Missings.Missing,Int}= missing,
    status::Int=1,
    element_id::Int=-1,
    resistance::Union{Missings.Missing,<:Number}= missing,
    reactance::Union{Missings.Missing,<:Number}= missing,
    susceptance::Union{Missings.Missing,<:Number}= missing,
    ang_min::Float64=-1.0,
    ang_max::Float64=1.0,
    p_to::Float64=0.0,
    p_from::Float64=0.0,
    q_to::Float64=0.0,
    q_from::Float64=0.0,
)
    ids = net.line[:element_id]
    if element_id == -1 || element_id in ids
        control = false
        if element_id in ids
            w = "Desired line id, $element_id, is already in use. "
            control = true
        end
        element_id = isempty(ids) ? 1 : maximum(ids) + 1
        if control
            warn(LOGGER, w * "Using id = $element_id instead.")
        end
    end
    push!(net.line, Any[
        ang_max,
        ang_min,
        from_bus,
        element_id,
        p_from,
        p_to,
        q_from,
        q_to,
        rate_a,
        rate_b,
        rate_c,
        reactance,
        resistance,
        status,
        susceptance,
        to_bus,
        transformer,
    ])
end

# We want to deal with piecewise linear costs and polynomial costs using multiple dispatch.
# For polynomial costs, the format of the cost curve consists of the sequence of coefficients
# of the polynomial, from the highest degree to the lowest degree
"""

    abstract type CostCurve end

Abstract type for cost curves. We have two types, one being PolynomialCost, the other PWLCost.
These types are introduced in order to deal with the possible formats in which cost curves can
be specified and to enforce type (and unit) consistency.
"""
abstract type CostCurve end

"""

PolynomialCost{T<:Number} <: CostCurve

# Field

- `coefficients::AbstractVector{T}`

The `coefficients` field is a vector that contains the numerical values of the coefficients of
the monomials of the polynomial cost function. If `n=length(coefficients)`, the cost is a
polynomial of degree n - 1. The first element of the array is the coefficient of the monomial
of degree n - 1, the last is the coefficient of the monomial of degree 0.

# NOTE:

Dimensions are not enforced as the coefficients have all different dimension (the element i
of a vector of length n has dimensions USD/(MWhr)^(n-i))
"""
struct PolynomialCost{T<:Number} <: CostCurve
    coefficients::AbstractVector{T}
end
coefficients(cost::PolynomialCost) = cost.coefficients
degree(cost::PolynomialCost) = length(cost.coefficients) - 1
n_cost(pol_cost::PolynomialCost) = length(pol_cost.coefficients)


"""

    is_convex(mw::AbstractVector{<:Number}, cost::AbstractVector{<:Number})

Determines if the data used for a PWL curve are defining a convex curve.

"""
function is_convex(mw::AbstractVector{<:Number}, cost::AbstractVector{<:Number})
    x = ustrip.(mw)
    y = ustrip.(cost)
    dx = x[2:end] - x[1:end-1]
    dy = y[2:end] - y[1:end-1]
    slopes = dy ./ dx
    convex = true
    if length(slopes) >= 2
        convex = prod((slopes[2:end] - slopes[1:end-1]) .>= 0)
    end
    return convex
end

"""

    PWLCost <: CostCurve

The cost function as a piecewise linear function of the MWhr.

# Fields

- `cost::AbstractVector`: the vector of the value of the cost function at the extrema
of the segments
- `mw::AbstractVector`: the vector of the extrema of the segments in which
the range for power generation is split.

# Constructor
For convenience, and to avoid errors due to the order, the default constructor uses kwargs:
```
    PWLCost(mw=mw, cost=cost)
```

# NOTE:
Units are explicitly given.

# TODO:
Make the dimensions and units parameters of the type.
"""
struct PWLCost <: CostCurve
    mw::AbstractVector{<:Units.PowerHour}
    cost::AbstractVector{<:Units.Currency}
# Inner constructor contains some built-in checks
    function PWLCost(;
        mw::AbstractVector{<:Units.PowerHour}=Units.PowerHour[],
        cost::AbstractVector{<:Units.Currency}=Units.Currency[],
    )
        if length(mw) != length(cost)
            error(LOGGER, "Malformed cost curve. Please check.")
        else
            if (eltype(mw) == asqtype(u"MW*hr")) && (eltype(cost) <: Units.Currency)
                if is_convex(mw, cost)
                    return new(mw, cost)
                else
                    error(LOGGER, "Cost curve is not convex. Please check.")
                end
            else
                error(LOGGER, "Bad units in cost curve. Please check.")
            end
        end
    end
end
costs(pwl_cost::PWLCost) = pwl_cost.cost
mws(pwl_cost::PWLCost) = pwl_cost.mw
n_segments(pwl_cost::PWLCost) = length(pwl_cost.cost) - 1
n_cost(pwl_cost::PWLCost) = length(pwl_cost.cost)

"""
    prices(pwl_cost::PWLCost)

This function returns the price for each block of the PWL curve, defined as
the slopes of the linear function for each block.
"""
function prices(pwl_cost::PWLCost)
    out = []
    mw = mws(pwl_cost)
    usd = costs(pwl_cost)
    for i in 1:n_segments(pwl_cost)
        push!(out, (usd[i+1] - usd[i])/(mw[i+1] - mw[i]))
    end
    return out
end

"""
    costcurve2pmc(c)

Convert the CostCurve types into PowerModels costs
"""
costcurve2pmc(c::CostCurve) = Float64[]
function costcurve2pmc(c::PWLCost)
    mw = ustrip(mws(c))
    cost = ustrip(costs(c))
    tmp = hcat(mw, cost)'
    return vec(tmp)
end
costcurve2pmc(c::PolynomialCost) = coefficients(c)

"""
    add_cost_gen!(
        net::Network;
        coeffs::Vector{<:Real}=Vector{Float64}(),
        gen_id::Union{Missings.Missing,Int}= missing,
        element_id::Int= -1,
        model::Int=2,
        ncost::Int=0,
    )

Add new generator cost to a Network `net`. If `element_id` is not specified, a reasonable
value will be adopted. This method applies to polynomial costs.
"""
function add_cost_gen!(
    net::Network;
    coeffs::PolynomialCost=PolynomialCost(Float64[]),
    gen_id::Union{Missings.Missing,Int}= missing,
    element_id::Int= -1,
)
    model = 2
    ncost = n_cost(coeffs)
    ids = net.cost_gen[:element_id]
    if element_id == -1 || element_id in ids
        control = false
        if element_id in ids && element_id != -1
            w = "Desired bus id, $element_id, is already in use. "
            control = true
        end
        element_id = isempty(ids) ? 1 : maximum(ids) + 1
        if control
            warn(LOGGER, w * "Using id = $element_id instead.")
        end
    end
    if ismissing(gen_id)
        push!(net.cost_gen, Any[coeffs, element_id, model, ncost])
    else
        gen_ids = net.gen[:element_id]
        if gen_id in gen_ids
            push!(net.cost_gen, Any[coeffs, element_id, model, ncost])
            net.gen[Array{Bool}(net.gen[:element_id] .== gen_id), :cost] = element_id
        else
            warn(LOGGER, "A generator with id $gen_id does not exist. Cost will not be created.")
        end
    end
end

"""
    add_cost_load!(
        net::Network;
        coeffs::Vector{<:Real}=Vector{Float64}(),
        load_id::Union{Missings.Missing,Int}= missing,
        element_id::Int= -1,
        model::Int=2,
        ncost::Int=0,
    )

Add new load cost to a Network `net`. If `element_id` is not specified, a reasonable
value will be adopted. This method applies to polynomial costs.
"""
function add_cost_load!(
    net::Network;
    coeffs::PolynomialCost=PolynomialCost(Float64[]),
    load_id::Union{Missings.Missing,Int}= missing,
    element_id::Int= -1,
)
    model = 2
    ncost = n_cost(coeffs)
    ids = net.cost_load[:element_id]
    if element_id == -1 || element_id in ids
        control = false
        if element_id in ids && element_id != -1
            w = "Desired bus id, $element_id, is already in use. "
            control = true
        end
        element_id = isempty(ids) ? 1 : maximum(ids) + 1
        if control
            warn(LOGGER, w * "Using id = $element_id instead.")
        end
    end
    if ismissing(load_id)
        push!(net.cost_load, Any[coeffs, element_id, model, ncost])
    else
        load_ids = net.gen[:element_id]
        if load_id in load_ids
            push!(net.cost_load, Any[coeffs, element_id, model, ncost])
            net.ps_load[
                Array{Bool}(net.ps_load[:element_id] .== load_id),
                :cost
            ] = element_id
        else
            warn(LOGGER, "A PS load with id $load_id does not exist. Cost will not be created.")
        end
    end
end


"""
    add_cost_gen!(
        net::Network;
        pwl_cost::PWLCost;
        gen_id::Union{NAtype,Int}= NA,
        element_id::Int= -1
        model::Int=1,
        ncost::Int=0,
    )

Add new generator cost to a Network `net`. If `element_id` is not specified, a reasonable
value will be adopted. This method applies to piecewise linear costs (1).
"""
function add_cost_gen!(
    net::Network,
    pwl_cost::PWLCost;
    gen_id::Union{Missings.Missing,Int}=missing,
    element_id::Int= -1,
)
    model = 1
    ncost = length(mws(pwl_cost))
    ids = net.cost_gen[:element_id]
    if element_id == -1 || element_id in ids
        control = false
        if element_id in ids && element_id != -1
            w = "Desired bus id, $element_id, is already in use. "
            control = true
        end
        element_id = isempty(ids) ? 1 : maximum(ids) + 1
        if control
            warn(LOGGER, w * "Using id = $element_id instead.")
        end
    end
    if ismissing(gen_id)
        push!(net.cost_gen, Any[pwl_cost, element_id, model, ncost])
    else
        gen_ids = net.gen[:element_id]
        if gen_id in gen_ids
            push!(net.cost_gen, Any[pwl_cost, element_id, model, ncost])
            net.gen[Array{Bool}(net.gen[:element_id] .== gen_id), :cost] = element_id
        else
            warn(LOGGER, "A generator with id $gen_id does not exist. Cost will not be created.")
        end
    end
end

"""
    add_cost_load!(
        net::Network;
        coeffs::Tuple{Vector}=([], []),
        load_id::Union{NAtype,Int}= NA,
        element_id::Int= -1,
        model::Int=1,
        ncost::Int=0,
    )

Add new load cost to a Network `net`. If `element_id` is not specified, a reasonable
value will be adopted. This method applies to piecewise linear costs (1).
"""
function add_cost_load!(
    net::Network,
    pwl_cost::PWLCost;
    load_id::Union{Missings.Missing, Int}=missing,
    element_id::Int=-1,
)
    model = 1
    ncost = length(mws(pwl_cost))
    ids = net.cost_load[:element_id]
    if element_id == -1 || element_id in ids
        control = false
        if element_id in ids && element_id != -1
            w = "Desired bus id, $element_id, is already in use. "
            control = true
        end
        element_id = isempty(ids) ? 1 : maximum(ids) + 1
        if control
            warn(LOGGER, w * "Using id = $element_id instead.")
        end
    end
    if ismissing(load_id)
        push!(net.cost_load, Any[pwl_cost, element_id, model, ncost])
    else
        load_ids = net.gen[:element_id]
        if load_id in load_ids
            push!(net.cost_load, Any[pwl_cost, element_id, model, ncost])
            net.ps_load[
                Array{Bool}(net.ps_load[:element_id] .== load_id),
                :cost
            ] = element_id
        else
            warn(LOGGER, "A PS load with id $load_id does not exist. Cost will not be created.")
        end
    end
end

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
        "dcline" => Dict() # We don't use this.
    )
    stripunits!(net)
    gen_dict = Dict()
    to_warn = Int[]
    for r in eachrow(gen(net))
        old = Dict()
        if in(keys(pmc(net)), "gen") && string(r[:element_id]) in keys(pmc(net)["gen"])
            old = pmc(net)["gen"][string(r[:element_id])]
        end
        there = !isempty(old)
        coeffs = there ? old["cost"] : [0.0]
        coeffs = !ismissing(r[:cost]) ? costcurve2pmc(cost_gen(net)[:coeffs][r[:cost]]) : coeffs
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
            "cost" => (!there || !ismissing(r[:cost]) ) ? coeffs : old["cost"], # Default or old value
            "model" => !there  ? -1 : old["model"], # Default or old value
            "ncost" => !there  ? 0 : old["ncost"], # Default or old value
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
        if !there && isa(coeffs, PolynomialCost)
            gen_dict["cost"] = coefficients(coeffs)
            gen_dict["model"] = 2
            gen_dict["ncost"] = n_cost(coeffs)
        elseif !there && isa(coeffs, PWLCost)
            mw_pts = ustrip.(mws(coeffs))
            cost_pts = ustrip.(costs(coeffs))
            mix = vec(hcat(mw_pts, cost_pts)') # This is the format used by MATPOWER and PowerModels
            gen_dict["cost"] = mix
            gen_dict["model"] = 1
            gen_dict["ncost"] = length(mw_pts)
        else
            push!(to_warn, r[:element_id])
        end
    end
    if !isempty(to_warn)
        warn(LOGGER, "Using default costs for generator cost ids $to_warn")
    end
    # Here we also need to add the price sensitive loads as generators
    last_gen = maximum(gen(net)[:element_id])
    to_warn = Int[]
    for r in eachrow(ps_load(net))
        old = Dict()
        if string(r[:element_id]) in keys(pmc(net)["gen"])
            old = pmc(net)["gen"][string(r[:element_id])]
        end
        there = !isempty(old)
        coeffs = there ? old["cost"] : [0.0]
        coeffs = !ismissing(r[:cost]) ? cost_load(net)[:coeffs][r[:cost]] : coeffs
        gen_dict[string(r[:element_id] + last_gen)] = Dict(
            "index" => (!there || !ismissing(r[:element_id])) ? r[:element_id] + last_gen : old["index"],
            "gen_bus" => (!there || !ismissing(r[:bus])) ? r[:bus] : old["gen_bus"],
            "pmax" => (!there || !ismissing(r[:load_max])) ? r[:load_max] : old["pmax"],
            "cost" => (!there || !ismissing(r[:cost]) ) ? costcurve2pmc(cost_load(net)[:coeffs][r[:cost]]) : old["cost"], # Default or old value
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
        if !there && isa(coeffs, PolynomialCost)
            # Load is negative generation. By serving load, the total cost receives a negative contribution
            # due to the amount payed by the utility.
            gen_dict["cost"] = -coefficients(coeffs) # The cost curve has to be flipped along the x axis.
            aux = gen_dict["pmin"]
            gen_dict["pmin"] = -gen_dict["pmax"] # The cost curve has to be flipped along the y axis.
            gen_dict["pmax"] = -aux
            gen_dict["model"] = 2
            gen_dict["ncost"] = n_cost(coeffs)
        elseif !there && isa(coeffs, PWLCost)
            mw_pts = -ustrip.(mws(coeffs)) # The mws have to be flipped.
            cost_pts = -ustrip.(costs(coeffs)) # The cost has to be flipped too. It's a cost for the utility, profit for generators.
            mix = vec(hcat(mw_pts, cost_pts)') # This is the format used by MATPOWER and PowerModels
            gen_dict["cost"] = mix
            gen_dict["model"] = 1
            gen_dict["ncost"] = length(mw_pts)
            gen_dict["pmin"] = -gen_dict["pmax"] # The cost curve has to be flipped along the y axis. Load is negative generation
            gen_dict["pmax"] = 0.0
        else
            push!(to_warn, r[:element_id])
        end
    end
    if !isempty(to_warn)
        warn(LOGGER, "Used default costs for ps load cost ids $to_warn.")
    end
    outnet["gen"] = gen_dict
    branch_dict = Dict()
    for r in eachrow(line(net))
        old = Dict()
        if in(keys(pmc(net)), "branch") && string(r[:element_id]) in keys(pmc(net)["branch"])
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
            "b_to" => (!there || !ismissing(r[:element_id])) ? r[:p_to] : old["b_to"],
            "b_fr" => (!there || !ismissing(r[:element_id])) ? r[:p_from] : old["b_fr"],
            "g_to" => (!there || !ismissing(r[:element_id])) ? r[:q_to] : old["g_to"],
            "g_fr" => (!there || !ismissing(r[:element_id])) ? r[:q_from] : old["g_fr"],
            "shift" => 0.0,
            "tap" => 1.0,
        )
    end
    outnet["branch"] = branch_dict
    bus_dict = Dict()
    for r in eachrow(bus(net))
        old = Dict()
        if in(keys(pmc(net)), "bus") && string(r[:element_id]) in keys(pmc(net)["bus"])
            old = pmc(net)["bus"][string(r[:element_id])]
        end
        there = !isempty(old)
        bus_dict[string(r[:element_id])] = Dict(
            "index" => (!there || !ismissing(r[:element_id])) ? r[:element_id] : old["index"],
            "bus_name" => (!there || !ismissing(r[:element_id])) ? r[:name] : old["bus_name"],
            "vm" => (!there || !ismissing(r[:element_id])) ? r[:voltage] : old["vm"],
            "vmax" => (!there || !ismissing(r[:element_id])) ? r[:volt_max] : old["vmax"],
            "vmin" => (!there || !ismissing(r[:element_id])) ? r[:volt_min] : old["vmin"],
            "bus_type" => (!there || !ismissing(r[:element_id])) ? r[:bus_type] : old["bus_type"],
            "base_kv" => (!there || !ismissing(r[:element_id])) ? r[:base_kv] : old["base_kv"],
            "zone" => (!there || !ismissing(r[:element_id])) ? r[:zone] : old["zone"],
            "bus_i" => (!there || !ismissing(r[:element_id])) ? r[:element_id] : old["bus_i"],
            "area" => (!there || !ismissing(r[:element_id])) ? r[:area] : old["area"],
            "va" => 0.0,
        )
    end
    outnet["bus"] = bus_dict
    load_dict = Dict()
    for r in eachrow(pi_load(net))
        old = Dict()
        if in(keys(pmc(net)), "load") && string(r[:element_id]) in keys(pmc(net)["load"])
            old = pmc(net)["load"][string(r[:element_id])]
        end
        there = !isempty(old)
        load_dict[string(r[:element_id])] = Dict(
            "index" => (!there || !ismissing(r[:element_id])) ? r[:element_id] : old["index"],
            "load_bus" => (!there || !ismissing(r[:element_id])) ? r[:bus] : old["load_bus"],
            "qd" => (!there || !ismissing(r[:element_id])) ? r[:load_q] : old["qd"],
            "pd" => (!there || !ismissing(r[:element_id])) ? r[:load_p] : old["pd"],
            "status" => (!there || !ismissing(r[:element_id])) ? r[:status] : old["status"],
        )
    end
    outnet["load"] = load_dict
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

### Getters ###
bus(net::Network) = net.bus
buses(net::Network) = bus(net)
gen(net::Network) = net.gen
gens(net::Network) = gen(net)
generator(net::Network) = gens(net)
generators(net::Network) = generator(net)
pi_load(net::Network) = net.pi_load
ps_load(net::Network) = net.ps_load
load(net::Network) = Dict("pi_load" => pi_load(net), "ps_load" => ps_load(net))
loads(net::Network) = load(net)
line(net::Network) = net.line
lines(net::Network) = line(net)
cost_gen(net::Network) = net.cost_gen
gen_cost(net::Network) = cost_gen(net)
cost_load(net::Network) = net.cost_load
load_cost(net::Network) = cost_load(net)
cost(net::Network) = Dict("cost_gen" => cost_gen(net), "cost_load" => cost_load(net))
costs(net::Network) = cost(net)
pmc(net::Network) = net.pmc
results(net::Network) = net.results

# GenericPowerModels
ACPPowerModel(net::Network) = ACCPowerModel(pmc(net))
APIACPPowerModel(net::Network) = APIACPPowerModel(pmc(net))
DCPPowerModel(net::Network) = DCPPowerModel(pmc(net))
SOCWRPowerModel(net::Network) = SOCWRPowerModel(pmc(net))
QCWRPowerModel(net::Network) = QCWRPowerModel(pmc(net))

# Run OPF
run_dc_opf(net::Network, solver) = run_dc_opf(pmc(net), solver)
run_dc_opf(net::Network) = run_dc_opf(pmc(net))
run_ac_opf(net::Network, solver) = run_ac_opf(pmc(net), solver)
run_ac_opf(net::Network) = run_ac_opf(pmc(net))
run_opf(net::Network, model::DataType, solver) = run_opf(pmc(net), model, solver)
run_opf(net::Network, model::DataType) = run_opf(pmc(net), model)

# Extract results
"""
    lmps(net::Network)

Return a DataFrame with the LMPs for all buses after an OPF run. The prices are extracted
from the Lagrange multipliers for the Kirchhoff's conservation law. Due to conventions
in PowerModels.jl, the LMP are the negative of these Lagrange multipliers.

"""
function lmps(net::Network)
    isempty(results(net)) && return DataFrame()
    lmps = Float64[]
    buses = Int64[] # Here assuming buses are just numbered. If we are to use general names
    # we should go to either strings or symbols
    for k in keys(results(net)["solution"]["bus"])
        push!(lmps, -results(net)["solution"]["bus"][k]["lam_kcl_r"]) # the relation between Lagrange
        # multipliers and LMPs requires a -1, here.
        push!(buses, parse(Int64, k))
    end
    # Let's sort it so it looks decent
    p = sortperm(buses)
    buses = buses[p]
    lmps = lmps[p]
    return DataFrame(:bus => buses, :lmp => lmps)
end

"""
    shadow_prices_lines(net::Network)

Return a DataFrame containing the shadow prices for all lines after an OPF run.
"""
function shadow_prices_lines(net::Network)
    isempty(results(net)) && return DataFrame()
    lower = Float64[]
    upper = Float64[]
    lines = Int64[]
    for k in keys(results(net)["solution"]["branch"])
        push!(lower, results(net)["solution"]["branch"][k]["mu_sm_to"])
        push!(upper, results(net)["solution"]["branch"][k]["mu_sm_fr"])
        push!(lines, parse(Int64, k))
    end
    # Let's sort it so it looks decent
    p = sortperm(lines)
    lines = lines[p]
    lower = lower[p]
    upper = upper[p]
    return DataFrame([lines, lower, upper], [:line, :lower, :upper])
end

"""
    res_gen(net::Network)

Return a DataFrame with the power supplied by each generator as a solution to an OPF.
"""
function res_gen(net::Network)
    isempty(results(net)) && return DataFrame()
    gen_p = Float64[]
    gen_q = Float64[]
    gens = []
    for k in keys(results(net)["solution"]["gen"])
        push!(gen_p, results(net)["solution"]["gen"][k]["pg"])
        push!(gen_q, results(net)["solution"]["gen"][k]["qg"])
        push!(gens, parse(Int64, k))
    end
    return DataFrame([gens, gen_p, gen_q], [:gen, :gen_p, :gen_q])
end

"""
    converged(net::Network)

Return `true` if an OPF was ran and successfully converged, `false` otherwise (inclunding if
no OPF was ran).
"""
function converged(net::Network)
    isempty(results(net)) && return false
    return results(net)["status"] == :LocalInfeasible ? false : true
end

"""
    max_load_percent!(net::Network, maxload::Real)

Scale line ratings such that the maximum load is equal to `maxload`% of the original value.
"""
function max_load_percent!(net::Network, maxload::Real)
    if maxload <= 0
        warn(LOGGER, "Maximum load percentage has to be strictly greater than zero.")
    else
        # Let's update both the dataframes and the pmc such that this function can be used
        # at any moment without the risk of changes not getting through.

        maxload /= 100 # % to fraction

        # First the dataframes:
        line(net)[:rate_a] *= maxload
        line(net)[:rate_b] *= maxload
        line(net)[:rate_c] *= maxload

        # Now the pmc:
        for k in keys(pmc(net)["branch"])
            pmc(net)["branch"][k]["rate_a"] *= maxload
            pmc(net)["branch"][k]["rate_b"] *= maxload
            pmc(net)["branch"][k]["rate_c"] *= maxload
        end
    end
end

"""
    max_load_percent!(pmc::Dict, maxload::Real)

Scale line ratings such that the maximum load is equal to `maxload`% of the original value.
"""
function max_load_percent!(pmc::Dict, maxload::Real)
    if maxload <= 0
        warn(LOGGER, "Maximum load percentage has to be strictly greater than zero.")
    else
        maxload /= 100 # % to fraction
        for k in keys(pmc["branch"])
            pmc["branch"][k]["rate_a"] *= maxload
            pmc["branch"][k]["rate_b"] *= maxload
            pmc["branch"][k]["rate_c"] *= maxload
        end
    end
end

"""
    infeasible(net::Network)

Check if the total demand can be satisfied by the generators. Returns `true` in case no
solution is possible.
"""
function infeasible(net::Network)
    tot_load = sum(pi_load(net)[:load_p])
    max_gen = sum(gen(net)[:p_max])
    return tot_load > max_gen
end
