# PowerModelsAnnex FrontEnd

## Motivation

Providing a user-friendly interface with the network format used by PowerModels. Networks are represented as DataFrame objects, allowing for easier visualisation and access. Simple functions for adding different types of elements are provided.

## Usage

All functions are built around the `Network` object. It stores both the dataframes for each component of the network as well as the dictionary used by PowerModels for computations. A `Network` can be created from scratch or from a Matpower case.
```julia
using PowerModelsAnnex
# Create empty Network
net = Network()
# Create Network from Matpower case
net = Network("data/case14.m")
```

Different elements of a `Network` can be inspected:
```julia
julia> PMA = PowerModelsAnnex
PowerModelsAnnex

julia> PMA.bus(net)
14×8 DataFrames.DataFrame
│ Row │ base_kv │ bus_type │ coords  │ element_id │ name         │ volt_max │ volt_min │ voltage │
├─────┼─────────┼──────────┼─────────┼────────────┼──────────────┼──────────┼──────────┼─────────┤
│ 1   │ 0.0     │ 1        │ missing │ 4          │ Bus 4     HV │ 1.06     │ 0.94     │ 1.019   │
│ 2   │ 0.0     │ 3        │ missing │ 1          │ Bus 1     HV │ 1.06     │ 0.94     │ 1.06    │
│ 3   │ 0.0     │ 1        │ missing │ 12         │ Bus 12    LV │ 1.06     │ 0.94     │ 1.055   │
│ 4   │ 0.0     │ 2        │ missing │ 2          │ Bus 2     HV │ 1.06     │ 0.94     │ 1.045   │
│ 5   │ 0.0     │ 2        │ missing │ 6          │ Bus 6     LV │ 1.06     │ 0.94     │ 1.07    │
│ 6   │ 0.0     │ 1        │ missing │ 11         │ Bus 11    LV │ 1.06     │ 0.94     │ 1.057   │
│ 7   │ 0.0     │ 1        │ missing │ 13         │ Bus 13    LV │ 1.06     │ 0.94     │ 1.05    │
│ 8   │ 0.0     │ 1        │ missing │ 5          │ Bus 5     HV │ 1.06     │ 0.94     │ 1.02    │
│ 9   │ 0.0     │ 1        │ missing │ 14         │ Bus 14    LV │ 1.06     │ 0.94     │ 1.036   │
│ 10  │ 0.0     │ 1        │ missing │ 7          │ Bus 7     ZV │ 1.06     │ 0.94     │ 1.062   │
│ 11  │ 0.0     │ 2        │ missing │ 8          │ Bus 8     TV │ 1.06     │ 0.94     │ 1.09    │
│ 12  │ 0.0     │ 1        │ missing │ 10         │ Bus 10    LV │ 1.06     │ 0.94     │ 1.051   │
│ 13  │ 0.0     │ 1        │ missing │ 9          │ Bus 9     LV │ 1.06     │ 0.94     │ 1.056   │
│ 14  │ 0.0     │ 2        │ missing │ 3          │ Bus 3     HV │ 1.06     │ 0.94     │ 1.01    │

julia> PMA.gen(net)
5×12 DataFrames.DataFrame
│ Row │ bus │ cost │ element_id │ gen_p │ gen_q  │ p_max │ p_min │ q_max │ q_min │ ramp │ startup_cost │ status │
├─────┼─────┼──────┼────────────┼───────┼────────┼───────┼───────┼───────┼───────┼──────┼──────────────┼────────┤
│ 1   │ 6   │ 1    │ 4          │ 0.0   │ 0.122  │ 1.0   │ 0.0   │ 0.24  │ -0.06 │ 0.0  │ 0.0          │ 1      │
│ 2   │ 1   │ 2    │ 1          │ 2.324 │ -0.169 │ 3.324 │ 0.0   │ 0.1   │ 0.0   │ 0.0  │ 0.0          │ 1      │
│ 3   │ 8   │ 3    │ 5          │ 0.0   │ 0.174  │ 1.0   │ 0.0   │ 0.24  │ -0.06 │ 0.0  │ 0.0          │ 1      │
│ 4   │ 2   │ 4    │ 2          │ 0.4   │ 0.424  │ 1.4   │ 0.0   │ 0.5   │ -0.4  │ 0.0  │ 0.0          │ 1      │
│ 5   │ 3   │ 5    │ 3          │ 0.0   │ 0.234  │ 1.0   │ 0.0   │ 0.4   │ 0.0   │ 0.0  │ 0.0          │ 1      │
# Load is divided in price sensitive load and price insensitive load.
julia> PMA.load(net)
Dict{String,DataFrames.DataFrame} with 2 entries:
  "ps_load" => 0×6 DataFrames.DataFrame…
  "pi_load" => 11×4 DataFrames.DataFrame…

julia> PMA.pi_load(net)
11×4 DataFrames.DataFrame
│ Row │ bus │ element_id │ load  │ status │
├─────┼─────┼────────────┼───────┼────────┤
│ 1   │ 4   │ 1          │ 0.478 │ 1      │
│ 2   │ 12  │ 2          │ 0.061 │ 1      │
│ 3   │ 2   │ 3          │ 0.217 │ 1      │
│ 4   │ 6   │ 4          │ 0.112 │ 1      │
│ 5   │ 11  │ 5          │ 0.035 │ 1      │
│ 6   │ 13  │ 6          │ 0.135 │ 1      │
│ 7   │ 5   │ 7          │ 0.076 │ 1      │
│ 8   │ 14  │ 8          │ 0.149 │ 1      │
│ 9   │ 10  │ 9          │ 0.09  │ 1      │
│ 10  │ 9   │ 10         │ 0.295 │ 1      │
│ 11  │ 3   │ 11         │ 0.942 │ 1      │
```

Elements can be added via `add` functions:
```julia
julia> PMA.add_bus!(net, name="MyBus", base_kv=110)

julia> PMA.bus(net)
15×8 DataFrames.DataFrame
│ Row │ base_kv │ bus_type │ coords  │ element_id │ name         │ volt_max │ volt_min │ voltage │
├─────┼─────────┼──────────┼─────────┼────────────┼──────────────┼──────────┼──────────┼─────────┤
│ 1   │ 0.0     │ 1        │ missing │ 4          │ Bus 4     HV │ 1.06     │ 0.94     │ 1.019   │
│ 2   │ 0.0     │ 3        │ missing │ 1          │ Bus 1     HV │ 1.06     │ 0.94     │ 1.06    │
│ 3   │ 0.0     │ 1        │ missing │ 12         │ Bus 12    LV │ 1.06     │ 0.94     │ 1.055   │
│ 4   │ 0.0     │ 2        │ missing │ 2          │ Bus 2     HV │ 1.06     │ 0.94     │ 1.045   │
│ 5   │ 0.0     │ 2        │ missing │ 6          │ Bus 6     LV │ 1.06     │ 0.94     │ 1.07    │
│ 6   │ 0.0     │ 1        │ missing │ 11         │ Bus 11    LV │ 1.06     │ 0.94     │ 1.057   │
│ 7   │ 0.0     │ 1        │ missing │ 13         │ Bus 13    LV │ 1.06     │ 0.94     │ 1.05    │
│ 8   │ 0.0     │ 1        │ missing │ 5          │ Bus 5     HV │ 1.06     │ 0.94     │ 1.02    │
│ 9   │ 0.0     │ 1        │ missing │ 14         │ Bus 14    LV │ 1.06     │ 0.94     │ 1.036   │
│ 10  │ 0.0     │ 1        │ missing │ 7          │ Bus 7     ZV │ 1.06     │ 0.94     │ 1.062   │
│ 11  │ 0.0     │ 2        │ missing │ 8          │ Bus 8     TV │ 1.06     │ 0.94     │ 1.09    │
│ 12  │ 0.0     │ 1        │ missing │ 10         │ Bus 10    LV │ 1.06     │ 0.94     │ 1.051   │
│ 13  │ 0.0     │ 1        │ missing │ 9          │ Bus 9     LV │ 1.06     │ 0.94     │ 1.056   │
│ 14  │ 0.0     │ 2        │ missing │ 3          │ Bus 3     HV │ 1.06     │ 0.94     │ 1.01    │
│ 15  │ 110     │ 0        │ missing │ 15         │ MyBus        │ missing  │ missing  │ missing │
```

It can be checked if the `Network` corresponds to the converged solution of an OPF run:
```julia
julia> converged(net)
false
```

`Network` objects can be exported to the dictionary format used by PowerModels:
```julia 
pmc = network2pmc(net)
```
