# PowerModelsAnnex JuMP Models

## Motivation

Integrating novel power model formulations in to PowerModels.jl can be time consuming and requires 
significant understanding of the PowerModels software design.  Consequently PowerModels.jl is not 
ideal for rapid prototyping of novel problem formulations.  To address this issue the files 
provided in this directory defined functions that build JuMP models from scratch, which are identical to those produced by 
PowerModels.jl, for a few of the most common formulations.  These functions leverage the 
data processing tools in PowerModels.jl and conform to the same naming conventions. 
Unit tests are used to ensure that the solutions produced by these from-scratch models match those of 
PowerModels.jl.


## Usage

The functions provided here designed to setup a JuMP model.  Notably, they are not concerned with
reading data files, solving the model, or reporting the solution.  These tasks are left to the user.

For example, a simple AC OPF work flow would be as follows,
```
using PowerModelsAnnex
using PowerModels
using JuMP
using Ipopt

ipopt = Ipopt.Optimizer
model = Model(with_optimizer(ipopt))

data = PowerModels.parse_file("case3.m")
build_ac_opf(data, model)
optimize!(model)
```

Once the model is solved the solution can be extracted as follows,
```
for (i, bus) in data["bus"]
    if bus["bus_type"] != 4
        println("$i, $(value(model[:t][bus["index"]])), $(value(model[:v][bus["index"]]))")
    end
end

for (i, gen) in data["gen"]
    if gen["gen_status"] != 0
        println("$i, $(value(model[:pg][gen["index"]])), $(value(model[:qg][gen["index"]]))")
    end
end
```
Note that all values are given in per unit and radians, as this is the internal PowerModels data format.

For additional examples of how to use these functions see the files in `PowerModelsAnnex.jl/test/model`.

