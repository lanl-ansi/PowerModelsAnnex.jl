using PowerModelsAnnex
using Unitful
using PowerSystemsUnits
PMA = PowerModelsAnnex

Unitful.register(PowerSystemsUnits)

import JuMP
import PowerModels
PMs = PowerModels

import InfrastructureModels

import Ipopt

import Memento
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModels), "error")

using Test

pms_path = joinpath(dirname(pathof(PowerModels)), "..")

# default setup for solvers
ipopt_solver = JuMP.with_optimizer(Ipopt.Optimizer, tol=1e-6, print_level=0)


# this will work because PowerModels is a dependency
case_files = Dict(
    "case3"      => "$(pms_path)/test/data/matpower/case3.m",
    "case5"      => "$(pms_path)/test/data/matpower/case5.m",
    "case5_asym" => "$(pms_path)/test/data/matpower/case5_asym.m",
    "case5_gap"  => "$(pms_path)/test/data/matpower/case5_gap.m",
    "case5_dc"   => "$(pms_path)/test/data/matpower/case5_dc.m",
    "case14"     => "$(pms_path)/test/data/matpower/case14.m",
    "case30"     => "$(pms_path)/test/data/matpower/case30.m"
)

@testset "PowerModelsAnnex" begin

    include("form/acr.jl")
    include("form/wr.jl")

    include("pglib/api.jl")
    include("pglib/sad.jl")

    include("model/pf.jl")
    include("model/opf.jl")

    include("frontend/frontend.jl")

end
