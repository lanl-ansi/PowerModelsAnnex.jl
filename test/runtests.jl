using PowerModelsAnnex
using Unitful
using PowerModelsAnnex.Units
PMA = PowerModelsAnnex

using Memento
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(PowerModels), "error")

using JuMP
using PowerModels
PMs = PowerModels

using Ipopt

using Base.Test

# default setup for solvers
ipopt_solver = IpoptSolver(tol=1e-6, print_level=0)

# this will work because PowerModels is a dependency
case_files = Dict(
    "case3"      => "../../PowerModels/test/data/matpower/case3.m",
    "case5"      => "../../PowerModels/test/data/matpower/case5.m",
    "case5_asym" => "../../PowerModels/test/data/matpower/case5_asym.m",
    "case5_gap"  => "../../PowerModels/test/data/matpower/case5_gap.m",
    "case5_dc"   => "../../PowerModels/test/data/matpower/case5_dc.m",
    "case14"     => "../../PowerModels/test/data/matpower/case14.m",
    "case30"     => "../../PowerModels/test/data/matpower/case30.m"
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
