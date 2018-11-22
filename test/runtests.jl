using PowerModelsAnnex
using Unitful
PMA = PowerModelsAnnex

using JuMP
using PowerModels
PMs = PowerModels

using InfrastructureModels

using Ipopt

using Memento
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(PowerModels), "error")

using Compat.Test

if VERSION < v"0.7.0-"
    pms_path = Pkg.dir("PowerModels")
end

if VERSION > v"0.7.0-"
    pms_path = joinpath(dirname(pathof(PowerModels)), "..")
end

# default setup for solvers
ipopt_solver = IpoptSolver(tol=1e-6, print_level=0)

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
