#### DC Optimal Power Flow ####

# This file provides a pedagogical example of modeling the DC Optimal Power
# Flow problem using the Julia Mathematical Programming package (JuMP) and the
# PowerModels package for data parsing.

# This file can be run by calling `include("dc-opf.jl")` from the Julia REPL or
# by calling `julia dc-opf.jl` in Julia v0.6.

# Developed by Line Roald (@lroald) and Carleton Coffrin (@ccoffrin)


###############################################################################
# 0. Initialization
###############################################################################

# Load Julia Packages
#--------------------
using PowerModels
using Ipopt
using JuMP


# Instancate a Solver
#--------------------

nlp_solver = IpoptSolver(print_level=0)
# note: print_level changes the amount of solver information printed to the terminal


# Load System Data
# ----------------

if VERSION < v"0.7.0"
    powermodels_path = Pkg.dir("PowerModels")
else
    powermodels_path = joinpath(dirname(pathof(PowerModels)), "..")
end

file_name = "$(powermodels_path)/test/data/matpower/case5.m"
# note: change this string to modify the network data that will be loaded

# load the data file
data = PowerModels.parse_file(file_name)

# use build_ref to filter out inactive components
ref = PowerModels.build_ref(data)[:nw][0]
# Note: ref contains all the relevant system parameters needed to build the OPF model
# When we introduce constraints and variable bounds below, we use the parameters in ref.


###############################################################################
# 1. Building the Optimal Power Flow Model
###############################################################################

# Initialize a JuMP Optimization Model
#-------------------------------------
model = Model(solver = nlp_solver)


# Add Optimization and State Variables
# ------------------------------------

# Add voltage angles va for each bus
@variable(model, va[i in keys(ref[:bus])])
# note: [i in keys(ref[:bus])] adds one `va` variable for each bus in the network

# Add active power generation variable pg for each generator (including limits)
@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])

# Add power flow variables p to represent the active power flow for each branch
@variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

# Build JuMP expressions for the value of p[(l,i,j)] and p[(l,j,i)] on the branches
p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))
# note: this is used to make the definition of nodal power balance simpler

# Add power flow variables p_dc to represent the active power flow for each HVDC line
@variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])


# Add Objective Function
# ----------------------

# index representing which side the HVDC line is starting
from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])

# Minimize the cost of active power generation and cost of HVDC line usage
# assumes costs are given as quadratic functions
@objective(model, Min,
    sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
    sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
)


# Add Constraints
# ---------------

# Fix the voltage angle to zero at the reference bus
for (i,bus) in ref[:ref_buses]
    @constraint(model, va[i] == 0)
end

# Nodal power balance constraints
for (i,bus) in ref[:bus]
    # Build a list of the loads and shunt elements connected to the bus i
    bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
    bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

    # Active power balance at node i
    @constraint(model,
        sum(p_expr[a] for a in ref[:bus_arcs][i]) +                  # sum of active power flow on lines from bus i +
        sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==     # sum of active power flow on HVDC lines from bus i =
        sum(pg[g] for g in ref[:bus_gens][i]) -                 # sum of active power generation at bus i -
        sum(load["pd"] for load in bus_loads) -                 # sum of active load consumption at bus i -
        sum(shunt["gs"] for shunt in bus_shunts)*1.0^2          # sum of active shunt element injections at bus i
    )
end

# Branch power flow physics and limit constraints
for (i,branch) in ref[:branch]
    # Build the from variable id of the i-th branch, which is a tuple given by (branch id, from bus, to bus)
    f_idx = (i, branch["f_bus"], branch["t_bus"])

    p_fr = p[f_idx]                     # p_fr is a reference to the optimization variable p[f_idx]

    va_fr = va[branch["f_bus"]]         # va_fr is a reference to the optimization variable va on the from side of the branch
    va_to = va[branch["t_bus"]]         # va_fr is a reference to the optimization variable va on the to side of the branch

    # Compute the branch parameters and transformer ratios from the data
    g, b = PowerModels.calc_branch_y(branch)

    # DC Power Flow Constraint
    @constraint(model, p_fr == -b*(va_fr - va_to))
    # note: that upper and lower limits on the power flow (i.e. p_fr) are not included here.
    #   these limits were already enforced for p (which is the same as p_fr) when
    #   the optimization variable p was defined (around line 65).


    # Voltage angle difference limit
    @constraint(model, va_fr - va_to <= branch["angmax"])
    @constraint(model, va_fr - va_to >= branch["angmin"])
end

# HVDC line constraints
for (i,dcline) in ref[:dcline]
    # Build the from variable id of the i-th HVDC line, which is a tuple given by (hvdc line id, from bus, to bus)
    f_idx = (i, dcline["f_bus"], dcline["t_bus"])
    # Build the to variable id of the i-th HVDC line, which is a tuple given by (hvdc line id, to bus, from bus)
    t_idx = (i, dcline["t_bus"], dcline["f_bus"])   # index of the ith HVDC line which is a tuple given by (line number, to bus, from bus)
    # note: it is necessary to distinguish between the from and to sides of a HVDC line due to power losses

    # Constraint defining the power flow and losses over the HVDC line
    @constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
end



###############################################################################
# 3. Solve the Optimal Power Flow Model and Review the Results
###############################################################################

# Solve the optimization problem
status = solve(model)

# Check the value of the objective function
cost = getobjectivevalue(model)
println("The cost of generation is $(cost).")

# Check the value of an optimization variable
# Example: Active power generated at generator 1
pg1 = getvalue(pg[1])
println("The active power generated at generator 1 is $(pg1*ref[:baseMVA]) MW.")
# note: the optimization model is in per unit, so the baseMVA value is used to restore the physical units
