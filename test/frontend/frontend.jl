using Missings

@testset "Costs" begin
    @testset "Polynomial Cost" begin
        pol = PMA.PolynomialCost([0.1, 1, 0.0])
        @test PMA.coefficients(pol) == [0.1, 1.0, 0.0]
        @test PMA.degree(pol) == 2
    end
    @testset "Convexity Checks" begin
        cost = Float64[0, 1, 3, 6]u"USD"
        mw = Float64[0, 1, 2, 3]u"MW*hr"
        @test PMA.is_convex(mw, cost)
    end
    @testset "PWL Cost" begin
        cost = Float64[0, 1, 3, 6] * 1u"USD"
        mw = Float64[0, 1, 2, 3] * 1u"MWh"
        pwl_cost = PMA.PWLCost(mw=mw, cost=cost)
        @test PMA.n_segments(pwl_cost) == 3
        @test PMA.costs(pwl_cost) == cost
        @test PMA.mws(pwl_cost) == mw
        @test length(PMA.prices(pwl_cost)) == 3
    end

    @testset "Conversion to PMC/MATPOWER" begin
        cost = Float64[0, 1, 3, 6] * 1u"USD"
        mw = Float64[0, 1, 2, 3] * 1u"MWh"
        pwl_cost = PMA.PWLCost(mw=mw, cost=cost)
        cc_pmc = PMA.costcurve2pmc(pwl_cost)
        @test isa(cc_pmc, Vector{Float64})
        @test cc_pmc == Float64[0, 0, 1, 1, 2, 3, 3, 6]
        pol = PMA.PolynomialCost([0.0, 1.0])
        cc_pmc = PMA.costcurve2pmc(pol)
        @test isa(cc_pmc, Vector{Float64})
        @test cc_pmc == [0.0, 1.0]
    end
end

@testset "PowerModels Network" begin
    @testset "Piecewise linear costs" begin
        net = Network(case_files["case14"])
        @test size(PMA.bus(net))[1] == 14
        add_bus!(net)
        add_gen!(net)
        add_ps_load!(net)
        add_pi_load!(net)
        add_load!(net)
        add_line!(net)
        example_pwl_cost = PMA.PWLCost(mw=[0, 1, 2] * 1.0u"MWh", cost=[0, 10, 20] * 1.0u"USD")
        ex_pwl_cost_2 = PMA.PWLCost(mw=[0, 2, 3] * 1.0u"MWh", cost=[0, 10, 20] * 1.0u"USD")
        add_cost_gen!(net, example_pwl_cost, gen_id = 6)
        add_cost_load!(net, ex_pwl_cost_2, load_id = 1)
        @test size(PMA.lines(net))[1] == 21
        @test size(PMA.buses(net))[1] == 15
        @test size(PMA.generators(net))[1] == 6
        @test size(PMA.loads(net)["ps_load"])[1] == 1
        @test size(PMA.loads(net)["pi_load"])[1] == 13
        @test size(PMA.cost_gen(net), 1) == 6
        @test size(PMA.load_cost(net), 1) == 1
        build_pmc!(net)
        @test length(PMA.pmc(net)["bus"]) == 15
        @test length(PMA.pmc(net)["gen"]) == 7 # 6 generators + 1 ps_load
        old_rate = PMA.pmc(net)["branch"]["1"]["rate_a"]
        PMA.max_load_percent!(PMA.pmc(net), 50)
        @test PMA.pmc(net)["branch"]["1"]["rate_a"] == 0.5 * old_rate
    end


    @testset "Polynomial costs" begin
        net = Network(case_files["case14"])
        @test size(PMA.bus(net))[1] == 14
        add_bus!(net)
        add_gen!(net)
        add_ps_load!(net)
        add_pi_load!(net)
        add_load!(net)
        add_line!(net)
        ex_pol_curve = PMA.PolynomialCost(Float64[0, 10, 0])
        add_cost_gen!(net, coeffs=ex_pol_curve, gen_id = 6)
        add_cost_load!(net, coeffs=ex_pol_curve, load_id = 1)
        @test size(PMA.lines(net))[1] == 21
        @test size(PMA.buses(net))[1] == 15
        @test size(PMA.generators(net))[1] == 6
        @test size(PMA.loads(net)["ps_load"])[1] == 1
        @test size(PMA.loads(net)["pi_load"])[1] == 13
        @test size(PMA.cost_gen(net), 1) == 6
        @test size(PMA.cost_load(net), 1) == 1
        build_pmc!(net)
        @test length(PMA.pmc(net)["bus"]) == 15
        @test length(PMA.pmc(net)["gen"]) == 7 # 6 generators + 1 ps_load
        old_rate = PMA.pmc(net)["branch"]["1"]["rate_a"]
        PMA.max_load_percent!(PMA.pmc(net), 50)
        @test PMA.pmc(net)["branch"]["1"]["rate_a"] == 0.5 * old_rate
    end

    net = Network(case_files["case14"])
    pm = PMA.network2pmc(net)
    @test isa(pm, Dict)
    for k in ["branch", "gen", "load", "shunt", "storage", "bus"]
        @test k in keys(pm)
    end
    PMA.applyunits!(net)
    @test isa(net.pi_load[:load_p], Array{Union{<:Unitful.Quantity, Missings.Missing}})
    PMA.stripunits!(net)
    @test isa(net.pi_load[:load_p], Array{Union{Missings.Missing,Float64}})
    @test_throws UndefVarError applyunits!(Dict{String, AbstractArray}())
    @test_throws UndefVarError stripunits!(Dict{String, AbstractArray}())
end
