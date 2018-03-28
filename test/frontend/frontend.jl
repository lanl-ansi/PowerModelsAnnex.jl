@testset "PowerModels FrontEnd" begin
    net = Network(case_files[5])
    @test size(FrontEnd.bus(net))[1] == 14
    add_bus!(net)
    add_gen!(net)
    add_ps_load!(net)
    add_pi_load!(net)
    add_load!(net)
    add_line!(net)
    add_cost_gen!(net)
    add_cost_load!(net)
    @test size(FrontEnd.lines(net))[1] == 21
    @test size(FrontEnd.buses(net))[1] == 15
    @test size(FrontEnd.generators(net))[1] == 6
    @test size(FrontEnd.loads(net)["ps_load"])[1] == 1
    @test size(FrontEnd.loads(net)["pi_load"])[1] == 13
    build_pmc!(net)
    @test length(FrontEnd.pmc(net)["bus"]) == 15
    old_rate = FrontEnd.pmc(net)["branch"]["1"]["rate_a"]
    FrontEnd.max_load_percent!(FrontEnd.pmc(net), 50)
    @test FrontEnd.pmc(net)["branch"]["1"]["rate_a"] == 0.5 * old_rate
end
