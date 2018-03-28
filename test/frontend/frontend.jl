@testset "PowerModels FrontEnd" begin
    net = Network(case_files[5])
    @test size(PMA.bus(net))[1] == 14
    add_bus!(net)
    add_gen!(net)
    add_ps_load!(net)
    add_pi_load!(net)
    add_load!(net)
    add_line!(net)
    add_cost_gen!(net)
    add_cost_load!(net)
    @test size(PMA.lines(net))[1] == 21
    @test size(PMA.buses(net))[1] == 15
    @test size(PMA.generators(net))[1] == 6
    @test size(PMA.loads(net)["ps_load"])[1] == 1
    @test size(PMA.loads(net)["pi_load"])[1] == 13
    build_pmc!(net)
    @test length(PMA.pmc(net)["bus"]) == 15
    old_rate = PMA.pmc(net)["branch"]["1"]["rate_a"]
    PMA.max_load_percent!(PMA.pmc(net), 50)
    @test PMA.pmc(net)["branch"]["1"]["rate_a"] == 0.5 * old_rate
end
