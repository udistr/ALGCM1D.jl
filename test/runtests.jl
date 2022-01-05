using ALGCM1D
using Test
using Dates


@testset "ALGCM1D.jl" begin
    # Write your tests here.
    stime=DateTime(0,7,1)
    ALGCM1D.run(stime)
end
