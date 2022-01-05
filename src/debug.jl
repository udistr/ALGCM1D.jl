using Dates

include("./src/ALGCM1D.jl")


stime=DateTime(0,7,1)
ALGCM1D.run(stime)
#include("constants.jl")
#include("holtslag_psis.jl")
#include("holtslag_psim.jl")
#include("VelocityScale.jl")

#VelocityScale.(1:100,100.,-0.01,0.01,300,0.06)
