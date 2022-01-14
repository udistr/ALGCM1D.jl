using Dates
using ALGCM1D
#include("./src/ALGCM1D.jl")

start_time=DateTime(0,7,1)
ALGCM1D.run(start_time)
#include("constants.jl")
#include("psim.jl")
#include("psis.jl")
#include("holtslag_psis.jl")
#include("holtslag_psim.jl")
#include("VelocityScale.jl")

#VelocityScale.(1:100,100.,-0.01,0.01,300,0.06)
#include("./surface.jl")
#TAU,LH,SH,EVAP,_,_,_,MOL=surface(0,5,25,26,0.01,0.02,1.2,25)##