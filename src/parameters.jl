#include("constants.jl")

# schemes
cpl=1;
domoist=1;

#boundary conditions
TAup=270; #top model temperature
UAup=(0.2/karman).*(log((ZAF[end]+ΔZA[end]/2)/1e-4)); #top model temperature
UOdown=0; #lower boundary model current

# numerics
Aimp=0.5; # Atmospheric implicitness
Oimp=0.5; # Oceanic implicitness



@with_kw struct snparameters{R}
  #soil numerical parameters and domain
  Zs::R = 2.; #soil's depth [m]
  Ns::Int = 201; #number of cells in soil colomn
  dz_s::R = Zs/(Ns-1); #soil cell size [m]
  z_s = Zs.-((1:Ns).-1).*dz_s; #soil depth at each cell [m]
  eps_h::R = 1e-5; # matric head tolerance
  eps_T::R = 1e-5; #soil temperature tolerance
end


soil_numerical_parameters=snparameters{Float64}()
