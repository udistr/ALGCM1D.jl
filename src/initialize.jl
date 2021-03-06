
############################################################################
# Initial arrays
############################################################################

global UA=zeros(vleva); global UA1=zeros(vleva); global UAst=zeros(vleva);
global qA=zeros(vleva); global qA1=zeros(vleva); global qAst=zeros(vleva);
global ΘA=zeros(vleva); global ΘA1=zeros(vleva); global ΘAst=zeros(vleva);
global ρA=zeros(vleva);
global TA=zeros(vleva);
global PA=zeros(vleva);

global qS=0.01

############################################################################
# Initial conditions
############################################################################
# atmospheric temperature in DegC
TA=28 .-0.007.*ZAC;
global PA0=101300;
PA=PA0.*exp.(-cumsum(ΔZA./(Rgas.*(TA.+d2k)./gravity_mks)));
ρA=PA./(Rgas.*(TA.+d2k));
# atmospheric potential temperature in kelvin
ΘA[:]=(TA.+d2k).*(PA[1]./PA).^(2/7);
ΘA0=ΘA[end].+1;
RH=0.8
qA=saltsat*cvapor_fac*exp.(-cvapor_exp./(TA.+d2k))./ρA.*RH;
#println("qA0")
#println(qA)

UA=(0.2/karman).*(log.(ZAC./1e-4));

############################################################################

global LH=zeros(1); #latent heat flux
global SH=zeros(1); #sensible heat flux
global TAU=zeros(1); # wind stress
global LW=zeros(1); # net longwave
global SW=zeros(1); # net shortwavwe
global Qnet=zeros(1); # net heat flux
global MOL=zeros(1); # monin obokuv length
global EVAP=zeros(1); # Evaporation
global KAm=zeros(vleva); # verical mixing of momntum
global KAt=zeros(vleva); # verical mixing of scalar
global γct=zeros(vleva); 
global γcq=zeros(vleva); 
global γcm=zeros(vleva);

global ustar=0.01;
FWflux=zeros(size(Qnet));
ρO=1000
q0=1.01*qA[1]
ΘA0=300



############################################################################
# SOIL
############################################################################



#initial conditions
T_initial = 273+28;
h_initial = -3;

global head=h_initial*ones(soil_numerical_parameters.Ns,1)
global TS=T_initial*ones(soil_numerical_parameters.Ns,1)
global q_l=zeros(soil_numerical_parameters.Ns,1)
global q_v=zeros(soil_numerical_parameters.Ns,1)
global rho_vs=zeros(soil_numerical_parameters.Ns,1)


Cp_soil=zeros(soil_numerical_parameters.Ns,1)

struct shydraulic{R}
  C::R
  K::R
  K_vh::R
  K_h::R
  K_lT::R
  K_vT::R
  K_T::R
end

theta_l, theta_v= hydraulic_variables(head, TS, constants, soil_parameters);
Cp_soil = constants.Cs.*(soil_parameters.theta_s.-theta_l).+constants.Cw.*theta_l+constants.Cv.*theta_v;