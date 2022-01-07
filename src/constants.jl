#constants


#Atmosphere
ρO0=1029; # kg/m^3ρO
cp=4184; #J/K/kg
cpa=1005; #J/K/kg
sb = 5.670e-8;
d2k=273.15;
Av=2500000;
gamma_blk=0.01;
cvapor_fac=640380;
cvapor_exp=5107.4;
cvapor_exp_ice=5897.8;
saltsat=0.980;
p0=101300; #reference pressure
gravity_mks=9.81;
humid_fac=0.606;
tref=25; #reference temperature
sref=25; #reference salinity
alpha=2e-4; #expansion coefficients for temperature
beta=7.4e-4; #expansion coefficients for salt
emissivity=0.97; # ocean emissivity
kappa=1e-5; # m^2/kg air absorption coefficient
Rgas=287.05;
karman=0.41;
albedo=0.04; #ocean albedo
SunConstant=1362
Rv=461 #gas constant for water vapor J kg-1 K-1
humid_fac=0.606;



#Land
using Parameters
@with_kw struct Const{R}
  Cs::R = 1.92e6; #soild heat capacity [J m^-3 K^-1]
  Cw::R = 4.18e6; # water heat capacity [J m^-3 K^-1]
  Cv::R = 1.8e6; # vapor heat capacity [J m^-3 K^-1]
  Ca::R = 1200; # heat capacity of air [J m^-3 K^-1]
  k::R = 0.41; # von-Karman constant [-]
  sigma::R = 5.6704e-8; # Stephan-Boltzman constant [W m^-2]
  M::R = 0.018015; #water molecular weight [kg mol^-1]
  g::R = 9.81; # gravity accelaration [m sec^-2]
  Rg::R = 8.314; # gas constant [J mol^-1 K^-1]
  gamma_o::R = 0.07189; # surface water tension [Kg sec^-2]
end

global constants=Const{Float64}()