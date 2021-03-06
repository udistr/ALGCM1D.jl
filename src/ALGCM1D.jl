module ALGCM1D

using Plots
using Plots.PlotMeasures
using LinearAlgebra
using Dates
using Dierckx
using Measures
using Parameters
using Statistics

include("constants.jl")
include("grid.jl")
include("parameters.jl")
include("time.jl")

include("heaviside.jl")
include("delta.jl")
include("angle_of_incidence.jl")
include("rad.jl")
include("psim.jl")
include("psis.jl")
include("surface.jl")
include("holtslag.jl")
include("turbulence.jl")
include("moist.jl")

include("atm_parameters.jl")
include("hydraulic_variables.jl")
include("soil_parameters.jl")
include("Temperature_dist.jl")
include("metric_head.jl")
include("soil.jl")

include("initialize.jl")
include("atmosphere.jl")

p1 = plot(
    1,
    #xlim = (0, 360*24*2),
    #ylim = (0, 50),
    marker = 2,
    legend=false
)
p2 = plot(
    1,
    #xlim = (0, 360*24*2),
    #ylim = (0.01, 0.015),
    marker = 2,
    legend=false
)

function init()
  include("initialize.jl")
end


function run(start_time)

  global model_time=start_time # model clock
  global plt
  global LH,SH,Qnet,SW,LW,TAU,EVAP # surface fluxes
  global ΘA,qA,ρA,PA,UA # atmosphere state variables
  global head,TS,theta_l,theta_v,rho_vs # soil state variables
  global KAm,KAt,γcq,γct,γcm # vertical mixing terms
   
  
  Q=zeros(n)
  sw=zeros(n)
  lw=zeros(n)
  lh=zeros(n)
  sh=zeros(n)

  for i=1:n

      ########################################################################
      # time
      ########################################################################

      EpochDays=Dates.date2epochdays(Date(model_time))
      TimeOfDay=(Dates.datetime2epochms(model_time)/1000-EpochDays*24*60*60)/3600;
      #julian day
      JulianDay=Dates.datetime2julian(model_time)-Dates.datetime2julian(DateTime(year(model_time)));

      println(i,"   ",model_time)
      ########################################################################
      # Radiation
      # Input: date, lat and lon
      # Output: longwave and shortwave radiation
      # Future input: cloud fraction and chemestry
      ########################################################################
      # Total preciptible water in cm
      PW=sum(qA ./ 1000 .* ρA .* ΔZA)*100 # (kg/kg) / (kg/m^3) * (kg/m^3) * m * (100 cm/m)= cm
      SW,LW=rad(TS[end],TAup,theta_l[end],PW,lat,lon,JulianDay,TimeOfDay)
      #println("SW")
      #println(SW)
      #println("---------")
    
      if rem(i,360)==0;
        #p1=scatter([i],[TS[end]])
        #display(p1)
        #p1=scatter!([Dates.format.(model_time,"dd:HH")],[theta_l[end]],xlabel="time [DAY:HOUR]",ylabel="PBL height",xrot=60)
        p1=scatter!([Dates.format.(model_time,"dd:HH")],[hpbl[1]],xlabel="time [DAY:HOUR]",ylabel="PBL height",xrot=60)
        display(p1)
        #p2=scatter!([Dates.format.(model_time,"dd:HH")],[hpbl],xlabel="time [DAY:HOUR]",ylabel="PBL height",xrot=60)
        #display(p2)
      end
      ########################################################################
      # Surface
      # Input from land: surface volumetric water content and temperature
      # Input from atmosphere: surface specific humidity, temperature and wind
      # Output: latent heat, sensible heat and wind stress.
      ########################################################################
      println("Qnet=",Qnet[1])
      #println("qA=",qA[1])


      TAU,LH,SH,EVAP=surface(0,UA[1],TS[end],ΘA[1],rho_vs[end]/ρA[1],qA[1],ρA[1],ZAC[1])
      Qnet=SW-LH-SH-LW
      println("SW=",SW[1])
      ########################################################################
      # turbulence
      # Input: temperature, specific humidity, air densiy, surface wind speed, 
      # latent heat, sensible heat, wind stress
      # Output: vertical diffusion for momentum, heat and water
      ########################################################################
      KAm,KAt,γcq,γct,γcm=turbulence(ΘA,qA,UA,UA*0,LH,SH,ustar,ρA,ZAF)

      ########################################################################
      # Soil
      # Input: fluxes
      ########################################################################
      TS,head,theta_l,theta_v,rho_vs=soil(TS,head,EVAP[1],ΔT)

      ########################################################################
      # Atmosphere
      # Input: fluxes, albedo
      ########################################################################
      UA,ΘA,qA=atmosphere(UA,ΘA,qA,TAU,EVAP,SH,KAm,KAt,γcq,γct,ΔZA)
      
      ########################################################################
      # Moist
      # Input: fluxes, albedo
      ########################################################################
      if domoist==1
        θA,qA=moist(ΘA,qA,ρA,PA)
      end
      TA=(ΘA).*(PA./PA0).^(2/7).-d2k;
      PA=PA0.*exp.(-cumsum(ΔZA./(Rgas.*(TA.+d2k)./gravity_mks)));
      ρA=PA./(Rgas.*(TA.+d2k));

      if rem(i,360000)==0;#3600*6/dt
        title = plot(title = model_time, grid = false, axis = false,
        showaxis = false, bottom_margin = -10pt,yaxis=nothing,titlefont = 100)
        #=
        println("UA")
        println(UA)
        println("Θv")
        println(Θv)
        println("qA")
        println(qA)
        println("KAm")
        println(KAm)
        =#
        f1=plot(ZAC,UA,ylabel=("UA"));
        f2=plot(ZAC,ΘA.-273.15,ylabel=("ΘA")); #ylim([19 21]);...
        f3=plot(ZAC,qA,ylabel=("qA"));
        f4=plot(ZAC,KAm,ylabel=("KAm"));
        f5=plot(soil_numerical_parameters.z_s,theta_l,ylabel=("theta_l"));
        f6=plot(soil_numerical_parameters.z_s,TS.-273.15,ylabel=("TS")); #ylim([19 21]);...
        plt=plot(title,f1,f2,f3,f4,f5,f6,layout=grid(7, 1, heights=[0.1 ,0.15, 0.15, 0.15, 0.15, 0.15, 0.15]),legend=false,titlefontsize=6,
        size = (700, 1200),margin=8mm)
        display(plt)
      end
      model_time=model_time+Second(ΔT)
      Q[i]=Qnet
      sw[i]=sum(SW)
      lw[i]=LW
      lh[i]=LH
      sh[i]=SH
  end

end

end

#if abspath(PROGRAM_FILE) == @__FILE__
#  main()
#end