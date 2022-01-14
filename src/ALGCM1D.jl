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
  global jd,t
  global plt
  global LH,SH,Qnet,SW,LW,TAU,EVAP
  
  Q=zeros(n)
  sw=zeros(n)
  lw=zeros(n)
  lh=zeros(n)
  sh=zeros(n)

  for i=1:n

      ########################################################################
      # time
      ########################################################################

      epochdays=Dates.date2epochdays(Date(model_time))
      t=(Dates.datetime2epochms(model_time)/1000-epochdays*24*60*60)/3600;
      #julian day
      jd=Dates.datetime2julian(model_time)-Dates.datetime2julian(DateTime(year(model_time)));

      println(i,"   ",model_time)
      ########################################################################
      # Radiation
      # Input: date, lat and lon
      # Output: longwave and shortwave radiation
      # Future input: cloud fraction and chemestry
      ########################################################################
      rad(model_time)
      #println("SW")
      #println(SW)
      #println("---------")
    
      if rem(i,360)==0;
        #p1=scatter([i],[TS[end]])
        #display(p1)
        p2=scatter!([Dates.format.(model_time,"dd:HH")],[hpbl],xlabel="time [DAY:HOUR]",ylabel="PBL height",xrot=60)
        display(p2)
      end
      ########################################################################
      # Surface
      # Input from land: surface volumetric water content and temperature
      # Input from atmosphere: surface specific humidity, temperature and wind
      # Output: latent heat, sensible heat and wind stress.
      ########################################################################
      TAU,LH,SH,EVAP=surface(0,UA[1],TS[end],ΘA[1],q_v[end],qA[1],ρA[1],ZAC[1])
      Qnet=SW-LH-SH-LW
      ########################################################################
      # turbulence
      # Input: temperature, specific humidity, air densiy, surface wind speed, 
      # latent heat, sensible heat, wind stress
      # Output: vertical diffusion for momentum, heat and water
      ########################################################################
      turbulence(model_time)

      ########################################################################
      # Soil
      # Input: fluxes
      ########################################################################
      soil(model_time)

      ########################################################################
      # Atmosphere
      # Input: fluxes, albedo
      ########################################################################
      atmosphere(model_time)    

      if rem(i,3600000)==0;#3600*6/dt
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
        f2=plot(ZAC,Θv.-273.15,ylabel=("Θv")); #ylim([19 21]);...
        f3=plot(ZAC,qA,ylabel=("q"));
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