function rad(TS,TAup,theta_l,PW,lat,lon,JulianDay,TimeOfDay)

  ########################################################################
  # radiation
  # currently only calculates albedo and net shortwave and longwave 
  # radiation
  ########################################################################
  AOI=angle_of_incidence(lat,lon,JulianDay,TimeOfDay); 
  AOI<0 ? AOI=0 : AOI

  #define albedo
  if theta_l<0.1
    albedo = 0.25;
  else
    if theta_l>0.25
      albedo = 0.1;
    else
      albedo = 0.35-theta_l;
    end
  end

  SW=(1-albedo).*(SunConstant.*AOI);
  # longwave radiation at the surface from Stensrud eq. 8.24-8.25
  ϵa=0.725+0.17*log10(sum(qA/1000*1.2*ΔZA)*100)
  LW=emissivity*sb*((TS)^4 - ϵa*(TAup)^4);

  return SW,LW
end