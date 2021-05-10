using Dates

function rad(mtime)

  global SW,LW,albedo

  ########################################################################
  # radiation
  ########################################################################
  AOI=angle_of_incidence(lat,lon,jd,t); 
  AOI<0 ? AOI=0 : AOI

  #define albedo
  if theta_l[end]<0.1
    albedo = 0.25;
  else
    if theta_l[end]>0.25
      albedo = 0.1;
    else
      albedo = 0.35-theta_l[end];
    end
  end

  SW=albedo.*(SunConstant.*AOI);
  ϵa=0.725+0.17*log10(sum(qA/1000*1.2*ΔZA)*100)
  LW=emissivity*sb*((TS[end]+d2k).^4 .- ϵa.*(TAup).^4);

end