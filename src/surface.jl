function surface(US,UA,TS,ΘA,qV,qA,ρA,Z)
  ########################################################################
  # surface
  ########################################################################
  
  #roughness length (stensrud, table 2.3)
  #Ice 10^-4
  #Grass (mown) 10^-2
  #Long grass, rocky ground 0.05
  #Pasture 0.20
  #Suburban housing 0.6
  #Forests, cities 1–5
  
  #Wikipedia
  # Mud flats, snow; no vegetation, no obstacles	0.005
  zm=0.005
  zh=0.005
  ze=0.005

  ψm=0
  ψs=0
  CD=karman^2/log(Z/zm)^2
  CH=sqrt(CD)*karman/log(Z/zh)
  CE=sqrt(CD)*karman/log(Z/ze)
  #kinematic surface flux
  ustar=sqrt(CD)*abs(UA-US)
  tstar=CH/sqrt(CD)*abs(UA-US)*(TS-ΘA)
  qstar=CE/sqrt(CD)*abs(UA-US)*(qV-qA)

  println("TS=",TS)
  println("ΘA=",ΘA)
  MOL=0;

  for i=1:5
    wθv=tstar*(1+humid_fac*qA)+humid_fac*ΘA*qstar
    MOL=-ΘA*ustar^3/(karman*gravity_mks*wθv)
    
    ψm1=max(min(psim(Z/MOL),1),-1)
    ψs1=max(min(psis(Z/MOL),1),-1)

    CD1=karman^2/(log(Z/zm)-ψm1)^2
    CH1=sqrt(CD1)*karman/(log(Z/zh)-ψs1)
    CE1=sqrt(CD1)*karman/(log(Z/ze)-ψs1)

    ∂CD∂ψm=(CD1-CD)/2
    ∂CH∂ψm=(CH1-CH)/2
    ∂CE∂ψm=(CE1-CE)/2

    CD=CD+∂CD∂ψm
    CH=CH+∂CH∂ψm
    CE=CE+∂CE∂ψm
    ψm=ψm1
    ψs=ψs1

    ustar=sqrt(CD)*abs(UA-US)
    tstar=sqrt(CH)*abs(UA-US)*(TS-ΘA)
    qstar=sqrt(CE)*abs(UA-US)*(qV-qA)

  end

  #rs = 10*exp(35.63*(0.15-theta_top)); # soil surface resistance to vapor flow
  #rv = 1/(u*constants.k^2)*(log((atm_parameters.z_ref-atm_parameters.d+atm_parameters.z_H)/atm_parameters.z_oH)+psi_H)*
  #    (log((atm_parameters.z_ref-atm_parameters.d+atm_parameters.z_m)/atm_parameters.z_om)+psi_m);
  #rH = rv;

  SH=ρA*CH*cp*abs(UA)*(TS-ΘA)
  LH=ρA*CE*Av*abs(UA)*(qV-qA)
  TAU=ρA*CD*abs(UA)*UA;

  EVAP=LH/Av;
  #Qnet=SW-LH-SH-LW;

  return TAU, LH, SH, EVAP, CD, CH, CE, MOL
end