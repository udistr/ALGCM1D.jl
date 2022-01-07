function surface(mtime)

  global SW,LW,LH,SH,E,TAU,Qnet,MOL
  ########################################################################
  # surface
  ########################################################################

  humid_fac=0.606;
  
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

  CD=karman^2/log(ZAC[1]/zm)^2
  CH=sqrt(CD)*karman/log(ZAC[1]/zh)
  CE=sqrt(CD)*karman/log(ZAC[1]/ze)

  #kinematic surface flux
  ustar=-sqrt(CD)*abs(UA[1])*UA[1]
  tstar=-sqrt(CH)*abs(UA[1])*(ΘA[1]-TS[end])
  qstar=-sqrt(CE)*abs(UA[1])*(qA[1]-q_v[end])

  L=0;

  for i=1:5
    wθv=tstar*(1+humid_fac*qA[1])+humid_fac*ΘA[1]*qstar
    L=-ΘA[1]*ustar^3/(karman*gravity_mks*wθv)
    
    #CD=karman^2/(log(ZAC[1]/zm)-holtslag_psim(ZAC[1]/L))
    #CH=CD/(log(ZAC[1]/zh)-holtslag_psis(ZAC[1]/L))
    #CE=CD/(log(ZAC[1]/ze)-holtslag_psis(ZAC[1]/L))

    CD=CD/(1+CD*(log(ZAC[1]/zm)-holtslag_psim(ZAC[1]/L)/karman))
    CH=CH/(1+CH*(log(ZAC[1]/zm)-holtslag_psis(ZAC[1]/L)/karman))
    CE=CE/(1+CE*(log(ZAC[1]/zm)-holtslag_psis(ZAC[1]/L)/karman))

    ustar=-sqrt(CD)*abs(UA[1])*UA[1]
    tstar=-sqrt(CH)*abs(UA[1])*(ΘA[1]-TS[end])
    qstar=-sqrt(CE)*abs(UA[1])*(qA[1]-q_v[end])

  end

  SH=-ρA[1]*CH*cp*abs(UA[1])*(ΘA[1]-TS[end])
  LH=-ρA[1]*CE*Av*abs(UA[1])*(qA[1]-q_v[end])
  TAU=ρA[1]*CD*abs(UA[1])*UA[1];

  MOL=L;
  E=-LH/Av;

  Qnet=SW-LH-SH-LW;

end