
function turbulence(mtime)

  global Θv,KAm,KAt,γcq,γct,γcm
  ########################################################################
  # Atmospheric vertical diffusion
  ########################################################################
  Θv=ΘA.*(1 .+humid_fac*qA);
  KAm,KAt,γcq,γct,γcm=holtslag(Θv,UA,UA.*0,LH,SH,ustar,ρA[1],ZAF);
  #println("----",maximum(γct),"----")


end