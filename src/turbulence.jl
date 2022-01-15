
function turbulence(ΘA,qA,UA,VA,LH,SH,ustar,ρA,ZAF)

  #global Θv,KAm,KAt,γcq,γct,γcm
  ########################################################################
  # Atmospheric vertical diffusion
  ########################################################################
  Θv=ΘA.*(1 .+humid_fac*qA);
  KAm,KAt,γcq,γct,γcm=holtslag(Θv,UA,UA.*0,LH[1],SH[1],ustar,ρA[1],ZAF);
  #println("----",maximum(γct),"----")

  return KAm,KAt,γcq,γct,γcm

end