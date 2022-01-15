function moist(ΘA,qA,ρA,PA)

  TA=(ΘA).*(PA./PA0).^(2/7).-d2k;
  qsat=saltsat*cvapor_fac*exp.(-cvapor_exp./(TA.+d2k))./ρA;
  dq=(qA.-qsat).*((qA.-qsat).>0);
  dq=min.(dq,qA);
  # change of temperature due to condensation
  dΘA1=Av./cpa.*(PA0./PA).^(2/7).*dq;
  # change of temperature due to change in staurated vapor pressure
  dΘA2=Av.^2 ./Rv ./cpa.*(PA0./PA).^(2*2/7).*qsat./ΘA.^2 .*dq;
  #println("dΘA")
  #println(dΘA1)
  #println("============")
  dΘA=dΘA1+dΘA2; #[J/kg][j-1 kg K]
  #println([i maximum(dTHETA) maximum(dq)])
  ind=qA.>qsat;
  qA[ind]=qA[ind]-dq[ind];
  ΘA[ind]=ΘA[ind]+dΘA[ind];
  
  return ΘA,qA
end