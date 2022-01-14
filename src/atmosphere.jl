function atmosphere(mtime)

  global Θv,TA,TS,ρA,ustar,LH,SH,E,LW,TAU,Qnet,SW
  global KAm,KAt
  global UA1,qA1,ΘA1
  global qA,ΘA,UA

  ########################################################################
  # Atmospheric explicit time stepping
  ########################################################################
  TAUS=sign(UA[1]);

  UAst[1]=UA[1]-TAU[1]*TAUS/(ΔZA*ρA[1])*ΔT-(1-Aimp)*(
      KAm[1].*(UA[1]-UA[2])/ΔZA^2*ΔT);
  UAst[end]=UA[end]-(1-Aimp)*(
      KAm[end]*(UA[end]-UAup)/ΔZA^2*ΔT);
  UAst[2:end-1]=UA[2:end-1]+(1-Aimp)*(
      KAm[1:end-2].*(UA[1:end-2]-UA[2:end-1])/ΔZA^2*ΔT-
      KAm[2:end-1].*(UA[2:end-1]-UA[3:end  ])/ΔZA^2*ΔT);

  qAst[1]=qA[1]+EVAP[1]/(ΔZA*ρO)*ΔT-
      (1-Aimp)*(KAt[1]*(qA[1]-qA[2])/ΔZA^2*ΔT+
      (-γcq[1].*KAt[1])/ΔZA*ΔT);
  qAst[end]=qA[end]-(1-Aimp)*(KAt[end]*(qA[end]-q0)/ΔZA^2*ΔT+
      (γcq[end-1].*KAt[end-1]-γcq[end].*KAt[end])/ΔZA*ΔT);
  qAst[2:end-1]=qA[2:end-1]+(1-Aimp)*(
      KAt[1:end-2].*(qA[1:end-2]-qA[2:end-1])/ΔZA^2*ΔT-
      KAt[2:end-1].*(qA[2:end-1]-qA[3:end  ])/ΔZA^2*ΔT+
      (γcq[1:end-2].*KAt[1:end-2]-γcq[2:end-1].*KAt[2:end-1])/ΔZA*ΔT);

  ΘAst[1]=ΘA[1]+(SH[1])/cpa/(ΔZA*ΘA[1])*ΔT-
        (1-Aimp)*(KAt[1]*(ΘA[1]-ΘA[2])/ΔZA^2*ΔT+
        (-γct[1].*KAt[1])/ΔZA*ΔT);
  ΘAst[end]=ΘA[end]-
        (1-Aimp)*(KAt[end]*(ΘA[end]-ΘA0)/ΔZA^2*ΔT+
        (γct[end-1].*KAt[end-1]-γct[end].*KAt[end])/ΔZA*ΔT);
  ΘAst[2:end-1]=ΘA[2:end-1]+(1-Aimp)*(
        KAt[1:end-2].*(ΘA[1:end-2]-ΘA[2:end-1])/ΔZA^2*ΔT-
        KAt[2:end-1].*(ΘA[2:end-1]-ΘA[3:end  ])/ΔZA^2*ΔT+
        (γct[1:end-2].*KAt[1:end-2]-γct[2:end-1].*KAt[2:end-1])/ΔZA*ΔT);

  if Aimp==0
      UA1=UAst; qA1=qAst; ΘA1=ΘAst;
  else
      du=-Aimp*KAm[1:end-1].*ΔT./ΔZA^2; dl=-Aimp*KAm[2:end].*ΔT./ΔZA^2;
      d=(1 .+Aimp.*ΔT./ΔZA^2 .*([KAm[1:end-1]; 0]+[0; KAm[2:end]]));
      A=inv(Tridiagonal(dl,d,du));
      UA1=A*UAst;
      du=-Aimp*KAt[1:end-1].*ΔT./ΔZA^2; dl=-Aimp*KAt[2:end].*ΔT./ΔZA^2;
      d=(1 .+Aimp.*ΔT./ΔZA^2 .*([KAt[1:end-1]; 0]+[0; KAt[2:end]]));
      A=inv(Tridiagonal(dl,d,du));
      qA1=A*qAst; ΘA1=A*ΘAst;
  end

  TA=(ΘA1).*(PA./PA[1]).^(2/7).-d2k;
  ρA=PA./(Rgas.*(TA.+d2k));

  if moist==1
      qsat=saltsat*cvapor_fac*exp.(-cvapor_exp./(TA.+d2k))./ρA;
      dq=(qA.-qsat).*((qA1.-qsat).>0);
      dq=min.(dq,qA1);
      # change of temperature due to condensation
      dΘA1=Av./cpa.*(PA[1]./PA).^(2/7).*dq;
      # change of temperature due to change in staurated vapor pressure
      dΘA2=Av.^2 ./Rv ./cpa.*(PA[1]./PA).^(2*2/7).*qsat./ΘA1.^2 .*dq;
      #println("dΘA")
      #println(dΘA1)
      #println("============")
      dΘA=dΘA1+dΘA2; #[J/kg][j-1 kg K]
      #println([i maximum(dTHETA) maximum(dq)])
      ind=qA1.>qsat;
      qA1[ind]=qA1[ind]-dq[ind];
      ΘA1[ind]=ΘA1[ind]+dΘA[ind];
      TA=(ΘA1).*(PA./PA[1]).^(2/7).-d2k;
  end

  qA=qA1; ΘA=ΘA1; UA=UA1;
end