function surface(mtime)

  global SW,LW,LH,SH,E,TAU,Qnet
  ########################################################################
  # surface
  ########################################################################

  CM=0.01
  CQ=0.01
  CD=0.01

  SH=-ρA[1]*CD*cp*(ΘA[1]-TS[end])
  LH=-ρA[1]*CD*Av*(qA[1]-q_v[end])
  ustar=0.01
  TAU=ρA[1].*ustar.^2 .*sign(ustar);

  E=-LH/Av;

  Qnet=SW-LH-SH-LW;

end