function psis(xi)
  if (xi>=0)
      ψs = -4.7*xi;
  else
      x=(1-15*xi)^(0.25)
      ψs=2*log((1+x^2)/2)
  end
  return ψs
  end