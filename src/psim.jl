function psim(xi)
if (xi>=0)
    ψm = -4.7*xi;
else
    x=(1-15*xi)^(0.25)
    ψm=2*log((1+x)/2)+log((1+x^2)/2)-2*atan(x)+pi/2
end
return ψm
end