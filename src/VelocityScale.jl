function VelocityScale(z,h,L,ustar,θv=nothing,wθv=nothing)
  xi=z/L
  psis=holtslag_psis(xi)
  psim=holtslag_psim(xi)
  if xi>0
      wt=ustar/psis
      wm=wt
  else
      if z/h<0.1
          wt=ustar/psis
          wm=ustar/psim
      else
          a=7.2
          xi1=0.1*h/L
          psis1=holtslag_psis(xi1)
          psim1=holtslag_psim(xi1)
          #wθv<0 ? println("wθv",wθv) : nothing
          #println("θv",θv)
          wstar=((gravity_mks/θv)*wθv*h)^(1/3)
          # this part is different than Holtslag et. al to ensure continueos 
          # velositiy scale at the surface and the mixed layer above.
          ct=(ustar/wstar)^3*(1/psis1^3-1)
          cm=(ustar/wstar)^3*(1/psim1^3-1)
          wm=(ustar^3+cm*wstar^3)^(1/3)
          wt=(ustar^3+ct*wstar^3)^(1/3)
          #Pr=max(min(psis1/psim1*abs(xi1)+a*karman*(0.1*(wstar/wm)),4),1)
          #wt=wm/Pr
      end
  end    
  return wt,wm
end