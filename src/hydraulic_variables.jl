function hydraulic_variables(h, T, constants, soil_parameters)

  d_gamma = -0.0001425 .- 4.76e-7 .*(T.-273);
  rho_w = 1000 .- 7.3e-3.*(T.-(273+4)).^2+3.79e-5.*(T.-(273+4)).^3;  #water density [kg m^-3]
  D_a = 2.12e-5.*(T./273).^2; # air diffusivity
  rho_sv = 1e-3.*exp.(19.84 .- 4975.9 ./T); # saturated vapor pressure [kg/m^3]
  drho_vs = 4.9759./(T.^2).*exp.(19.84 .-4975.9 ./T); #saturated vapor pressure derivative [kg m^-3 K^-1]

  Se = 1. ./(1 .+(soil_parameters.alpha.*abs.(h)).^soil_parameters.n).^soil_parameters.m;
  Se[h.>0] .= 1;
  #Se(Se<1e-3) = 1e-3;
  theta_l = soil_parameters.theta_r.+(soil_parameters.theta_s-soil_parameters.theta_r).*Se;
  theta_a = soil_parameters.theta_s.-theta_l;
  Hr = exp.(h.*constants.M.*constants.g./(constants.Rg.*T));
  Hr[h.>0] .= 1;
  rho_vs = Hr.*rho_sv;
  theta_v = rho_vs.*theta_a./rho_w;
      
  Cm = -(soil_parameters.alpha.*abs.(h)).^soil_parameters.n.*soil_parameters.m.*soil_parameters.n.*
      (soil_parameters.theta_s.-soil_parameters.theta_r).*((soil_parameters.alpha.*abs.(h)).^soil_parameters.n.+1).^(.-soil_parameters.m.-1)./h;
  Cm[h.>=0] .= 1e-13;
  Cv = rho_sv./rho_w.*(soil_parameters.theta_s.-soil_parameters.theta_r).*(constants.g.*constants.M./(constants.Rg.*T).*
      Hr.*(1 .-Se).-Cm.*Hr./(soil_parameters.theta_s.-soil_parameters.theta_r));    
  Cv[h.>0] .= 0;


  tau = theta_a.^(7/3)./soil_parameters.theta_s.^2;
  D = tau.*theta_a.*D_a;
  eta = 9.5 .+ 3 .*theta_l./soil_parameters.theta_s.-8.5.*exp.(-((1 .+2.6 ./sqrt.(soil_parameters.fc)).*(theta_l./soil_parameters.theta_s)).^4);
    
  #hydraulic parameters
  tau = theta_a.^(7/3)./soil_parameters.theta_s.^2;
  D = tau.*theta_a.*D_a;
  eta = 9.5 .+ 3 .* theta_l./soil_parameters.theta_s .- 8.5 .* exp.(-((1+2.6/sqrt(soil_parameters.fc)).*(theta_l./soil_parameters.theta_s)).^4);
      
  soil_hydraulicC = Cm+Cv;
  soil_hydraulicK = soil_parameters.Ks.*sqrt.(Se).*(1 .-(1 .-Se.^(1/soil_parameters.m)).^soil_parameters.m).^2;
  soil_hydraulicK_vh = D./(constants.Rg.*T.*rho_w).*rho_sv.*constants.M.*constants.g.*Hr;
  soil_hydraulicK_h = soil_hydraulicK+soil_hydraulicK_vh;
      
  soil_hydraulicK_lT = soil_hydraulicK.*h.*d_gamma./(soil_parameters.G_wT.*constants.gamma_o);
  soil_hydraulicK_vT = D.*eta.*Hr.*drho_vs./rho_w;
  soil_hydraulicK_T = soil_hydraulicK_lT + soil_hydraulicK_vT;

  soil_hydraulic=shydraulic(soil_hydraulicC,soil_hydraulicK,soil_hydraulicK_vh,soil_hydraulicK_h,
                    soil_hydraulicK_lT,soil_hydraulicK_vT,soil_hydraulicK_T)

  rho_vs = rho_vs[end];

  return theta_l, theta_v, soil_hydraulic, rho_vs

end
