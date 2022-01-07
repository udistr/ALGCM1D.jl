function Temperature_dist(T,theta_l,theta_v,q_l,q_v,rho_vs,soil_parameters,
  Ta, Hr_atm, St,u, psi_H, psi_m, soil_numerical_parameters,constants,atm_parameters,dt)

theta_top = theta_l[end]; # top layer water content
rs = 10*exp(35.63*(0.15-theta_top)); # soil surface resistance to vapor flow
rv = 1/(u*constants.k^2)*(log((atm_parameters.z_ref-atm_parameters.d+atm_parameters.z_H)/atm_parameters.z_oH)+psi_H)*
    (log((atm_parameters.z_ref-atm_parameters.d+atm_parameters.z_m)/atm_parameters.z_om)+psi_m);
rH = rv;

e_a = 0.611*exp(17.27*(Ta-273.15)/(Ta-35.85))*Hr_atm; # atm. vapor pressure
rho_sa = 1e-3.*exp(19.84-4975.9./Ta);
rho_va = rho_sa*Hr_atm;
eps_a = 0.7+5.95e-5*e_a*exp(1500/Ta); # atm. emmisivity;
eps_s = min(0.9+0.18*theta_top,1); #soil emmisivity

Cp = constants.Cs.*(soil_parameters.theta_s.-theta_l)+constants.Cw.*theta_l.+constants.Cv.*theta_v;
lambda = soil_parameters.b1.+soil_parameters.b2.*theta_l+soil_parameters.b3.*sqrt.(theta_l);
rho_w = 1000 .-7.3e-3.*(T.-(273+4)).^2 .+ 3.79e-5.*(T.-(273+4)).^3;  #water density [kg m^-3]
Lw = 2.501e6 .- 2369.2.*(T .- 273); # latent heat [J kg^-1]
Lo = Lw.*rho_w; # vol. latent heat [J m^-3]
L = 2.260e6; # latent heat for vaporazation [J kg^-1]

T_old = T;
O=zeros(2)
O[1] = 1;
LE=0
while O[1]>soil_numerical_parameters.eps_T 
    temp_T_2 = T;
    temp_T = (T+T_old)./2;
    for i = 2:soil_numerical_parameters.Ns-1
        T[i] = T_old[i]+dt/(soil_numerical_parameters.dz_s^2*Cp[i])*(lambda[i]*(temp_T[i+1]-2*temp_T[i]+temp_T[i-1])+
            (lambda[i+1]-lambda[i])*(temp_T[i+1]-temp_T[i])-(constants.Cw*q_l[i]+constants.Cv*q_v[i])*soil_numerical_parameters.dz_s*(temp_T[i+1]-temp_T[i]));
    end
    #set bottom boundary
    T[1] = T_old[1];
    #set top boundary
    Rn = (1-albedo)*St+eps_s*eps_a*constants.sigma*Ta^4-
        eps_s*constants.sigma*temp_T[end]^4;
    H = constants.Ca*(temp_T[end]-Ta)/rH;
    LE = L.*(rho_vs[end].-rho_va)./(rv.+rs);
    G = Rn.-H.-LE;
    T[end] = (G+T[end-1]*lambda[end]/soil_numerical_parameters.dz_s+Lo[end]*q_v[end])/(lambda[end]/soil_numerical_parameters.dz_s-constants.Cv*q_v[end]-constants.Cw*q_l[end]);
    
    rho_w = 1000 .- 7.3e-3.*(T.-(273+4)).^2 .+ 3.79e-5.*(T.-(273+4)).^3;  #water density [kg m^-3]
    Lw = 2.501e6.-2369.2.*(T.-273); # latent heat [J kg^-1]
    Lo = Lw.*rho_w; # vol. latent heat [J m^-3]
    
    O[2] = mean(abs.(temp_T_2.-T)./temp_T);
    
    if O[2]<O[1]
        O[1] = O[2];
    else
        dt = dt/2;
        T = T_old;
        O[1] = 1;
    end
    if dt==0
        dt = 1;
    end

    #println(dt)
end

E = LE/(L*rho_w[end]);
d
return T,E,dt

end
