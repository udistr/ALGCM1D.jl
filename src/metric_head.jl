function metric_head(h,T, q_o,constants, soil_parameters, soil_numerical_parameters, dt)

# initializing
h_old = h;
theta_l, theta_v, soil_hydraulic = hydraulic_variables(h, T, constants, soil_parameters);

O=zeros(2)
O[1] = 1; 

while O[1]>soil_numerical_parameters.eps_h
    temp_h_2 = h;
    temp_h = (h+h_old)./2;
    
    for i = 2:soil_numerical_parameters.Ns-1
        h[i] = h_old[i] + dt/(soil_numerical_parameters.dz_s^2*soil_hydraulic.C[i])*
            (soil_hydraulic.K_h[i]*(temp_h[i+1]-2*temp_h[i]+temp_h[i-1]) +
            (soil_hydraulic.K_h[i+1]-soil_hydraulic.K_h[i])*(temp_h[i+1]-temp_h[i]) +
            soil_hydraulic.K_T[i]*(T[i+1]-2*T[i]+T[i-1]) +
            (soil_hydraulic.K_T[i+1]-soil_hydraulic.K_T[i])*(T[i+1]-T[i]) +
            (soil_hydraulic.K[i+1]-soil_hydraulic.K[i])*soil_numerical_parameters.dz_s);
    end
    h[1] = h_old[1]+dt/(soil_numerical_parameters.dz_s*soil_hydraulic.C[1])*(soil_hydraulic.K[2]-soil_hydraulic.K[1]);
    h[end] = h_old[end]-dt/(soil_hydraulic.C[end]*soil_numerical_parameters.dz_s^2)*(soil_hydraulic.K_h[end]*(temp_h[end]-temp_h[end-1])+
        soil_hydraulic.K_T[end]*(T[end]-T[end-1])+(soil_hydraulic.K[end]+q_o)*soil_numerical_parameters.dz_s);
        
    
    h[h.<soil_parameters.h_max] .= soil_parameters.h_max;
    
    theta_l, theta_v, soil_hydraulic, rho_vs = hydraulic_variables(h, T, constants, soil_parameters);
    
    O[2] = mean(abs.(temp_h_2-h)./abs.(temp_h));
    
    if O[2]<O[1]
        O[1] = O[2];
    else
        dt = dt/2;
        h = h_old;
        theta_l, theta_v, soil_hydraulic = hydraulic_variables(h, T, constants, soil_parameters);        
        O[1] = 1;
    end
    #println(dt)
end

q_l[1] = -soil_hydraulic.K[1];
q_v[1] = 0;

for i = 2:soil_numerical_parameters.Ns   
    q_l[i] = -soil_hydraulic.K_lT[i]*(T[i]-T[i-1])/soil_numerical_parameters.dz_s-
        soil_hydraulic.K[i]*((h[i]-h[i-1])/soil_numerical_parameters.dz_s+1);    
    q_v[i] = -soil_hydraulic.K_vT[i]*(T[i]-T[i-1])/soil_numerical_parameters.dz_s-
        soil_hydraulic.K_vh[i]*((h[i]-h[i-1])/soil_numerical_parameters.dz_s);    
end


return h, theta_l, theta_v, q_l, q_v, rho_vs, dt

end
