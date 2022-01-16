function soil(TS,head,EVAP,dt)

  #global head,TS,soil_parameters,theta_l,theta_v,EVAP,rho_w,q_l,q_v

  rho_w = 1000 .-7.3e-3.*(TS.-(273+4)).^2 .+ 3.79e-5.*(TS.-(273+4)).^3;  #water density [kg m^-3]
  q_o = EVAP[1]./rho_w[end];

  head1, theta_l1, theta_v1, q_l1, q_v1, rho_vs1, dt2 = metric_head(head, 
                TS,q_o,constants,soil_parameters,soil_numerical_parameters,dt);  

  TS1,dt2 = Temperature_dist(TS,theta_l1,theta_v1,q_l1,q_v1,rho_w,Qnet[1],soil_parameters,soil_numerical_parameters,constants,dt2);
 # need to change qA to relative humidity
 if dt2<dt
  dt2_2 = dt2;
      while dt2_2<dt
        head1,theta_l1,theta_v1,q_l1,q_v1,rho_vs1,dt2 = metric_head(head1, 
          TS1,q_o,constants,soil_parameters,soil_numerical_parameters,dt);  

        TS1,dt2 = Temperature_dist(TS1,theta_l1,theta_v1,q_l1,q_v1,rho_vs1,Qnet,soil_parameters,soil_numerical_parameters,constants,dt2);
          
          dt2_2 = dt2_2+dt2;
          #println(dt2)

      end
  end     

TS=TS1
head=head1
theta_l=theta_l1
theta_v=theta_v1
rho_vs=rho_vs1

return TS,head,theta_l,theta_v,rho_vs

end