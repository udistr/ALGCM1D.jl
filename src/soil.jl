function soil(mtime)

  global ΘA,qA,SW,UA,head,TS,soil_parameters,theta_l,theta_v,MOL
  #set atmosphric values
  println(MOL)
  psi_m = holtslag_psim(ZAC[1]/MOL);
  psi_H = holtslag_psis(ZAC[1]/MOL);

  head1, theta_l1, theta_v1, q_l1, q_v1, rho_vs1, dt2 = metric_head(head, 
                TS,E,constants,soil_parameters,soil_numerical_parameters,dt);  

  TS1,E1,dt2 = Temperature_dist(TS,theta_l1,theta_v1,q_l1,q_v1,rho_vs1, soil_parameters,
    ΘA[1],qA[1],SW,UA[1],psi_H, psi_m,soil_numerical_parameters,constants,atm_parameters,dt2);

if dt2<dt
   dt2_2 = dt2;
      while dt2_2<dt
        head1,theta_l1,theta_v1,q_l1,q_v1,rho_vs1,dt2 = metric_head(head1, 
          TS1,E,constants,soil_parameters,soil_numerical_parameters,dt);  

        TS1,E1,dt2 = Temperature_dist(TS1,theta_l1,theta_v1,q_l1,q_v1,rho_vs1,soil_parameters,ΘA[1],qA[1],SW,UA[1],psi_m,soil_numerical_parameters,constants,atm_parameters,dt2);
          
          dt2_2 = dt2_2+dt2;
          println(dt2)

      end
end     

TS=TS1
head=head1
theta_l=theta_l1
theta_v=theta_v1

end