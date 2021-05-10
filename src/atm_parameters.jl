@with_kw struct aparameters{R}
  z_ref = 2;
  z_H = 1e-3;
  z_oH = 1e-3;
  z_m = 1e-3;
  z_om = 1e-3;
  d = 0;
end

global atm_parameters=aparameters{Float64}()