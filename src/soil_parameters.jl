@with_kw struct sparameters{R}
  #thermal constants
  b1::R = 0.243;
  b2::R = 0.393;
  b3::R = 1.534;
  #retention - VG model
  alpha::R = 2;
  n::R = 1.5;
  m::R = 1-1/n;
  theta_s::R = 0.4;
  theta_r::R = 0.01;
  Ks::R = 1/(24*60*60);
  #matric head top constrain
  h_max::R = -1e20;
  #additional parameters
  G_wT::R = 7; # gain factor [-] taken from Jiangbo et al. (2017)
  fc::R = 0.02; #taken from Jiangbo et al. (2017)
end

global soil_parameters=sparameters{Float64}()


